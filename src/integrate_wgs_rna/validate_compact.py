from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import zarr

HAP_LABELS = ("HAP1", "HAP2")


@dataclass
class CompactValidationConfig:
    compact_contig_dir: Path
    read_coordinate_mode: str = "personal_realigned"
    validation_level: str = "standard"
    output_json: Path | None = None


def run_compact_validation(config: CompactValidationConfig) -> dict[str, Any]:
    root = config.compact_contig_dir
    metadata = json.loads((root/'contig_metadata.json').read_text())
    if metadata.get('wgs_payload')=='pileup_feature':
        from .pileup_features.export import validate_contig
        result=validate_contig(root)
        report=dict(status='pass',wgs_payload='pileup_feature',by_haplotype=result)
        if config.output_json:config.output_json.write_text(json.dumps(report))
        return report
    if config.read_coordinate_mode not in {"personal_realigned", "reference_projected"}:
        raise ValueError("read_coordinate_mode must be personal_realigned or reference_projected")
    if config.validation_level not in {"fast", "standard", "full"}:
        raise ValueError("validation_level must be fast, standard, or full")
    required = [root / "contig_metadata.json", root / "windows.tsv", root / "hap_sequence.zarr", root / "depth.zarr"]
    metadata = json.loads((root / "contig_metadata.json").read_text(encoding="utf-8"))
    if config.read_coordinate_mode == "personal_realigned":
        required.extend([root / "readset.HAP1.zarr", root / "readset.HAP2.zarr"])
    else:
        required.append(root / "readset.zarr")
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("missing compact outputs: " + ", ".join(missing))

    hap_root = zarr.open_group(str(root / "hap_sequence.zarr"), mode="r")
    depth_root = zarr.open_group(str(root / "depth.zarr"), mode="r")
    supplementary_path = root / "supplementary_depth.zarr"
    configured_supplementary = tuple(metadata.get("filter_policy", {}).get("supplementary_depth_segment_types", ()))
    if configured_supplementary and not supplementary_path.exists():
        raise RuntimeError("supplementary_depth.zarr is required by filter_policy")
    supplementary_root = zarr.open_group(str(supplementary_path), mode="r") if supplementary_path.exists() else None

    hap_results = {}
    for hap in HAP_LABELS:
        if hap not in hap_root or hap not in depth_root:
            raise RuntimeError(f"missing {hap} in hap_sequence/depth zarr")
        hap_group = hap_root[hap]
        depth_group = depth_root[hap]
        for name in ("bases_packed", "offsets", "lengths", "raw_start", "actual_length"):
            if name not in hap_group:
                raise RuntimeError(f"missing hap_sequence/{hap}/{name}")
        for name in ("depth", "offsets", "lengths"):
            if name not in depth_group:
                raise RuntimeError(f"missing depth/{hap}/{name}")
        seq_offsets = hap_group["offsets"][:]
        seq_lengths = hap_group["lengths"][:]
        depth_offsets = depth_group["offsets"][:]
        depth_lengths = depth_group["lengths"][:]
        if depth_group["depth"].dtype != np.dtype("uint32"):
            raise RuntimeError(f"depth/{hap}/depth dtype is {depth_group['depth'].dtype}, expected uint32")
        if len(seq_offsets) != len(depth_offsets):
            raise RuntimeError(f"{hap} sequence/depth offset length mismatch")
        if len(seq_lengths) != len(depth_lengths):
            raise RuntimeError(f"{hap} sequence/depth window length mismatch")
        if len(seq_offsets) != len(seq_lengths) + 1:
            raise RuntimeError(f"{hap} invalid sequence offsets length")
        if len(depth_offsets) != len(depth_lengths) + 1:
            raise RuntimeError(f"{hap} invalid depth offsets length")
        if int(seq_offsets[0]) != 0 or int(depth_offsets[0]) != 0:
            raise RuntimeError(f"{hap} offsets must start at 0")
        if np.any(np.diff(seq_offsets) < 0) or np.any(np.diff(depth_offsets) < 0):
            raise RuntimeError(f"{hap} offsets are not monotonic")
        if np.any(seq_lengths < 0) or np.any(depth_lengths < 0):
            raise RuntimeError(f"{hap} negative window length")
        if np.any(seq_lengths != depth_lengths):
            raise RuntimeError(f"{hap} sequence/depth lengths differ")
        if np.any(depth_offsets[1:] - depth_offsets[:-1] != depth_lengths):
            raise RuntimeError(f"{hap} depth offset deltas differ from lengths")
        if supplementary_root is not None:
            if hap not in supplementary_root:
                raise RuntimeError(f"missing {hap} in supplementary_depth zarr")
            supplementary_group = supplementary_root[hap]
            for name in ("depth", "offsets", "lengths"):
                if name not in supplementary_group:
                    raise RuntimeError(f"missing supplementary_depth/{hap}/{name}")
            supplementary_offsets = supplementary_group["offsets"][:]
            supplementary_lengths = supplementary_group["lengths"][:]
            if supplementary_group["depth"].dtype != np.dtype("uint32"):
                raise RuntimeError(f"supplementary_depth/{hap}/depth dtype is {supplementary_group['depth'].dtype}, expected uint32")
            if len(supplementary_offsets) != len(seq_offsets) or len(supplementary_lengths) != len(seq_lengths):
                raise RuntimeError(f"{hap} supplementary depth window shape mismatch")
            if int(supplementary_offsets[0]) != 0 or np.any(np.diff(supplementary_offsets) < 0):
                raise RuntimeError(f"{hap} supplementary depth offsets are invalid")
            if np.any(supplementary_lengths != seq_lengths) or np.any(supplementary_offsets[1:] - supplementary_offsets[:-1] != supplementary_lengths):
                raise RuntimeError(f"{hap} supplementary depth lengths differ from sequence")

        readset_path = root / f"readset.{hap}.zarr" if config.read_coordinate_mode == "personal_realigned" else root / "readset.zarr"
        readset = zarr.open_group(str(readset_path), mode="r")
        segments = readset["segments"][hap]
        if len(segments["window_read_offsets"]) != len(seq_lengths) + 1:
            raise RuntimeError(f"{hap} segment window_read_offsets length mismatch")
        hap_results[hap] = {
            "window_count": int(len(seq_lengths)),
            "sequence_packed_bytes": int(hap_group["bases_packed"].shape[0]),
            "depth_bases": int(depth_group["depth"].shape[0]),
            "segment_count": int(segments["read_id"].shape[0]),
            "readset": str(readset_path),
            "supplementary_depth_bases": int(supplementary_root[hap]["depth"].shape[0]) if supplementary_root is not None else None,
        }
    payload = {
        "schema": "mfm_integrated_compact_v2_validation",
        "status": "pass",
        "compact_contig_dir": str(root),
        "read_coordinate_mode": config.read_coordinate_mode,
        "validation_level": config.validation_level,
        "haplotypes": hap_results,
        "config": {key: str(value) if isinstance(value, Path) else value for key, value in asdict(config).items()},
    }
    if config.output_json is not None:
        config.output_json.parent.mkdir(parents=True, exist_ok=True)
        config.output_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload
