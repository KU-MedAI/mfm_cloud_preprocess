#!/usr/bin/env python3
"""Stage 02 RNA 입력 계약 검사와 변경 감지용 fingerprint 생성."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path


def fail(message: str) -> None:
    raise SystemExit(f"[FATAL] {message}")


def sha256(path: Path, chunk: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def stat_record(path: Path, *, content_hash: bool = False) -> dict:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing or empty input: {path}")
    st = path.stat()
    out = {
        "path": str(path.resolve()), "size": st.st_size, "mtime_ns": st.st_mtime_ns,
        "ctime_ns": st.st_ctime_ns, "inode": st.st_ino,
    }
    if content_hash:
        out["sha256"] = sha256(path)
    return out


def alignment_index(path: Path) -> Path:
    candidates = [Path(str(path) + ".bai"), path.with_suffix(".bai")]
    if path.suffix.lower() == ".cram":
        candidates = [Path(str(path) + ".crai"), path.with_suffix(".crai")]
    for candidate in candidates:
        if candidate.is_file() and candidate.stat().st_size:
            return candidate
    fail(f"alignment index missing: {path}")
    raise AssertionError


def run_text(command: list[str]) -> str:
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode:
        fail(f"command failed ({' '.join(command)}): {result.stderr.strip()}")
    return result.stdout


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--expected-samples", default="", help="쉼표 구분 sample/public alias")
    parser.add_argument("--rna", required=True, type=Path)
    parser.add_argument("--decode-reference", type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--hap-root", type=Path)
    parser.add_argument("--skip-hap-assets", action="store_true")
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--samtools", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--hap-read-mode", required=True)
    parser.add_argument("--min-paired-fraction", type=float, default=0.90)
    parser.add_argument(
        "--ratio-qc-policy", choices=("record_only", "strict"),
        default=os.environ.get("QC_RATIO_POLICY", "record_only"),
        help="record_only records paired fraction without gating; strict reproduces the legacy cutoff",
    )
    parser.add_argument("--require-sm-match", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--sm-placeholders", default="sm1")
    parser.add_argument("--quickcheck-json", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    contigs = [token for token in args.contigs.replace(",", " ").split() if token]
    if not contigs:
        fail("empty RNA contig list")
    if not 0 <= args.min_paired_fraction <= 1:
        fail("--min-paired-fraction must be in [0,1]")

    rna = args.rna.resolve()
    index = alignment_index(rna)
    reference_args: list[str] = []
    if rna.suffix.lower() == ".cram":
        if not args.decode_reference:
            fail("RNA CRAM requires --decode-reference")
        reference = args.decode_reference.resolve()
        stat_record(reference)
        stat_record(Path(str(reference) + ".fai"), content_hash=True)
        reference_args = ["-T", str(reference)]

    quickcheck_reused = False
    if args.quickcheck_json:
        from handoff.input_check import reuse
        try:
            reuse(rna,args.quickcheck_json,args.sample)
            quickcheck_reused = True
        except (ValueError,OSError) as exc:
            # 재개 시 입력을 교체했거나 캐시가 없으면 현재 파일을 다시 검사한다.
            # 아래 전체 flagstat/SM/contig QC는 항상 수행하여 검증을 생략하지 않는다.
            import sys
            print(f"[RNA] quickcheck 캐시 무효 → 현재 입력 재검사: {exc}",file=sys.stderr)
    if not quickcheck_reused:
        quick = subprocess.run(
            [args.samtools, "quickcheck", "-v", str(rna)], text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        if quick.returncode:
            fail(f"samtools quickcheck failed: {(quick.stdout + quick.stderr).strip()}")

    header = run_text([args.samtools, "view", "-H", *reference_args, str(rna)])
    sms: set[str] = set()
    sq: dict[str, dict[str, str]] = {}
    sort_order = ""
    for line in header.splitlines():
        fields = line.split("\t")
        if line.startswith("@HD"):
            sort_order = next((f[3:] for f in fields[1:] if f.startswith("SO:")), "")
        elif line.startswith("@RG"):
            sms.update(f[3:] for f in fields[1:] if f.startswith("SM:"))
        elif line.startswith("@SQ"):
            values = {f[:2]: f[3:] for f in fields[1:] if len(f) > 3 and f[2] == ":"}
            if values.get("SN"):
                sq[values["SN"]] = values
    missing_contigs = sorted(set(contigs) - set(sq))
    if missing_contigs:
        fail(f"RNA header missing contigs: {','.join(missing_contigs)}")

    expected = {x for x in (args.expected_samples or args.sample).split(",") if x}
    placeholders = {x for x in args.sm_placeholders.split(",") if x}
    informative_sms = sms - placeholders
    sm_status = "absent" if not sms else "placeholder" if not informative_sms else "match"
    if informative_sms and not informative_sms.issubset(expected):
        sm_status = "mismatch"
        if args.require_sm_match:
            fail(f"RNA @RG SM mismatch: expected={sorted(expected)} observed={sorted(sms)}")

    # view -c와 view -c -f1은 파일을 두 번 전수 scan한다. flagstat JSON은 같은
    # SAM flag 정의의 total/paired를 한 번의 scan에서 함께 제공한다.
    flagstat_ref: list[str] = []
    if rna.suffix.lower() == ".cram":
        flagstat_ref = ["--input-fmt-option", f"reference={args.decode_reference.resolve()}"]
    flagstat = json.loads(run_text([
        args.samtools, "flagstat", "-@", str(args.threads), "-O", "json",
        *flagstat_ref, str(rna),
    ]))
    categories = [value for value in flagstat.values() if isinstance(value, dict)]
    total = sum(int(value.get("total", 0)) for value in categories)
    primary = sum(int(value.get("primary", 0)) for value in categories)
    paired = sum(int(value.get("paired in sequencing", 0)) for value in categories)
    if total <= 0:
        fail("RNA alignment has no records")
    if primary <= 0:
        fail("RNA alignment has no primary records")
    # RNA BAM에는 multi-mapper의 secondary/supplementary segment가 추가된다. 이들은
    # library가 single-end라는 뜻이 아니므로, paired-end library 계약의 분모는
    # 전체 alignment segment가 아니라 primary sequencing record여야 한다.
    # samtools flagstat의 paired in sequencing은 primary read 기준 집계다.
    if paired > primary:
        fail(f"RNA paired record count exceeds primary record count: paired={paired} primary={primary}")
    paired_fraction = paired / primary
    below_reference = paired_fraction < args.min_paired_fraction
    if args.ratio_qc_policy == "strict" and below_reference:
        fail(f"RNA paired fraction below threshold: {paired_fraction:.6f} < {args.min_paired_fraction:.6f}")

    assets: list[dict] = []
    if not args.skip_hap_assets:
        if not args.hap_root:
            fail("--hap-root is required unless --skip-hap-assets is used")
        for hap in ("HAP1", "HAP2"):
            for contig in contigs:
                base = args.hap_root / hap
                assets.extend([
                    stat_record(base / f"{contig}.fa"),
                    stat_record(base / f"{contig}.fa.fai", content_hash=True),
                    stat_record(base / f"{contig}.chain", content_hash=True),
                ])

    contract = {
        "schema": "mfm_rna_preflight_contract_v2", "sample": args.sample,
        "hap_read_mode": args.hap_read_mode, "contigs": contigs,
        "rna": stat_record(rna), "rna_index": stat_record(index, content_hash=True),
        "gtf": stat_record(args.gtf.resolve(), content_hash=True), "hap_assets": assets,
        "hap_assets_deferred": bool(args.skip_hap_assets),
        "samtools": run_text([args.samtools, "--version"]).splitlines()[0],
        "checks": {
            "quickcheck": "pass", "sort_order": sort_order,
            "read_group_samples": sorted(sms), "sm_status": sm_status,
            "expected_samples": sorted(expected), "header_contigs": len(sq),
            "total_records": total, "primary_records": primary,
            "paired_records": paired, "paired_fraction_denominator": "primary_records",
            "paired_fraction": paired_fraction,
            "min_paired_fraction": args.min_paired_fraction,
            "paired_fraction_reference": args.min_paired_fraction,
            "paired_fraction_below_reference": below_reference,
            "paired_fraction_qc_action": args.ratio_qc_policy,
        },
    }
    encoded = json.dumps(contract, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    contract["contract_sha256"] = hashlib.sha256(encoded).hexdigest()
    contract["status"] = "pass"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(f".{args.output.name}.part.{os.getpid()}")
    temporary.write_text(json.dumps(contract, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, args.output)
    print(json.dumps({"sample": args.sample, "status": "pass", "contract_sha256": contract["contract_sha256"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
