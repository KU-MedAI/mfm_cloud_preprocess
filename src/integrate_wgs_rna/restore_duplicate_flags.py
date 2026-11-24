"""BWA-MEM2 SAM stream에서 운반용 OD tag를 source duplicate flag로 복원한다."""
from __future__ import annotations

import argparse
import json
import os
import sys

import pysam


def restore_stream_pysam(mode: str) -> dict[str, int]:
    counts = {"records": 0, "od_records": 0, "duplicate_flags_restored": 0}
    with pysam.AlignmentFile("-", "r") as source:
        with pysam.AlignmentFile("-", "wb", header=source.header) as destination:
            for read in source.fetch(until_eof=True):
                counts["records"] += 1
                carried = read.has_tag("OD") and int(read.get_tag("OD")) == 1
                if carried:
                    counts["od_records"] += 1
                    if mode == "preserve_source":
                        read.is_duplicate = True
                        counts["duplicate_flags_restored"] += 1
                    read.set_tag("OD", None)
                destination.write(read)
    return counts


def restore_stream_sam(mode: str) -> dict[str, int]:
    """BWA의 SAM 바이트를 유지하고 OD=1 레코드만 수정한다."""
    counts = {"records": 0, "od_records": 0, "duplicate_flags_restored": 0}
    source, destination = sys.stdin.buffer, sys.stdout.buffer
    header = True
    for line in source:
        if header and line.startswith(b"@"):
            destination.write(line)
            continue
        header = False
        counts["records"] += 1
        if b"\tOD:" not in line:
            destination.write(line)
            continue
        ending = b"\r\n" if line.endswith(b"\r\n") else b"\n" if line.endswith(b"\n") else b""
        fields = line[:-len(ending)].split(b"\t") if ending else line.split(b"\t")
        if len(fields) < 11:
            raise ValueError("SAM 레코드의 필수 열이 11개보다 적습니다.")
        kept = fields[:11]
        carried = False
        for field in fields[11:]:
            if field.startswith(b"OD:") and int(field.split(b":", 2)[2]) == 1:
                carried = True
            else:
                kept.append(field)
        if carried:
            counts["od_records"] += 1
            if mode == "preserve_source":
                kept[1] = str(int(kept[1]) | 0x400).encode()
                counts["duplicate_flags_restored"] += 1
        destination.write(b"\t".join(kept) + ending)
    destination.flush()
    return counts


def restore_stream(mode: str) -> dict[str, int]:
    engine = os.environ.get("MFM_DUP_RESTORE_ENGINE", "sam")
    if engine == "sam":
        return restore_stream_sam(mode)
    if engine == "pysam":
        return restore_stream_pysam(mode)
    raise ValueError("MFM_DUP_RESTORE_ENGINE은 sam 또는 pysam이어야 합니다.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("preserve_source", "none"), required=True)
    args = parser.parse_args()
    counts = restore_stream(args.mode)
    print(json.dumps(counts, sort_keys=True), file=sys.stderr)


if __name__ == "__main__":
    main()
