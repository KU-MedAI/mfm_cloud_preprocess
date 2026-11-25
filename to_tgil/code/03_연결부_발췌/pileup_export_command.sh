# 검토용 원문 발췌임. 단독 실행 모듈이 아니며 import/호출 문맥은 전체 전달본을 사용함.

# 원본: 10_compact_zarr.sh:16-24
if sample_pileup_feature "$S"; then
  export PYTHONPATH="$PWD:$PWD/src:$PWD/vendor/MFM-DNARNAFM:${PYTHONPATH:-}"
  "$PYTHON" -m integrate_wgs_rna.pileup_features.export \
    --root "$(sample_final_root "$S")/$S" --alignment "$(sample_haplotagged_alignment "$S")" \
    --reference "$(sample_sequence_reference "$S")" --alignment-reference "$(sample_decode_reference "$S")" \
    --variants "$(sample_persample_vcf_dir "$S")/$S.seq.vcf.gz" --variant-sample "$(sample_source_sample_id "$S")" \
    --read-type "$(sample_read_type "$S")" --min-mapq "$MIN_MAPQ" --contigs "$(chroms_space "$CHROMS")"
  exit 0
fi
