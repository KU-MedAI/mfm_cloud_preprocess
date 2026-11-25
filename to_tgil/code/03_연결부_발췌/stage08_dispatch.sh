# 검토용 원문 발췌임. 단독 실행 모듈이 아니며 import/호출 문맥은 전체 전달본을 사용함.

# 원본: run_wgs_sample.sh:72-76
if sample_hap_realign_required "$S"; then
  pipeline_run_stage 08 08_realign_wgs.sh
else
  pipeline_skip_stage 08 "$(sample_output_mode "$S"): 개인 HAP WGS 재정렬 비적용"
fi
