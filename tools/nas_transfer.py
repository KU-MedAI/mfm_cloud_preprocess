#!/usr/bin/env python3
"""원본 서버의 수동 NAS 이관 명령 진입점임. conda jj_mfm에서 실행함."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src'),str(ROOT/'vendor/MFM-DNARNAFM')]
from handoff.nas_transfer import main
if __name__=='__main__':main()
