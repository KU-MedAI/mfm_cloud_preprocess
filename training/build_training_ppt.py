#!/usr/bin/env python3
"""MFM v2 전처리 파이프라인 교육용 PPTX를 생성한다."""
from __future__ import annotations

import csv
import math
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_CONNECTOR, MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches, Pt


ROOT = Path(__file__).resolve().parent
ASSETS = ROOT / "assets"
OUT = ROOT / "전처리_파이프라인_교육자료.pptx"
NOTES_OUT = ROOT / "교육자료_슬라이드원고.md"
MAP_OUT = ROOT / "교육자료_근거맵.tsv"

SW, SH = 13.333, 7.5
FONT = "Noto Sans CJK KR"
MONO = "DejaVu Sans Mono"

NAVY = "0D1B2A"
NAVY2 = "16283D"
BLUE = "2F6BFF"
CYAN = "13B8A6"
ORANGE = "F59E0B"
PURPLE = "7C5CE7"
RED = "D94F4F"
GREEN = "258B68"
INK = "1D2733"
MUTED = "65758B"
LIGHT = "F4F7FB"
WHITE = "FFFFFF"
LINE = "D9E2EC"
PALE_BLUE = "EAF0FF"
PALE_CYAN = "E5F7F4"
PALE_ORANGE = "FFF3D9"
PALE_RED = "FCEBEC"
PALE_PURPLE = "F0ECFF"


def rgb(value: str) -> RGBColor:
    return RGBColor.from_string(value)


prs = Presentation()
prs.slide_width = Inches(SW)
prs.slide_height = Inches(SH)
blank = prs.slide_layouts[6]
notes_rows: list[dict[str, object]] = []


def set_bg(slide, color: str = WHITE) -> None:
    fill = slide.background.fill
    fill.solid()
    fill.fore_color.rgb = rgb(color)


def add_text(slide, text: str, x: float, y: float, w: float, h: float, *,
             size: float = 18, color: str = INK, bold: bool = False,
             align: PP_ALIGN = PP_ALIGN.LEFT, font: str = FONT,
             valign: MSO_ANCHOR = MSO_ANCHOR.TOP, margin: float = 0.06,
             line_spacing: float = 1.05):
    box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    frame = box.text_frame
    frame.clear()
    frame.word_wrap = True
    frame.margin_left = frame.margin_right = Inches(margin)
    frame.margin_top = frame.margin_bottom = Inches(margin)
    frame.vertical_anchor = valign
    para = frame.paragraphs[0]
    para.text = text
    para.alignment = align
    para.line_spacing = line_spacing
    for run in para.runs:
        run.font.name = font
        run.font.size = Pt(size)
        run.font.bold = bold
        run.font.color.rgb = rgb(color)
    return box


def add_rich_text(slide, parts: list[tuple[str, str, bool]], x: float, y: float, w: float,
                  h: float, *, size: float = 18, align: PP_ALIGN = PP_ALIGN.LEFT):
    box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    frame = box.text_frame
    frame.clear(); frame.word_wrap = True
    frame.margin_left = frame.margin_right = Inches(0.04)
    para = frame.paragraphs[0]; para.alignment = align
    for text, color, bold in parts:
        run = para.add_run(); run.text = text
        run.font.name = FONT; run.font.size = Pt(size); run.font.bold = bold
        run.font.color.rgb = rgb(color)
    return box


def add_bullets(slide, bullets: list[str], x: float, y: float, w: float, h: float,
                *, size: float = 18, color: str = INK, accent: str = BLUE,
                gap: float = 8, numbered: bool = False):
    box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    frame = box.text_frame
    frame.clear(); frame.word_wrap = True
    frame.margin_left = frame.margin_right = Inches(0.04)
    frame.margin_top = frame.margin_bottom = Inches(0.02)
    for idx, item in enumerate(bullets):
        para = frame.paragraphs[0] if idx == 0 else frame.add_paragraph()
        marker = f"{idx + 1}." if numbered else "●"
        marker_run = para.add_run(); marker_run.text = marker + " "
        marker_run.font.name = FONT; marker_run.font.size = Pt(size - 1)
        marker_run.font.bold = True; marker_run.font.color.rgb = rgb(accent)
        text_run = para.add_run(); text_run.text = item
        text_run.font.name = FONT; text_run.font.size = Pt(size)
        text_run.font.color.rgb = rgb(color)
        para.space_after = Pt(gap)
        para.line_spacing = 1.08
    return box


def rect(slide, x: float, y: float, w: float, h: float, *, fill: str = WHITE,
         line: str = LINE, radius: bool = True, width: float = 1.0):
    shape_type = MSO_SHAPE.ROUNDED_RECTANGLE if radius else MSO_SHAPE.RECTANGLE
    shape = slide.shapes.add_shape(shape_type, Inches(x), Inches(y), Inches(w), Inches(h))
    shape.fill.solid(); shape.fill.fore_color.rgb = rgb(fill)
    shape.line.color.rgb = rgb(line); shape.line.width = Pt(width)
    return shape


def pill(slide, text: str, x: float, y: float, w: float, *, fill: str, color: str = WHITE):
    shape = rect(slide, x, y, w, 0.36, fill=fill, line=fill)
    add_text(slide, text, x, y + 0.015, w, 0.29, size=10.5, color=color,
             bold=True, align=PP_ALIGN.CENTER, valign=MSO_ANCHOR.MIDDLE)
    return shape


def card(slide, title: str, body: str, x: float, y: float, w: float, h: float, *,
         accent: str = BLUE, fill: str = WHITE, title_size: float = 15,
         body_size: float = 12.5):
    rect(slide, x, y, w, h, fill=fill, line=LINE)
    rect(slide, x, y, 0.08, h, fill=accent, line=accent, radius=False)
    add_text(slide, title, x + 0.22, y + 0.18, w - 0.38, 0.38,
             size=title_size, color=INK, bold=True)
    add_text(slide, body, x + 0.22, y + 0.68, w - 0.4, h - 0.83,
             size=body_size, color=MUTED)


def io_step_card(slide, stage: str, title: str, input_text: str, output_text: str,
                 x: float, y: float, w: float, h: float, *, accent: str = BLUE,
                 fill: str = WHITE, tag: str = ""):
    # 단계별 입력·출력을 좁은 교육용 카드에 동일한 규칙으로 표시한다.
    rect(slide, x, y, w, h, fill=fill, line=LINE)
    pill(slide, stage, x + 0.16, y + 0.14, 0.58, fill=accent)
    title_w = w - (1.85 if tag else 1.05)
    add_text(slide, title, x + 0.86, y + 0.15, title_w, 0.34,
             size=13.2, color=INK, bold=True, valign=MSO_ANCHOR.MIDDLE)
    if tag:
        pill(slide, tag, x + w - 0.9, y + 0.14, 0.72, fill=ORANGE)
    add_text(slide, "INPUT", x + 0.18, y + 0.62, 0.58, 0.22,
             size=8.5, color=accent, bold=True)
    add_text(slide, input_text, x + 0.82, y + 0.55, w - 1.0, 0.48,
             size=10.1, color=INK, line_spacing=0.98)
    add_text(slide, "OUTPUT", x + 0.18, y + 1.20, 0.62, 0.22,
             size=8.5, color=accent, bold=True)
    add_text(slide, output_text, x + 0.82, y + 1.09, w - 1.0, h - 1.2,
             size=10.1, color=MUTED, line_spacing=0.98)


def header(slide, title: str, section: str, number: int, *, accent: str = BLUE,
           subtitle: str | None = None):
    pill(slide, section, 0.55, 0.35, 1.35, fill=accent)
    add_text(slide, title, 0.55, 0.87, 12.1, 0.58, size=27, color=NAVY, bold=True)
    if subtitle:
        add_text(slide, subtitle, 0.58, 1.43, 11.9, 0.4, size=13, color=MUTED)
    add_text(slide, f"{number:02d}", 12.25, 0.37, 0.52, 0.28, size=10.5,
             color=MUTED, bold=True, align=PP_ALIGN.RIGHT)


def footer(slide, sources: list[str], number: int) -> None:
    slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(0.55), Inches(7.08),
                               Inches(12.78), Inches(7.08)).line.color.rgb = rgb(LINE)
    source = "근거: " + " · ".join(sources)
    add_text(slide, source, 0.58, 7.12, 11.75, 0.2, size=7.2, color=MUTED)
    add_text(slide, str(number), 12.4, 7.11, 0.35, 0.2, size=8, color=MUTED,
             align=PP_ALIGN.RIGHT)


def add_notes(slide, number: int, title: str, talking: list[str], sources: list[str]) -> None:
    note_text = title + "\n\n" + "\n".join(f"- {item}" for item in talking)
    note_text += "\n\n근거\n" + "\n".join(f"- {item}" for item in sources)
    try:
        slide.notes_slide.notes_text_frame.text = note_text
    except Exception:
        pass
    notes_rows.append({"number": number, "title": title, "talking": talking, "sources": sources})


def new_slide(title: str, section: str, *, accent: str = BLUE,
              subtitle: str | None = None, sources: list[str] | None = None,
              talking: list[str] | None = None):
    slide = prs.slides.add_slide(blank)
    set_bg(slide)
    number = len(prs.slides)
    header(slide, title, section, number, accent=accent, subtitle=subtitle)
    footer(slide, sources or [], number)
    add_notes(slide, number, title, talking or [], sources or [])
    return slide


def section_slide(title: str, subtitle: str, section_no: str, accent: str, bullets: list[str],
                  sources: list[str]):
    slide = prs.slides.add_slide(blank); set_bg(slide, NAVY)
    number = len(prs.slides)
    pill(slide, section_no, 0.75, 0.7, 1.25, fill=accent)
    add_text(slide, title, 0.75, 1.45, 11.7, 0.9, size=36, color=WHITE, bold=True)
    add_text(slide, subtitle, 0.78, 2.48, 11.3, 0.62, size=18, color="B9C7D8")
    for idx, item in enumerate(bullets):
        y = 3.55 + idx * 0.72
        rect(slide, 0.8, y, 0.35, 0.35, fill=accent, line=accent)
        add_text(slide, str(idx + 1), 0.8, y + 0.01, 0.35, 0.27, size=11,
                 color=WHITE, bold=True, align=PP_ALIGN.CENTER)
        add_text(slide, item, 1.38, y - 0.01, 10.7, 0.42, size=17, color=WHITE)
    add_text(slide, f"{number:02d}", 12.1, 6.78, 0.55, 0.28, size=11,
             color="8293A8", bold=True, align=PP_ALIGN.RIGHT)
    add_notes(slide, number, title, bullets, sources)
    return slide


def arrow(slide, x1: float, y1: float, x2: float, y2: float, color: str = MUTED, width: float = 2):
    conn = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(x1), Inches(y1),
                                      Inches(x2), Inches(y2))
    conn.line.color.rgb = rgb(color); conn.line.width = Pt(width)
    conn.line.end_arrowhead = True
    return conn


def flow(slide, labels: list[str], x: float, y: float, total_w: float, *,
         colors: list[str] | None = None, box_h: float = 0.82, font_size: float = 13):
    gap = 0.28
    box_w = (total_w - gap * (len(labels) - 1)) / len(labels)
    centers = []
    for idx, label in enumerate(labels):
        bx = x + idx * (box_w + gap)
        fill = (colors or [PALE_BLUE] * len(labels))[idx]
        accent = BLUE if fill == PALE_BLUE else CYAN if fill == PALE_CYAN else ORANGE if fill == PALE_ORANGE else PURPLE
        rect(slide, bx, y, box_w, box_h, fill=fill, line=accent)
        add_text(slide, label, bx + 0.06, y + 0.08, box_w - 0.12, box_h - 0.16,
                 size=font_size, color=INK, bold=True, align=PP_ALIGN.CENTER,
                 valign=MSO_ANCHOR.MIDDLE)
        centers.append((bx + box_w / 2, y + box_h / 2))
        if idx:
            arrow(slide, bx - gap + 0.02, y + box_h / 2, bx - 0.03, y + box_h / 2, MUTED, 1.4)
    return centers


def add_table(slide, headers: list[str], rows: list[list[str]], x: float, y: float,
              w: float, h: float, *, col_widths: list[float] | None = None,
              font_size: float = 11.5, accent: str = BLUE):
    table = slide.shapes.add_table(len(rows) + 1, len(headers), Inches(x), Inches(y),
                                   Inches(w), Inches(h)).table
    if col_widths:
        total = sum(col_widths)
        for idx, weight in enumerate(col_widths):
            table.columns[idx].width = Inches(w * weight / total)
    for c, value in enumerate(headers):
        cell = table.cell(0, c); cell.text = value
        cell.fill.solid(); cell.fill.fore_color.rgb = rgb(accent)
        cell.margin_left = cell.margin_right = Inches(0.05)
        for para in cell.text_frame.paragraphs:
            para.alignment = PP_ALIGN.CENTER
            for run in para.runs:
                run.font.name = FONT; run.font.size = Pt(font_size); run.font.bold = True
                run.font.color.rgb = rgb(WHITE)
    for r_idx, row in enumerate(rows, start=1):
        for c, value in enumerate(row):
            cell = table.cell(r_idx, c); cell.text = value
            cell.fill.solid(); cell.fill.fore_color.rgb = rgb(WHITE if r_idx % 2 else LIGHT)
            cell.margin_left = cell.margin_right = Inches(0.06)
            cell.margin_top = cell.margin_bottom = Inches(0.03)
            for para in cell.text_frame.paragraphs:
                para.alignment = PP_ALIGN.LEFT if c else PP_ALIGN.CENTER
                for run in para.runs:
                    run.font.name = FONT; run.font.size = Pt(font_size)
                    run.font.color.rgb = rgb(INK); run.font.bold = c == 0
    return table


def add_code(slide, code: str, x: float, y: float, w: float, h: float, *, size: float = 11.2):
    rect(slide, x, y, w, h, fill=NAVY2, line=NAVY2)
    add_text(slide, code, x + 0.18, y + 0.15, w - 0.36, h - 0.3,
             size=size, color="E8EFF7", font=MONO)


def add_image(slide, path: Path, x: float, y: float, w: float, h: float):
    if not path.is_file():
        card(slide, "그림 파일 누락", str(path), x, y, w, h, accent=RED, fill=PALE_RED)
        return
    from PIL import Image
    with Image.open(path) as image:
        iw, ih = image.size
    ratio = min(w / iw, h / ih)
    pw, ph = iw * ratio, ih * ratio
    slide.shapes.add_picture(str(path), Inches(x + (w - pw) / 2), Inches(y + (h - ph) / 2),
                             width=Inches(pw), height=Inches(ph))


# 01 표지
slide = prs.slides.add_slide(blank); set_bg(slide, NAVY)
rect(slide, 0, 0, 0.16, SH, fill=CYAN, line=CYAN, radius=False)
pill(slide, "MFM · preprocessing v2", 0.85, 0.75, 2.45, fill=CYAN)
add_text(slide, "개인 HAP 기반\nWGS·RNA 전처리 파이프라인", 0.85, 1.55, 8.5, 1.7,
         size=38, color=WHITE, bold=True)
add_text(slide, "입력 계약부터 HAP 재정렬·Compact Zarr·ASE QC·외부 전달까지", 0.9, 3.55, 9.8, 0.55,
         size=18, color="C4D2E3")
card(slide, "현재 production 기준", "3,212명 · WGS–RNA pair 449명 · 37 shards\nshort: BWA-MEM2 · long: minimap2 · RNA: STAR",
     8.95, 4.55, 3.55, 1.4, accent=ORANGE, fill=NAVY2, title_size=14, body_size=13)
add_text(slide, "2026-09-28 · 코드/문서 동기화본", 0.9, 6.75, 5.5, 0.3, size=11, color="8497AC")
add_notes(slide, 1, "표지", ["이 자료는 현재 v2 코드 동작을 기준으로 한다.", "과거 stage 번호는 호환성 설명에만 사용한다."],
          ["README.md", "GUIDE.md", "manifests/portable_summary.json"])

# 02 학습 목표
slide = new_slide("오늘 무엇을 이해해야 하는가", "INTRO", accent=CYAN,
                  subtitle="도구 이름 암기보다 ‘왜 이 분기가 있고 무엇을 검증하는가’에 집중",
                  sources=["GUIDE.md", "PIPELINE_DIAGRAM.md"],
                  talking=["입력·reference·VCF 계약을 이해한다.", "short/long/RNA 정책 차이를 설명할 수 있다.", "실패 sample을 격리하고 재실행할 수 있다."])
cards = [
    ("입력 계약", "CRAM/BAM/uBAM/FASTQ를 canonical BAM 하나로 표준화", BLUE),
    ("개인 HAP", "phased VCF→HAP1/HAP2 FASTA→read 재정렬", PURPLE),
    ("모델 입력", "window·readset·depth를 Compact Zarr로 고정", CYAN),
    ("운영 안전", "QC fail-fast + sample fail-continuation + provenance", ORANGE),
]
for idx,(t,b,c) in enumerate(cards):
    card(slide,t,b,0.7+(idx%2)*6.1,2.05+(idx//2)*1.65,5.65,1.32,accent=c,fill=WHITE,title_size=16,body_size=13)
add_text(slide,"핵심 질문",0.72,5.62,1.35,0.35,size=13,color=CYAN,bold=True)
add_text(slide,"“이 read는 어떤 reference·어떤 filter·어떤 HAP 규칙을 거쳐 모델에 들어왔는가?”",
         2.05,5.53,10.1,0.55,size=18,color=NAVY,bold=True)

# 03 데이터 구성
slide = new_slide("현재 cohort와 데이터 역할", "INTRO", accent=CYAN,
                  subtitle="PacBio 우선 선택 + 1KGP Illumina 보완, 개인 중복 제거",
                  sources=["manifests/portable_summary.json", "cohorts/pacbio_1kgp/README.md"],
                  talking=["최종 고유 개인은 3,212명이다.", "449명은 WGS 완료 후 RNA branch까지 진행한다.", "HPRC R3와 ONT는 현재 production 선택 source가 아니다."])
metrics=[("3,212","고유 개인",BLUE),("2,918","Illumina 30×",BLUE),("294","PacBio HiFi",PURPLE),("449","WGS–RNA pair",CYAN)]
for idx,(num,label,c) in enumerate(metrics):
    x=0.7+idx*3.12
    rect(slide,x,2.0,2.72,1.35,fill=WHITE,line=c,width=1.6)
    add_text(slide,num,x+0.1,2.2,2.52,0.5,size=28,color=c,bold=True,align=PP_ALIGN.CENTER)
    add_text(slide,label,x+0.1,2.8,2.52,0.3,size=12.5,color=MUTED,align=PP_ALIGN.CENTER)
add_table(slide,["data group","인원","처리 범위"],[["WGS-only","2,763","01→12→23 검증→24 cleanup"],["WGS–RNA pair","449","WGS 완료 후 13→26"]],
          1.25,4.05,10.8,1.55,col_widths=[2,1,4],font_size=13,accent=CYAN)
add_text(slide,"현재 manifest: Illumina 2,918 + HiFi 294 · CLR/ONT 0",3.15,6.15,7.1,0.35,size=14,color=MUTED,align=PP_ALIGN.CENTER)

# 04 전체 그림
slide = new_slide("전체 파이프라인 한 장으로 보기", "INTRO", accent=CYAN,
                  subtitle="직선 경로가 아니라 platform·VCF·optional branch가 있는 조건부 workflow",
                  sources=["PIPELINE_DIAGRAM.md", "run_wgs_sample.sh", "run_rna_sample.sh"],
                  talking=["01~12가 WGS 본선이다.", "13~26은 WGS–RNA pair만 실행한다.", "Clair3와 de novo는 기본 경로가 아닌 opt-in이다."])
flow(slide,["01 BAM","02 Preflight","03 VCF","06 HP/PS","09 FASTA","10 Realign","11 Windows","12 Zarr"],
     0.55,2.05,12.2,colors=[PALE_BLUE,PALE_RED,PALE_BLUE,PALE_PURPLE,PALE_PURPLE,PALE_ORANGE,PALE_BLUE,PALE_CYAN],box_h=0.83,font_size=9.6)
flow(slide,["13 FASTQ","14 Lift","15 Index","16 STAR","17 MarkDup","18 Coverage","19 QC"],
     0.7,3.75,11.9,colors=[PALE_BLUE,PALE_PURPLE,PALE_PURPLE,PALE_PURPLE,PALE_ORANGE,PALE_CYAN,PALE_CYAN],box_h=0.78,font_size=9.5)
flow(slide,["20 HAP variants","21 Counts","22 CNV","23 ASE","24 Metadata","25 Validate","26 Cleanup","27 Catalog"],
     0.7,5.25,11.9,colors=[PALE_ORANGE,PALE_PURPLE,PALE_ORANGE,PALE_CYAN,PALE_CYAN,PALE_RED,PALE_ORANGE,PALE_CYAN],box_h=0.78,font_size=8.8)

# 05 단계 번호
slide = new_slide("현재 실행 번호와 과거 이름을 구분한다", "INTRO", accent=CYAN,
                  sources=["README.md", "MIGRATION.md"],
                  talking=["실행은 01~20 파일명을 기준으로 한다.", "JSON schema와 환경변수에는 stage22/stage70 이름이 호환성 때문에 남아 있다."])
flow(slide,["01 입력 BAM","02 Preflight","03 VCF","06 Haplotag","09 HAP FASTA","10 재정렬","11 Windows","12 Zarr"],
     0.65,2.15,12.0,colors=[PALE_BLUE,PALE_RED,PALE_BLUE,PALE_PURPLE,PALE_PURPLE,PALE_ORANGE,PALE_CYAN,PALE_CYAN],box_h=0.88,font_size=10.8)
card(slide,"과거 Stage 22", "현재 파일: 10_realign_wgs.sh\nPython/JSON의 stage22 명칭은 호환 유지",0.95,3.75,5.35,1.42,accent=ORANGE,fill=PALE_ORANGE)
card(slide,"과거 Stage 70", "현재 파일: 12_compact_zarr.sh\nlayout/schema 명칭 일부는 호환 유지",7.0,3.75,5.35,1.42,accent=CYAN,fill=PALE_CYAN)
add_text(slide,"원칙: 발표·운영 문서에서는 현재 파일명을 먼저 쓰고, 과거 번호는 괄호로 병기",
         1.35,5.75,10.7,0.52,size=16,color=NAVY,bold=True,align=PP_ALIGN.CENTER)

# 06 실행 구조
slide = new_slide("한 파일에서 시작하지만 내부는 모듈식이다", "INTRO", accent=CYAN,
                  sources=["run_pipeline.sh", "config_v2.sh", "scripts/run_integrated_preprocess.sh"],
                  talking=["run_pipeline.sh가 대표 진입점이다.", "config_v2.sh가 manifest→환경변수→기본값 우선순위를 해석한다.", "실제 Python 기능은 integrate_wgs_rna CLI로 모인다."])
flow(slide,["run_pipeline.sh","runner / sample wrapper","01~20 stage scripts","Python CLI modules"],1.0,2.25,11.25,
     colors=[PALE_BLUE,PALE_ORANGE,PALE_PURPLE,PALE_CYAN],box_h=1.02,font_size=13)
card(slide,"설정 우선순위","1) manifest sample 행\n2) 환경변수 override\n3) config 기본값",1.25,4.05,3.3,1.55,accent=BLUE,fill=PALE_BLUE)
card(slide,"왜 모듈화했나","플랫폼별 분기·검증·재시작을 분리하고 sample 실패를 격리",5.0,4.05,3.3,1.55,accent=PURPLE,fill=PALE_PURPLE)
card(slide,"결과 경로","output/{wgs|wgs_rna_pair}/<project>/<sample>",8.75,4.05,3.3,1.55,accent=CYAN,fill=PALE_CYAN)

# 07 section
section_slide("WGS: 입력에서 개인 HAP 모델 입력까지","reference 계약과 HAP routing이 결과 해석을 결정한다","PART 1",BLUE,
              ["입력 형식을 canonical BAM으로 통일","22개 autosome preflight 후 VCF·HP/PS 확정","개인 HAP에 재정렬하고 Compact Zarr로 저장"],
              ["GUIDE.md","PIPELINE_DIAGRAM.md","단계별_문헌근거.md"])

# 08 input contract
slide = new_slide("네 가지 입력 형식, 하나의 canonical BAM", "WGS", accent=BLUE,
                  sources=["01_prepare_analysis_bam.sh", "src/integrate_wgs_rna/input_manifest.py", "prepare_alignment.py"],
                  talking=["aligned CRAM/BAM은 원 alignment를 재현 가능한 BAM으로 표준화한다.", "uBAM/FASTQ는 analysis reference에 새로 정렬한다.", "aligned와 unaligned 입력은 한 sample에서 섞지 않는다."])
inputs=[("CRAM","decode FASTA 필수","aligned"),("BAM","header·sort·index 검증","aligned"),("uBAM","SEQ/QUAL 추출 후 정렬","unaligned"),("FASTQ","R1/R2 또는 long-read","unaligned")]
for idx,(t,b,state) in enumerate(inputs):
    x=0.7+idx*3.1
    card(slide,t,f"{state}\n{b}",x,2.0,2.7,1.6,accent=BLUE if state=="aligned" else PURPLE,fill=WHITE,title_size=18,body_size=12.5)
arrow(slide,2.05,4.1,11.1,4.1,BLUE,2.2)
pill(slide,"canonical analysis BAM",4.75,4.66,3.85,fill=BLUE)
add_text(slide,"공통 검증: @SQ · @RG/SM · coordinate sort · index · base quality · record/tag 보존",
         1.15,5.55,11.0,0.55,size=15,color=NAVY,bold=True,align=PP_ALIGN.CENTER)

# 09 WGS diagram
slide = new_slide("WGS 본선: 조건 분기까지 포함한 실제 흐름", "WGS", accent=BLUE,
                  sources=["run_wgs_sample.sh", "PIPELINE.md"],
                  talking=["public VCF가 정상이라면 Clair3를 실행하지 않는다.", "de novo assembly도 명시적 opt-in이다.", "정상 경로는 01→02→03→06→09→10→11→12다."])
flow(slide,["01 Canonical BAM","02 Preflight","03 VCF","06 Haplotag","09 HAP FASTA","10 Realign","11 Windows","12 Compact"],
     0.55,2.0,12.2,colors=[PALE_BLUE,PALE_RED,PALE_BLUE,PALE_PURPLE,PALE_PURPLE,PALE_ORANGE,PALE_BLUE,PALE_CYAN],box_h=0.9,font_size=9.4)
card(slide,"04–05 Clair3", "public VCF가 없는 승인된 long-read 비교 실행만",1.25,4.0,4.75,1.25,accent=ORANGE,fill=PALE_ORANGE,title_size=16,body_size=12)
card(slide,"07–08 de novo", "명시적 opt-in · Stage 09 자동 승격 없음",7.3,4.0,4.75,1.25,accent=ORANGE,fill=PALE_ORANGE,title_size=16,body_size=12)
add_text(slide,"검증 실패는 fallback이 아니라 해당 sample 실패로 기록하고 다음 sample을 진행",1.2,5.85,10.9,0.46,size=15,color=RED,bold=True,align=PP_ALIGN.CENTER)

# WGS 단계별 입력·출력 01~06
slide = new_slide("WGS 01–06: reference 좌표에서 HAP tag까지", "WGS I/O", accent=BLUE,
                  subtitle="각 단계의 출력이 다음 단계에서 어떤 역할을 하는가",
                  sources=["01_prepare_analysis_bam.sh", "02_preflight.sh", "03_prepare_vcf.sh",
                           "04_call_longread_clair3.sh", "05_clair3_qc.sh", "06_haplotag.sh"],
                  talking=["Stage 01 BAM은 개인 유전체 정렬본이 아니다.",
                           "seq VCF와 tag VCF는 목적이 다르다.",
                           "Stage 06에서 read에 HP/PS가 붙지만 좌표는 아직 공통 reference다."])
wgs_front = [
    ("01", "Canonical BAM", "CRAM/BAM 또는 uBAM/FASTQ\nmanifest · analysis reference",
     "sample.analysis.bam + index\ninput contract · QC", BLUE, ""),
    ("02", "Preflight", "canonical BAM · reference\nVCF/model · manifest",
     "preflight.tsv/json\nPASS일 때만 다음 단계", RED, ""),
    ("03", "Sample VCF", "1KGP panel 또는\npublic phased VCF",
     "S.seq.vcf.gz\nS.tag.vcf.gz", PURPLE, ""),
    ("04", "Clair3 call", "long-read BAM · model\ndecode reference",
     "self-called phased VCF\nClair3 HP/PS BAM", ORANGE, "선택"),
    ("05", "Clair3 QC", "Stage 04 VCF · BAM\nmodel · reference",
     "clair3_output_qc.tsv\nQC metadata JSON", ORANGE, "선택"),
    ("06", "Haplotag", "canonical BAM · tag VCF\ndecode reference",
     "HP/PS-tagged BAM + index\nsummary · metadata", CYAN, ""),
]
for idx, item in enumerate(wgs_front):
    stage, title, inp, out, color, tag = item
    io_step_card(slide, stage, title, inp, out,
                 0.58 + (idx % 3) * 4.12, 1.93 + (idx // 3) * 2.12,
                 3.88, 1.86, accent=color, fill=WHITE, tag=tag)
add_text(slide, "01·06의 BAM은 공통 reference 좌표 · 개인 HAP 좌표 재정렬은 Stage 10",
         1.0, 6.28, 11.3, 0.37, size=14.2, color=RED, bold=True,
         align=PP_ALIGN.CENTER)

# WGS 단계별 입력·출력 07~12
slide = new_slide("WGS 07–12: 개인 HAP 좌표와 모델 입력 생성", "WGS I/O", accent=PURPLE,
                  subtitle="개인 서열 생성과 read 재정렬은 서로 다른 단계다",
                  sources=["07_denovo_assembly.sh", "08_anchor_assembly.sh", "09_personal_fasta.sh",
                           "10_realign_wgs.sh", "11_windows.sh", "12_compact_zarr.sh"],
                  talking=["de novo 경로는 현재 Stage 09로 자동 연결되지 않는다.",
                           "Stage 09는 서열을 만들고 Stage 10은 read를 그 서열에 재정렬한다.",
                           "Stage 12 Compact Zarr가 WGS 모델의 주 입력이다."])
wgs_back = [
    ("07", "de novo assembly", "long-read BAM · reference", "GFA/contig FASTA\nassembly QC", ORANGE, "선택"),
    ("08", "Anchoring QC", "assembly FASTA\nanchor reference", "PAF · anchor summary\n자동 HAP 승격 없음", ORANGE, "선택"),
    ("09", "Personal FASTA", "sequence reference\nphased seq VCF", "HAP1/HAP2 FASTA\nFAI · chain · QC", PURPLE, ""),
    ("10", "HAP 재정렬", "HP/PS BAM\nHAP1/HAP2 FASTA", "HAP별 CRAM + index\ntarget.fa · realign QC", PURPLE, ""),
    ("11", "Windows", "HAP FASTA · chain\nStage 10 QC · blacklist", "windows.tsv\nblacklist · metadata", BLUE, ""),
    ("12", "Compact Zarr", "HAP CRAM · FASTA\nwindows · source BAM", "sequence/readset/depth Zarr\nmanifest · full trace", CYAN, ""),
]
for idx, item in enumerate(wgs_back):
    stage, title, inp, out, color, tag = item
    io_step_card(slide, stage, title, inp, out,
                 0.58 + (idx % 3) * 4.12, 1.93 + (idx // 3) * 2.12,
                 3.88, 1.86, accent=color, fill=WHITE, tag=tag)
add_text(slide, "09 = 개인 서열 생성 · 10 = 개인 좌표 재정렬 · 11 = window 계약 · 12 = 최종 WGS 모델 입력",
         1.0, 6.28, 11.3, 0.37, size=14.2, color=NAVY, bold=True,
         align=PP_ALIGN.CENTER)

# 10 aligner comparison
slide = new_slide("Short와 long read는 왜 다른 aligner를 쓰는가", "WGS", accent=BLUE,
                  sources=["prepare_alignment.py", "personal_realign.py", "단계별_문헌근거.md R1/R3/R7"],
                  talking=["BWA-MEM2는 Illumina paired-end production workflow와 호환된다.", "minimap2 preset은 long-read 오류 구조에 맞춘다.", "도구 차이는 기술 차이를 숨기지 않고 metadata로 기록한다."])
add_table(slide,["구분","Short read","Long read"],[
    ["기본 도구","BWA-MEM2 mem","minimap2"],
    ["입력 특성","짧은 paired-end, 낮은 오류율","긴 read, platform별 오류 profile"],
    ["preset/options","-p -C -Y","map-hifi / map-pb / map-ont"],
    ["핵심 보존","pair, HP/PS/XU/OD","긴 indel·split alignment"],
    ["대체 사용","minimap2 -x sr는 명시 opt-in","BWA 사용 안 함"],
],0.8,2.0,11.75,3.65,col_widths=[1.4,3,3],font_size=12.3,accent=BLUE)
add_text(slide,"같은 HAP 목적이라도 read technology에 맞는 aligner를 사용하고, downstream schema를 통일한다.",
         1.0,6.05,11.3,0.48,size=16,color=NAVY,bold=True,align=PP_ALIGN.CENTER)

# 11 preflight
slide = new_slide("Stage 02: 잘못된 입력은 계산 전에 멈춘다", "WGS", accent=BLUE,
                  sources=["02_preflight.sh", "tests/test_preflight_22contig.sh", "LIBRARY_PROVENANCE.md"],
                  talking=["chr1 한 개가 아니라 22개 autosome을 검사한다.", "reference·VCF·CRAM M5/contig/sample 중 하나라도 불일치하면 해당 sample만 중단한다.", "batch 전체는 계속 진행한다."])
checks=[("Reference","FASTA/FAI · contig 길이 · M5",BLUE),("Alignment","CRAM decode · @SQ · sample · index",PURPLE),("VCF","donor · REF · phase · contig",CYAN),("Library","project PCR provenance",ORANGE)]
for idx,(t,b,c) in enumerate(checks):
    card(slide,t,b,0.72+idx*3.08,2.0,2.7,1.35,accent=c,fill=WHITE,title_size=16,body_size=12)
flow(slide,["검증 PASS","Stage 03 이후 진행"],2.0,4.25,4.25,colors=[PALE_CYAN,PALE_BLUE],box_h=0.85,font_size=14)
flow(slide,["검증 FAIL","status.tsv 기록","다음 sample"],7.05,4.25,4.25,colors=[PALE_RED,PALE_ORANGE,PALE_CYAN],box_h=0.85,font_size=12)
add_text(slide,"기준을 낮춰 통과시키지 않는다 — 문제 sample을 격리해 전체 cohort 오염을 막는다.",
         1.05,5.8,11.2,0.48,size=16,color=RED,bold=True,align=PP_ALIGN.CENTER)

# 12 VCF decision
slide = new_slide("VCF 경로는 ‘있으면 사용’이 아니라 계약으로 선택", "WGS", accent=BLUE,
                  sources=["03_prepare_vcf.sh", "public_vcf.py", "config_v2.sh"],
                  talking=["Illumina는 1KGP phased panel을 쓴다.", "long read는 manifest의 제공 public VCF가 기본이다.", "Clair3는 VCF가 없고 명시적으로 승인한 재현/비교에서만 사용한다."])
card(slide,"Illumina", "1KGP 2022 phased SNV+indel panel\nseq VCF + tag VCF",0.8,2.0,3.55,1.6,accent=BLUE,fill=PALE_BLUE,title_size=18)
card(slide,"Long read production", "project/donor public phased VCF\ndonor·REF·phase 전수 검증",4.9,2.0,3.55,1.6,accent=PURPLE,fill=PALE_PURPLE,title_size=18)
card(slide,"Clair3 opt-in", "제공 VCF 없음 + 승인된 model\n검증 실패 자동 fallback 금지",9.0,2.0,3.55,1.6,accent=ORANGE,fill=PALE_ORANGE,title_size=18)
flow(slide,["S.seq.vcf.gz","HAP FASTA"],1.25,4.55,4.4,colors=[PALE_BLUE,PALE_PURPLE],box_h=0.85,font_size=14)
flow(slide,["S.tag/public VCF","WhatsHap HP/PS"],7.0,4.55,4.9,colors=[PALE_CYAN,PALE_PURPLE],box_h=0.85,font_size=14)

# 13 haplotag
slide = new_slide("HP/PS는 무엇을 말해 주는가", "WGS", accent=BLUE,
                  sources=["06_haplotag.sh", "WhatsHap", "단계별_문헌근거.md R5/R15"],
                  talking=["HP는 read가 지지하는 haplotype 번호다.", "PS는 phase block을 식별한다.", "HP가 없다고 read가 나쁜 것은 아니며, marker 정보가 부족할 수 있다."])
card(slide,"HP=1", "HAP1 allele을 지지하는 read/fragment",0.85,2.05,3.45,1.35,accent=BLUE,fill=PALE_BLUE,title_size=19)
card(slide,"HP=2", "HAP2 allele을 지지하는 read/fragment",4.95,2.05,3.45,1.35,accent=PURPLE,fill=PALE_PURPLE,title_size=19)
card(slide,"HP=0 / missing", "정보 부족·동점·marker 비중첩\n정책에 따라 양쪽/제외",9.05,2.05,3.45,1.35,accent=ORANGE,fill=PALE_ORANGE,title_size=19)
add_table(slide,["tag","의미","주의"],[["HP","haplotype label","reference coordinate 자체가 아님"],["PS","phase-set block","block 간 방향은 독립 가능"],["XU","원 source unmapped 표시","재정렬 후 mapped 여부와 구분"],["OD","source duplicate 운반","FASTQ 변환에서 0x400 복원용"]],
          1.15,4.05,11.0,2.05,col_widths=[1,3,4],font_size=11.8,accent=PURPLE)

# 14 personal fasta
slide = new_slide("개인 HAP FASTA와 chain: 좌표계가 달라진다", "WGS", accent=BLUE,
                  sources=["09_personal_fasta.sh", "personal_fasta.py", "bcftools consensus"],
                  talking=["phased genotype의 HAP1/HAP2 allele을 각각 consensus에 적용한다.", "indel 때문에 HAP 좌표는 GRCh38 좌표와 같지 않다.", "chain은 지금 projection하지 않고 후속 분석을 위해 검증·보존한다."])
flow(slide,["GRCh38 autosomes","phased seq VCF","bcftools consensus"],0.8,2.0,5.55,
     colors=[PALE_BLUE,PALE_PURPLE,PALE_ORANGE],box_h=0.86,font_size=12)
arrow(slide,6.55,2.43,7.35,2.43,PURPLE,2.0)
card(slide,"HAP1 FASTA", "chr1–chr22 + FAI + chain",7.55,1.72,4.65,1.25,accent=BLUE,fill=PALE_BLUE,title_size=17)
card(slide,"HAP2 FASTA", "chr1–chr22 + FAI + chain",7.55,3.18,4.65,1.25,accent=PURPLE,fill=PALE_PURPLE,title_size=17)
add_text(slide,"Stage 10 결과 POS/CIGAR = HAP 좌표",1.0,4.85,5.25,0.45,size=17,color=RED,bold=True,align=PP_ALIGN.CENTER)
add_text(slide,"chain 방향·contig 길이·MD5를 QC에 연결",6.65,4.85,5.5,0.45,size=17,color=GREEN,bold=True,align=PP_ALIGN.CENTER)
add_text(slide,"GRCh38 좌표가 필요한 후속 분석만 chain liftover를 수행한다.",1.45,5.75,10.4,0.45,size=16,color=NAVY,bold=True,align=PP_ALIGN.CENTER)

# 15 routing
slide = new_slide("Stage 10 이전: 어떤 read를 어느 HAP로 보낼까", "WGS", accent=BLUE,
                  sources=["10_realign_wgs.sh", "personal_realign.py", "cli.py --hap-read-mode"],
                  talking=["default all_to_both는 coverage 보존용이다.", "hp_routed_with_unknown은 기존 HP를 우선한다.", "hap_specific은 HAP-consistent read만 남기는 별도 분석 옵션이다."])
add_table(slide,["mode","HP=1","HP=2","HP 없음/0","목적"],[
    ["all_to_both","양쪽","양쪽","양쪽","coverage 보존 기본"],
    ["hp_routed_with_unknown","HAP1","HAP2","양쪽 또는 drop","HP provenance 활용"],
    ["hap_specific","HAP1만","HAP2만","제외","HAP-specific 분석"],
],0.7,2.0,12.0,2.35,col_widths=[2.1,1.2,1.2,2.3,3],font_size=12,accent=BLUE)
card(slide,"중요", "unknown을 양 HAP에 복제한 결과를 HAP-specific count로 이중 집계하면 안 된다.",1.05,4.9,11.25,1.05,accent=RED,fill=PALE_RED,title_size=15,body_size=14)
add_text(slide,"routing mode는 output 폴더명과 metadata에 포함되어 서로 섞이지 않는다.",1.25,6.23,10.8,0.35,size=14,color=MUTED,align=PP_ALIGN.CENTER)

# 16 realign
slide = new_slide("Stage 10: 같은 read를 개인 HAP 좌표에서 다시 해석", "WGS", accent=BLUE,
                  sources=["10_realign_wgs.sh", "personal_realign.py", "HIGH_PARALLEL_POLICY.md"],
                  talking=["read routing 후 HAP1/HAP2 reference에 독립 정렬한다.", "short는 paired FASTQ와 HP/PS/XU/OD를 보존한다.", "HAP1/HAP2를 병렬 실행해 결과 의미는 바꾸지 않고 시간을 줄인다."])
flow(slide,["canonical BAM","route / collate","FASTQ 또는 long read","HAP1 정렬"],0.7,2.05,8.15,
     colors=[PALE_BLUE,PALE_ORANGE,PALE_BLUE,PALE_PURPLE],box_h=0.9,font_size=12)
card(slide,"HAP2 정렬", "동일 routing input\n별도 HAP2 FASTA",9.45,1.85,2.65,1.3,accent=PURPLE,fill=PALE_PURPLE,title_size=15)
arrow(slide,8.95,2.5,9.4,2.5,PURPLE,1.8)
add_table(slide,["확인 항목","왜 필요한가"],[["input/output record·pair 수","FASTQ 변환 누락·오-pair 탐지"],["mapped fraction·NM·MAPQ·AS","개인 reference 정렬 품질"],["HP/PS/XU/OD","provenance·duplicate 보존"],["chain fingerprint","출력 HAP 좌표의 reference 연결"]],
          1.05,3.7,11.2,2.35,col_widths=[2.7,6],font_size=12,accent=PURPLE)

# 17 filters
slide = new_slide("Short와 long read의 filter가 다른 이유", "WGS", accent=BLUE,
                  sources=["12_compact_zarr.sh", "full_preprocess.py", "단계별_문헌근거.md R11/R14/R16"],
                  talking=["short-read supplementary는 split/ambiguous evidence로 모델 입력에서 제외한다.", "long-read supplementary는 SV/complex alignment evidence이므로 별도 channel에 보존한다.", "depth는 두 플랫폼 모두 primary 중심으로 계산한다."])
add_table(slide,["항목","Short read base","Long read SNV/HAP","Long read SV evidence"],[
    ["secondary","제외","제외","QC/near-best만"],
    ["supplementary","제외","primary depth에서 제외","별도 evidence 보존"],
    ["duplicate","제외","일괄 제외 안 함","source 비율 QC"],
    ["depth","primary·non-dup","primary","supp depth 별도 기록"],
],0.65,2.0,12.05,3.15,col_widths=[1.6,2.7,2.7,3.1],font_size=11.7,accent=BLUE)
card(slide,"해석", "같은 flag라도 기술과 분석 목적이 다르다. 보존 여부와 depth 기여 여부를 분리한다.",1.0,5.55,11.3,0.95,accent=CYAN,fill=PALE_CYAN,title_size=14,body_size=14)

# 18 duplicate/PCR
slide = new_slide("PCR 사용 여부와 duplicate 처리는 같은 질문이 아니다", "WGS", accent=BLUE,
                  sources=["LIBRARY_PROVENANCE.md", "library_provenance.py", "prepare_alignment.py"],
                  talking=["1KGP는 PCR-free지만 Picard duplicate flag가 존재한다.", "PacBio에는 blanket MarkDuplicates를 적용하지 않는다.", "GEUVADIS RNA는 10-cycle PCR이며 STAR 뒤 mark only다."])
add_table(slide,["자료원","PCR provenance","duplicate 정책"],[
    ["1KGP Illumina 2,918","pcr_free_confirmed","upstream Picard 0x400 보존"],
    ["HPRC R1/HGSVC 63","pcr_free_protocol","blanket marking 없음"],
    ["HPRC R2 232","pcr_free_inferred","불확실성 유지·blanket 없음"],
    ["GEUVADIS RNA 449","pcr_amplified_confirmed · 10 cycles","STAR 뒤 Picard mark only"],
],0.7,2.0,11.95,3.35,col_widths=[2.8,3.3,4],font_size=11.8,accent=ORANGE)
add_text(slide,"REMOVE_DUPLICATES=false: record는 보존하고 downstream read-set에서 명시적으로 제외",
         1.0,5.8,11.3,0.5,size=16,color=ORANGE,bold=True,align=PP_ALIGN.CENTER)

# 19 compact
slide = new_slide("Stage 12: 모델이 읽을 구조로 압축·검증", "WGS", accent=BLUE,
                  sources=["12_compact_zarr.sh", "full_preprocess.py", "direct_compact_writer.py"],
                  talking=["1 Mb window index를 기준으로 HAP sequence, readset, depth를 묶는다.", "layout 이름에 routing/filter 정책이 포함된다.", "full trace는 source read·CIGAR·depth·metadata corruption을 탐지한다."])
tree = """compact/<LAYOUT>__<HAP_READ_MODE>/<contig>/
├── hap_sequence.zarr
├── readset.HAP1.zarr
├── readset.HAP2.zarr
├── depth.zarr
├── segment / read metadata
└── QC + manifest hashes"""
add_code(slide,tree,0.8,2.0,5.6,3.65,size=12.5)
card(slide,"모델 입력 단위", "1 Mb window\nHAP sequence + reads + per-base depth",7.05,2.0,5.1,1.45,accent=CYAN,fill=PALE_CYAN,title_size=17,body_size=14)
card(slide,"검증 최적화", "65,536-record prefetch\nSQLite 8,192-row batch\n결과 hash 동일",7.05,3.8,5.1,1.45,accent=GREEN,fill=WHITE,title_size=17,body_size=14)
add_text(slide,"최종 compact 검증 성공 뒤에만 재생성 가능한 work 중간물을 cleanup",
         1.25,6.0,10.9,0.4,size=15,color=NAVY,bold=True,align=PP_ALIGN.CENTER)

# 20 section RNA
section_slide("RNA: WGS 개인 유전체 위에서 ASE를 계산","RNA는 독립 파이프라인이 아니라 검증된 WGS HAP 산출물의 후속 branch","PART 2",CYAN,
              ["GTF를 HAP 좌표로 lift","STAR 재정렬 후 duplicate mark","allele count·CNV·ASE QC를 함께 기록"],
              ["RNASEQ_ALLELESEQ2_REVIEW.md","run_rna_sample.sh","13_rna_to_fastq.py","21_allelic_counts.py"])

# 21 RNA diagram
slide = new_slide("RNA 13~26 단계의 전체 흐름", "RNA", accent=CYAN,
                  sources=["run_rna_sample.sh", "PIPELINE.md"],
                  talking=["WGS HAP FASTA와 chain이 선행 조건이다.", "HAP1/HAP2 STAR 정렬과 coverage를 분리한다.", "AlleleSeq2 원칙을 참고하되 현재 구현 차이를 QC에 기록한다."])
flow(slide,["13 FASTQ","14 GTF lift","15 STAR index","16 STAR align","17 MarkDup","18 Coverage","19 RNA QC"],
     0.6,2.05,12.0,colors=[PALE_BLUE,PALE_PURPLE,PALE_PURPLE,PALE_PURPLE,PALE_ORANGE,PALE_CYAN,PALE_CYAN],box_h=0.88,font_size=9.4)
flow(slide,["20 HAP variants","21 Allele count","22 CNV","23 ASE QC","24 Metadata","25 Validate","26 Cleanup"],
     0.6,4.15,12.0,colors=[PALE_ORANGE,PALE_PURPLE,PALE_ORANGE,PALE_CYAN,PALE_CYAN,PALE_RED,PALE_ORANGE],box_h=0.88,font_size=9.2)
add_text(slide,"27_make_catalog.py는 sample 실행 밖에서 cohort 단위로 수동 집계",2.0,5.75,9.4,0.48,size=15,color=NAVY,bold=True,align=PP_ALIGN.CENTER)

# RNA 단계별 입력·출력
slide = new_slide("RNA 13–26: 개인 좌표 signal에서 ASE label·최종 검증까지", "RNA I/O", accent=CYAN,
                  subtitle="각 번호가 runner에서 독립 stage로 기록되며 HAP 작업은 단계 안에서 병렬 실행된다",
                  sources=["13_rna_to_fastq.py", "14_lift_annotation.py", "15_star_index.py", "16_star_align.py",
                           "17_mark_duplicates.py", "18_coverage_zarr.py", "19_rna_realign_qc.py",
                           "20_hap_variants.py", "21_allelic_counts.py", "22_cnv_filter.py", "23_ase_qc.py",
                           "24_stamp_rna_metadata.py", "25_validate_sample_package.sh", "26_cleanup_intermediates.sh"],
                  talking=["RNA는 WGS HAP FASTA·chain·개인 WGS 정렬을 선행 입력으로 사용한다.",
                           "HAP별 coverage는 총 발현이며 Stage 21 AS 비교가 allele count를 만든다.",
                           "Stage 23이 학습 가능한 ASE label과 sample 적격성을 결정한다."])
add_table(slide, ["단계", "입력", "주요 출력", "다음 사용처"], [
    ["13", "RNA BAM/CRAM · preflight JSON", "paired FASTQ · singleton/pair QC", "HAP별 공통 read"],
    ["14", "GENCODE GTF · HAP chain", "annotation.gtf · exons TSV · lift QC", "STAR junction · gene 구간"],
    ["15–16", "HAP FASTA · lifted GTF · FASTQ", "STAR index · HAP별 정렬 BAM", "duplicate mark"],
    ["17", "STAR BAM", "marked BAM · duplicate metrics", "coverage · allele count"],
    ["18–19", "marked BAM · WGS windows", "전장/window coverage · RNA QC", "모델 signal · provenance"],
    ["20–21", "VCF/WGS evidence · 양 HAP RNA BAM", "validated variants · fragment allele counts", "ASE 통계"],
    ["22–24", "WGS depth · counts/sidecar", "CNV · ASE QC · sample metadata", "training gate"],
    ["25–26", "전체 package · QC · hash", "final validation · 조건부 cleanup", "완료 sample"],
], 0.55, 1.94, 12.22, 4.55, col_widths=[1.35, 3.25, 4.0, 2.5],
   font_size=8.9, accent=CYAN)
add_text(slide, "주의: HAP1/HAP2 coverage만 비교하면 ASE가 아니다 — 같은 fragment의 AS를 Stage 21에서 비교",
         0.85, 6.61, 11.6, 0.3, size=12.2, color=RED, bold=True,
         align=PP_ALIGN.CENTER)

# 22 GTF lift
slide = new_slide("Stage 14: annotation도 HAP 좌표에 맞춘다", "RNA", accent=CYAN,
                  sources=["14_lift_annotation.py", "RNASEQ_ALLELESEQ2_REVIEW.md"],
                  talking=["GENCODE v29 feature를 GRCh38→HAP chain으로 옮긴다.", "입력/출력 feature 수와 snap/drop 비율을 JSON으로 남긴다.", "chain 누락은 조용히 건너뛰지 않고 실패한다."])
flow(slide,["GENCODE v29 GTF","GRCh38→HAP chain","HAP1/HAP2 GTF"],1.15,2.25,11.0,
     colors=[PALE_BLUE,PALE_ORANGE,PALE_CYAN],box_h=1.0,font_size=15)
add_table(slide,["QC","의미"],[["input_features","원 annotation feature 수"],["output_features","HAP 좌표로 살아남은 feature"],["dropped","chain 밖/변환 불가"],["snap_ratio","경계 보정이 필요했던 비율"],["exon_count","coverage annotation sanity check"]],
          2.0,3.85,9.4,2.25,col_widths=[2.3,5.5],font_size=11.5,accent=CYAN)

# 23 STAR
slide = new_slide("Stage 13–18: FASTQ 복원부터 HAP별 coverage까지", "RNA", accent=CYAN,
                  sources=["13_rna_to_fastq.py", "15_star_index.py", "16_star_align.py", "17_mark_duplicates.py", "18_coverage_zarr.py", "TOPMed RNA pipeline"],
                  talking=["HAP별 STAR index와 lifted GTF를 사용한다.", "R1/R2 pair 수와 singleton discard를 먼저 검증한다.", "Picard는 REMOVE_DUPLICATES=false로 flag만 표시한다."])
flow(slide,["RNA BAM","name-collated FASTQ","STAR HAP1/HAP2","Picard mark","coverage Zarr"],0.65,2.05,12.05,
     colors=[PALE_BLUE,PALE_ORANGE,PALE_PURPLE,PALE_ORANGE,PALE_CYAN],box_h=0.92,font_size=11.5)
add_table(slide,["TOPMed와 동일","개인 HAP 목적상 다른 점"],[
    ["STAR 2-pass·splice-aware 정렬","HAP1/HAP2 index를 각각 생성"],
    ["MarkDuplicates로 표시","두 HAP alignment를 후속 비교"],
    ["junction/mismatch/multimap QC","lift·AS·CNV·ASE provenance 추가"],
],1.05,3.85,11.2,2.15,col_widths=[4,5],font_size=12,accent=CYAN)

# 24 AS assignment
slide = new_slide("Stage 21: AS는 절대 cutoff보다 HAP 간 상대 비교", "RNA", accent=CYAN,
                  sources=["21_allelic_counts.py", "RNASEQ_ALLELESEQ2_REVIEW.md", "단계별_문헌근거.md"],
                  talking=["같은 fragment의 HAP1/HAP2 alignment score를 비교한다.", "동점·비정보성 fragment는 ambiguous로 남긴다.", "split-half와 both-haps 산출물로 민감도를 본다."])
card(slide,"HAP1 승", "AS(H1) > AS(H2)\n→ allele 1 evidence",0.85,2.0,3.45,1.45,accent=BLUE,fill=PALE_BLUE,title_size=18)
card(slide,"HAP2 승", "AS(H2) > AS(H1)\n→ allele 2 evidence",4.95,2.0,3.45,1.45,accent=PURPLE,fill=PALE_PURPLE,title_size=18)
card(slide,"동점/저정보", "같은 score 또는 유효 비교 불가\n→ ambiguous",9.05,2.0,3.45,1.45,accent=ORANGE,fill=PALE_ORANGE,title_size=18)
add_table(slide,["sidecar","용도"],[[".split_half.tsv","query name hash 기반 재현성"],[".both_haps.tsv","양쪽에 정렬된 fragment만 민감도 분석"],[".meta.json","filter·score·입출력 count provenance"]],
          1.45,4.25,10.45,1.85,col_widths=[3,5],font_size=12,accent=PURPLE)

# 25 CNV/ASE
slide = new_slide("Stage 20–24: 검증된 variant와 CNV를 반영해 ASE label을 만든다", "RNA", accent=CYAN,
                  sources=["20_hap_variants.py", "21_allelic_counts.py", "22_cnv_filter.py", "23_ase_qc.py", "24_stamp_rna_metadata.py", "AlleleSeq2"],
                  talking=["WGS HAP depth ratio로 copy-number imbalance 후보를 거른다.", "기본 허용 구간은 0.5≤rd≤1.5다.", "split-half·single-HAP sensitivity·overdispersion을 QC한다."])
flow(slide,["allelic counts","WGS HAP depth","CNV rd filter","beta-binomial ASE QC","metadata stamp"],0.65,2.0,12.05,
     colors=[PALE_BLUE,PALE_PURPLE,PALE_ORANGE,PALE_CYAN,PALE_CYAN],box_h=0.92,font_size=11.2)
metrics=[("0.5–1.5","기본 rd 허용",ORANGE),("split-half","재현성",CYAN),("both-haps","민감도",PURPLE),("metadata","정책 고정",BLUE)]
for idx,(num,label,c) in enumerate(metrics):
    x=0.8+idx*3.05
    rect(slide,x,4.05,2.65,1.25,fill=WHITE,line=c,width=1.5)
    add_text(slide,num,x+0.08,4.24,2.49,0.42,size=20,color=c,bold=True,align=PP_ALIGN.CENTER)
    add_text(slide,label,x+0.08,4.78,2.49,0.25,size=11.5,color=MUTED,align=PP_ALIGN.CENTER)
add_text(slide,"CNV filter는 생물학적 allelic imbalance와 copy-number artifact를 구분하기 위한 보호막",
         1.0,5.8,11.3,0.48,size=15,color=NAVY,bold=True,align=PP_ALIGN.CENTER)

# 26 RNA QC
slide = new_slide("RNA QC는 ‘몇 개 남았는지’를 단계별로 설명한다", "RNA", accent=CYAN,
                  sources=["13_rna_to_fastq.py", "14_lift_annotation.py", "16_star_align.py", "17_mark_duplicates.py", "21_allelic_counts.py"],
                  talking=["FASTQ pair/singleton, duplicate, GTF lift, STAR, allele assignment를 각각 기록한다.", "단일 최종 count만 보면 어느 단계에서 손실됐는지 알 수 없다."])
items=[("FASTQ","R1/R2 pair · singleton discard",BLUE),("STAR","mapped · unique · multimapped",PURPLE),("Duplicate","marked count · REMOVE=false",ORANGE),("GTF lift","input/output · drop · snap",CYAN),("ASE","H1/H2 · tie · ambiguous",GREEN)]
for idx,(t,b,c) in enumerate(items):
    y=1.95+idx*0.86
    rect(slide,0.9,y,2.0,0.58,fill=c,line=c)
    add_text(slide,t,0.95,y+0.08,1.9,0.3,size=13,color=WHITE,bold=True,align=PP_ALIGN.CENTER)
    add_text(slide,b,3.25,y+0.04,8.7,0.42,size=15,color=INK)
    arrow(slide,2.93,y+0.29,3.18,y+0.29,c,1.5)

# 27 section ops
section_slide("운영: 한 번에 돌리되 sample 단위로 안전하게","병렬도·저장공간·shard·복구를 데이터 계약의 일부로 본다","PART 3",ORANGE,
              ["대표 entrypoint와 output 경로","실패 격리·재시작·resource guard","37개 shard 전달과 검증"],
              ["PORTABLE_GUIDE.md","HIGH_PARALLEL_POLICY.md","SHARDED_HANDOFF.md"])

# 운영 단계별 입력·출력
slide = new_slide("Stage 23–27: 최종 검증·cleanup·cohort 집계", "OPS I/O", accent=ORANGE,
                  subtitle="WGS/RNA 검증 PASS 뒤에만 cleanup · Catalog는 cohort 완료 뒤 수동 실행",
                  sources=["23_validate_sample_package.sh", "24_cleanup_intermediates.sh", "25_validate_sample_package.sh", "26_cleanup_intermediates.sh", "27_make_catalog.py"],
                  talking=["Stage 24/26은 새 생물학적 데이터를 만들지 않고 재생성 가능한 중간물을 제거한다.",
                           "최종 Compact/RNA/QC/metadata와 원본 입력은 보존한다.",
                           "Stage 27은 sample별 실제 경로와 학습 적격성을 cohort TSV로 집계한다."])
io_step_card(slide, "23/25", "Final validation", "WGS: compact/trace/QC\nRNA: coverage/counts/ASE/metadata",
             "sample_package_validation.json\n불완전 단계·파일 명시", 0.55, 2.02, 3.9, 2.05,
             accent=RED, fill=PALE_RED, tag="필수")
io_step_card(slide, "24/26", "Cleanup", "최신 PASS 보고서·hash\n보존 목록·삭제 대상 재확인",
             "cleanup_intermediates.tsv\n재생성 가능한 중간물 제거", 4.72, 2.02, 3.9, 2.05,
             accent=ORANGE, fill=PALE_ORANGE, tag="조건")
io_step_card(slide, "27", "Cohort Catalog", "cohort sample roots · sample.json\ncompact/trace/RNA · 선택적 split/meta",
             "sample_catalog_v2.<cohort>.tsv\n경로·완료·학습 적격성", 8.9, 2.02, 3.9, 2.05,
             accent=GREEN, fill=PALE_CYAN, tag="수동")
flow(slide, ["WGS-only: 12", "23 validate", "24 cleanup", "27 catalog"],
     0.8, 4.78, 5.85, colors=[PALE_CYAN, PALE_ORANGE, PALE_BLUE, PALE_CYAN],
     box_h=0.72, font_size=10)
flow(slide, ["WGS–RNA: 12", "13–24 RNA", "25 validate", "26 cleanup", "27 catalog"],
     6.8, 4.78, 5.85, colors=[PALE_CYAN, PALE_PURPLE, PALE_RED, PALE_ORANGE, PALE_CYAN],
     box_h=0.72, font_size=10)
add_text(slide, "삭제 대상보다 먼저 readiness guard를 검사하며, catalog는 존재하는 경로만 기록한다.",
         1.05, 6.08, 11.2, 0.4, size=14.2, color=NAVY, bold=True,
         align=PP_ALIGN.CENTER)

# 28 command
slide = new_slide("대표 실행 명령은 하나", "OPS", accent=ORANGE,
                  sources=["run_pipeline.sh", "PORTABLE_GUIDE.md"],
                  talking=["manifest, sample list, input/output root를 명시한다.", "output root를 생략하면 v2/output이 기본이다.", "대량 실행은 wgs-phased-list와 shard wrapper를 권장한다."])
add_code(slide,"""./run_pipeline.sh sample-list \\
  --manifest manifests/preprocess_manifest.tsv \\
  --sample-list sample_lists/all_pair_first_3212.txt \\
  --input-root /data/mfm_input \\
  --output-root /data/mfm_output \\
  --sample-count all""",0.8,2.0,11.75,2.5,size=14)
card(slide,"기본 output 분리", "output/wgs/<project>/<sample>\noutput/wgs_rna_pair/<project>/<sample>",1.0,5.0,5.25,1.1,accent=BLUE,fill=PALE_BLUE,title_size=14,body_size=13)
card(slide,"대량 권장", "wgs-phased-list 또는 run_shard.sh\n성공 sample만 다음 phase로 전달",7.05,5.0,5.25,1.1,accent=ORANGE,fill=PALE_ORANGE,title_size=14,body_size=13)

# 29 fail continuation
slide = new_slide("실패는 즉시 기록하되 batch는 멈추지 않는다", "OPS", accent=ORANGE,
                  sources=["runner_local.sh", "run_wgs_phased.sh", "HIGH_PARALLEL_POLICY.md"],
                  talking=["sample 내부에서는 fail-fast다.", "runner는 실패 sample을 status.tsv에 남기고 다음 sample로 이동한다.", "phased scheduler는 성공 sample만 다음 phase에 전달한다."])
flow(slide,["Sample A PASS","Sample B FAIL","Sample C PASS"],0.9,2.1,11.55,
     colors=[PALE_CYAN,PALE_RED,PALE_CYAN],box_h=0.92,font_size=15)
arrow(slide,2.8,3.55,2.8,4.25,GREEN,2)
arrow(slide,6.65,3.55,6.65,4.25,RED,2)
arrow(slide,10.55,3.55,10.55,4.25,GREEN,2)
card(slide,"다음 phase", "A, C만 진행",1.35,4.4,3.0,1.05,accent=GREEN,fill=PALE_CYAN,title_size=15,body_size=14)
card(slide,"상태 기록", "B: stage·error·resource",5.15,4.4,3.0,1.05,accent=RED,fill=PALE_RED,title_size=15,body_size=14)
card(slide,"다음 sample", "queue는 계속",8.95,4.4,3.0,1.05,accent=ORANGE,fill=PALE_ORANGE,title_size=15,body_size=14)
add_text(slide,"전체 중단이 아니라 ‘문제 sample 격리’가 production 기본 계약",
         1.55,6.0,10.2,0.4,size=17,color=NAVY,bold=True,align=PP_ALIGN.CENTER)

# 30 parallel
slide = new_slide("72-core 정책: sample 수 × 내부 thread를 함께 제한", "OPS", accent=ORANGE,
                  sources=["HIGH_PARALLEL_POLICY.md", "run_wgs_phased.sh", "config_v2.sh"],
                  talking=["예시 예산은 01–10과 11은 4 samples, 12 이후 고부하 단계는 2 samples다.", "RNA는 12 samples×6 cores 예시다.", "실제 값은 CPU뿐 아니라 RAM과 scratch guard가 동적으로 하향한다."])
add_table(slide,["phase","동시 sample","sample CPU","전체"],[
    ["WGS 01–10","4","18","≤72"],
    ["WGS 11","4","18","≤72"],
    ["WGS 12/final","2","36","≤72"],
    ["RNA","12","6","≤72"],
],0.85,2.0,6.0,3.15,col_widths=[2.8,1.5,1.5,1.2],font_size=12.5,accent=ORANGE)
card(slide,"Resource guard", "WGS 64–96 GiB RAM/sample\n임시공간 400 GiB/sample\nRNA 48 GiB·80 GiB/sample",7.35,2.0,4.95,1.75,accent=RED,fill=PALE_RED,title_size=17,body_size=13)
card(slide,"Local scratch", "Stage 10 route/sort와 Stage 12 temp만 local\n최종 output·QC는 OUTPUT_ROOT",7.35,4.2,4.95,1.55,accent=CYAN,fill=PALE_CYAN,title_size=17,body_size=13)

# 31 shard
slide = new_slide("외부 전달은 37개 stratified shard로 나눈다", "OPS", accent=ORANGE,
                  sources=["SHARDED_HANDOFF.md", "scripts/build_stratified_shards.py"],
                  talking=["목표 8.5 TiB, hard cap 10 TiB다.", "WGS-only/pair와 project를 각 shard에 층화한다.", "현재 예상 output은 307.309 TiB이며 pilot 전까지 provisional이다."])
metrics=[("37","shards",ORANGE),("≈87","samples/shard",BLUE),("8.5 TiB","목표",CYAN),("10 TiB","hard cap",RED)]
for idx,(num,label,c) in enumerate(metrics):
    x=0.75+idx*3.1
    rect(slide,x,2.0,2.7,1.3,fill=WHITE,line=c,width=1.6)
    add_text(slide,num,x+0.08,2.18,2.54,0.48,size=24,color=c,bold=True,align=PP_ALIGN.CENTER)
    add_text(slide,label,x+0.08,2.76,2.54,0.25,size=11.5,color=MUTED,align=PP_ALIGN.CENTER)
add_table(slide,["층화 기준","보존 조건"],[["data_group","WGS-only 2,763 + pair 449"],["platform","Illumina / HiFi 비율"],["project","1KGP / HPRC / HGSVC 분포"],["estimated output","모든 shard ≤10 TiB"]],
          1.3,4.0,10.75,2.1,col_widths=[3,5],font_size=12,accent=ORANGE)

# 32 validation
slide = new_slide("검증은 결과 파일 존재 확인보다 훨씬 깊다", "OPS", accent=ORANGE,
                  sources=["tests/", "HIGH_PARALLEL_POLICY.md", "RESULTS_2026-09-18.md"],
                  talking=["sequence, CIGAR, depth, metadata corruption을 의도적으로 넣어 탐지한다.", "최적화 전후 output count/hash/QC를 비교한다.", "실제 pilot의 CPU/RSS/time으로 ETA를 재보정한다."])
levels=[("L1","입력 계약","reference·VCF·sample·checksum",BLUE),("L2","단계 QC","count·tag·mapping·lift",PURPLE),("L3","Compact full trace","sequence·CIGAR·depth·metadata",CYAN),("L4","Pilot 동등성","hash·QC·resource·timing",ORANGE)]
for idx,(level,title,body,c) in enumerate(levels):
    y=1.95+idx*1.08
    pill(slide,level,0.85,y+0.1,0.72,fill=c)
    add_text(slide,title,1.85,y,2.45,0.42,size=16,color=INK,bold=True)
    add_text(slide,body,4.45,y,7.55,0.42,size=14,color=MUTED)
    if idx < len(levels)-1: arrow(slide,1.2,y+0.55,1.2,y+1.02,c,1.5)
add_text(slide,"최근 trace prefetch: 3개 pilot에서 7.7–124.1×, QC/manifest hash 동일",
         1.0,6.25,11.35,0.42,size=15,color=GREEN,bold=True,align=PP_ALIGN.CENTER)

# 33 output for model designers
slide = new_slide("모델 설계자가 실제로 받아야 하는 것", "OUTPUT", accent=GREEN,
                  sources=["output/README.md", "GUIDE.md", "LIBRARY_PROVENANCE.md"],
                  talking=["compact layout 이름만 보고 routing/filter를 알 수 있어야 한다.", "HAP 좌표 reference와 chain fingerprint를 함께 보존한다.", "unknown을 양쪽에 넣은 layout을 HAP-specific label로 오해하면 안 된다."])
cards=[
    ("HAP sequence","개인 HAP1/HAP2 1 Mb window",BLUE),
    ("Readset","SEQ·CIGAR·MAPQ·tags·segments",PURPLE),
    ("Depth","primary depth + long supp evidence",CYAN),
    ("Provenance","reference·VCF·chain·PCR·filter",ORANGE),
    ("RNA pair","coverage·allelic counts·ASE QC",GREEN),
    ("QC","count/hash/corruption/resource JSON",RED),
]
for idx,(t,b,c) in enumerate(cards):
    card(slide,t,b,0.7+(idx%3)*4.1,1.95+(idx//3)*1.65,3.7,1.3,accent=c,fill=WHITE,title_size=15,body_size=12)
card(slide,"사용 전 체크", "layout_base · hap_read_mode · filter_preset · coordinate system · QC PASS",
     1.15,5.45,11.0,0.9,accent=GREEN,fill=PALE_CYAN,title_size=14,body_size=14)

# 34 final checklist
slide = new_slide("운영 체크리스트와 다음 읽을 문서", "WRAP-UP", accent=GREEN,
                  sources=["README.md", "GUIDE.md", "PORTABLE_GUIDE.md", "단계별_문헌근거.md"],
                  talking=["manifest와 exact reference를 먼저 고정한다.", "QC를 통과한 sample만 다음 단계로 보낸다.", "모델 입력은 layout metadata를 확인한 뒤 사용한다."])
add_bullets(slide,[
    "manifest sample·project·platform·PCR provenance 확인",
    "reference/VCF/alignment 22-autosome preflight PASS",
    "HAP routing·filter preset·좌표계가 목적과 일치",
    "Compact full validation 및 RNA ASE QC PASS",
    "resource·scratch·output 용량 확인 후 batch 시작",
],0.8,2.0,6.0,3.8,size=15.5,accent=GREEN,gap=11)
card(slide,"운영", "README.md\nGUIDE.md\nPORTABLE_GUIDE.md",7.35,2.0,2.25,1.8,accent=BLUE,fill=PALE_BLUE,title_size=16,body_size=13)
card(slide,"근거", "단계별_문헌근거.md\nLIBRARY_PROVENANCE.md",9.95,2.0,2.25,1.8,accent=PURPLE,fill=PALE_PURPLE,title_size=16,body_size=13)
card(slide,"성능·전달", "HIGH_PARALLEL_POLICY.md\nSHARDED_HANDOFF.md",7.35,4.3,4.85,1.35,accent=ORANGE,fill=PALE_ORANGE,title_size=16,body_size=13)
add_text(slide,"최종 원칙: 검증은 강하게, 실패 범위는 sample로 작게, provenance는 끝까지 유지",
         1.05,6.25,11.2,0.43,size=17,color=NAVY,bold=True,align=PP_ALIGN.CENTER)


def write_notes() -> None:
    lines = ["# 전처리 파이프라인 교육자료 슬라이드 원고", "",
             "이 문서는 PPT의 발표자 노트와 근거 파일을 텍스트로도 보존한다.", ""]
    for row in notes_rows:
        lines.extend([f"## {int(row['number']):02d}. {row['title']}", "", "### 설명"])
        lines.extend(f"- {item}" for item in row["talking"])
        lines.extend(["", "### 근거"])
        lines.extend(f"- `{item}`" for item in row["sources"])
        lines.append("")
    NOTES_OUT.write_text("\n".join(lines), encoding="utf-8")
    with MAP_OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["slide", "title", "source"])
        for row in notes_rows:
            for source in row["sources"]:
                writer.writerow([row["number"], row["title"], source])


prs.core_properties.title = "개인 HAP 기반 WGS·RNA 전처리 파이프라인"
prs.core_properties.subject = "MFM preprocessing v2 교육자료"
prs.core_properties.author = "MFM preprocessing team"
prs.core_properties.keywords = "WGS, RNA, haplotype, BWA-MEM2, minimap2, STAR, Zarr"
prs.core_properties.comments = "2026-10-02 v2 코드·문서 기준 자동 생성"
prs.save(OUT)
write_notes()
print(f"pptx={OUT} slides={len(prs.slides)}")
print(f"notes={NOTES_OUT}")
print(f"source_map={MAP_OUT}")
