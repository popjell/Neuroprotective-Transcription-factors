"""Export the short HTML results report as a three-page PDF."""

from html import escape
from pathlib import Path
import sys

from lxml import html
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph, Table, TableStyle
from PIL import Image


ROOT = Path(__file__).resolve().parent
RANK = sys.argv[1] if len(sys.argv) > 1 else "rank2"
assert RANK in ("rank1", "rank2")
SOURCE = ROOT / ("scRNAseq_rank1_results.html" if RANK == "rank1" else "scRNAseq_key_results.html")
OUTPUT = ROOT / f"output/pdf/neuron_{RANK}_results.pdf"
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
page = html.parse(str(SOURCE))
figures = page.xpath("//figure")
assert len(figures) == 8
score_tables = page.xpath("//table[@class='singscore-table']")

W, H = A4
LEFT = 44
RIGHT = W - 44
CONTENT = RIGHT - LEFT
INK = colors.HexColor("#24333F")
MUTED = colors.HexColor("#52616D")
RULE = colors.HexColor("#CBD3D8")
pdfmetrics.registerFont(TTFont("LiberationSans", "/usr/share/fonts/liberation/LiberationSans-Regular.ttf"))
pdfmetrics.registerFont(TTFont("LiberationSans-Bold", "/usr/share/fonts/liberation/LiberationSans-Bold.ttf"))

body_style = ParagraphStyle("body", fontName="LiberationSans", fontSize=9.5,
                            leading=14, textColor=INK)
caption_style = ParagraphStyle("caption", parent=body_style, fontSize=8.7,
                               leading=12, textColor=MUTED)
small_style = ParagraphStyle("small", parent=body_style, fontSize=8.4,
                             leading=11.7, textColor=MUTED)
header_style = ParagraphStyle("table_header", parent=small_style,
                              fontName="LiberationSans-Bold", textColor=INK)
cell_style = ParagraphStyle("table_cell", parent=body_style, fontSize=8.8,
                            leading=11.5)


def clean(value):
    return " ".join(value.replace("−", "-").replace("×", "x")
                    .replace("·", "-").split())


def node_text(node):
    return clean(node.text_content())


def paragraph(pdf, value, y, style=body_style, width=CONTENT):
    item = Paragraph(escape(clean(value)), style)
    _, height = item.wrap(width, H)
    item.drawOn(pdf, LEFT, y - height)
    return y - height


def heading(pdf, value, y):
    pdf.setFont("LiberationSans-Bold", 13)
    pdf.setFillColor(INK)
    pdf.drawString(LEFT, y, clean(value))
    pdf.setStrokeColor(RULE)
    pdf.line(LEFT, y - 6, RIGHT, y - 6)
    return y - 22


def figure(pdf, node, y, width):
    path = ROOT / node.xpath("./img/@src")[0]
    with Image.open(path) as image:
        image_width, image_height = image.size
    height = width * image_height / image_width
    x = LEFT + (CONTENT - width) / 2
    pdf.drawImage(str(path), x, y - height, width=width, height=height,
                  preserveAspectRatio=True, mask="auto")
    y -= height + 7
    return paragraph(pdf, node_text(node.xpath("./figcaption")[0]), y,
                     caption_style)


def footer(pdf, number):
    pdf.setStrokeColor(RULE)
    pdf.line(LEFT, 37, RIGHT, 37)
    pdf.setFont("LiberationSans", 8)
    pdf.setFillColor(MUTED)
    pdf.drawString(LEFT, 24, "Schirmer et al., 2019  |  Trobisch et al., 2022")
    pdf.drawRightString(RIGHT, 24, str(number))


def score_table(pdf, node, y):
    rows = []
    for row in node.xpath(".//tr"):
        style = header_style if not rows else small_style
        rows.append([Paragraph(escape(node_text(cell)), style)
                     for cell in row.xpath("./th|./td")])
    table = Table(rows, colWidths=[139, 76, 42, 65, 47, 47, 91], hAlign="LEFT")
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#EEF2F4")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LINEBELOW", (0, 0), (-1, -1), 0.3, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    _, height = table.wrap(CONTENT, H)
    table.drawOn(pdf, LEFT, y - height)
    return y - height


pdf = canvas.Canvas(str(OUTPUT), pagesize=A4)
pdf.setTitle(node_text(page.xpath("//h1")[0]))
pdf.setAuthor("Neuronal MS analysis")

# Page 1: numerical results and donor-level score.
y = H - 49
title = node_text(page.xpath("//h1")[0])
pdf.setFont("LiberationSans-Bold", 19)
pdf.setFillColor(INK)
assert stringWidth(title, "LiberationSans-Bold", 19) < CONTENT
pdf.drawString(LEFT, y, title)
y -= 23
y = paragraph(pdf, node_text(page.xpath("//p[@class='subtitle']")[0]), y,
              small_style) - 21
y = heading(pdf, "Gene-set tests", y)

rows = []
for row in page.xpath("(//table)[1]//tr"):
    style = header_style if not rows else cell_style
    rows.append([Paragraph(escape(node_text(cell)), style)
                 for cell in row.xpath("./th|./td")])
table = Table(rows, colWidths=[151, 57, 65, 89, 72, 73], hAlign="LEFT")
table.setStyle(TableStyle([
    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#EEF2F4")),
    ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ("LINEBELOW", (0, 0), (-1, -1), 0.35, RULE),
    ("LEFTPADDING", (0, 0), (-1, -1), 5),
    ("RIGHTPADDING", (0, 0), (-1, -1), 5),
    ("TOPPADDING", (0, 0), (-1, -1), 6),
    ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
]))
_, table_height = table.wrap(CONTENT, H)
table.drawOn(pdf, LEFT, y - table_height)
y -= table_height + 15
summary = page.xpath("//table/following-sibling::p[1]")[0]
y = paragraph(pdf, node_text(summary), y) - 24
y = heading(pdf, f"Donor-level Rank {RANK[-1]} scores", y)
y = figure(pdf, figures[0], y, 456)
assert y > 56, f"Page 1 content reached footer: {y}"
footer(pdf, 1)
pdf.showPage()

# Page 2: both annotation maps.
y = heading(pdf, "Cell annotations", H - 56)
y = figure(pdf, figures[1], y, 445) - 18
y = figure(pdf, figures[2], y, 445)
assert y > 56, f"Page 2 content reached footer: {y}"
footer(pdf, 2)
pdf.showPage()

# Page 3: both volcano plots and methods.
y = heading(pdf, "Gene-level differential expression", H - 56)
y = figure(pdf, figures[3], y, 410) - 10
y = figure(pdf, figures[4], y, 410) - 10
methods = node_text(page.xpath("//p[@class='methods']")[0])
methods = methods.split("Detailed report")[0].strip()
assert y > 55, f"Page 3 content reached footer: {y}"
footer(pdf, 3)
pdf.showPage()
y = heading(pdf, "Singscore across all neuronal subtypes", H - 56)
note = "One point per donor-subtype profile. Trobisch includes cortex and cerebellum; tests adjust for region. Panel y-axes vary by subtype."
y = paragraph(pdf, note, y, small_style) - 12
y = figure(pdf, figures[5], y, 470)
assert y > 55, f"Page 4 content reached footer: {y}"
footer(pdf, 4)
pdf.showPage()

y = heading(pdf, "Trobisch: other neuronal subtypes", H - 56)
y = figure(pdf, figures[6], y, 470)
assert y > 55, f"Page 5 content reached footer: {y}"
footer(pdf, 5)
pdf.showPage()

y = heading(pdf, f"Trobisch: Rank {RANK[-1]} singscore tests", H - 56)
y = score_table(pdf, score_tables[0], y)
y -= 16
y = paragraph(pdf, "q values are corrected within singscore across testable Rank 1 and Rank 2 subtype tests.", y, small_style)
assert y > 55, f"Page 6 content reached footer: {y}"
footer(pdf, 6)
pdf.showPage()

y = heading(pdf, f"Schirmer: Rank {RANK[-1]} singscore results", H - 56)
y = figure(pdf, figures[7], y, 470)
assert y > 55, f"Page 7 content reached footer: {y}"
footer(pdf, 7)
pdf.showPage()

y = heading(pdf, f"Schirmer: Rank {RANK[-1]} singscore tests", H - 56)
y = score_table(pdf, score_tables[1], y) - 20
y = heading(pdf, "Methods", y)
y = paragraph(pdf, methods, y, small_style)
assert y > 55, f"Page 8 content reached footer: {y}"
footer(pdf, 8)
pdf.save()
print(OUTPUT)
