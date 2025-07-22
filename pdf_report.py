import csv
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4

def generate_report(csv_file, pdf_file):
    rows = []
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    c = canvas.Canvas(pdf_file, pagesize=A4)
    y = A4[1] - 40
    c.drawString(40, y, "Nesting Report")
    y -= 20
    if rows:
        fill = float(rows[0].get('fill_pct',0))
        waste = float(rows[0].get('waste_mm2',0))
        c.drawString(40, y, f"Fill: {fill:.2f}%  Waste: {waste:.2f} mm2")
        y -= 20
    for r in rows:
        txt = f"run {r['run']} part {r['part']} x={r['x_mm']} y={r['y_mm']} angle={r['angle']}"
        c.drawString(40, y, txt)
        y -= 12
        if y < 40:
            c.showPage()
            y = A4[1] - 40
    c.save()

if __name__ == '__main__':
    import sys
    if len(sys.argv)!=3:
        print('usage: pdf_report.py lay.csv report.pdf')
        sys.exit(1)
    generate_report(sys.argv[1], sys.argv[2])
