import argparse
import csv
from pathlib import Path
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.lib.utils import ImageReader
import matplotlib.pyplot as plt
import ezdxf
import io

def _plot_layout(dxf_file, fill, waste):
    doc = ezdxf.readfile(dxf_file)
    msp = doc.modelspace()
    fig, ax = plt.subplots()
    for e in msp:
        if e.dxftype() == 'LWPOLYLINE':
            pts = [(p[0], p[1]) for p in e.get_points('xy')]
            if len(pts) < 2:
                continue
            xs, ys = zip(*pts)
            ax.plot(xs, ys)
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set_title(f'Fill {fill:.1f}%  Waste {waste:.1f}mm2')
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    plt.close(fig)
    buf.seek(0)
    return buf

def generate_report(csv_file, pdf_file, dxf_file=None):
    rows = []
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    c = canvas.Canvas(pdf_file, pagesize=A4)

    # Summary table
    y = A4[1] - 40
    c.drawString(40, y, "Nesting Report")
    y -= 20
    summary = {}
    for r in rows:
        k = (r['run'], r['strategy'])
        if k not in summary:
            summary[k] = r
    c.drawString(40, y, "Run  Strategy  Fill%  Waste")
    y -= 15
    best = None
    for k, r in sorted(summary.items()):
        fill = float(r.get('fill_pct',0))
        waste = float(r.get('waste_mm2',0))
        c.drawString(40, y, f"{r['run']}    {r['strategy']}    {fill:.1f}    {waste:.1f}")
        if best is None or waste < float(best.get('waste_mm2',1e9)):
            best = r
        y -= 12
        if y < 60:
            c.showPage()
            y = A4[1] - 40

    # Layout page
    if dxf_file:
        buf = _plot_layout(dxf_file, float(best.get('fill_pct',0)), float(best.get('waste_mm2',0)))
        img = ImageReader(buf)
        c.showPage()
        c.drawImage(img, 40, 80, width=A4[0]-80, height=A4[1]-160)
        buf.close()

    c.save()

def parse_args(argv=None):
    parser = argparse.ArgumentParser(description='Generate PDF report')
    parser.add_argument('csv')
    parser.add_argument('pdf')
    parser.add_argument('--dxf', default=None)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    dxf_file = args.dxf
    generate_report(args.csv, args.pdf, dxf_file)


if __name__ == '__main__':
    raise SystemExit(main())
