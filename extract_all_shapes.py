import sys
import json
import math
import ezdxf
from ezdxf.units import InsertUnits, conversion_factor

def dist(p, q):
    return math.hypot(p[0] - q[0], p[1] - q[1])

def extract_all_details(dxf_file: str, tol: float = 0.1):
    doc = ezdxf.readfile(dxf_file)
    insunits = doc.header.get("$INSUNITS", 0)
    try:
        factor = conversion_factor(InsertUnits(insunits), InsertUnits.Millimeters)
    except Exception:
        factor = 1.0

    msp = doc.modelspace()
    # --- Группируем все полигоны в детали ---
    # Этот простой вариант: каждая LWPOLYLINE или POLYLINE считается отдельной деталью.
    # Для реальных DXF лучше доработать группировку!
    details = []
    for e in msp:
        t = e.dxftype()
        path = []
        if t == "LINE":
            p1 = e.dxf.start
            p2 = e.dxf.end
            path = [[p1.x * factor, p1.y * factor], [p2.x * factor, p2.y * factor]]
        elif t == "LWPOLYLINE":
            for x, y in e.get_points("xy"):
                path.append([x * factor, y * factor])
        elif t == "POLYLINE":
            for v in e.vertices:
                x, y = v.dxf.location.x, v.dxf.location.y
                path.append([x * factor, y * factor])
        elif t == "SPLINE":
            for v in e.flattening(tol * 0.1 / factor):
                path.append([v.x * factor, v.y * factor])
        elif t in ("ARC", "CIRCLE"):
            for v in e.flattening(tol * 0.1 / factor):
                path.append([v.x * factor, v.y * factor])
        else:
            continue

        if not path:
            continue

        if dist(path[0], path[-1]) <= tol:
            path[-1] = path[0]

        # Можно доработать: искать отверстия, группировать по слоям и т.д.
        details.append([path])  # Каждый path — отдельная деталь (outer + [holes])

    # --- Сдвигаем каждую деталь в (0,0) ---
    for detail_paths in details:
        all_points = [pt for path in detail_paths for pt in path]
        min_x = min(pt[0] for pt in all_points)
        min_y = min(pt[1] for pt in all_points)
        for path in detail_paths:
            for i, (x, y) in enumerate(path):
                path[i] = [x - min_x, y - min_y]

    return details

def main(argv=None):
    argv = argv or sys.argv[1:]
    if len(argv) != 2:
        print("Usage: python extract_multi_details.py <input.dxf> <output.json>")
        return 1
    in_dxf, out_json = argv
    details = extract_all_details(in_dxf, tol=0.01)
    with open(out_json, "w") as f:
        json.dump(details, f, indent=0)
    print(f"Saved {len(details)} details to {out_json}")

if __name__ == "__main__":
    sys.exit(main())
