import sys
import json
import math
import ezdxf
from ezdxf.units import InsertUnits, conversion_factor


def dist(p, q):
    return math.hypot(p[0]-q[0], p[1]-q[1])


def extract_all_paths(dxf_file: str, tol: float = 0.2):
    doc = ezdxf.readfile(dxf_file)
    insunits = doc.header.get("$INSUNITS", 0)
    try:
        factor = conversion_factor(InsertUnits(insunits), InsertUnits.Millimeters)
    except Exception:
        factor = 1.0

    msp = doc.modelspace()
    paths = []
    autoclosed = 0
    min_x = float("inf")
    min_y = float("inf")

    for e in msp:
        t = e.dxftype()
        path = []
        if t == "LINE":
            p1 = e.dxf.start
            p2 = e.dxf.end
            path = [(p1.x * factor, p1.y * factor), (p2.x * factor, p2.y * factor)]
        elif t in ("LWPOLYLINE", "POLYLINE"):
            for x, y in e.get_points("xy"):
                path.append((x * factor, y * factor))
        elif t == "SPLINE":
            for v in e.flattening(tol * 0.5 / factor):
                path.append((v.x * factor, v.y * factor))
        elif t in ("ARC", "CIRCLE"):
            for v in e.flattening(tol * 0.5 / factor):
                path.append((v.x * factor, v.y * factor))
        else:
            continue

        if len(path) == 0:
            continue

        if dist(path[0], path[-1]) <= tol:
            path[-1] = path[0]
            autoclosed += 1

        for x, y in path:
            if x < min_x:
                min_x = x
            if y < min_y:
                min_y = y
        paths.append(path)

    for path in paths:
        for i, (x, y) in enumerate(path):
            path[i] = [x - min_x, y - min_y]

    return paths, autoclosed


def main(argv=None):
    argv = argv or sys.argv[1:]
    if len(argv) != 2:
        print("Usage: python extract_all_shapes.py <input.dxf> <output.json>")
        return 1
    in_dxf, out_json = argv
    paths, auto = extract_all_paths(in_dxf, tol=0.2)
    with open(out_json, "w") as f:
        json.dump(paths, f, indent=0)
    print(f"Saved {len(paths)} paths (autoclosed {auto}) to {out_json}")


if __name__ == "__main__":
    sys.exit(main())
