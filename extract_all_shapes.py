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
    # --- Собираем все замкнутые контуры ---
    all_paths = []
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

        # Только замкнутые контуры
        if len(path) > 2 and dist(path[0], path[-1]) <= tol:
            all_paths.append(path)

    # --- Группируем по вложенности: внешний контур + отверстия ---
    def point_in_polygon(pt, poly):
        # Четное-нечетное правило
        x, y = pt
        inside = False
        n = len(poly)
        for i in range(n):
            x0, y0 = poly[i]
            x1, y1 = poly[(i+1)%n]
            if ((y0 > y) != (y1 > y)):
                xinters = (x1 - x0) * (y - y0) / (y1 - y0 + 1e-12) + x0
                if x < xinters:
                    inside = not inside
        return inside

    # Для каждого контура ищем, сколько других его содержат
    parents = [None] * len(all_paths)
    for i, path in enumerate(all_paths):
        for j, other in enumerate(all_paths):
            if i == j:
                continue
            # Проверяем, лежит ли точка path внутри other
            if point_in_polygon(path[0], other):
                # Если уже есть родитель, выбираем самый "маленький" (наиболее вложенный)
                if parents[i] is None or len(other) < len(all_paths[parents[i]]):
                    parents[i] = j

    # Строим дерево: внешний контур (parent=None) + все вложенные в него (holes)
    details = []
    used = set()
    for i, parent in enumerate(parents):
        if parent is not None:
            continue  # не внешний контур
        # Собираем все, у кого parent == i (отверстия)
        holes = [all_paths[j] for j, p in enumerate(parents) if p == i]
        detail = [all_paths[i]] + holes
        details.append(detail)
        used.add(i)
        used.update(j for j, p in enumerate(parents) if p == i)

    # Если остались неиспользованные (например, вложенность >2), добавляем их как отдельные детали
    for i in range(len(all_paths)):
        if i not in used:
            details.append([all_paths[i]])

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
    if not (2 <= len(argv) <= 3):
        print("Usage: python extract_multi_details.py <input.dxf> <output.json> [repeat=N]")
        return 1
    in_dxf, out_json = argv[:2]
    repeat = 1
    if len(argv) == 3:
        arg = argv[2]
        if arg.startswith("repeat="):
            try:
                repeat = int(arg.split("=", 1)[1])
            except Exception:
                print("Invalid repeat argument, must be repeat=N")
                return 1
    details = extract_all_details(in_dxf, tol=0.01)
    if repeat > 1:
        details = details * repeat
    with open(out_json, "w") as f:
        json.dump(details, f, indent=0)
    print(f"Saved {len(details)} details to {out_json} (repeat={repeat})")

if __name__ == "__main__":
    sys.exit(main())
