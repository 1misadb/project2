import argparse
import json
import math
import sys
import ezdxf
from ezdxf.units import InsertUnits, conversion_factor
from config_util import load_config

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
    # --- Собираем все примитивы как отдельные сегменты ---
    segments = []
    for e in msp:
        t = e.dxftype()
        if t == "LINE":
            p1 = e.dxf.start
            p2 = e.dxf.end
            segments.append([[p1.x * factor, p1.y * factor], [p2.x * factor, p2.y * factor]])
        elif t == "LWPOLYLINE":
            pts = [[x * factor, y * factor] for x, y in e.get_points("xy")]
            if len(pts) > 1:
                closed = dist(pts[0], pts[-1]) <= tol or getattr(e, 'closed', False)
                for i in range(len(pts) - 1):
                    segments.append([pts[i], pts[i+1]])
                if closed:
                    segments.append([pts[-1], pts[0]])
        elif t == "POLYLINE":
            pts = [[v.dxf.location.x * factor, v.dxf.location.y * factor] for v in e.vertices]
            if len(pts) > 1:
                closed = dist(pts[0], pts[-1]) <= tol or getattr(e, 'is_closed', False)
                for i in range(len(pts) - 1):
                    segments.append([pts[i], pts[i+1]])
                if closed:
                    segments.append([pts[-1], pts[0]])
        elif t == "SPLINE":
            pts = [[v.x * factor, v.y * factor] for v in e.flattening(tol * 0.1 / factor)]
            for i in range(len(pts) - 1):
                segments.append([pts[i], pts[i+1]])
        elif t in ("ARC", "CIRCLE"):
            pts = [[v.x * factor, v.y * factor] for v in e.flattening(tol * 0.1 / factor)]
            for i in range(len(pts) - 1):
                segments.append([pts[i], pts[i+1]])

    # --- Автоматическое объединение сегментов в замкнутые контуры ---
    def join_segments(segments, tol=0.1):
        paths = []
        used = [False] * len(segments)
        for i, seg in enumerate(segments):
            if used[i]:
                continue
            path = seg[:]
            used[i] = True
            changed = True
            while changed:
                changed = False
                for j, seg2 in enumerate(segments):
                    if used[j]:
                        continue
                    # Совпадает конец path и начало seg2
                    if dist(path[-1], seg2[0]) <= tol:
                        path.extend(seg2[1:])
                        used[j] = True
                        changed = True
                        break
                    # Совпадает конец path и конец seg2
                    if dist(path[-1], seg2[-1]) <= tol:
                        path.extend(reversed(seg2[:-1]))
                        used[j] = True
                        changed = True
                        break
                    # Совпадает начало path и конец seg2
                    if dist(path[0], seg2[-1]) <= tol:
                        path = seg2[:-1] + path
                        used[j] = True
                        changed = True
                        break
                    # Совпадает начало path и начало seg2
                    if dist(path[0], seg2[0]) <= tol:
                        path = list(reversed(seg2[1:])) + path
                        used[j] = True
                        changed = True
                        break
            # Если замкнуто
            if len(path) > 2 and dist(path[0], path[-1]) <= tol:
                path[-1] = path[0]
                paths.append(path)
        return paths

    all_paths = join_segments(segments, tol)

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

def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Extract details from DXF")
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("repeat", nargs="?", default="repeat=1")
    parser.add_argument("--config", default="config.yaml")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    cfg = load_config(args.config)
    in_dxf = args.input
    out_json = args.output
    repeat = 1
    if isinstance(args.repeat, str) and args.repeat.startswith("repeat="):
        try:
            repeat = int(args.repeat.split("=", 1)[1])
        except Exception:
            print("Invalid repeat argument, must be repeat=N")
            return 1
    tol = cfg.get("tolerance", 0.05)
    details = extract_all_details(in_dxf, tol=tol)
    if repeat > 1:
        details = details * repeat
    
    # === БЛОК ДЛЯ ВЫВОДА ГАБАРИТОВ КАЖДОЙ ДЕТАЛИ ===
    print("\nГабариты деталей:")
    for idx, detail_paths in enumerate(details):
        xs = []
        ys = []
        for path in detail_paths:
            xs.extend([pt[0] for pt in path])
            ys.extend([pt[1] for pt in path])
        min_x = min(xs)
        max_x = max(xs)
        min_y = min(ys)
        max_y = max(ys)
        width = max_x - min_x
        height = max_y - min_y
        print(f"Деталь {idx+1}: X от {min_x:.2f} до {max_x:.2f} (ширина {width:.2f} мм), "
            f"Y от {min_y:.2f} до {max_y:.2f} (высота {height:.2f} мм)")
    with open(out_json, "w") as f:
        json.dump(details, f, indent=0)
    print(f"Saved {len(details)} details to {out_json} (repeat={repeat})")

if __name__ == "__main__":
    raise SystemExit(main())
