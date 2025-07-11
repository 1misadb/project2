import ezdxf
import json
import math
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union

def arc_to_polyline(entity, segments=32):
    center = entity.dxf.center
    radius = entity.dxf.radius
    start_angle = math.radians(entity.dxf.start_angle)
    end_angle = math.radians(entity.dxf.end_angle)
    if end_angle < start_angle:
        end_angle += 2 * math.pi
    angle_range = end_angle - start_angle
    points = []
    for i in range(segments + 1):
        angle = start_angle + (angle_range * i) / segments
        x = center.x + radius * math.cos(angle)
        y = center.y + radius * math.sin(angle)
        points.append([float(x), float(y)])
    return points

def entity_to_polyline(entity, spline_tol=0.5):
    if entity.dxftype() == 'LINE':
        s, e = entity.dxf.start, entity.dxf.end
        return [[float(s.x), float(s.y)], [float(e.x), float(e.y)]]
    elif entity.dxftype() == 'LWPOLYLINE':
        return [[float(x), float(y)] for x, y, *rest in entity]
    elif entity.dxftype() == 'POLYLINE':
        return [[float(v.dxf.location.x), float(v.dxf.location.y)] for v in entity.vertices]
    elif entity.dxftype() == 'SPLINE':
        return [[float(x), float(y)] for x, y, *_ in entity.flattening(spline_tol)]
    elif entity.dxftype() == 'ARC':
        return arc_to_polyline(entity)
    elif entity.dxftype() == 'CIRCLE':
        class Dummy: pass
        dummy = Dummy()
        dummy.dxf = entity.dxf
        dummy.dxf.start_angle = 0
        dummy.dxf.end_angle = 360
        return arc_to_polyline(dummy)
    return None

def autoclose_path(path, tol=0.2):
    """Если концы близко — замыкаем путь"""
    if len(path) > 2:
        x0, y0 = path[0]
        x1, y1 = path[-1]
        dist = ((x0 - x1)**2 + (y0 - y1)**2) ** 0.5
        if dist < tol:
            path[-1] = [x0, y0]
    return path

def merge_paths_into_shapes(all_paths):
    # Переводим пути в shapely-полигоны (игнорируем короткие)
    polys = []
    for path in all_paths:
        if len(path) < 3:
            continue
        try:
            poly = Polygon(path)
            if poly.is_valid:
                polys.append(poly)
        except Exception as e:
            print('Ошибка poly:', e)
    if not polys:
        print("No valid polygons found!")
        return []

    united = unary_union(polys)

    out = []
    if isinstance(united, Polygon):
        one = [list(map(list, united.exterior.coords))]
        for hole in united.interiors:
            one.append(list(map(list, hole.coords)))
        out.append(one)
    elif isinstance(united, MultiPolygon):
        for poly in united.geoms:
            one = [list(map(list, poly.exterior.coords))]
            for hole in poly.interiors:
                one.append(list(map(list, hole.coords)))
            out.append(one)
    else:
        print("Unknown geometry:", type(united))
    return out

def extract_all_shapes(dxf_path, out_json="parts.json", spline_tol=0.5):
    doc = ezdxf.readfile(dxf_path)
    INSUNITS = doc.header.get('$INSUNITS', 0)
    SCALE = {1: 25.4, 2: 25.4*12, 3: 25.4*12*5280, 4: 1,
             5: 10, 6: 1000}.get(INSUNITS, 1)
    msp = doc.modelspace()
    all_paths = []
    for entity in msp:
        poly = entity_to_polyline(entity, spline_tol=spline_tol)
        if poly and len(poly) > 1:
            scaled = [[x * SCALE, y * SCALE] for x, y in poly]
            closed = autoclose_path(scaled, tol=0.2)
            all_paths.append(closed)
    if not all_paths:
        print("No paths found!")
        return
    min_x = min(pt[0] for path in all_paths for pt in path)
    min_y = min(pt[1] for path in all_paths for pt in path)
    shifted_paths = []
    for path in all_paths:
        shifted = [[x - min_x, y - min_y] for x, y in path]
        shifted_paths.append(shifted)
    # Просто сохраняй каждый контур как отдельную деталь
    with open(out_json, "w") as f:
        json.dump([ [path] for path in shifted_paths ], f)
    print(f"Extracted {len(shifted_paths)} polylines as shapes from {dxf_path} → {out_json}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("dxf_file")
    parser.add_argument("-o", "--output", default="parts.json")
    args = parser.parse_args()
    extract_all_shapes(args.dxf_file, out_json=args.output)
