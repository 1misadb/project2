import ezdxf
import json
import math

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
    """Returns list of [ [x, y], ... ] for the entity, or None if not supported."""
    if entity.dxftype() == 'LINE':
        s, e = entity.dxf.start, entity.dxf.end
        return [[float(s.x), float(s.y)], [float(e.x), float(e.y)]]
    elif entity.dxftype() in ('LWPOLYLINE', 'POLYLINE'):
        return [[float(x), float(y)] for x, y, *_ in entity.get_points()]
    elif entity.dxftype() == 'SPLINE':
        return [[float(x), float(y)] for x, y, *_ in entity.approximate(spline_tol)]
    elif entity.dxftype() == 'ARC':
        return arc_to_polyline(entity)
    elif entity.dxftype() == 'CIRCLE':
        # Тот же подход, что и для ARC, только полный круг
        class Dummy: pass
        dummy = Dummy()
        dummy.dxf = entity.dxf
        dummy.dxf.start_angle = 0
        dummy.dxf.end_angle = 360
        return arc_to_polyline(dummy)
    return None

def extract_all_paths(dxf_path, out_json="paths.json", spline_tol=0.5):
    doc = ezdxf.readfile(dxf_path)
    msp = doc.modelspace()
    all_paths = []
    for entity in msp:
        poly = entity_to_polyline(entity, spline_tol=spline_tol)
        if poly and len(poly) > 1:
            all_paths.append(poly)
    with open(out_json, "w") as f:
        json.dump(all_paths, f)
    print(f"Extracted {len(all_paths)} polylines from {dxf_path} → {out_json}")

if __name__ == "__main__":
    # Использование: python thisfile.py input.dxf
    import sys
    if len(sys.argv) < 2:
        print("Usage: python extract_polylines.py input.dxf [output.json]")
        sys.exit(1)
    dxf_in = sys.argv[1]
    json_out = sys.argv[2] if len(sys.argv) > 2 else "paths.json"
    extract_all_paths(dxf_in, out_json=json_out)
