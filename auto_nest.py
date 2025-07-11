import subprocess
import sys
import os

def usage():
    print("Пример использования:")
    print("  python auto_nest.py -s 2000x2000 [-i N] file1.dxf file2.dxf ...")
    sys.exit(1)

if __name__ == "__main__":
    if "-s" not in sys.argv or len(sys.argv) < 4:
        usage()

    s_idx = sys.argv.index("-s")
    try:
        sheet = sys.argv[s_idx + 1]

        # Проверим, есть ли -i N
        if "-i" in sys.argv:
            i_idx = sys.argv.index("-i")
            iterations = sys.argv[i_idx + 1]
            # DXF файлы = всё после sheet, без -i и N
            dxf_files = sys.argv[s_idx + 2:i_idx] + sys.argv[i_idx + 2:]
        else:
            iterations = None
            dxf_files = sys.argv[s_idx + 2:]

        if not dxf_files:
            usage()

    except Exception:
        usage()

    json_files = []
    for dxf in dxf_files:
        json_out = os.path.splitext(dxf)[0] + ".json"
        print(f"[PY] DXF → JSON: {dxf} -> {json_out}")

        # ✅ Путь к твоему Python39
        python39_path = r"C:\Users\User\AppData\Local\Programs\Python\Python39\python.exe"

        # ✅ Вызов конвертера
        subprocess.check_call([python39_path, "extract_all_shapes.py", dxf, json_out])

        json_files.append(json_out)

    # Формируем команду nest.exe
    nest_cmd = ["nest.exe", "-s", sheet]

    if iterations is not None:
        nest_cmd += ["-i", iterations]

    nest_cmd += ["-o", "lay.csv"] + json_files

    print("[PY] RUN:", " ".join(nest_cmd))
    subprocess.check_call(nest_cmd)
    print("[PY] Success!")
