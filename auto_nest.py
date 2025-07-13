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

        iterations = None
        nums = []
        dxf_files = []

        i = s_idx + 2
        while i < len(sys.argv):
            if sys.argv[i] == "-i":
                iterations = sys.argv[i + 1]
                i += 2
            elif sys.argv[i] == "-n":
                i += 1
                while i < len(sys.argv):
                    arg = sys.argv[i]
                    if arg.startswith("-"):
                        break
                    try:
                        nums.append(int(arg))
                    except ValueError:
                        break
                    i += 1
            else:
                dxf_files.append(sys.argv[i])
                i += 1

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
    if nums:
        nest_cmd += ["-n"] + list(map(str, nums))
    nest_cmd += ["-o", "lay.csv"] + json_files

    print("[PY] RUN:", " ".join(nest_cmd))
    subprocess.check_call(nest_cmd)
    print("[PY] Success!")
