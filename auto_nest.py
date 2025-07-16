import subprocess
import sys
import os
sys.stdout.reconfigure(encoding='utf-8')

def usage():
    print("Пример использования:")
    print("  python auto_nest.py -s 2000x2000 [-i N] [--флаги] part1.dxf:3 part2.dxf:1 ...")
    sys.exit(1)

if __name__ == "__main__":
    if "-s" not in sys.argv or len(sys.argv) < 4:
        usage()

    # --- добавляем список известных флагов, которые прокидываем в nest.exe ---
    allowed_flags = {"-v", "--verbose", "-j", "--join-segments"}  # можно дополнять

    s_idx = sys.argv.index("-s")
    try:
        sheet = sys.argv[s_idx + 1]

        iterations = None
        nums = []
        dxf_files = []  # список кортежей (имя, количество)
        nest_flags = [] # сюда кидаем "-v", "-j" и любые другие поддерживаемые

        i = s_idx + 2
        while i < len(sys.argv):
            arg = sys.argv[i]
            if arg == "-i":
                iterations = sys.argv[i + 1]
                i += 2
            elif arg == "-n":
                i += 1
                while i < len(sys.argv):
                    arg2 = sys.argv[i]
                    if arg2.startswith("-"):
                        break
                    try:
                        nums.append(int(arg2))
                    except ValueError:
                        break
                    i += 1
            elif arg in allowed_flags:
                nest_flags.append(arg)
                i += 1
            else:
                # Файл
                if ":" in arg:
                    fname, nstr = arg.rsplit(":", 1)
                    try:
                        n = int(nstr)
                    except Exception:
                        print(f"[ERROR] Некорректный формат количества: {arg}")
                        usage()
                    dxf_files.append((fname, n))
                else:
                    dxf_files.append((arg, 1))
                i += 1

        if not dxf_files:
            usage()

    except Exception:
        usage()

    json_files = []
    python_exec = r"C:\Users\User\AppData\Local\Programs\Python\Python39\python.exe"
    for dxf, repeat in dxf_files:
        json_out = os.path.splitext(dxf)[0] + ".json"
        print(f"[PY] DXF → JSON: {dxf} -> {json_out} (repeat={repeat})")
        cmd = [python_exec, "extract_all_shapes.py", dxf, json_out]
        if repeat > 1:
            cmd.append(f"repeat={repeat}")
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Failed to convert {dxf}: {e}")
            sys.exit(1)
        json_files.append(json_out)

    # Формируем команду nest.exe
    nest_binary = "nest.exe" if os.name == "nt" else "./nest"
    nest_cmd = [nest_binary, "-s", sheet]

    if iterations is not None:
        nest_cmd += ["-i", iterations]
    if nums:
        nest_cmd += ["-n"] + list(map(str, nums))

    nest_cmd += nest_flags  # <-- передаем все поддерживаемые флаги
    nest_cmd += ["-o", "lay.csv"] + json_files

    print("[PY] RUN:", " ".join(nest_cmd))
    subprocess.check_call(nest_cmd)
    print("[PY] Success!")
