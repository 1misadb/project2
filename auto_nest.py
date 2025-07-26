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
    allowed_flags = {"-v", "--verbose", "-j", "--join-segments",
                     "--pop", "--gen", "--polish", "--fill-gaps"}

    s_idx = sys.argv.index("-s")
    try:
        sheet = sys.argv[s_idx + 1]

        iterations = None
        nums = []
        dxf_files = []  # список кортежей (имя, количество)
        nest_flags = [] # сюда кидаем "-v", "-j" и любые другие поддерживаемые
        strategies = ["area"]
        runs = 1
        fill_gaps = False

        i = s_idx + 2
        while i < len(sys.argv):
            arg = sys.argv[i]
            if arg == "-i":
                iterations = sys.argv[i + 1]
                i += 2
            elif arg == "--strategies":
                strategies = sys.argv[i+1].split(',')
                i += 2
            elif arg == "--strategy":
                strategies = [sys.argv[i+1]]
                i += 2
            elif arg == "--runs":
                runs = int(sys.argv[i+1])
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
                if arg == "--fill-gaps":
                    fill_gaps = True
                    i += 1
                elif arg.startswith('--') and arg not in ("-v", "--verbose", "-j", "--join-segments"):
                    if i + 1 < len(sys.argv):
                        nest_flags.append(sys.argv[i+1])
                        i += 2
                    else:
                        usage()
                else:
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
    python_exec = sys.executable
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

    nest_binary = "nest.exe" if os.name == "nt" else "./nest"

    tasks = []
    for strat in strategies:
        for r in range(runs):
            out_csv = f"lay_{r}_{strat}.csv"
            out_dxf = f"layout_{r}_{strat}.dxf"
            cmd = [nest_binary, "-s", sheet]
            if iterations is not None:
                cmd += ["-i", iterations]
            if nums:
                cmd += ["-n"] + list(map(str, nums))
            cmd += nest_flags
            if fill_gaps:
                cmd += ["--fill-gaps"]
            cmd += ["--strategy", strat, "--run", str(r), "--dxf", out_dxf, "-o", out_csv]
            cmd += json_files
            tasks.append((cmd, out_csv, out_dxf, strat, r))

    from multiprocessing import Pool
    def _run(t):
        cmd=t[0]
        print("[PY] RUN:", " ".join(cmd))
        subprocess.check_call(cmd)
        return t[1], t[2], t[3], t[4]

    with Pool() as pool:
        results = pool.map(_run, tasks)

    # Aggregate lay.csv
    with open('lay.csv','w') as fout:
        header_written=False
        for csvf, dxff, strat, r in results:
            with open(csvf) as fin:
                if not header_written:
                    fout.write(fin.readline())
                    header_written=True
                else:
                    fin.readline()
                fout.writelines(fin.readlines())

    # choose best by minimal waste
    import csv
    best = None
    with open('lay.csv') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if best is None or float(row['waste_mm2']) < float(best['waste_mm2']):
                best = row
    best_dxf = f"layout_{best['run']}_{best['strategy']}.dxf"
    if os.path.exists(best_dxf):
        os.replace(best_dxf, 'layout.dxf')
    try:
        import pdf_report
        pdf_report.generate_report('lay.csv', 'nest_report.pdf', 'layout.dxf')
        print("[PY] PDF saved to nest_report.pdf")
    except Exception as e:
        print("[PY] PDF generation failed", e)
    print("[PY] Done")