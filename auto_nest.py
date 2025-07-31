"""Helper script for running the nesting binary.

The binary path is hard-coded as ``Release/nest.exe`` relative to this file.
"""

import argparse
import os
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path
from tqdm import tqdm  # прогресс-бар

sys.stdout.reconfigure(encoding='utf-8')


def _run(task):
    cmd = task[0]
    print("[PY] RUN:", " ".join(cmd))
    subprocess.check_call(cmd)
    return task[1], task[2], task[3], task[4]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Auto nesting pipeline")
    parser.add_argument("files", nargs="*", help="DXF files with optional :N count")
    parser.add_argument("-s", "--sheet", help="Sheet size WxH")
    parser.add_argument("-i", "--iterations", type=int)
    parser.add_argument("--strategy", default=None)
    parser.add_argument("--runs", type=int)
    parser.add_argument("-n", "--num", nargs="*")
    parser.add_argument("--fill-gaps", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-j", "--join-segments", action="store_true")
    parser.add_argument("--pop", type=int)
    parser.add_argument("--gen", type=int)
    parser.add_argument("--polish", type=int)
    parser.add_argument("--rot", type=int, help="rotation step (degrees)")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    sheet = args.sheet
    if not sheet or len(args.files) == 0:
        print("usage: auto_nest.py -s WxH part1.dxf[:N] ...")
        return 1

    iterations = args.iterations if args.iterations is not None else 1
    strategies = [args.strategy] if args.strategy else ["area"]
    runs = args.runs if args.runs is not None else 1
    nums = [int(n) for n in (args.num or [])]

    nest_flags = []
    if args.verbose:
        nest_flags.append("-v")
    if args.join_segments:
        nest_flags.append("-j")
    if args.fill_gaps:
        nest_flags.append("--fill-gaps")
    if args.pop is not None:
        nest_flags += ["--pop", str(args.pop)]
    if args.gen is not None:
        nest_flags += ["--gen", str(args.gen)]
    if args.rot is not None:
        nest_flags += ["--rot", str(args.rot)]
    if args.polish is not None:
        polish = args.polish
        nest_flags += ["--polish", str(polish)]

    dxf_files = []
    for item in args.files:
        if ":" in item:
            fname, n = item.rsplit(":", 1)
            dxf_files.append((fname, int(n)))
        else:
            dxf_files.append((item, 1))

    python_exec = sys.executable
    json_files = []
    for dxf, repeat in dxf_files:
        json_out = Path(dxf).with_suffix(".json")
        cmd = [python_exec, "extract_all_shapes.py", dxf, str(json_out)]
        if repeat > 1:
            cmd.append(f"repeat={repeat}")
        subprocess.check_call(cmd)

        if not json_out.exists() or json_out.stat().st_size == 0:
            raise RuntimeError(f"JSON output {json_out} missing or empty")
        with open(json_out, "r", encoding="utf-8") as chk:
            first = chk.read(1)
            if first != "[":
                raise RuntimeError(f"Invalid JSON produced: {json_out}")
        print(f"[PY] JSON generated: {json_out}")
        json_files.append(str(json_out))

    nest_binary = os.path.abspath("Release/nest.exe")
    tasks = []
    for strat in strategies:
        for r in range(runs):
            out_csv = f"lay_{r}_{strat}.csv"
            out_dxf = f"layout_{r}_{strat}.dxf"
            cmd = [nest_binary, "-s", sheet]
            if iterations:
                cmd += ["--iter", str(iterations)]
            if nums:
                cmd += ["-n"] + [str(n) for n in nums]
            cmd += nest_flags
            cmd += ["--strategy", strat, "--dxf", out_dxf, "-o", out_csv]
            cmd += json_files
            tasks.append((cmd, out_csv, out_dxf, strat, r))

    # Прогресс-бар по всем задачам
    total = len(tasks)
    with Pool() as pool:
        results = list(tqdm(pool.imap_unordered(_run, tasks), total=total, desc="Nesting progress", ncols=80))

    with open("lay.csv", "w") as fout:
        header_written = False
        for csvf, dxff, strat, r in results:
            with open(csvf) as fin:
                if not header_written:
                    fout.write(fin.readline())
                    header_written = True
                else:
                    fin.readline()
                fout.writelines(fin.readlines())

    import csv
    best = None
    with open("lay.csv") as f:
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
