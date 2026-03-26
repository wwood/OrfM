import os

ORFM_C = os.path.expanduser("~/git/OrfM/orfm")
ORFM_RS = "/tmp/orfm-target/release/orfm"

# Number of random sequences and their lengths for benchmarking
NUM_SEQS = 1_000_000
SEQ_LEN = 150
REPLICATES = range(1, 4)

# Input types and their file extensions
INPUT_TYPES = {
    "fasta_unwrapped": ".fasta",
    "fasta_wrapped": ".fasta",
    "fasta_gzipped": ".fasta.gz",
    "fastq": ".fastq",
    "fastq_gzipped": ".fastq.gz",
}

# getorf (EMBOSS) sequence format prefix for each input type; None = not supported
GETORF_FORMAT = {
    "fasta_unwrapped": "fasta",
    "fasta_wrapped": "fasta",
    "fasta_gzipped": None,
    "fastq": "fastq",
    "fastq_gzipped": None,
}

DATA_DIR = "benchmark/data"

rule all:
    input:
        "benchmark/results.tsv",
        "benchmark/correctness.txt"

rule build_rust:
    output:
        ORFM_RS
    shell:
        "CARGO_HOME=/container_home/.cargo cargo build --release"

rule generate_fasta_unwrapped:
    output:
        f"{DATA_DIR}/random_seqs_fasta_unwrapped_{{rep}}.fasta"
    params:
        num_seqs=NUM_SEQS,
        seq_len=SEQ_LEN
    shell:
        """
        python3 -c "
import random
random.seed({wildcards.rep})
bases = 'ACGT'
with open('{output}', 'w') as f:
    for i in range({params.num_seqs}):
        seq = ''.join(random.choices(bases, k={params.seq_len}))
        f.write(f'>seq_{{i}}\\n{{seq}}\\n')
        "
        """

rule generate_fasta_wrapped:
    output:
        f"{DATA_DIR}/random_seqs_fasta_wrapped_{{rep}}.fasta"
    params:
        num_seqs=NUM_SEQS,
        seq_len=SEQ_LEN,
        wrap_width=60
    shell:
        """
        python3 -c "
import random, textwrap
random.seed({wildcards.rep})
bases = 'ACGT'
with open('{output}', 'w') as f:
    for i in range({params.num_seqs}):
        seq = ''.join(random.choices(bases, k={params.seq_len}))
        f.write(f'>seq_{{i}}\\n')
        for line in textwrap.wrap(seq, {params.wrap_width}):
            f.write(line + '\\n')
        "
        """

rule generate_fasta_gzipped:
    input:
        f"{DATA_DIR}/random_seqs_fasta_unwrapped_{{rep}}.fasta"
    output:
        f"{DATA_DIR}/random_seqs_fasta_gzipped_{{rep}}.fasta.gz"
    shell:
        "gzip -c {input} > {output}"

rule generate_fastq:
    output:
        f"{DATA_DIR}/random_seqs_fastq_{{rep}}.fastq"
    params:
        num_seqs=NUM_SEQS,
        seq_len=SEQ_LEN
    shell:
        """
        python3 -c "
import random
random.seed({wildcards.rep})
bases = 'ACGT'
with open('{output}', 'w') as f:
    for i in range({params.num_seqs}):
        seq = ''.join(random.choices(bases, k={params.seq_len}))
        qual = 'I' * {params.seq_len}
        f.write(f'@seq_{{i}}\\n{{seq}}\\n+\\n{{qual}}\\n')
        "
        """

rule generate_fastq_gzipped:
    input:
        f"{DATA_DIR}/random_seqs_fastq_{{rep}}.fastq"
    output:
        f"{DATA_DIR}/random_seqs_fastq_gzipped_{{rep}}.fastq.gz"
    shell:
        "gzip -c {input} > {output}"


def get_input_file(wildcards):
    itype = wildcards.itype
    rep = wildcards.rep
    ext = INPUT_TYPES[itype]
    return f"{DATA_DIR}/random_seqs_{itype}_{rep}{ext}"


rule benchmark_orfm_c:
    input:
        seqs=get_input_file,
        bin=ORFM_C
    output:
        time="benchmark/time_c_{itype}_{rep}.txt",
        result="benchmark/output_c_{itype}_{rep}.fasta"
    run:
        import subprocess, time, resource
        start = time.monotonic()
        with open(output.result, 'w') as out:
            subprocess.run([input.bin, input.seqs], stdout=out, check=True)
        elapsed = time.monotonic() - start
        rusage = resource.getrusage(resource.RUSAGE_CHILDREN)
        with open(output.time, 'w') as f:
            f.write(f"wall_clock_s\t{elapsed:.3f}\nmax_rss_kb\t{rusage.ru_maxrss}\n")

rule benchmark_orfm_rs:
    input:
        seqs=get_input_file,
        bin=ORFM_RS
    output:
        time="benchmark/time_rs_{itype}_{rep}.txt",
        result="benchmark/output_rs_{itype}_{rep}.fasta"
    run:
        import subprocess, time, resource
        start = time.monotonic()
        with open(output.result, 'w') as out:
            subprocess.run([input.bin, input.seqs], stdout=out, check=True)
        elapsed = time.monotonic() - start
        rusage = resource.getrusage(resource.RUSAGE_CHILDREN)
        with open(output.time, 'w') as f:
            f.write(f"wall_clock_s\t{elapsed:.3f}\nmax_rss_kb\t{rusage.ru_maxrss}\n")

rule benchmark_getorf:
    input:
        seqs=get_input_file,
    output:
        time="benchmark/time_getorf_{itype}_{rep}.txt",
        result="benchmark/output_getorf_{itype}_{rep}.fasta"
    run:
        import subprocess, time, resource
        fmt = GETORF_FORMAT[wildcards.itype]
        if fmt is None:
            # getorf does not support this format (e.g. gzipped)
            with open(output.time, 'w') as f:
                f.write("wall_clock_s\tN/A\nmax_rss_kb\tN/A\n")
            open(output.result, 'w').close()
        else:
            start = time.monotonic()
            subprocess.run(
                ["getorf",
                 "-sequence", f"{fmt}::{input.seqs}",
                 "-outseq", f"fasta::{output.result}",
                 "-minsize", "96",
                 "-find", "0"],
                check=True, capture_output=True
            )
            elapsed = time.monotonic() - start
            rusage = resource.getrusage(resource.RUSAGE_CHILDREN)
            with open(output.time, 'w') as f:
                f.write(f"wall_clock_s\t{elapsed:.3f}\nmax_rss_kb\t{rusage.ru_maxrss}\n")

rule check_correctness:
    input:
        c=expand("benchmark/output_c_{itype}_{rep}.fasta", itype=INPUT_TYPES, rep=REPLICATES),
        rs=expand("benchmark/output_rs_{itype}_{rep}.fasta", itype=INPUT_TYPES, rep=REPLICATES)
    output:
        "benchmark/correctness.txt"
    run:
        import subprocess
        all_ok = True
        with open(output[0], 'w') as f:
            for itype in INPUT_TYPES:
                for rep in REPLICATES:
                    c_file = f"benchmark/output_c_{itype}_{rep}.fasta"
                    rs_file = f"benchmark/output_rs_{itype}_{rep}.fasta"
                    result = subprocess.run(["diff", c_file, rs_file], capture_output=True, text=True)
                    if result.returncode == 0:
                        f.write(f"{itype} replicate {rep}: PASS (outputs identical)\n")
                    else:
                        f.write(f"{itype} replicate {rep}: FAIL\n")
                        f.write(result.stdout[:500] + "\n")
                        all_ok = False
            if all_ok:
                f.write("\nAll replicates produce identical output.\n")

rule collect_results:
    input:
        c_times=expand("benchmark/time_c_{itype}_{rep}.txt", itype=INPUT_TYPES, rep=REPLICATES),
        rs_times=expand("benchmark/time_rs_{itype}_{rep}.txt", itype=INPUT_TYPES, rep=REPLICATES),
        go_times=expand("benchmark/time_getorf_{itype}_{rep}.txt", itype=INPUT_TYPES, rep=REPLICATES),
    output:
        "benchmark/results.tsv"
    run:
        with open(output[0], 'w') as out:
            out.write("tool\tinput_type\treplicate\twall_clock_s\tmax_rss_kb\n")
            for itype in INPUT_TYPES:
                for rep in REPLICATES:
                    for tool, prefix in [("orfm_c", "c"), ("orfm_rs", "rs"), ("getorf", "getorf")]:
                        timefile = f"benchmark/time_{prefix}_{itype}_{rep}.txt"
                        vals = {}
                        with open(timefile) as f:
                            for line in f:
                                key, val = line.strip().split("\t")
                                vals[key] = val
                        wall = vals.get("wall_clock_s", "NA")
                        rss = vals.get("max_rss_kb", "NA")
                        out.write(f"{tool}\t{itype}\t{rep}\t{wall}\t{rss}\n")

        # Print summary table (median of replicates)
        import statistics

        def median_or_na(vals):
            nums = [float(v) for v in vals if v != "N/A"]
            return statistics.median(nums) if nums else None

        data = {}  # (tool, itype) -> {wall: [...], rss: [...]}
        for itype in INPUT_TYPES:
            for rep in REPLICATES:
                for tool, prefix in [("C", "c"), ("Rust", "rs"), ("getorf", "getorf")]:
                    timefile = f"benchmark/time_{prefix}_{itype}_{rep}.txt"
                    vals = {}
                    with open(timefile) as f:
                        for line in f:
                            key, val = line.strip().split("\t")
                            vals[key] = val
                    key = (tool, itype)
                    if key not in data:
                        data[key] = {"wall": [], "rss": []}
                    data[key]["wall"].append(vals.get("wall_clock_s", "N/A"))
                    data[key]["rss"].append(vals.get("max_rss_kb", "N/A"))

        # Box-drawing table
        # Columns: Input | C s | C MB | Rust s | Rust MB | getorf s | getorf MB | Rust/C
        max_name = max(len(t) for t in INPUT_TYPES)
        col_w  = max(max_name, 5)   # input name
        ts_w   = 6   # time (s)
        mb_w   = 6   # RAM (MB)
        rat_w  = 6   # ratio

        def fmt_s(v):
            return f"{v:.2f}" if v is not None else "N/A"

        def fmt_mb(v):
            return f"{v/1024:.0f}" if v is not None else "N/A"

        def fmt_rat(c, r):
            return f"{r/c:.2f}x" if (c is not None and r is not None) else "N/A"

        # Build column specs: (header, width)
        cols = [
            ("Input",    col_w),
            ("C (s)",    ts_w),
            ("C (MB)",   mb_w),
            ("Rs (s)",   ts_w),
            ("Rs (MB)",  mb_w),
            ("go (s)",   ts_w),
            ("go (MB)",  mb_w),
            ("Rs/C",     rat_w),
        ]

        def hline(left, mid, right, fill="─"):
            parts = [fill * (w + 2) for _, w in cols]
            return left + mid.join(parts) + right

        lines = []
        lines.append(hline("┌", "┬", "┐"))
        header = "│" + "│".join(f" {h:^{w}} " for h, w in cols) + "│"
        lines.append(header)
        lines.append(hline("├", "┼", "┤"))

        first = True
        for itype in INPUT_TYPES:
            if not first:
                lines.append(hline("├", "┼", "┤"))
            first = False

            c_w_val  = median_or_na(data[("C",      itype)]["wall"])
            c_mb_val = median_or_na(data[("C",      itype)]["rss"])
            r_w_val  = median_or_na(data[("Rust",   itype)]["wall"])
            r_mb_val = median_or_na(data[("Rust",   itype)]["rss"])
            g_w_val  = median_or_na(data[("getorf", itype)]["wall"])
            g_mb_val = median_or_na(data[("getorf", itype)]["rss"])

            cells = [
                f" {itype:<{col_w}} ",
                f" {fmt_s(c_w_val):>{ts_w}} ",
                f" {fmt_mb(c_mb_val):>{mb_w}} ",
                f" {fmt_s(r_w_val):>{ts_w}} ",
                f" {fmt_mb(r_mb_val):>{mb_w}} ",
                f" {fmt_s(g_w_val):>{ts_w}} ",
                f" {fmt_mb(g_mb_val):>{mb_w}} ",
                f" {fmt_rat(c_w_val, r_w_val):>{rat_w}} ",
            ]
            lines.append("│" + "│".join(cells) + "│")

        lines.append(hline("└", "┴", "┘"))
        print("\n".join(lines))
