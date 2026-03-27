import os

ORFM_C = os.path.expanduser("~/git/OrfM/orfm")
ORFM_RS = "/tmp/orfm-target/release/orfm"

# Number of random sequences and their lengths for benchmarking
NUM_SEQS = 1_000_000
SEQ_LEN = 150
LONG_NUM_SEQS = 100
LONG_SEQ_LEN = 1_000_000
REPLICATES = range(1, 4)

# Input types and their file extensions
INPUT_TYPES = {
    "fasta_unwrapped": ".fasta",
    "fasta_wrapped": ".fasta",
    "fasta_gzipped": ".fasta.gz",
    "fastq": ".fastq",
    "fastq_gzipped": ".fastq.gz",
}
# Long-sequence benchmark types (stored in LONG_DATA_DIR = /tmp)
LONG_INPUT_TYPES = {
    "long_fasta": ".fasta",
    "long_fasta_iupac": ".fasta",
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
LONG_DATA_DIR = "/tmp/orfm_long_data"

ruleorder: generate_long_fasta_iupac > generate_long_fasta

rule all:
    input:
        "benchmark/results.tsv",
        "benchmark/results_long.tsv",
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

rule generate_long_fasta:
    output:
        f"{LONG_DATA_DIR}/random_seqs_long_fasta_{{rep}}.fasta"
    params:
        num_seqs=LONG_NUM_SEQS,
        seq_len=LONG_SEQ_LEN
    shell:
        """
        mkdir -p {LONG_DATA_DIR}
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

rule generate_long_fasta_iupac:
    wildcard_constraints:
        rep=r"\d+"
    output:
        f"{LONG_DATA_DIR}/random_seqs_long_fasta_iupac_{{rep}}.fasta"
    params:
        num_seqs=LONG_NUM_SEQS,
        seq_len=LONG_SEQ_LEN
    shell:
        """
        mkdir -p {LONG_DATA_DIR}
        python3 -c "
import random
random.seed({wildcards.rep})
bases = 'ACGT'
iupac_extra = 'RYSWKMBDHVN'
with open('{output}', 'w') as f:
    for i in range({params.num_seqs}):
        seq = list(''.join(random.choices(bases, k={params.seq_len})))
        for pos in random.sample(range({params.seq_len}), int({params.seq_len} * 0.04)):
            seq[pos] = random.choice(iupac_extra)
        f.write(f'>seq_{{i}}\\n' + ''.join(seq) + '\\n')
        "
        """


def get_long_input_file(wildcards):
    itype = wildcards.itype
    rep = wildcards.rep
    ext = LONG_INPUT_TYPES[itype]
    return f"{LONG_DATA_DIR}/random_seqs_{itype}_{rep}{ext}"


rule benchmark_long_orfm_c:
    input:
        seqs=get_long_input_file,
        bin=ORFM_C
    output:
        time="benchmark/time_long_c_{itype}_{rep}.txt",
        result="/tmp/orfm_long_output/output_long_c_{itype}_{rep}.fasta"
    wildcard_constraints:
        itype="|".join(LONG_INPUT_TYPES.keys()),
        rep=r"\d+"
    run:
        import subprocess, time, resource
        os.makedirs("/tmp/orfm_long_output", exist_ok=True)
        start = time.monotonic()
        with open(output.result, 'w') as out:
            subprocess.run([input.bin, input.seqs], stdout=out, check=True)
        elapsed = time.monotonic() - start
        rusage = resource.getrusage(resource.RUSAGE_CHILDREN)
        with open(output.time, 'w') as f:
            f.write(f"wall_clock_s\t{elapsed:.3f}\nmax_rss_kb\t{rusage.ru_maxrss}\n")


rule benchmark_long_orfm_rs:
    input:
        seqs=get_long_input_file,
        bin=ORFM_RS
    output:
        time="benchmark/time_long_rs_{itype}_{rep}.txt",
        result="/tmp/orfm_long_output/output_long_rs_{itype}_{rep}.fasta"
    wildcard_constraints:
        itype="|".join(LONG_INPUT_TYPES.keys()),
        rep=r"\d+"
    run:
        import subprocess, time, resource
        os.makedirs("/tmp/orfm_long_output", exist_ok=True)
        start = time.monotonic()
        with open(output.result, 'w') as out:
            subprocess.run([input.bin, input.seqs], stdout=out, check=True)
        elapsed = time.monotonic() - start
        rusage = resource.getrusage(resource.RUSAGE_CHILDREN)
        with open(output.time, 'w') as f:
            f.write(f"wall_clock_s\t{elapsed:.3f}\nmax_rss_kb\t{rusage.ru_maxrss}\n")


rule collect_long_results:
    input:
        c_times=expand("benchmark/time_long_c_{itype}_{rep}.txt", itype=LONG_INPUT_TYPES, rep=REPLICATES),
        rs_times=expand("benchmark/time_long_rs_{itype}_{rep}.txt", itype=LONG_INPUT_TYPES, rep=REPLICATES),
    output:
        "benchmark/results_long.tsv"
    run:
        import statistics

        def read_time(prefix, itype, rep):
            vals = {}
            with open(f"benchmark/time_long_{prefix}_{itype}_{rep}.txt") as f:
                for line in f:
                    k, v = line.strip().split("\t")
                    vals[k] = v
            return vals.get("wall_clock_s", "N/A"), vals.get("max_rss_kb", "N/A")

        with open(output[0], 'w') as out:
            out.write("tool\tinput_type\treplicate\twall_clock_s\tmax_rss_kb\n")
            for itype in LONG_INPUT_TYPES:
                for rep in REPLICATES:
                    for tool, prefix in [("orfm_c", "c"), ("orfm_rs", "rs")]:
                        wall, rss = read_time(prefix, itype, rep)
                        out.write(f"{tool}\t{itype}\t{rep}\t{wall}\t{rss}\n")

        def median_or_na(vals):
            nums = [float(v) for v in vals if v != "N/A"]
            return statistics.median(nums) if nums else None

        data = {}
        for itype in LONG_INPUT_TYPES:
            for rep in REPLICATES:
                for tool, prefix in [("C", "c"), ("Rs", "rs")]:
                    wall, rss = read_time(prefix, itype, rep)
                    key = (tool, itype)
                    if key not in data:
                        data[key] = {"wall": []}
                    data[key]["wall"].append(wall)

        col_w = max(max(len(t) for t in LONG_INPUT_TYPES), 5)
        ts_w, rat_w = 8, 6

        def fmt_s(v):      return f"{v:.2f}" if v is not None else "N/A"
        def fmt_rat(c, r): return f"{r/c:.2f}x" if (c and r) else "N/A"

        cols = [("Input", col_w), ("C (s)", ts_w), ("Rs (s)", ts_w), ("Rs/C", rat_w)]

        def hline(l, m, r):
            return l + m.join("─" * (w + 2) for _, w in cols) + r

        lines = [hline("┌", "┬", "┐")]
        lines.append("│" + "│".join(f" {h:^{w}} " for h, w in cols) + "│")
        lines.append(hline("├", "┼", "┤"))
        for i, itype in enumerate(LONG_INPUT_TYPES):
            if i > 0:
                lines.append(hline("├", "┼", "┤"))
            c_w  = median_or_na(data[("C",  itype)]["wall"])
            r_w  = median_or_na(data[("Rs", itype)]["wall"])
            cells = [
                f" {itype:<{col_w}} ",
                f" {fmt_s(c_w):>{ts_w}} ",
                f" {fmt_s(r_w):>{ts_w}} ",
                f" {fmt_rat(c_w, r_w):>{rat_w}} ",
            ]
            lines.append("│" + "│".join(cells) + "│")
        lines.append(hline("└", "┴", "┘"))
        print("\n".join(lines))


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
        import statistics

        def read_time(prefix, itype, rep):
            vals = {}
            with open(f"benchmark/time_{prefix}_{itype}_{rep}.txt") as f:
                for line in f:
                    k, v = line.strip().split("\t")
                    vals[k] = v
            return vals.get("wall_clock_s", "N/A"), vals.get("max_rss_kb", "N/A")

        with open(output[0], 'w') as out:
            out.write("tool\tinput_type\treplicate\twall_clock_s\tmax_rss_kb\n")
            for itype in INPUT_TYPES:
                for rep in REPLICATES:
                    for tool, prefix in [("orfm_c", "c"), ("orfm_rs", "rs"), ("getorf", "getorf")]:
                        wall, rss = read_time(prefix, itype, rep)
                        out.write(f"{tool}\t{itype}\t{rep}\t{wall}\t{rss}\n")

        def median_or_na(vals):
            nums = [float(v) for v in vals if v != "N/A"]
            return statistics.median(nums) if nums else None

        data = {}
        for itype in INPUT_TYPES:
            for rep in REPLICATES:
                for tool, prefix in [("C", "c"), ("Rs", "rs"), ("getorf", "getorf")]:
                    wall, rss = read_time(prefix, itype, rep)
                    key = (tool, itype)
                    if key not in data:
                        data[key] = {"wall": [], "rss": []}
                    data[key]["wall"].append(wall)
                    data[key]["rss"].append(rss)

        def fmt_s(v):   return f"{v:.2f}" if v is not None else "N/A"
        def fmt_mb(v):  return f"{v/1024:.0f}" if v is not None else "N/A"
        def fmt_rat(c, r): return f"{r/c:.2f}x" if (c and r) else "N/A"

        def print_table(type_list, include_getorf):
            col_w = max(max(len(t) for t in type_list), 5)
            ts_w, rat_w = 6, 6
            cols = [("Input", col_w), ("C (s)", ts_w), ("Rs (s)", ts_w)]
            if include_getorf:
                cols.append(("go (s)", ts_w))
            cols.append(("Rs/C", rat_w))

            def hline(l, m, r):
                return l + m.join("─" * (w + 2) for _, w in cols) + r

            lines = [hline("┌", "┬", "┐")]
            lines.append("│" + "│".join(f" {h:^{w}} " for h, w in cols) + "│")
            lines.append(hline("├", "┼", "┤"))
            for i, itype in enumerate(type_list):
                if i > 0:
                    lines.append(hline("├", "┼", "┤"))
                c_w  = median_or_na(data[("C",  itype)]["wall"])
                r_w  = median_or_na(data[("Rs", itype)]["wall"])
                g_w  = median_or_na(data[("getorf", itype)]["wall"]) if include_getorf else None
                cells = [
                    f" {itype:<{col_w}} ",
                    f" {fmt_s(c_w):>{ts_w}} ",
                    f" {fmt_s(r_w):>{ts_w}} ",
                ]
                if include_getorf:
                    cells.append(f" {fmt_s(g_w):>{ts_w}} ")
                cells.append(f" {fmt_rat(c_w, r_w):>{rat_w}} ")
                lines.append("│" + "│".join(cells) + "│")
            lines.append(hline("└", "┴", "┘"))
            print("\n".join(lines))

        print_table(list(INPUT_TYPES.keys()), include_getorf=True)
