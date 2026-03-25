import os

ORFM_C = os.path.expanduser("~/git/OrfM/orfm")
ORFM_RS = "target/release/orfm"

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

rule all:
    input:
        "benchmark/results.tsv",
        "benchmark/correctness.txt"

rule build_rust:
    output:
        ORFM_RS
    shell:
        "cargo build --release"

rule generate_fasta_unwrapped:
    output:
        "benchmark/random_seqs_fasta_unwrapped_{rep}.fasta"
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
        "benchmark/random_seqs_fasta_wrapped_{rep}.fasta"
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
        "benchmark/random_seqs_fasta_unwrapped_{rep}.fasta"
    output:
        "benchmark/random_seqs_fasta_gzipped_{rep}.fasta.gz"
    shell:
        "gzip -c {input} > {output}"

rule generate_fastq:
    output:
        "benchmark/random_seqs_fastq_{rep}.fastq"
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
        "benchmark/random_seqs_fastq_{rep}.fastq"
    output:
        "benchmark/random_seqs_fastq_gzipped_{rep}.fastq.gz"
    shell:
        "gzip -c {input} > {output}"


def get_input_file(wildcards):
    itype = wildcards.itype
    rep = wildcards.rep
    ext = INPUT_TYPES[itype]
    return f"benchmark/random_seqs_{itype}_{rep}{ext}"


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
        rs_times=expand("benchmark/time_rs_{itype}_{rep}.txt", itype=INPUT_TYPES, rep=REPLICATES)
    output:
        "benchmark/results.tsv"
    run:
        with open(output[0], 'w') as out:
            out.write("tool\tinput_type\treplicate\twall_clock_s\tmax_rss_kb\n")
            for itype in INPUT_TYPES:
                for rep in REPLICATES:
                    for tool, prefix in [("orfm_c", "c"), ("orfm_rs", "rs")]:
                        timefile = f"benchmark/time_{prefix}_{itype}_{rep}.txt"
                        vals = {}
                        with open(timefile) as f:
                            for line in f:
                                key, val = line.strip().split("\t")
                                vals[key] = val
                        wall = vals.get("wall_clock_s", "NA")
                        rss = vals.get("max_rss_kb", "NA")
                        out.write(f"{tool}\t{itype}\t{rep}\t{wall}\t{rss}\n")
        # Print summary
        with open(output[0]) as f:
            print(f.read())
