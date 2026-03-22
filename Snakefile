import os

ORFM_C = os.path.expanduser("~/git/OrfM/orfm")
ORFM_RS = "target/release/orfm"

# Number of random sequences and their lengths for benchmarking
NUM_SEQS = 1_000_000
SEQ_LEN = 150
REPLICATES = range(1, 4)

rule all:
    input:
        "benchmark/results.tsv",
        "benchmark/correctness.txt"

rule build_rust:
    output:
        ORFM_RS
    shell:
        "cargo build --release"

rule generate_random_sequences:
    output:
        "benchmark/random_seqs_{rep}.fasta"
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

rule benchmark_orfm_c:
    input:
        seqs="benchmark/random_seqs_{rep}.fasta",
        bin=ORFM_C
    output:
        time="benchmark/time_c_{rep}.txt",
        result="benchmark/output_c_{rep}.fasta"
    shell:
        """
        /usr/bin/time -v {input.bin} {input.seqs} > {output.result} 2> {output.time} || true
        """

rule benchmark_orfm_rs:
    input:
        seqs="benchmark/random_seqs_{rep}.fasta",
        bin=ORFM_RS
    output:
        time="benchmark/time_rs_{rep}.txt",
        result="benchmark/output_rs_{rep}.fasta"
    shell:
        """
        /usr/bin/time -v {input.bin} {input.seqs} > {output.result} 2> {output.time} || true
        """

rule check_correctness:
    input:
        c=expand("benchmark/output_c_{rep}.fasta", rep=REPLICATES),
        rs=expand("benchmark/output_rs_{rep}.fasta", rep=REPLICATES)
    output:
        "benchmark/correctness.txt"
    run:
        all_ok = True
        with open(output[0], 'w') as f:
            for rep in REPLICATES:
                c_file = f"benchmark/output_c_{rep}.fasta"
                rs_file = f"benchmark/output_rs_{rep}.fasta"
                import subprocess
                result = subprocess.run(["diff", c_file, rs_file], capture_output=True, text=True)
                if result.returncode == 0:
                    f.write(f"Replicate {rep}: PASS (outputs identical)\n")
                else:
                    f.write(f"Replicate {rep}: FAIL\n")
                    f.write(result.stdout[:500] + "\n")
                    all_ok = False
            if all_ok:
                f.write("\nAll replicates produce identical output.\n")

rule collect_results:
    input:
        c_times=expand("benchmark/time_c_{rep}.txt", rep=REPLICATES),
        rs_times=expand("benchmark/time_rs_{rep}.txt", rep=REPLICATES)
    output:
        "benchmark/results.tsv"
    run:
        import re
        with open(output[0], 'w') as out:
            out.write("tool\treplicate\twall_clock_s\tmax_rss_kb\n")
            for rep in REPLICATES:
                for tool, prefix in [("orfm_c", "c"), ("orfm_rs", "rs")]:
                    timefile = f"benchmark/time_{prefix}_{rep}.txt"
                    wall = rss = "NA"
                    with open(timefile) as f:
                        for line in f:
                            m = re.search(r"Elapsed \(wall clock\) time.*: (.+)", line)
                            if m:
                                # Parse h:mm:ss or m:ss.ss format
                                parts = m.group(1).strip().split(":")
                                if len(parts) == 3:
                                    wall = str(int(parts[0])*3600 + int(parts[1])*60 + float(parts[2]))
                                elif len(parts) == 2:
                                    wall = str(int(parts[0])*60 + float(parts[1]))
                                else:
                                    wall = parts[0]
                            m = re.search(r"Maximum resident set size.*: (\d+)", line)
                            if m:
                                rss = m.group(1)
                    out.write(f"{tool}\t{rep}\t{wall}\t{rss}\n")
        # Print summary
        with open(output[0]) as f:
            print(f.read())
