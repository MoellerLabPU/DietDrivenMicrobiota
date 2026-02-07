"""
Snakemake workflow: FastQC (raw) -> fastp -> FastQC (trimmed) -> MultiQC
for isolate Bacteroides muris reads
"""

READS_DIR = "/scratch/gpfs/AMOELLER/sidd/diet_manip/isolate_Bacteroides_muris/reads"
OUT_DIR   = "/scratch/gpfs/AMOELLER/sidd/diet_manip/isolate_Bacteroides_muris/qc"

import glob, os

# Discover all sample basenames (everything before _R1/_R2)
SAMPLES = sorted(set(
    os.path.basename(f).replace("_R1.fastq.gz", "").replace("_R2.fastq.gz", "")
    for f in glob.glob(os.path.join(READS_DIR, "*.fastq.gz"))
))

rule all:
    input:
        os.path.join(OUT_DIR, "multiqc_raw", "multiqc_report.html"),
        os.path.join(OUT_DIR, "multiqc_trimmed", "multiqc_report.html"),


# ---------------------------------------------------------------------------
#  FastQC on raw reads
# ---------------------------------------------------------------------------
rule fastqc_raw:
    input:
        r1 = os.path.join(READS_DIR, "{sample}_R1.fastq.gz"),
        r2 = os.path.join(READS_DIR, "{sample}_R2.fastq.gz"),
    output:
        html1 = os.path.join(OUT_DIR, "fastqc_raw", "{sample}_R1_fastqc.html"),
        zip1  = os.path.join(OUT_DIR, "fastqc_raw", "{sample}_R1_fastqc.zip"),
        html2 = os.path.join(OUT_DIR, "fastqc_raw", "{sample}_R2_fastqc.html"),
        zip2  = os.path.join(OUT_DIR, "fastqc_raw", "{sample}_R2_fastqc.zip"),
    params:
        outdir = os.path.join(OUT_DIR, "fastqc_raw"),
    threads: 8
    conda:
        "diet_manip"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2}
        """


# ---------------------------------------------------------------------------
#  fastp — adapter trimming & quality filtering
# ---------------------------------------------------------------------------
rule fastp:
    input:
        r1 = os.path.join(READS_DIR, "{sample}_R1.fastq.gz"),
        r2 = os.path.join(READS_DIR, "{sample}_R2.fastq.gz"),
    output:
        r1 = os.path.join(OUT_DIR, "fastp", "{sample}_R1.fastq.gz"),
        r2 = os.path.join(OUT_DIR, "fastp", "{sample}_R2.fastq.gz")
    threads: 8
    conda:
        "diet_manip"
    shell:
        """
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --thread {threads} \
            --detect_adapter_for_pe \
            --allow_gap_overlap_trimming
        """


# ---------------------------------------------------------------------------
#  FastQC on trimmed reads
# ---------------------------------------------------------------------------
rule fastqc_trimmed:
    input:
        r1 = os.path.join(OUT_DIR, "fastp", "{sample}_R1.fastq.gz"),
        r2 = os.path.join(OUT_DIR, "fastp", "{sample}_R2.fastq.gz"),
    output:
        html1 = os.path.join(OUT_DIR, "fastqc_trimmed", "{sample}_R1_fastqc.html"),
        zip1  = os.path.join(OUT_DIR, "fastqc_trimmed", "{sample}_R1_fastqc.zip"),
        html2 = os.path.join(OUT_DIR, "fastqc_trimmed", "{sample}_R2_fastqc.html"),
        zip2  = os.path.join(OUT_DIR, "fastqc_trimmed", "{sample}_R2_fastqc.zip"),
    params:
        outdir = os.path.join(OUT_DIR, "fastqc_trimmed"),
    threads: 8
    conda:
        "diet_manip"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2}
        """


# ---------------------------------------------------------------------------
#  MultiQC on raw reads (before trimming)
# ---------------------------------------------------------------------------
rule multiqc_raw:
    input:
        expand(os.path.join(OUT_DIR, "fastqc_raw", "{sample}_R{read}_fastqc.zip"),
               sample=SAMPLES, read=["1", "2"]),
    output:
        os.path.join(OUT_DIR, "multiqc_raw", "multiqc_report.html"),
    params:
        indir  = os.path.join(OUT_DIR, "fastqc_raw"),
        outdir = os.path.join(OUT_DIR, "multiqc_raw"),
    conda:
        "diet_manip"
    shell:
        """
        multiqc {params.indir} -o {params.outdir} --force
        """


# ---------------------------------------------------------------------------
#  MultiQC on trimmed reads (after fastp + FastQC)
# ---------------------------------------------------------------------------
rule multiqc_trimmed:
    input:
        fastqc = expand(os.path.join(OUT_DIR, "fastqc_trimmed", "{sample}_R{read}_fastqc.zip"),
                        sample=SAMPLES, read=["1", "2"])
    output:
        os.path.join(OUT_DIR, "multiqc_trimmed", "multiqc_report.html"),
    params:
        indir  = os.path.join(OUT_DIR, "fastqc_trimmed"),
        outdir = os.path.join(OUT_DIR, "multiqc_trimmed"),
    conda:
        "diet_manip"
    shell:
        """
        multiqc {params.indir} -o {params.outdir} --force
        """
