"""
Snakemake workflow: Bowtie2 index + align trimmed reads + SAM→BAM→sort→index
for isolate Bacteroides muris
"""

import glob, os

READS_DIR = "/scratch/gpfs/AMOELLER/sidd/diet_manip/isolate_Bacteroides_muris/qc/fastp"
OUT_DIR   = "/scratch/gpfs/AMOELLER/sidd/diet_manip/isolate_Bacteroides_muris/alignment"
REF       = "/scratch/gpfs/AMOELLER/sidd/diet_manip/isolate_Bacteroides_muris/redo_representative_megamag.fa"
IDX_BASE  = os.path.join(OUT_DIR, "index", "megamag")

# Discover sample basenames from trimmed reads
SAMPLES = sorted(set(
    os.path.basename(f).replace("_R1.fastq.gz", "").replace("_R2.fastq.gz", "")
    for f in glob.glob(os.path.join(READS_DIR, "*.fastq.gz"))
))

rule all:
    input:
        expand(os.path.join(OUT_DIR, "sorted", "{sample}.sorted.bam.bai"), sample=SAMPLES),


# ---------------------------------------------------------------------------
#  Bowtie2 index
# ---------------------------------------------------------------------------
rule bowtie2_index:
    input:
        REF,
    output:
        multiext(IDX_BASE, ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    params:
        idx_base = IDX_BASE,
    threads: 8
    conda:
        "diet_manip"
    shell:
        """
        mkdir -p $(dirname {params.idx_base})
        bowtie2-build --threads {threads} {input} {params.idx_base}
        """


# ---------------------------------------------------------------------------
#  Bowtie2 align
# ---------------------------------------------------------------------------
rule bowtie2_align:
    input:
        r1  = os.path.join(READS_DIR, "{sample}_R1.fastq.gz"),
        r2  = os.path.join(READS_DIR, "{sample}_R2.fastq.gz"),
        idx = multiext(IDX_BASE, ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    output:
        sam = temp(os.path.join(OUT_DIR, "sam", "{sample}.sam")),
    params:
        idx_base = IDX_BASE,
    threads: 8
    conda:
        "diet_manip"
    shell:
        """
        mkdir -p $(dirname {output.sam})
        bowtie2 -x {params.idx_base} \
            -1 {input.r1} -2 {input.r2} \
            --threads {threads} \
            -S {output.sam}
        """


# ---------------------------------------------------------------------------
#  SAM → BAM
# ---------------------------------------------------------------------------
rule sam_to_bam:
    input:
        os.path.join(OUT_DIR, "sam", "{sample}.sam"),
    output:
        bam = os.path.join(OUT_DIR, "bam", "{sample}.bam"),
    threads: 8
    conda:
        "diet_manip"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        samtools view --threads {threads} -bS {input} -o {output.bam}
        """


# ---------------------------------------------------------------------------
#  Sort BAM
# ---------------------------------------------------------------------------
rule sort_bam:
    input:
        os.path.join(OUT_DIR, "bam", "{sample}.bam"),
    output:
        os.path.join(OUT_DIR, "sorted", "{sample}.sorted.bam"),
    threads: 8
    conda:
        "diet_manip"
    shell:
        """
        mkdir -p $(dirname {output})
        samtools sort --threads {threads} -o {output} {input}
        """


# ---------------------------------------------------------------------------
#  Index sorted BAM
# ---------------------------------------------------------------------------
rule index_bam:
    input:
        os.path.join(OUT_DIR, "sorted", "{sample}.sorted.bam"),
    output:
        os.path.join(OUT_DIR, "sorted", "{sample}.sorted.bam.bai"),
    threads: 4
    conda:
        "diet_manip"
    shell:
        """
        samtools index --threads {threads} {input}
        """
