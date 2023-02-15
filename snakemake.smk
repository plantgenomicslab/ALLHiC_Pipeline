#!snakemake
# conda activate allhic
# snakemake -s allhic.snakefile --cores 10 -p --dry-run --latency-wait 60
# snakemake -s allhic.snakefile --cores 10 -p --latency-wait 60 --notemp 

import os

os.makedirs("samples/trimmed/", exist_ok=True)
os.makedirs("samples/bam/", exist_ok=True)
os.makedirs("samples/sam/", exist_ok=True)
os.makedirs("samples/output/", exist_ok=True)

REFERENCE = "reference/genome.fa"
SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fastq.gz")
CTG_TABLE = "reference/Allele.ctg.table"
PARTITION = "17"

MAX_PARTITION = []
for i in range(1, int(PARTITION)+1):
    MAX_PARTITION.append(i)

SPLITS = ["{0:03}".format(i) for i in range(1, 11)]

ENZYME_SITE = "MBOI"
ENZYME_SEQ = ""
if ENZYME_SITE == "Arima" or ENZYME_SITE == "GATCGATC_GANTGATC_GANTANTC_GATCANTC":
    ENZYME_SEQ = "GATCGATC_GANTGATC_GANTANTC_GATCANTC"
elif ENZYME_SITE == "MBOI" or ENZYME_SITE == "GATC":
    ENZYME_SEQ = "GATC"
elif ENZYME_SITE == "HINDIII" or ENZYME_SITE == "AAGCTT":
    ENZYME_SEQ = "AAGCTT"
else:
    raise ValueError("ValueError: enzyme_sites parameter is an invalid value.")

rule all:
    input:
        index = "reference/genome.fa.index",
        fwd = expand("samples/trimmed/{split}.{sample}_R1.fastq.gz", sample=SAMPLES, split=SPLITS),
        rev = expand("samples/trimmed/{split}.{sample}_R2.fastq.gz", sample=SAMPLES, split=SPLITS),
        sam = expand("samples/sam/{split}.{sample}.sam", split=SPLITS, sample=SAMPLES),
        tours = expand("samples/bam/prune.counts_{enzyme_seq}.{parition}g{max_partition}.tour", enzyme_seq = ENZYME_SEQ, parition = PARTITION, max_partition = MAX_PARTITION),
        clm = "samples/bam/prune.clm"

#rule index:
#       message: "~-~ Trimming fastq files... ~-~"
#    threads:
#    params:
#    input:
#    output:
#    shell:

rule refIndex:
    message: "~-~ Trimming fastq files... ~-~"
    params: REFERENCE=REFERENCE
    output:
        index = "reference/genome.fa.index",
        fai = "reference/genome.fa.fai"
    shell: """
        chromap -i -r ${input.ref} -o ${output.index}
        samtools faidx ${input.ref}
    """
rule readTrim:
    message: "~-~ Trimming fastq files... ~-~"
    threads: 10
    input:
        fwd = "samples/raw/{sample}_R1.fastq.gz",
        rev = "samples/raw/{sample}_R2.fastq.gz"
    params:
        fwd = "samples/trimmed/{sample}_R1.fastq.gz",
        rev = "samples/trimmed/{sample}_R2.fastq.gz"
    output:
        fwd = "samples/trimmed/{split}.{sample}_R1.fastq.gz",
        rev = "samples/trimmed/{split}.{sample}_R2.fastq.gz"
    shell: """
        fastp -i {input.fwd} -I {input.rev} \
        -o {params.fwd} -O {params.rev}  \
        -h samples/trimmed/{wildcards.sample}.html -j samples/trimmed/{wildcards.sample}.json \
        -R {wildcards.sample} -w {threads} \
        --split {threads} \
        --split_prefix_digits 3 \
        --average_qual 20 \
        --length_limit 400 \
        --cut_mean_quality 10
    """

rule readMapping:
    message: "~-~ Trimming fastq files... ~-~"
    threads: 10
    params: REFERENCE=REFERENCE
    input:
        index = "reference/genome.fa.index",
        fwd = "samples/trimmed/{split}.{sample}_R1.fastq.gz",
        rev = "samples/trimmed/{split}.{sample}_R2.fastq.gz"
    output:
        sam = "samples/sam/{split}.{sample}.sam"
    shell: """
        chromap -r ${input.index} -r -1 ${input.fwd} -2 ${input.rev} --SAM -o ${output.sam} -t ${threads}
    """

rule PreprocessSAMs:
    message: "~-~ Trimming fastq files... ~-~"
    threads: 10
    params:
        REFERENCE=REFERENCE,
        ENZYME_SITE=ENZYME_SITE,
        fai = "reference/genome.fa.fai",
        sampaired = "samples/sam/{split}.{sample}.REduced.paired_only.bam",
        clean_sam = "samples/sam/{split}.{sample}.clean.sam"
    input:
        sam = "samples/sam/{split}.{sample}.sam"
    output:
        clean_bam = "samples/bam/{split}.{sample}.clean.bam"

    shell: """
        PreprocessSAMs.pl ${input.sam} ${params.REFERENCE} ${params.ENZYME_SITE}
        filterBAM_forHiC.pl ${params.sampaired} ${params.clean_sam}
        samtools view -bt ${params.fai} ${params.clean_sam} > ${output.clean_bam}
    """
rule mergeBAMs:
    message: "~-~ Trimming fastq files... ~-~"
    threads: 10
    input:
        clean_bam = expand("samples/bam/{split}.{sample}.clean.bam", split=SPLITS, sample=SAMPLES)
    output:
        merge_bam = "samples/bam/sample.clean.bam"

    shell: """
        samtools merge ${output.merge_bam} ${input.clean_bam}
    """

if (os.path.exists(CTG_TABLE)):
    rule Prune:
        input:
            CTG_TABLE=CTG_TABLE,
            REFERENCE=REFERENCE,
            merge_bam = "samples/bam/sample.clean.bam"
        output:
            prune_bam = "samples/bam/prune.bam"
        shell: """
        ALLHiC_prune -i ${CTG_TABLE} -b ${input.merge_bam} -r ${REFERENCE}
        """
# If indexdir do not exists => build the index from fasta
else:
    rule PruneLink:
        input:
            merge_bam = "samples/bam/sample.clean.bam"
        output:
            prune_bam = "samples/bam/prune.bam"
        shell: """
        ln -s ${input.merge_bam} ${output.prune_bam}
        """
rule Partitioning:
    message: "~-~ Trimming fastq files... ~-~"
    threads: 10
    params:
        REFERENCE=REFERENCE,
        ENZYME_SITE=ENZYME_SITE,
        ENZYME_SEQ=ENZYME_SEQ,
        PARTITION = PARTITION,
        MAX_PARTITION =MAX_PARTITION
    input:
        prune_bam = "samples/bam/prune.bam"
    output:
        expand("samples/bam/prune.counts_{ENZYME_SEQ}.{PARTITION}g{MAX_PARTITION}.txt", ENZYME_SEQ = ENZYME_SEQ, PARTITION = PARTITION, MAX_PARTITION = MAX_PARTITION),
        "samples/bam/prune.clm"

    shell: """
        ALLHiC_partition -b ${input.prune_bam} -r ${params.REFERENCE} -e ${params.ENZYME_SITE} -k ${params.PARTITION}
    """

rule Optimize:
    message: "~-~ Optimizing... ~-~"
    threads: 10
    params:
        ENZYME_SEQ=ENZYME_SEQ,
        PARTITION = PARTITION,
        MAX_PARTITION =MAX_PARTITION
    input:
        partitions = expand("samples/bam/prune.counts_{ENZYME_SEQ}.{PARTITION}g{MAX_PARTITION}.txt", ENZYME_SEQ = ENZYME_SEQ, PARTITION = PARTITION, MAX_PARTITION = MAX_PARTITION),
        clm = "samples/bam/prune.clm"
    output:
        "samples/bam/prune.counts_{ENZYME_SEQ}.{PARTITION}g{MAX_PARTITION}.tour"
    shell:
        "allhic optimize {input.partitions} {input.clm}"