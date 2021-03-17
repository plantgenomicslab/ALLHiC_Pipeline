import os
import subprocess
from bin import collectGZ as cgz
from bin import genSplitString as genSS
configfile: "./config.json"
os.makedirs("./output/", exist_ok=True)
os.makedirs("./output/benchmarks/", exist_ok=True)
os.makedirs("./output/cluster/logs", exist_ok=True)

CWD = os.getcwd()
REF_DIR = config["__default__"]["datadir"]
if REF_DIR[-1] == "/":
	REF_DIR = REF_DIR[:-1]
SNAKE_DIR = config["__default__"]["snakedir"]
if SNAKE_DIR[-1] == "/":
	SNAKE_DIR = SNAKE_DIR[:-1]
PAIRS_R1, PAIRS_R2 = cgz.getGZPairs(cgz.rename_gz(cgz.collectGZ(REF_DIR + "/")))
PAIR_NAMES_1 = []
PAIR_NAMES_2 = []
FQ_DIR = ""
for i in PAIRS_R1[0].split("/")[0:-1]:
	FQ_DIR = FQ_DIR +i + "/"
for i in range(0,len(PAIRS_R1)):
	title = PAIRS_R1[i].split("/")[-1]
	title = title.replace(".fq.gz","")
	part1,part2 = cgz.get_parts(title)
	PAIR_NAMES_1.append(part1)
	PAIR_NAMES_2.append(part2)
PARTS = genSS.partSuffixes(str(config["__default__"]["split"]))
GROUP_LIST = []
for i in range(1, int(config["__default__"]["groups"])+1):
	GROUP_LIST.append(i)

ENZYME_SITE = config["__default__"]["enzyme_sites"]
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
		expand("{ref_dir}/{fasta}.amb", ref_dir = REF_DIR, fasta = config["__default__"]["fasta"]),
		expand("{ref_dir}/{fasta}.ann", ref_dir = REF_DIR, fasta = config["__default__"]["fasta"]),
		expand("{ref_dir}/{fasta}.bwt", ref_dir = REF_DIR, fasta = config["__default__"]["fasta"]),
		expand("{ref_dir}/{fasta}.pac", ref_dir = REF_DIR, fasta = config["__default__"]["fasta"]),
		expand("{ref_dir}/{fasta}.sa", ref_dir = REF_DIR, fasta = config["__default__"]["fasta"]),
		expand("{ref_dir}/{fasta}.fai", ref_dir = REF_DIR, fasta = config["__default__"]["fasta"]),
		#expand("output/trimmed_reads/{prefix}r1{suffix}.fq.gz", prefix = PAIR_NAMES_1, suffix = PAIR_NAMES_2)
		expand("output/ALLHiC/sorted_collection.clean.fixed.counts_{enzyme_seq}.{groups}g{group_section}.tour", enzyme_seq = ENZYME_SEQ, groups = config["__default__"]["groups"], group_section = GROUP_LIST),
		"output/graphs/500K_all_chrs.pdf"

rule trim_galore:
	message: "~-~ Trimming fastq files... ~-~"
	benchmark: "output/benchmarks/trim_galore.prefix={prefix}_r1-2_{suffix}.bm"
	#log: "output/cluster/logs/trim_galore_file={fastq}.log"
	threads: config["trim_galore"]["threads"]
	input:
		fwd = FQ_DIR + "{prefix}r1{suffix}.fq.gz",
		rev = FQ_DIR + "{prefix}r2{suffix}.fq.gz"
	output:
		"output/trimmed_reads/{prefix}r1{suffix}_val_1.fq.gz",
		"output/trimmed_reads/{prefix}r2{suffix}_val_2.fq.gz"
	shell: "trim_galore --paired --fastqc -j {threads} --gzip --max_n 50 --three_prime_clip_R1 10 --three_prime_clip_R2 10 -o output/trimmed_reads/ {input.fwd} {input.rev} || true" 

rule bwa_index:
	message: "~-~ Indexing fasta file with bwa... ~-~"
	benchmark: "output/benchmarks/bwa_index_file={fasta}.bm"
	#log: "output/cluster/logs/bwa_index_file={fasta}.log"
	input:
		fasta = "ref/{fasta}"
	output:
		"ref/{fasta}.amb",
		"ref/{fasta}.ann",
		"ref/{fasta}.bwt",
		"ref/{fasta}.pac",
		"ref/{fasta}.sa"
	shell:
		"bwa index -a bwtsw {input.fasta}"

rule samtools_index:
	message: "~-~ Indexing fasta file with samtools... ~-~"
	benchmark: "output/benchmarks/samtools_index_file={fasta}.bm"
	#log: "output/cluster/logs/samtools_index_file={fasta}.log"
	input:
		fasta = "ref/{fasta}"
	output:
		"ref/{fasta}.fai"
	shell:
		"samtools faidx {input.fasta}"

rule split:
	message: "~-~ Splitting fastq files... ~-~"
	threads: config["split"]["threads"]
	benchmark: "output/benchmarks/splitting_file={prefix}_{strand}.bm"
	#log: "output/cluster/logs/splitting_file={pair}_{strand}_{parts}.log"
	input:
		"output/trimmed_reads/{prefix}_val_{strand}.fq.gz",
		expand("ref/{fasta}.amb", fasta = config["__default__"]["fasta"]),
		expand("ref/{fasta}.fai", fasta = config["__default__"]["fasta"])
	output:
		expand("output/trimmed_reads/{{prefix}}_val_{{strand}}.part_{parts}.fq.gz", parts = PARTS)
	params:
		split = config["__default__"]["split"]
	shell: "seqkit split2 -j {threads} -p {params.split} {input[0]} --out-dir output/trimmed_reads"

rule bwa_align:
	message: "~-~ Aligning trimmed reads... ~-~"
	threads: config["bwa_align"]["threads"]
	benchmark: "output/benchmarks/bwa_algn_file={prefix}_{strand}.part_{parts}.bm"
	#log: "output/cluster/logs/bwa_algn_file={pair}{strand}_{parts}.log"
	input:
		"output/trimmed_reads/{prefix}_val_{strand}.part_{parts}.fq.gz"
	output:
		"output/bwa_algn/{prefix}_val_{strand}.part_{parts}.sai"
	params:
		refdir = config["__default__"]["datadir"],
		fasta = config["__default__"]["fasta"]
	shell: "mkdir bwa_algn/ || true;\
			bwa aln -t {threads} {params.refdir}{params.fasta} {input[0]} > {output[0]}"

rule bwa_sampe:
	message: "~-~ Bwa sampe... ~-~"
	benchmark: "output/benchmarks/bwa_sampe_pre={prefix}r1-2{suffix}.part_{parts}.bm"
	#log: "output/cluster/logs/bwa_sampe_pre={pair}_{parts}.log"
	input:
		"output/bwa_algn/{prefix}r1{suffix}_val_1.part_{parts}.sai",
		"output/bwa_algn/{prefix}r2{suffix}_val_2.part_{parts}.sai",
		"output/trimmed_reads/{prefix}r1{suffix}_val_1.part_{parts}.fq.gz",
		"output/trimmed_reads/{prefix}r2{suffix}_val_2.part_{parts}.fq.gz"
	output:
		"output/bwa_algn/{prefix}r1-2{suffix}.part_{parts}.sam"
	params:
		fasta = config["__default__"]["fasta"],
		datadir = config["__default__"]["datadir"]
	shell: "bwa sampe ref/{params.fasta} {input[0]} {input[1]} {input[2]} {input[3]} > {output[0]}"

rule preprocess_sams:
	message: "~-~ Preparing SAM files... ~-~"
	benchmark: "output/benchmarks/preprocessSAMs_pre={prefix}.bm"
	#log: "output/cluster/logs/preprocessSAMs_pre={pair}.part_{part}.log"
	input:
		"output/bwa_algn/{prefix}.sam"
	output:
		"output/bwa_algn/{prefix}.REduced.paired_only.bam"
	params:
		fasta = config["__default__"]["fasta"],
		snakedir = SNAKE_DIR,
		cwd = CWD,
		refdir = REF_DIR,
		enzyme = config["__default__"]["enzyme_sites"]
	shell: "cd {params.snakedir}/bin/ALLHiC/scripts/;\
			perl PreprocessSAMs.pl {params.cwd}/{input[0]} {params.cwd}/{params.refdir}/{params.fasta} {params.enzyme};"

## SORT BEFORE AND AFTER MERGE ; MERGE SPLIT FILES, NOT ALL FILES INTO ONE
#rule merge_bam:
#	message: "~-~ Merging bam files... ~-~"
#	threads: config["threads"]
#	benchmark: "output/benchmarks/merge_bams.pre_{prefix}.bm"
#	input:
#		expand("output/bwa_algn/{prefix}r1-2{suffix}.part_{parts}.REduced.paired_only.bam", prefix = PAIR_NAMES_1, suffix = PAIR_NAMES_2, parts = PARTS)
#	output:
#		"output/filterSAM/{prefix}.merged_sorted.bam"
#	shell: "mkdir output/filterSAM/ || true;\
#			cd output/filterSAM/;\
#			mv ../bwa_algn/{wildcards.prefix}.REduced* . || true;\
#			for file in {wildcards.prefix}.*.REduced.paired_only.bam; do sambamba sort -t {threads} $file -o sorted.$file; done;\
#			sambamba merge -t {threads} {wildcards.prefix}.merged.bam sorted.{wildcards.prefix}.*.REduced.paired_only.bam;\
#			sambamba sort -t {threads} -o {wildcards.prefix}.merged_sorted.bam {wildcards.prefix}.merged.bam;\
#			cd ../../" 

#rule filterBAM:
#	message: "~-~ Filtering BAM for ALLHiC ~-~"
#	benchmark: "output/benchmarks/filterBAM.pre_{prefix}.bm"
#	#log: "output/cluster/logs/filterBAM_pre={pair}.log"
#	input:
#		"output/filterSAM/{prefix}.merged_sorted.bam"
#	output:
#		"output/filterSAM/{prefix}.clean.sam",
#		"output/ALLHiC/{prefix}.clean.bam"
#	params:
#		fasta = config["__default__"]["fasta"],
#		snakedir = SNAKE_DIR
#	shell: "perl {params.snakedir}/bin/ALLHiC/scripts/filterBAM_forHiC.pl {input[0]} {output[0]};\
#			mkdir output/ALLHiC/ || true;\
#			samtools view -bt ref/{params.fasta}.fai {output[0]} > {output[1]}"

rule filterBAM:
	message: "~-~ Filtering BAM for ALLHiC ~-~"
	benchmark: "output/benchmarks/filterBAM.pre_{prefix}.part_{parts}.bm"
	input:
		"output/bwa_algn/{prefix}.part_{parts}.REduced.paired_only.bam"
	output:
		"output/filterSAM/{prefix}.part_{parts}.clean.sam",
		"output/ALLHiC/{prefix}.part_{parts}.clean.bam"
	params:
		fasta = config["__default__"]["fasta"],
		snakedir = SNAKE_DIR
	shell: "perl {params.snakedir}/bin/ALLHiC/scripts/filterBAM_forHiC.pl {input[0]} {output[0]};\
			mkdir output/ALLHiC/ || true;\
			samtools view -bt ref/{params.fasta}.fai {output[0]} > {output[1]}"

rule merge_clean:
	message: "~-~ Merging Cleaned BAM Files... ~-~"
	benchmark: "output/benchmarks/merged_clean.bm"
	threads: config["merge_clean"]["threads"]
	input:
		expand("output/ALLHiC/{prefix}r1-2{suffix}.part_{parts}.clean.bam", prefix = PAIR_NAMES_1, suffix = PAIR_NAMES_2, parts = PARTS)
	output:
		"output/ALLHiC/sorted_collection.clean.bam"
	shell: "cd output/ALLHiC/;\
			for file in *.part_*.clean.bam; do sambamba sort -p -t {threads} $file -o sorted.$file; done;\
			sambamba merge -p -t {threads} collection.clean.bam sorted.*.part_*.clean.bam;\
			sambamba sort -p -t {threads} -o sorted_collection.clean.bam collection.clean.bam;\
			cd ../../"

rule allhic_corrector:
	message: "~-~ Correcting ALLHiC... ~-~"
	benchmark: "output/benchmarks/allhic_corrector.bm"
	threads: config["allhic_corrector"]["threads"]
	input:
		"output/ALLHiC/sorted_collection.clean.bam"
	output:
		"output/ref_fixed/seq.corrected.fasta"
	params:
		snakedir = SNAKE_DIR,
		refdir = REF_DIR,
		fasta = config["__default__"]["fasta"]
	shell: "mkdir output/ref_fixed/ || true;\
			sambamba index -t {threads} {input[0]} {input[0]}.bai;\
			python {params.snakedir}/bin/ALLHiC/bin/ALLHiC_corrector -t {threads} -m {input[0]} -r {params.refdir}/{params.fasta} -o {output[0]};"

## BWA IDEX  output/ALLHiC/seq.corrected.fasta
rule bwa_index_fixed:
	message: "Indexing fixed fasta file with bwa... ~-~"
	benchmark: "output/benchmarks/bwa_index.fixed.bm"
	input:
		"output/ref_fixed/seq.corrected.fasta"
	output:
		"output/ref_fixed/seq.corrected.fasta.amb",
		"output/ref_fixed/seq.corrected.fasta.ann",
		"output/ref_fixed/seq.corrected.fasta.bwt",
		"output/ref_fixed/seq.corrected.fasta.pac",
		"output/ref_fixed/seq.corrected.fasta.sa"
	shell: "bwa index -a bwtsw {input[0]}"

## Samtools index
rule samtools_index_fixed:
	message: "Indexing fixed fasta file with samtools... ~-~"
	benchmark: "output/benchmarks/samtools_index.fixed.bm"
	input:
		"output/ref_fixed/seq.corrected.fasta"
	output:
		"output/ref_fixed/seq.corrected.fasta.fai"
	shell: "samtools faidx {input[0]}"

## BWA ALN fastq
rule bwa_align_fixed:
	message: "Aligning trimmed reads with fixed fasta file... ~-~"
	benchmark: "output/benchmarks/bwa_algn.{prefix}_{strand}.part_{parts}.fixed.bm"
	threads: config["bwa_align_fixed"]["threads"]
	input:
		trims = "output/trimmed_reads/{prefix}_val_{strand}.part_{parts}.fq.gz",
		fasta = "output/ref_fixed/seq.corrected.fasta",
		index = "output/ref_fixed/seq.corrected.fasta.amb"
	output:
		result = "output/bwa_algn_fixed/{prefix}_val_{strand}.part_{parts}.fixed.sai"
	shell: "mkdir bwa_algn_fixed/ || true;\
			bwa aln -t {threads} {input.fasta} {input.trims} > {output.result}"

## BWA sampe
rule bwa_sampe_fixed:
	message: "Bwa sampe fixed... ~-~"
	benchmark: "output/benchmarks/bwa_sampe_pre={prefix}r1-2{suffix}.part_{parts}.fixed.bm"
	input:
		map1 = "output/bwa_algn_fixed/{prefix}r1{suffix}_val_1.part_{parts}.fixed.sai",
		map2 = "output/bwa_algn_fixed/{prefix}r2{suffix}_val_2.part_{parts}.fixed.sai",
		trims1 = "output/trimmed_reads/{prefix}r1{suffix}_val_1.part_{parts}.fq.gz",
		trims2 = "output/trimmed_reads/{prefix}r2{suffix}_val_2.part_{parts}.fq.gz",
		fasta = "output/ref_fixed/seq.corrected.fasta",
		index = "output/ref_fixed/seq.corrected.fasta.amb"
	output:
		result = "output/bwa_algn_fixed/{prefix}r1-2{suffix}.part_{parts}.fixed.sam"
	shell: "bwa sampe {input.fasta} {input.map1} {input.map2} {input.trims1} {input.trims2} > {output.result}"

## Preprocess sams
rule preprocess_sams_fixed:
	message: "Preparing fixed SAM files... ~-~"
	benchmark: "output/benchmarks/preprocessSAMs_pre={prefix}.fixed.bm"
	#log: "output/cluster/logs/preprocessSAMs_pre={pair}.part_{part}.log"
	input:
		sam = "output/bwa_algn_fixed/{prefix}.fixed.sam",
		fasta = "output/ref_fixed/seq.corrected.fasta"
	output:
		result = "output/bwa_algn_fixed/{prefix}.fixed.REduced.paired_only.bam"
	params:
		snakedir = SNAKE_DIR,
		cwd = CWD,
		enzyme = config["__default__"]["enzyme_sites"]
	shell: "cd {params.snakedir}/bin/ALLHiC/scripts/;\
			perl PreprocessSAMs.pl {params.cwd}/{input.sam} {params.cwd}/{input.fasta} {params.enzyme};"

## Merge bam
#rule merge_bam_fixed:
#	message: "Merging fixed bam files... ~-~"
#	benchmark: "output/benchmarks/merge_bams.pre_{prefix}.bm"
#	threads: config["threads"]
#	input:
#		expand("output/bwa_algn_fixed/{prefix}r1-2{suffix}.part_{parts}.fixed.REduced.paired_only.bam", prefix = PAIR_NAMES_1, suffix = PAIR_NAMES_2, parts = PARTS)
#	output:
#		"output/filterSAM_fixed/{prefix}.merged_sorted.fixed.bam"
#	shell: "mkdir filterSAM_fixed;\
#			cd output/filterSAM_fixed/;\
#			mv ../bwa_algn_fixed/{wildcards.prefix}.fixed.REduced* . || true;\
#			for file in {wildcards.prefix}.*.fixed.REduced.paired_only.bam; do sambamba sort -t {threads} $file -o sorted.$file; done;\
#			sambamba merge -t {threads} {wildcards.prefix}.merged.fixed.bam sorted.{wildcards.prefix}.*.fixed.REduced.paired_only.bam;\
#			sambamba sort -t {threads} -o {wildcards.prefix}.merged_sorted.fixed.bam {wildcards.prefix}.merged.fixed.bam;\
#			cd ../../" 

## Filter bam
rule filterBAM_fixed:
	message: "~-~ Filtering fixed BAM for ALLHiC ~-~"
	benchmark: "output/benchmarks/filterBAM.pre_{prefix}.part_{parts}.fixed.bm"
	input:
		in_bam = "output/bwa_algn_fixed/{prefix}.part_{parts}.fixed.REduced.paired_only.bam",
		fasta = "output/ref_fixed/seq.corrected.fasta.fai"
	output:
		out_sam = "output/filterSAM_fixed/{prefix}.part_{parts}.clean.fixed.sam",
		out_bam = "output/ALLHiC/{prefix}.part_{parts}.clean.fixed.bam"
	params:
		fasta = config["__default__"]["fasta"],
		snakedir = SNAKE_DIR
	shell: "perl {params.snakedir}/bin/ALLHiC/scripts/filterBAM_forHiC.pl {input.in_bam} {output.out_sam};\
			samtools view -bt {input.fasta} {output.out_sam} > {output.out_bam}"

## Merge clean
rule merge_fixed_bam:
	message: "~-~ Merging fixed clean bam files... ~-~"
	benchmark: "output/benchmarks/merged_clean.fixed.bm"
	threads: config["merge_fixed_bam"]["threads"]
	input:
		expand("output/ALLHiC/{prefix}r1-2{suffix}.part_{parts}.clean.fixed.bam", prefix = PAIR_NAMES_1, suffix = PAIR_NAMES_2, parts = PARTS)
	output:
		"output/ALLHiC/sorted_collection.clean.fixed.bam"
	shell: "cd output/ALLHiC/;\
			for file in *.part_*.clean.fixed.bam; do sambamba sort -t {threads} $file -o sorted.$file; done;\
			sambamba merge -t {threads} collection.clean.fixed.bam sorted.*.part_*.clean.fixed.bam;\
			sambamba sort -t {threads} -o sorted_collection.clean.fixed.bam collection.clean.fixed.bam;\
			cd ../../"

rule partition:
	message: "~-~ Partitioning... ~-~"
	benchmark: "output/benchmarks/ALLHiC_partition.bm"
	input:
		clean_bam = "output/ALLHiC/sorted_collection.clean.fixed.bam",
		fasta = "output/ref_fixed/seq.corrected.fasta"
	output:
		expand("output/ALLHiC/sorted_collection.clean.fixed.counts_{enzyme_seq}.{groups}g{group_selection}.txt", enzyme_seq = ENZYME_SEQ, groups = config["__default__"]["groups"], group_selection = GROUP_LIST),
		"output/ALLHiC/sorted_collection.clean.fixed.clm"
	params:
		snakedir = SNAKE_DIR,
		cwd = CWD,
		refdir = REF_DIR,
		groups = config["__default__"]["groups"],
		ezs = config["__default__"]["enzyme_sites"]
	shell: "cd {params.snakedir}/bin/ALLHiC/bin/;\
			ALLHiC_partition -b {params.cwd}/{input.clean_bam} -r {params.cwd}/{input.fasta} -e {params.ezs} -k {params.groups};\
			cd {params.cwd}"

rule optimize:
	message: "~-~ Optimizing... ~-~"
	benchmark: "output/benchmarks/allhic_optimizing.enzyme_{enzyme}.groups_{groups}.part_{group_section}.bm"
	#log: "output/cluster/logs/allhic_optimizing_pre={fastq}.log"
	input:
		group_section = "output/ALLHiC/sorted_collection.clean.fixed.counts_{enzyme}.{groups}g{group_section}.txt",
		clm = "output/ALLHiC/sorted_collection.clean.fixed.clm"
	output:
		"output/ALLHiC/sorted_collection.clean.fixed.counts_{enzyme}.{groups}g{group_section}.tour"
	params:
		snakedir = SNAKE_DIR
	shell: "{params.snakedir}/bin/ALLHiC/bin/allhic optimize {input.group_section} {input.clm}"

### Command must be called in same directory as 
rule build:
	message: "~-~ Building ALLHiC file... ~-~"
	benchmark: "output/benchmarks/allhic_build.bm"
	#log: "output/cluster/logs/allhic_build.log"
	input:
		fixed_fasta = "output/ref_fixed/seq.corrected.fasta",
		ref = expand("output/ALLHiC/sorted_collection.clean.fixed.counts_{enzyme}.{groups}g{group_section}.tour", enzyme = ENZYME_SEQ, groups = config["__default__"]["groups"], group_section = GROUP_LIST)
	output:
		"output/ALLHiC/groups.agp",
		"output/ALLHiC/groups.asm.fasta"
	params:
		snakedir = SNAKE_DIR,
		cwd = CWD
	shell: "cd output/ALLHiC/;\
			perl {params.snakedir}/bin/ALLHiC/bin/ALLHiC_build {params.cwd}/{input.fixed_fasta};\
			cd {params.cwd}"

rule chrn_list:
	message: "~-~ Generating chromosome list... ~-~"
	benchmark: "output/benchmarks/chrnList.bm"
	#log: "output/cluster/logs/chrnList.log"
	input:
		"output/ALLHiC/groups.asm.fasta"
	output:
		"output/ALLHiC/chrn.list",
		"output/ALLHiC/len.txt"
	params:
		snakedir = SNAKE_DIR
	shell: "perl {params.snakedir}/bin/getFaLen.pl -i {input[0]} -o output/ALLHiC/len.txt;\
			grep '.fixed.counts.' output/ALLHiC/len.txt > {output[0]}"

rule plot:
	message: "~-~ Plotting heatmap... ~-~"
	benchmark: "output/benchmarks/plotting_pre.bm"
	#log: "output/cluster/logs/plotting_pre.log"
	input:
		"output/ALLHiC/sorted_collection.clean.fixed.bam",
		"output/ALLHiC/groups.agp",
		"output/ALLHiC/chrn.list"
	output:
		"output/graphs/500K_all_chrs.pdf"
	params:
		snakedir = SNAKE_DIR
	shell: "perl {params.snakedir}/bin/ALLHiC/bin/ALLHiC_plot {input[0]} {input[1]} {input[2]} 500k pdf;\
			mkdir output/graphs || true;\
			mv 500* output/graphs"
