import os
import subprocess
from bin import collectGZ as cgz
from bin import genSplitString as genSS
configfile: "./config.json"
os.makedirs("./output/", exist_ok=True)
os.makedirs("./output/benchmarks/", exist_ok=True)
os.makedirs("./output/cluster/logs", exist_ok=True)

CWD = os.getcwd()
REF_DIR = config["datadir"]
if REF_DIR[-1] == "/":
	REF_DIR = REF_DIR[:-1]
SNAKE_DIR = config["snakedir"]
if SNAKE_DIR[-1] == "/":
	SNAKE_DIR = SNAKE_DIR[:-1]
#print(cgz.collectGZ(REF_DIR + "/"))
PAIRS_R1, PAIRS_R2 = cgz.getGZPairs(cgz.collectGZ(REF_DIR + "/"))
PAIR_NAMES_1 = []
PAIR_NAMES_2 = []
DIFFS = []
FQ_DIR = ""
for i in PAIRS_R1[0].split("/")[0:-1]:
	FQ_DIR = FQ_DIR +i + "/"
for i in range(0,len(PAIRS_R1)):
	title = PAIRS_R1[i].split("/")[-1]
	title = title.replace(".fastq.gz","")
	title = title.replace(".fq.gz","")
	part1,part2,diff = cgz.get_parts(title)
	PAIR_NAMES_1.append(part1)
	PAIR_NAMES_2.append(part2)
	DIFFS.append(diff[0])
PARTS = genSS.partSuffixes(str(config["split"]))
GROUP_LIST = []
for i in range(1, int(config["groups"])+1):
	GROUP_LIST.append(i)
STRANDS = ["1","2"]
rule all:
	input:
		expand("{ref_dir}/{fasta}.amb", ref_dir = REF_DIR, fasta = config["fasta"]),
		expand("{ref_dir}/{fasta}.ann", ref_dir = REF_DIR, fasta = config["fasta"]),
		expand("{ref_dir}/{fasta}.bwt", ref_dir = REF_DIR, fasta = config["fasta"]),
		expand("{ref_dir}/{fasta}.pac", ref_dir = REF_DIR, fasta = config["fasta"]),
		expand("{ref_dir}/{fasta}.sa", ref_dir = REF_DIR, fasta = config["fasta"]),
		expand("{ref_dir}/{fasta}.fai", ref_dir = REF_DIR, fasta = config["fasta"]),
		expand("output/trimmed_reads/{part1}" + DIFFS[0] + "1{part2}_val_1.fq.gz", part1 = PAIR_NAMES_1, diff = DIFFS, part2 = PAIR_NAMES_2),
		expand("output/trimmed_reads/{part1}{diff}{strand}{part2}_val_{strand}.part_{parts}.fq.gz", part1 = PAIR_NAMES_1, diff = DIFFS, part2 = PAIR_NAMES_2, strand = STRANDS, parts = PARTS),
		expand("output/bwa_algn/{part1}_strand={diff}{strand}_{part2}_{strand}_part_{parts}.sai", part1 = PAIR_NAMES_1, diff = DIFFS, part2 = PAIR_NAMES_2, strand = STRANDS, parts = PARTS),
		expand("output/bwa_algn/p1={part1}_strand={diff}{strand1}-{strand2}_p2={part2}.part_{part}.sam", part1 = PAIR_NAMES_1, diff = DIFFS, strand1 = "1", strand2 = "2", part2 = PAIR_NAMES_2, part = PARTS),
		expand("output/filterSAM/p1={part1}_strand={diff}1-2_p2={part2}.merged.bam", part1 = PAIR_NAMES_1, diff = DIFFS, part2 = PAIR_NAMES_2),
		#expand("output/ALLHiC/p1={part1}_strand={diff}1-2_p2={part2}.clean.bam", part1 = PAIR_NAMES_1, diff = DIFFS, part2 = PAIR_NAMES_2),
		#"output/ALLHiC/collection.clean.bam",
		expand("output/ALLHiC/sorted_collection.clean.counts_{enzyme}.{groups}g1.txt", enzyme = config["enzyme_sites"], groups = config["groups"]),
		expand("output/ALLHiC/sorted_collection.clean.counts_{enzyme}.{groups}g{group_section}.tour", enzyme = config["enzyme_sites"], groups = config["groups"], group_section = GROUP_LIST)
		#"output/graphs/500K_all_chrs.pdf"

rule trim_galore:
	message: "~-~ Trimming fastq files... ~-~"
	benchmark: "output/benchmarks/trim_galore_file={part1}{diff}{part2}.bm"
	#log: "output/cluster/logs/trim_galore_file={fastq}.log"
	threads: config["threads"]
	input:
		fwd = FQ_DIR + "{part1}" + DIFFS[0] + "1{part2}.fastq.gz",
		rev = FQ_DIR + "{part1}" + DIFFS[0] + "2{part2}.fastq.gz"
	output:
		"output/trimmed_reads/{part1}{diff}1{part2}_val_1.fq.gz",
		"output/trimmed_reads/{part1}{diff}2{part2}_val_2.fq.gz"
	shell: "trim_galore --paired -j {threads} --gzip --max_n 50 --three_prime_clip_R1 10 --three_prime_clip_R2 10 -o output/trimmed_reads/ {input.fwd} {input.rev} || true" 

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
	threads: config["threads"]
	benchmark: "output/benchmarks/splitting_file={part1}{diff}{strand}{part2}_{strand}_parts={parts}.bm"
	#log: "output/cluster/logs/splitting_file={pair}_{strand}_{parts}.log"
	input:
		"output/trimmed_reads/{part1}{diff}{strand}{part2}_val_{strand}.fq.gz",
		expand("ref/{fasta}.amb", fasta = config["fasta"]),
		expand("ref/{fasta}.fai", fasta = config["fasta"])
	output:
		"output/trimmed_reads/{part1}{diff}{strand}{part2}_val_{strand}.part_{parts}.fq.gz"
	params:
		split = config["split"]
	shell: "seqkit split2 -j {threads} -p {params.split} {input[0]} --out-dir output/trimmed_reads"

rule bwa_align:
	message: "~-~ Aligning trimmed reads... ~-~"
	threads: config["threads"]
	benchmark: "output/benchmarks/bwa_algn_file={part1}{diff}{strand}{part2}_{parts}.bm"
	#log: "output/cluster/logs/bwa_algn_file={pair}{strand}_{parts}.log"
	input:
		"output/trimmed_reads/{part1}{diff}{strand}{part2}_val_{strand}.part_{parts}.fq.gz"
	output:
		"output/bwa_algn/{part1}_strand={diff}{strand}_{part2}_{strand}_part_{parts}.sai"
	params:
		refdir = config["datadir"],
		fasta = config["fasta"]
	shell: "mkdir bwa_algn/ || true;\
			bwa aln -t {threads} {params.refdir}{params.fasta} {input[0]} > {output[0]}"

rule bwa_sampe:
	message: "~-~ Bwa sampe... ~-~"
	benchmark: "output/benchmarks/bwa_sampe_pre={part1}{diff}{strand1}{strand2}{part2}_{parts}.bm"
	#log: "output/cluster/logs/bwa_sampe_pre={pair}_{parts}.log"
	input:
		"output/bwa_algn/{part1}_strand={diff}1_{part2}_1_part_{parts}.sai",
		"output/bwa_algn/{part1}_strand={diff}2_{part2}_2_part_{parts}.sai",
		"output/trimmed_reads/{part1}{diff}1{part2}_val_1.part_{parts}.fq.gz",
		"output/trimmed_reads/{part1}{diff}2{part2}_val_2.part_{parts}.fq.gz"
	output:
		"output/bwa_algn/{prefix}.sam"
	params:
		fasta = config["fasta"],
		datadir = config["datadir"],
		trims = config["trims"]
	shell: "bwa sampe ref/{params.fasta} {input[0]} {input[1]} {input[2]} {input[3]} > {output[0]}"

rule preprocess_sams:
	message: "~-~ Preparing SAM files... ~-~"
	benchmark: "output/benchmarks/preprocessSAMs_pre={part1}{diff}{part2}.part_{part}.bm"
	#log: "output/cluster/logs/preprocessSAMs_pre={pair}.part_{part}.log"
	input:
		"output/bwa_algn/{prefix}.sam"
	output:
		"output/bwa_algn/{prefix}.REduced.paired_only.bam"
	params:
		fasta = config["fasta"],
		snakedir = SNAKE_DIR,
		cwd = CWD,
		refdir = REF_DIR,
		enzyme = config["enzyme_sites"]
	shell: "cd {params.snakedir}/bin/ALLHiC/scripts/;\
			perl PreprocessSAMs.pl {params.cwd}/{input[0]} {params.cwd}/{params.refdir}/{params.fasta} {params.enzyme};"
			#cd {params.cwd};\
			#mkdir output/filterSAM/ || true;\
			#mv output/bwa_algn/*REduced* output/filterSAM/ || true"

## SORT BEFORE AND AFTER MERGE ; MERGE SPLIT FILES, NOT ALL FILES INTO ONE
rule merge_bam:
	message: "~-~ Merging bam files... ~-~"
	threads: config["threads"]
	benchmark: "output/benchmarks/merge_bams.pre_{prefix}.bm"
	output:
		"output/filterSAM/{prefix}.merged_sorted.bam"
	shell: "mkdir output/filterSAM/ || true;\
			cd output/filterSAM/;\
			mv ../bwa_algn/*.REduced* . || true;\
			for file in {wildcards.prefix}.*.REduced.paired_only.bam; do samtools sort -@ {threads} $file -o sorted.$file; done;\
			samtools merge -f -u {wildcards.prefix}.merged.bam sorted.{wildcards.prefix}.*.REduced.paired_only.bam -@ {threads};\
			samtools sort -@ {threads} {wildcards.prefix}.merged.bam -o {wildcards.prefix}.merged_sorted.bam;\
			cd ../../" 

rule filterBAM:
	message: "~-~ Filtering BAM for ALLHiC ~-~"
	benchmark: "output/benchmarks/filterBAM.pre_{prefix}.bm"
	#log: "output/cluster/logs/filterBAM_pre={pair}.log"
	input:
		"output/filterSAM/{prefix}.merged_sorted.bam"
	output:
		"output/filterSAM/{prefix}.clean.sam",
		"output/ALLHiC/{prefix}.clean.bam"
	params:
		fasta = config["fasta"],
		snakedir = SNAKE_DIR
	shell: "perl {params.snakedir}/bin/ALLHiC/scripts/filterBAM_forHiC.pl {input[0]} {output[0]};\
			mkdir output/ALLHiC/ || true;\
			samtools view -bt ref/{params.fasta}.fai {output[0]} > {output[1]}"

## MERGE ALONE FIRST (IF ERROR, SORT) ; MERGE ALL FILES INTO ONE
rule merge_clean: ### ADD
	message: "~-~ Merging Cleaned BAM Files... ~-~"
	benchmark: "output/benchmarks/merged_clean.bm"
	threads: config["threads"]
	output:
		"output/ALLHiC/sorted_collection.clean.bam"
	shell: "cd output/ALLHiC/;\
			for file in *.clean.bam; do samtools sort -@ {threads} $file -o sorted.$file; done;\
			samtools merge -f -u collection.clean.bam sorted.*.clean.bam -@ {threads};\
			samtools sort -@ {threads} collection.clean.bam -o sorted_collection.clean.bam;\
			cd ../../"

rule partition:
	message: "~-~ Partitioning... ~-~"
	benchmark: "output/benchmarks/ALLHiC_partition.enzyme={enzyme}.groups={groups}.bm"
	input:
		"output/ALLHiC/sorted_collection.clean.bam"
	output:
		"output/ALLHiC/sorted_collection.clean.counts_{enzyme}.{groups}g1.txt"
		#"output/ALLHiC/p1={part1}_strand={diff}1-2_p2={part2}.clean.clm"
	params:
		fasta = config["fasta"],
		snakedir = SNAKE_DIR,
		cwd = CWD,
		refdir = REF_DIR,
		groups = config["groups"],
		ezs = config["enzyme_sites"]
	shell: "cd {params.snakedir}/bin/ALLHiC/bin/;\
			ALLHiC_partition -b {params.cwd}/{input[0]} -r {params.cwd}/{params.refdir}/{params.fasta} -e {params.ezs} -k {params.groups};\
			cd {params.cwd}"

rule optimize:
	message: "~-~ Optimizing... ~-~"
	benchmark: "output/benchmarks/allhic_optimizing.enzyme_{enzyme}.groups_{groups}_sec_{group_section}.bm"
	#log: "output/cluster/logs/allhic_optimizing_pre={fastq}.log"
	input:
		group_section = "output/ALLHiC/sorted_collection.clean.counts_{enzyme}.{groups}g{group_section}.txt",
		clm = "output/ALLHiC/sorted_collection.clean.clm"
	output:
		"output/ALLHiC/sorted_collection.clean.counts_{enzyme}.{groups}g{group_section}.tour"
	params:
		snakedir = SNAKE_DIR
	shell: "{params.snakedir}/bin/ALLHiC/bin/allhic optimize {input.group_section} {input.clm}"

rule build:
	message: "~-~ Building ALLHiC file... ~-~"
	benchmark: "output/benchmarks/allhic_build.bm"
	#log: "output/cluster/logs/allhic_build.log"
	#input:
	#	expand("output/ALLHiC/{fastq}.clean.counts_AAGCTT.16g1.tour", fastq=PAIR_NAMES)
	output:
		"output/ALLHiC/groups.agp",
		"output/ALLHiC/groups.asm.fasta"
	params:
		fasta = config["fasta"],
		refdir = REF_DIR,
		snakedir = SNAKE_DIR,
		cwd = CWD
	shell: "cd output/ALLHiC/;\
			perl {params.snakedir}/bin/ALLHiC/bin/ALLHiC_build {params.cwd}/{params.refdir}/{params.fasta};\
			cd {params.cwd}"

rule chrn_list:
	message: "~-~ Generating chromosome list... ~-~"
	benchmark: "output/benchmarks/chrnList.bm"
	#log: "output/cluster/logs/chrnList.log"
	input:
		"output/ALLHiC/groups.asm.fasta"
	output:
		"output/ALLHiC/chrn.list"
		#"output/ALLHiC/len.txt"
	params:
		snakedir = SNAKE_DIR
	shell: "perl {params.snakedir}/bin/getFaLen.pl -i {input[0]} -o output/ALLHiC/len.txt;\
			grep 'merge.clean.counts' output/ALLHiC/len.txt > {output[0]}"

rule plot:
	message: "~-~ Plotting heatmap... ~-~"
	benchmark: "output/benchmarks/plotting_pre.bm"
	#log: "output/cluster/logs/plotting_pre.log"
	input:
		"output/ALLHiC/sorted_collection.clean.bam",
		"output/ALLHiC/groups.agp",
		"output/ALLHiC/chrn.list"
	output:
		"output/graphs/500K_all_chrs.pdf"
	params:
		snakedir = SNAKE_DIR
	shell: "perl {params.snakedir}/bin/ALLHiC/bin/ALLHiC_plot {input[0]} {input[1]} {input[2]} 500k pdf;\
			mkdir output/graphs || true;\
			mv 500* output/graphs"
