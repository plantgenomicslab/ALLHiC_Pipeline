import os
from bin import collectGZ as cgz
from bin import genSplitString as genSS
configfile: "./config.json"
os.makedirs("./output/", exist_ok=True)
os.makedirs("./output/benchmarks/", exist_ok=True)
#os.makedirs("./output/cluster/logs", exist_ok=True)

CWD = os.getcwd()
REF_DIR = config["datadir"]
if REF_DIR[-1] == "/":
	REF_DIR = REF_DIR[:-1]
SNAKE_DIR = config["snakedir"]
if SNAKE_DIR[-1] == "/":
	SNAKE_DIR = SNAKE_DIR[:-1]
PAIRS = cgz.getGZPairs(cgz.collectGZ(REF_DIR + "/"))
PAIR_NAMES = []
for i in range(0,len(PAIRS)):
	PAIRS[i] = PAIRS[i].replace("ref/","")
	PAIR_NAMES.append(PAIRS[i].split("/")[-1])
PARTS = genSS.partSuffixes(str(config["split"]))
STRANDS = ["1","2"]

rule all:
	input:
		expand("{ref_dir}/{fasta}.amb", ref_dir = REF_DIR, fasta = config["fasta"]),
		expand("{ref_dir}/{fasta}.ann", ref_dir = REF_DIR, fasta = config["fasta"]),
		expand("{ref_dir}/{fasta}.bwt", ref_dir = REF_DIR, fasta = config["fasta"]),
		expand("{ref_dir}/{fasta}.pac", ref_dir = REF_DIR, fasta = config["fasta"]),
		expand("{ref_dir}/{fasta}.sa", ref_dir = REF_DIR, fasta = config["fasta"]),
		expand("{ref_dir}/{fasta}.fai", ref_dir = REF_DIR, fasta = config["fasta"]),
		expand("output/bwa_algn/{pair}_{strand}_part_{parts}.sai", pair = PAIR_NAMES, strand = STRANDS, parts = PARTS),
		expand("output/bwa_algn/{pair}_{part}.sam", pair = PAIR_NAMES, part = PARTS),
		expand("output/filterSAM/{pair}_{part}.REduced.paired_only.bam", pair = PAIR_NAMES, part = PARTS),
		expand("output/filterSAM/{pair}.REduced.paired_only.bam.merged", pair = PAIR_NAMES),
		expand("output/filterSAM/{pair}.clean.bam", pair = PAIR_NAMES),
		"output/graphs/500K_all_chrs.pdf"

rule trim_galore:
	message: "~-~ Trimming fastq files... ~-~"
	benchmark: "output/benchmarks/trim_galore_file={fastq}.bm"
	#log: "output/cluster/logs/trim_galore_file={fastq}.log"
	threads: config["threads"]
	input:
		fwd = expand("ref/{pre}1.fastq.gz", pre = PAIRS),
		rev = expand("ref/{pre}2.fastq.gz", pre = PAIRS)
	output:
		"output/trimmed_reads/{fastq}1_val_1.fq.gz",
		"output/trimmed_reads/{fastq}2_val_2.fq.gz"
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
	benchmark: "output/benchmarks/splitting_file={pair}_{strand}_{parts}.bm"
	#log: "output/cluster/logs/splitting_file={pair}_{strand}_{parts}.log"
	input:
		"output/trimmed_reads/{pair}{strand}_val_{strand}.fq.gz",
		expand("ref/{fasta}.amb", fasta = config["fasta"]),
		expand("ref/{fasta}.fai", fasta = config["fasta"])
	output:
		"output/trimmed_reads/{pair}{strand}_val_{strand}.part_{parts}.fq.gz"
	params:
		split = config["split"]
	shell: "seqkit split2 -j {threads} -p {params.split} {input[0]} --out-dir output/trimmed_reads"

rule bwa_align:
	message: "~-~ Aligning trimmed reads... ~-~"
	threads: config["threads"]
	benchmark: "output/benchmarks/bwa_algn_file={pair}{strand}_{parts}.bm"
	#log: "output/cluster/logs/bwa_algn_file={pair}{strand}_{parts}.log"
	input:
		"output/trimmed_reads/{pair}{strand}_val_{strand}.part_{parts}.fq.gz"
	output:
		"output/bwa_algn/{pair}_{strand}_part_{parts}.sai"
	shell: "mkdir bwa_algn/ || true;\
			bwa aln -t {threads} ref/canu.contigs.fasta {input[0]} > {output[0]}"

rule bwa_sampe:
	message: "~-~ Bwa sampe... ~-~"
	benchmark: "output/benchmarks/bwa_sampe_pre={pair}_{parts}.bm"
	#log: "output/cluster/logs/bwa_sampe_pre={pair}_{parts}.log"
	input:
		"output/bwa_algn/{pair}_1_part_{parts}.sai",
		"output/bwa_algn/{pair}_2_part_{parts}.sai",
		"output/trimmed_reads/{pair}1_val_1.part_{parts}.fq.gz",
		"output/trimmed_reads/{pair}2_val_2.part_{parts}.fq.gz"
	output:
		"output/bwa_algn/{pair}_{parts}.sam"
	params:
		fasta = config["fasta"],
		datadir = config["datadir"],
		trims = config["trims"]
	shell: "bwa sampe ref/{params.fasta} {input[0]} {input[1]} {input[2]} {input[3]} > {output[0]}"

rule preprocess_sams:
	message: "~-~ Preparing SAM files... ~-~"
	benchmark: "output/benchmarks/preprocessSAMs_pre={pair}.part_{part}.bm"
	#log: "output/cluster/logs/preprocessSAMs_pre={pair}.part_{part}.log"
	input:
		"output/bwa_algn/{pair}_{part}.sam"
	output:
		"output/filterSAM/{pair}_{part}.REduced.paired_only.bam"
	params:
		fasta = config["fasta"],
		snakedir = SNAKE_DIR,
		cwd = CWD,
		refdir = REF_DIR
	shell: "cd {params.snakedir}/bin/ALLHiC/scripts/;\
			perl PreprocessSAMs.pl {params.cwd}/{input[0]} {params.cwd}/{params.refdir}/{params.fasta} MBOI;\
			cd {params.cwd};\
			mkdir output/filterSAM/ || true;\
			mv output/bwa_algn/*REduced* output/filterSAM/ || true"

rule merge_bam:
	message: "~-~ Merging bam files... ~-~"
	threads: config["threads"]
	benchmark: "output/benchmarks/merge_bams.file={pair}.bm"
	input:
		expand("output/filterSAM/{pair}_{part}.REduced.paired_only.bam", pair = PAIR_NAMES, part = PARTS)
	output:
		"output/filterSAM/{pair}.REduced.paired_only.bam.merged"
	shell: "samtools merge -u {output[0]} {input[0]} -@ {threads}" 

rule filterBAM:
	message: "~-~ Filtering BAM for ALLHiC ~-~"
	benchmark: "output/benchmarks/filterBAM_pre={pair}.bm"
	#log: "output/cluster/logs/filterBAM_pre={pair}.log"
	input:
		"output/filterSAM/{pair}.REduced.paired_only.bam.merged"
	output:
		"output/filterSAM/{pair}.clean.sam",
		"output/filterSAM/{pair}.clean.bam"
	params:
		fasta = config["fasta"],
		snakedir = SNAKE_DIR
	shell: "perl {params.snakedir}/bin/ALLHiC/scripts/filterBAM_forHiC.pl {input[0]} {output[0]};\
			samtools view -bt ref/{params.fasta}.fai {output[0]} > {output[1]}"

#rule prune:
#	message: "~-~ Performing ALLHiC_prune... ~-~"
#	benchmark: "output/benchmarks/ALLHiC_prune.bm"
#	input:
#		# allele table?
#		# filterBAM output
#	output:
#		"output/prune/prunning.bam"
#	shell: "bin/ALLHiC/bin/ALLHiC_prune -i [input table] -b [input filterBAM] -r [draft fasta];\
#			mkdir output/prune;\
#			mv prunning.bam output/prune"
#
#####################################################
##
## - What to use for enzyme site [-e]?
## - What value to use for # of groups [-k]?
##
#####################################################
rule partition:
	message: "~-~ Partitioning... ~-~"
	benchmark: "output/benchmarks/ALLHiC_partition_pre={fastq}.bm"
	#log: "output/cluster/logs/ALLHiC_partition_pre={fastq}.log"
	input:
		"output/filterSAM/{fastq}.clean.bam"
	output:
		"output/ALLHiC/{fastq}.clean.counts_AAGCTT.16g1.txt",
		"output/ALLHiC/{fastq}.clean.clm"
	params:
		fasta = config["fasta"],
		snakedir = SNAKE_DIR,
		cwd = CWD,
		refdir = REF_DIR,
		groups = config["groups"],
		ezs = config["enzyme_sites"]
	shell: "mkdir output/ALLHiC/ || true;\
			cp {input[0]} output/ALLHiC/;\
			cd {params.snakedir}/bin/ALLHiC/bin/;\
			ALLHiC_partition -b {params.cwd}/output/ALLHiC/{wildcards.fastq}.clean.bam -r {params.cwd}/{params.refdir}/{params.fasta} -e {params.ezs} -k {params.groups};\
			cd {params.cwd}"

#####################################################
##
## - Note: clusters and counts can be generated using
##		"allhic extract"
##
#####################################################
#rule rescue:
#	message: "~-~ Performing rescue... ~-~"
#	benchmark: "output/benchmarks/ALLHiC_rescue.bm"
#	input:
#		# Partition output
#		# clusters
#		# counts
#	output:
#		# To be determined.
#	shell: "bin/ALLHiC/bin/ALLHiC_rescue -b [input] -r [draft fasta] -c [clusters???] -i [counts???];\
#			mkdir output/rescue/;\
#			mv [output] output/rescue/"
#
rule optimize:
	message: "~-~ Optimizing... ~-~"
	benchmark: "output/benchmarks/allhic_optimizing_pre={fastq}.bm"
	#log: "output/cluster/logs/allhic_optimizing_pre={fastq}.log"
	input:
		group = "output/ALLHiC/{fastq}.clean.counts_AAGCTT.16g1.txt",
		clm = "output/ALLHiC/{fastq}.clean.clm"
	output:
		"output/ALLHiC/{fastq}.clean.counts_AAGCTT.16g1.tour"
	params:
		snakedir = SNAKE_DIR
	shell: "{params.snakedir}/bin/ALLHiC/bin/allhic optimize {input.group} {input.clm}"

rule build:
	message: "~-~ Building ALLHiC file... ~-~"
	benchmark: "output/benchmarks/allhic_build.bm"
	#log: "output/cluster/logs/allhic_build.log"
	input:
		expand("output/ALLHiC/{fastq}.clean.counts_AAGCTT.16g1.tour", fastq=PAIR_NAMES)
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
			grep 'clean.counts' output/ALLHiC/len.txt > {output[0]}"

rule plot:
	message: "~-~ Plotting heatmap... ~-~"
	benchmark: "output/benchmarks/plotting_pre.bm"
	#log: "output/cluster/logs/plotting_pre.log"
	input:
		expand("output/filterSAM/{fastq}.clean.bam", fastq=PAIR_NAMES),
		"output/ALLHiC/groups.agp",
		"output/ALLHiC/chrn.list"
	output:
		"output/graphs/500K_all_chrs.pdf"
	params:
		snakedir = SNAKE_DIR
	shell: "perl {params.snakedir}/bin/ALLHiC/bin/ALLHiC_plot {input[0]} {input[1]} {input[2]} 500k pdf;\
			mkdir output/graphs || true;\
			mv 500* output/graphs"
