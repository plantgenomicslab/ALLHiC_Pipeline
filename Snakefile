import os
from bin import collectGZ as cgz
configfile: "./config.json"
os.makedirs("./output/", exist_ok=True)
os.makedirs("./output/benchmarks/", exist_ok=True)
os.makedirs("./output/cluster/logs", exist_ok=True)

REF_DIR = config["datadir"]
pairs = cgz.getGZPairs(cgz.collectGZ(REF_DIR))
pair_names = []
for i in range(0,len(pairs)):
	pairs[i] = pairs[i].replace("ref/","")
	pair_names.append(pairs[i].split("/")[-1])

rule all:
	input:
		#expand("{datadir}{fasta}.amb", datadir = config["datadir"], fasta = config["fasta"]),
		#expand("{datadir}{fasta}.ann", datadir = config["datadir"], fasta = config["fasta"]),
		#expand("{datadir}{fasta}.bwt", datadir = config["datadir"], fasta = config["fasta"]),
		#expand("{datadir}{fasta}.pac", datadir = config["datadir"], fasta = config["fasta"]),
		#expand("{datadir}{fasta}.sa", datadir = config["datadir"], fasta = config["fasta"]),
		#expand("{datadir}{fasta}.fai", datadir = config["datadir"], fasta = config["fasta"]),
		#expand("output/bwa_algn/{pair}.bwa_algn.sam", pair = pair_names)
		#expand("output/bwa_algn/sample.bwa_algn.part_001.sam")
		#expand("output/bwa_algn/sample.bwa_algn.sam")
		#expand("output/trims/{fastq}1_trimmed.fq.gz", fastq = pairs),
		#expand("output/trims/{fastq}2_trimmed.fq.gz", fastq = pairs)
		#expand("output/filterSAM/{fastq}.REduced.paired_only.bam", fastq=pair_names)
		#expand("output/filterSAM/{fastq}.clean.sam", fastq=pair_names)
		expand("output/ALLHiC/{fastq}.clean.counts_AAGCTT.16g1.tour", fastq=pair_names),
		"output/ALLHiC/groups.agp",
		#expand("500K_{fastq}.clean.counts_AAGCTT.16g1.pdf", fastq = pair_names)
		"output/graphs/500K_all_chrs.pdf"

rule trim_galore:
	message: "~-~ Trimming fastq files... ~-~"
	benchmark: "output/benchmarks/trim_galore_file={fastq}"
	threads: config["threads"]
	input:
		fwd = expand("ref/{pair}1.fastq.gz", pair = pairs),
		rev = expand("ref/{pair}2.fastq.gz", pair = pairs)
	output:
		"output/trimmed_reads/{fastq}1_val_1.fq.gz",
		"output/trimmed_reads/{fastq}2_val_2.fq.gz"
	shell: "trim_galore --paired -j {threads} --gzip --max_n 50 --three_prime_clip_R1 10 --three_prime_clip_R2 10 -o output/trimmed_reads/ {input.fwd} {input.rev} || true" 

rule bwa_index:
	message: "~-~ Indexing fasta file with bwa... ~-~"
	benchmark: "output/benchmarks/bwa_index_file={fasta}"
	input:
		fasta = "{fasta}"
	output:
		"{fasta}.amb",
		"{fasta}.ann",
		"{fasta}.bwt",
		"{fasta}.pac",
		"{fasta}.sa"
	shell:
		"bwa index -a bwtsw {input.fasta}"

rule samtools_index:
	message: "~-~ Indexing fasta file with samtools... ~-~"
	benchmark: "output/benchmarks/samtools_index_file={fasta}"
	input:
		fasta = "{fasta}"
	output:
		"{fasta}.fai"
	shell:
		"samtools faidx {input.fasta}"

#rule bwa_align:
#	message: "~-~ Aligning trimmed reads... ~-~"
#	benchmark: "output/benchmarks/bwa_algn"
#	log: "output/clusters/logs/bwa_align.log"
#	output:
#		"output/bwa_algn/hic_r1.part_001.fastq.gz",
#		"output/bwa_algn/hic_r2.part_001.fastq.gz"
#	params:
#		split = config["split"],
#		fasta = config["fasta"],
#		datadir = config["datadir"],
#		trims = config["trims"]
#	run:
#		part_decision = ""
#		shell("seqkit split2 -p {params.split} {params.trims}/hic_r1.fastq.gz --out-dir output/bwa_algn/")
#		shell("seqkit split2 -p {params.split} {params.trims}/hic_r2.fastq.gz --out-dir output/bwa_algn/")
#		for i in range(1,params.split+1):
#			part_decision = ".part_00" + str(i)
#			R1 = "output/bwa_algn/sample_R1" + part_decision + ".sai"
#			R2 = "output/bwa_algn/sample_R2" + part_decision + ".sai"
#			shell("bwa aln -t {threads} {params.datadir}{params.fasta} output/bwa_algn/hic_r1" + part_decision + ".fastq.gz > " + R1)
#			shell("bwa aln -t {threads} {params.datadir}{params.fasta} output/bwa_algn/hic_r2" + part_decision + ".fastq.gz > " + R2)
#			#shell("bwa sampe {params.datadir}{params.fasta} " + R1 + " " + R2 + " {params.trims}hic_R1.fastq.gz {params.trims}hic_R2.fastq.gz > {output[0]}")

#rule bwa_sampe:
#	message: "~-~ Bwa sampe... ~-~"
#	benchmark: "output/benchmarks/bwa_sampe.bm"
#	input:
#		"output/bwa_algn/hic_r1.part_001.fastq.gz",
#		"output/bwa_algn/hic_r2.part_001.fastq.gz"
#	output:
#		"output/bwa_algn/sample.bwa_algn.part_001.sam"
#	params:
#		split = config["split"],
#		fasta = config["fasta"],
#		datadir = config["datadir"],
#		trims = config["trims"]
#	run:
#		for i in range(1,params.split+1):
#			part_decision = ".part_00" + str(i)
#			R1 = "output/bwa_algn/sample_R1" + part_decision + ".sai"
#			R2 = "output/bwa_algn/sample_R2" + part_decision + ".sai"
#			shell("bwa sampe {params.datadir}{params.fasta} " + R1 + " " + R2 + " {params.trims}hic_r1.fastq.gz {params.trims}hic_r2.fastq.gz > output/bwa_algn/sample.bwa_algn.part_00" + str(i) + ".sam")
rule bwa_align:
	message: "~-~ Aligning trimmed reads... ~-~"
	benchmark: "output/benchmarks/bwa_algn_pre={fastq}.bm"
	log: "output/clusters/logs/bwa_align_pre={fastq}.log"
	input:
		"output/trimmed_reads/{fastq}1_val_1.fq.gz",
		"output/trimmed_reads/{fastq}2_val_2.fq.gz"
	output:
		"output/bwa_algn/{fastq}_R1.sai",
		"output/bwa_algn/{fastq}_R2.sai"
	params:
		fasta = config["fasta"],
		datadir = config["datadir"],
		trims = config["trims"]
	run:
		os.makedirs("output/bwa_algn/", exist_ok=True)
		shell("bwa aln -t {threads} {params.datadir}{params.fasta} {input[0]} > {output[0]}")
		shell("bwa aln -t {threads} {params.datadir}{params.fasta} {input[1]} > {output[1]}")

rule bwa_sampe:
	message: "~-~ Bwa sampe... ~-~"
	benchmark: "output/benchmarks/bwa_sampe_pre={fastq}.bm"
	input:
		"output/bwa_algn/{fastq}_R1.sai",
		"output/bwa_algn/{fastq}_R2.sai"
	output:
		"output/bwa_algn/{fastq}.sam"
	params:
		fasta = config["fasta"],
		datadir = config["datadir"],
		trims = config["trims"]
	run:
		shell("bwa sampe {params.datadir}{params.fasta} {input[0]} {input[1]} output/trimmed_reads/hic_r1_val_1.fq.gz output/trimmed_reads/hic_r2_val_2.fq.gz > {output[0]}")

rule preprocess_sams:
	message: "~-~ Preparing SAM files... ~-~"
	benchmark: "output/benchmarks/preprocessSAMs_pre={fastq}.bm"
	input:
		"output/bwa_algn/{fastq}.sam"
	output:
		"output/filterSAM/{fastq}.REduced.paired_only.bam"
	params:
		fasta = config["fasta"]
	shell: "cd bin/ALLHiC/scripts/;\
			perl PreprocessSAMs.pl ../../../{input[0]} ../../../ref/{params.fasta} MBOI;\
			mkdir ../../../output/filterSAM/ || true;\
			mv ../../../output/bwa_algn/*REduced* ../../../output/filterSAM/;\
			cd ../../../"

rule filterBAM:
	message: "~-~ Filtering BAM for ALLHiC ~-~"
	benchmark: "output/benchmarks/filterBAM_pre={fastq}.bm"
	input:
		"output/filterSAM/{fastq}.REduced.paired_only.bam"
	output:
		"output/filterSAM/{fastq}.clean.sam",
		"output/filterSAM/{fastq}.clean.bam"
	params:
		fasta = config["fasta"]
	shell: "perl bin/ALLHiC/scripts/filterBAM_forHiC.pl {input[0]} {output[0]};\
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
	input:
		"output/filterSAM/{fastq}.clean.bam"
	output:
		"output/ALLHiC/{fastq}.clean.counts_AAGCTT.16g1.txt",
		"output/ALLHiC/{fastq}.clean.clm"
	params:
		fasta = config["fasta"]
	shell: "mkdir output/ALLHiC/ || true;\
			cp {input[0]} output/ALLHiC/;\
			cd bin/ALLHiC/bin/;\
			ALLHiC_partition -b ../../../output/ALLHiC/{wildcards.fastq}.clean.bam -r ../../../ref/{params.fasta} -e AAGCTT -k 16;\
			cd ../../.."

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
	input:
		group = "output/ALLHiC/{fastq}.clean.counts_AAGCTT.16g1.txt",
		clm = "output/ALLHiC/{fastq}.clean.clm"
	output:
		"output/ALLHiC/{fastq}.clean.counts_AAGCTT.16g1.tour"
	shell: "bin/ALLHiC/bin/allhic optimize {input.group} {input.clm}"

rule build:
	message: "~-~ Building ALLHiC file... ~-~"
	benchmark: "output/benchmarks/allhic_build.bm"
	input:
		expand("output/ALLHiC/{fastq}.clean.counts_AAGCTT.16g1.tour", fastq=pair_names)
	output:
		"output/ALLHiC/groups.agp",
		"output/ALLHiC/groups.asm.fasta"
	params:
		fasta = config["fasta"]
	shell: "cd output/ALLHiC/;\
			perl ../../bin/ALLHiC/bin/ALLHiC_build ../../ref/{params.fasta};\
			cd ../../"

rule chrn_list:
	message: "~-~ Generating chromosome list... ~-~"
	benchmark: "output/benchmarks/chrnList.bm"
	log: "output/clusters/logs/chrnList.log"
	input:
		"output/ALLHiC/groups.asm.fasta"
	output:
		"output/ALLHiC/chrn.list"
		#"output/ALLHiC/len.txt"
	shell: "perl bin/getFaLen.pl -i {input[0]} -o output/ALLHiC/len.txt;\
			grep 'clean.counts' output/ALLHiC/len.txt > {output[0]}"

rule plot:
	message: "~-~ Plotting heatmap... ~-~"
	benchmark: "output/benchmarks/plotting_pre.bm"
	input:
		expand("output/filterSAM/{fastq}.clean.bam", fastq=pair_names),
		"output/ALLHiC/groups.agp",
		"output/ALLHiC/chrn.list"
	output:
		"output/graphs/500K_all_chrs.pdf"
	shell: "perl bin/ALLHiC/bin/ALLHiC_plot {input[0]} {input[1]} {input[2]} 500k pdf;\
			mkdir output/graphs || true;\
			mv 500* output/graphs"
