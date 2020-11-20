import os
from bin import collectGZ as cgz
configfile: "./config.json"
os.makedirs("./output/", exist_ok=True)
os.makedirs("./output/benchmarks/", exist_ok=True)
os.makedirs("./output/cluster/logs", exist_ok=True)

REF_DIR = config["datadir"]
pairs = cgz.getGZPairs(cgz.collectGZ(REF_DIR))
for i in range(0,len(pairs)):
	pairs[i] = pairs[i].replace("ref/","")

rule all:
	input:
		#expand("{datadir}{fasta}.amb", datadir = config["datadir"], fasta = config["fasta"]),
		#expand("{datadir}{fasta}.ann", datadir = config["datadir"], fasta = config["fasta"]),
		#expand("{datadir}{fasta}.bwt", datadir = config["datadir"], fasta = config["fasta"]),
		#expand("{datadir}{fasta}.pac", datadir = config["datadir"], fasta = config["fasta"]),
		#expand("{datadir}{fasta}.sa", datadir = config["datadir"], fasta = config["fasta"]),
		#expand("{datadir}{fasta}.fai", datadir = config["datadir"], fasta = config["fasta"]),
		#expand("output/bwa_algn/sample.bwa_algn.part_001.sam")
		#expand("output/bwa_algn/sample.bwa_algn.sam")
		expand("output/trims/{fastq}1_val_1.fq.gz", fastq = pairs),
		expand("output/trims/{fastq}2_val_2.fq.gz", fastq = pairs)

rule trim_galore:
	message: "~-~ Trimming fastq files... ~-~"
	benchmark: "output/benchmarks/trim_galore_file={fastq}"
	threads: config["threads"]
	input:
		fwd = "ref/{fastq}1.fastq.gz",
		rev = "ref/{fastq}2.fastq.gz"
	output:
		"output/trims/{fastq}1_val_1.fq.gz",
		"output/trims/{fastq}2_val_2.fq.gz"
	shell: "trim_galore -j {threads} --paired --gzip --max_n 50 --three_prime_clip_R1 10 --three_prime_clip_R2 10 -o output/trims/ {input.fwd} {input.rev}" 

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
	benchmark: "output/benchmarks/bwa_algn"
	log: "output/clusters/logs/bwa_align.log"
	input:
		"{fastq1}_val_1.fq.gz",
		"{fastq2}_val_2.fq.gz"
	output:
		"output/bwa_algn/sample_R1.sai",
		"output/bwa_algn/sample_R2.sai"
	params:
		fasta = config["fasta"],
		datadir = config["datadir"],
		trims = config["trims"]
	run:
		shell("bwa aln -t {threads} {params.datadir}{params.fasta} {params.trims}hic_r1.fastq.gz > {output[0]}")
		shell("bwa aln -t {threads} {params.datadir}{params.fasta} {params.trims}hic_r2.fastq.gz > {output[1]}")

rule bwa_sampe:
	message: "~-~ Bwa sampe... ~-~"
	benchmark: "output/benchmarks/bwa_sampe.bm"
	input:
		"output/bwa_algn/sample_R1.sai",
		"output/bwa_algn/sample_R2.sai"
	output:
		"output/bwa_algn/sample.bwa_algn.sam"
	params:
		fasta = config["fasta"],
		datadir = config["datadir"],
		trims = config["trims"]
	run:
		shell("bwa sampe {params.datadir}{params.fasta} {input[0]} {input[1]} {params.trims}hic_r1.fastq.gz {params.trims}hic_r2.fastq.gz > {output[0]}")

