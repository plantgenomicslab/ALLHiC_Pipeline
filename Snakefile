import os
configfile: "./config.json"
os.makedirs("./output/", exist_ok=True)
os.makedirs("./output/benchmarks/", exist_ok=True)
os.makedirs("./output/cluster/logs", exist_ok=True)

rule all:
	input:
		expand("{datadir}{fasta}.amb", datadir = config["datadir"], fasta = config["fasta"]),
		expand("{datadir}{fasta}.ann", datadir = config["datadir"], fasta = config["fasta"]),
		expand("{datadir}{fasta}.bwt", datadir = config["datadir"], fasta = config["fasta"]),
		expand("{datadir}{fasta}.pac", datadir = config["datadir"], fasta = config["fasta"]),
		expand("{datadir}{fasta}.sa", datadir = config["datadir"], fasta = config["fasta"]),
		expand("{datadir}{fasta}.fai", datadir = config["datadir"], fasta = config["fasta"])

rule bwa_index:
	message: "~-~ Indexing fasta file with bwa... ~-~"
	benchmark: "./output/benchmarks/bwa_index_file={fasta}"
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
	benchmark: "./output/benchmarks/samtools_index_file={fasta}"
	input:
		fasta = "{fasta}"
	output:
		"{fasta}.fai"
	shell:
		"samtools faidx {input.fasta}"

rule bwa_align:
	message: "~-~ Aligning trimmed reads... ~-~"
	benchmark: "./output/benchmarks"
	output:
	params:
		split = config["split"]
		fasta = config["fasta"]
		datadir = config["datadir"]
	run:
		part_decision = ""
		if(int("{params.split}") > 1):
			shell("seqkit split2 -p {params.split} --out-dir output/bwa_alignment/")
		for i in range(1:int("{params.split}")):
			if(int("{params.split}") > 1):
				part_decision = ".part_00" + str(i)
			shell("bwa aln -t {threads} {params.datadir}/{params.fasta} hic_r1" + part_decision + ".fastq.gz > sample_R1" + part_decision + ".sai")
			shell("bwa aln -t {threads} {params.datadir}/{params.fasta} hic_r2" + part_decision + ".fastq.gz > sample_R2" + part_decision + ".sai")
