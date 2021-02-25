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
PAIRS_R1, PAIRS_R2 = cgz.getGZPairs(cgz.rename_gz(cgz.collectGZ(REF_DIR + "/")))
PAIR_NAMES_1 = []
PAIR_NAMES_2 = []
FQ_DIR = ""
for i in PAIRS_R1[0].split("/")[0:-1]:
	FQ_DIR = FQ_DIR +i + "/"
for i in range(0,len(PAIRS_R1)):
	title = PAIRS_R1[i].split("/")[-1]
	title = title.replace(".fastq.gz","")
	title = title.replace(".fq.gz","")
	part1,part2 = cgz.get_parts(title)
	PAIR_NAMES_1.append(part1)
	PAIR_NAMES_2.append(part2)
PARTS = genSS.partSuffixes(str(config["split"]))
GROUP_LIST = []
for i in range(1, int(config["groups"])+1):
	GROUP_LIST.append(i)

rule all:
	input:
#		expand("{ref_dir}/{fasta}.amb", ref_dir = REF_DIR, fasta = config["fasta"]),
#		expand("{ref_dir}/{fasta}.ann", ref_dir = REF_DIR, fasta = config["fasta"]),
#		expand("{ref_dir}/{fasta}.bwt", ref_dir = REF_DIR, fasta = config["fasta"]),
#		expand("{ref_dir}/{fasta}.pac", ref_dir = REF_DIR, fasta = config["fasta"]),
#		expand("{ref_dir}/{fasta}.sa", ref_dir = REF_DIR, fasta = config["fasta"]),
#		expand("{ref_dir}/{fasta}.fai", ref_dir = REF_DIR, fasta = config["fasta"]),
#		expand("output/ALLHiC/sorted_collection.clean.counts_{enzyme}.{groups}g{group_section}.tour", enzyme = config["enzyme_sites"], groups = config["groups"], group_section = GROUP_LIST)
#		"output/graphs/500K_all_chrs.pdf"

rule trim_galore:
	message: "~-~ Trimming fastq files... ~-~"
	benchmark: "output/benchmarks/trim_galore.prefix={prefix}_r1-2_{suffix}.bm"
	#log: "output/cluster/logs/trim_galore_file={fastq}.log"
	threads: config["threads"]
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
	threads: config["threads"]
	benchmark: "output/benchmarks/splitting_file={prefix}.part_{parts}.bm"
	#log: "output/cluster/logs/splitting_file={pair}_{strand}_{parts}.log"
	input:
		"output/trimmed_reads/{prefix}.fq.gz",
		expand("ref/{fasta}.amb", fasta = config["fasta"]),
		expand("ref/{fasta}.fai", fasta = config["fasta"])
	output:
		"output/trimmed_reads/{prefix}.part_{parts}.fq.gz"
	params:
		split = config["split"]
	shell: "seqkit split2 -j {threads} -p {params.split} {input[0]} --out-dir output/trimmed_reads"

rule bwa_align:
	message: "~-~ Aligning trimmed reads... ~-~"
	threads: config["threads"]
	benchmark: "output/benchmarks/bwa_algn_file={prefix}.part_{parts}.bm"
	#log: "output/cluster/logs/bwa_algn_file={pair}{strand}_{parts}.log"
	input:
		"output/trimmed_reads/{prefix}.part_{parts}.fq.gz"
	output:
		"output/bwa_algn/{prefix}.part_{parts}.sai"
	params:
		refdir = config["datadir"],
		fasta = config["fasta"]
	shell: "mkdir bwa_algn/ || true;\
			bwa aln -t {threads} {params.refdir}{params.fasta} {input[0]} > {output[0]}"

rule bwa_sampe:
	message: "~-~ Bwa sampe... ~-~"
	benchmark: "output/benchmarks/bwa_sampe_pre={prefix}r1-2{suffix}.part_{parts}.bm"
	#log: "output/cluster/logs/bwa_sampe_pre={pair}_{parts}.log"
	input:
		"output/bwa_algn/{prefix}r1{suffix}_1.part_{parts}.sai",
		"output/bwa_algn/{prefix}r2{suffix}_2.part_{parts}.sai",
		"output/trimmed_reads/{prefix}r1{suffix}_val_1.part_{parts}.fq.gz",
		"output/trimmed_reads/{prefix}r2{suffix}_val_2.part_{parts}.fq.gz"
	output:
		"output/bwa_algn/{prefix}r1-2{suffix}.part_{parts}.sam"
	params:
		fasta = config["fasta"],
		datadir = config["datadir"],
		trims = config["trims"]
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
			mv ../bwa_algn/{wildcards.prefix}.REduced* . || true;\
			for file in {wildcards.prefix}.*.REduced.paired_only.bam; do sambamba sort -t {threads} $file -o sorted.$file; done;\
			sambamba merge -t {threads} {wildcards.prefix}.merged.bam sorted.{wildcards.prefix}.*.REduced.paired_only.bam;\
			sambamba sort -t {threads} -o {wildcards.prefix}.merged_sorted.bam {wildcards.prefix}.merged.bam;\
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

rule merge_clean:
	message: "~-~ Merging Cleaned BAM Files... ~-~"
	benchmark: "output/benchmarks/merged_clean.bm"
	threads: config["threads"]
	output:
		"output/ALLHiC/sorted_collection.clean.bam"
	shell: "cd output/ALLHiC/;\
			for file in *.clean.bam; do sambamba sort -t {threads} $file -o sorted.$file; done;\
			sambamba merge -t {threads} collection.clean.bam sorted.*.clean.bam;\
			sambamba sort -t {threads} -o sorted_collection.clean.beam collection.clean.bam;\
			cd ../../"

rule allhic_corrector:
	message: "~-~ Correcting ALLHiC... ~-~"
	benchmark: "output/benchmarks/allhic_corrector.bm"
	threads: config["threads"]
	input:
		"output/ALLHiC/sorted_collection.clean.bam"
	output:
		"output/ALLHiC/seq.corrected.fasta"
	params:
		snakedir = SNAKEDIR,
		refdir = REFDIR,
		fasta = config["fasta"]
	shell: "python {params.snakedir}/bin/ALLHiC/bin/ALLHiC_corrector -m {input[0]} -r {params.refdir}/{params.fasta} -o {output[0]}"

rule partition:
	message: "~-~ Partitioning... ~-~"
	benchmark: "output/benchmarks/ALLHiC_partition.enzyme={enzyme}.groups={groups}.bm"
	input:
		"output/ALLHiC/sorted_collection.clean.bam",
		"output/ALLHiC/seq.corrected.fasta"
	output:
		"output/ALLHiC/sorted_collection.clean.counts_{enzyme}.{groups}g1.txt"
		#"output/ALLHiC/p1={part1}_strand={diff}1-2_p2={part2}.clean.clm"
	params:
		snakedir = SNAKE_DIR,
		cwd = CWD,
		refdir = REF_DIR,
		groups = config["groups"],
		ezs = config["enzyme_sites"]
	shell: "cd {params.snakedir}/bin/ALLHiC/bin/;\
			ALLHiC_partition -b {params.cwd}/{input[0]} -r {params.cwd}/output/ALLHiC/seq.corrected.fasta -e {params.ezs} -k {params.groups};\
			cd {params.cwd}"

rule optimize:
	message: "~-~ Optimizing... ~-~"
	benchmark: "output/benchmarks/allhic_optimizing.enzyme_{enzyme}.groups_{groups}.part_{group_section}.bm"
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
	input:
		"output/ALLHiC/seq.corrected.fasta"
	output:
		"output/ALLHiC/groups.agp",
		"output/ALLHiC/groups.asm.fasta"
	params:
		fasta = config["fasta"],
		refdir = REF_DIR,
		snakedir = SNAKE_DIR,
		cwd = CWD
	shell: "cd output/ALLHiC/;\
			perl {params.snakedir}/bin/ALLHiC/bin/ALLHiC_build {input[0]};\
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
