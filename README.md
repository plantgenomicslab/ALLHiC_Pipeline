# ALLHiC_Pipeline
Pipeline for Straightforward Usage of ALLHic

### PREREQUISITES
1. ALLHiC Software (https://github.com/tangerzhang/ALLHiC)
2. getFalen.pl (obtained from tangerzhang's my_script repository on GitHub) (https://github.com/tangerzhang/my_script/blob/master/getFaLen.pl)
3. samtools
4. seqkit
5. trim_galore
6. snakemake
7. bwa
8. Anaconda (optional, necessary for simple environment)
*Note: 3-7 included in env.yml if installed using Anaconda.

### Environment
The pipeline may be activated wherever the user wishes, but the snakefile needs to have the appropriate subfiles and subdirectories accompanying it in the same directory.
1. The ALLHiC software file must be downloaded whole (such as with git clone) and placed into the bin subdirectory.
2. getFalen.pl must be downloaded individually and placed in the bin subdirectory.
3. The rest of the software may be downloaded individually or through Anaconda.
	 conda env create -f env.yml
	 *Note: Default environment name is "allhic." Name can be modified in env.yml file, simply change the text after "name:"
   
### Usage
(A) Without Cluster
1. Modify the config.json file to meet your needs and environment.
	* datadir: Path to directory that contains all your reference files, such as fasta and fastq files.
	* snakedir: Path to the directory where the ALLHiC pipeline snakefile is contained (should also contain the accompanied bin subdirectory).
	* fasta: Name of the fasta file to use.
	* groups:
	* enzyme_sites:
	* split:
	* threads:
2. Run snakemake
	* Example Command: snakemake -s {path to snakefile} --cores {# of cores}

(B) With Cluster
1. Same as step 1 from (A).
2. Modify the cluster.json file to meet your needs.
3. Modify the run.sh file to meet your needs.
4. Run run.sh
	* sh run.sh
