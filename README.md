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
> Note: 3-7 included in env.yml if installed using Anaconda.
___
### Environment
The pipeline may be activated wherever the user wishes, but the snakefile needs to have the appropriate subfiles and subdirectories accompanying it in the same directory.
1. The ALLHiC software file must be downloaded whole (such as with git clone) and placed into the bin subdirectory.  
`git clone https://github.com/tangerzhang/ALLHiC.git`
2. getFalen.pl must be downloaded individually and placed in the bin subdirectory.  
`wget https://raw.githubusercontent.com/tangerzhang/my_script/master/getFaLen.pl`
3. The rest of the software may be downloaded individually or through Anaconda.  
`conda env create -f env/env.yml`
> Note: Default environment name is "allhic." Name can be modified in env/env.yml file, simply change the text after "name:"
___   
### Usage
(A) Without Cluster
1. Modify the config.json file to meet your needs and environment.
2. Run snakemake  
`snakemake -s {path to snakefile} --cores {# of cores}`

(B) With Cluster
1. Same as step 1 from (A).
2. Modify the cluster.json file to meet your needs.
3. Modify the run.sh file to meet your needs.
4. Run run.sh  
`sh run.sh`
___
### Options Description  
There are a few files that the pipeline works with that need to be modified to fit your desired requirements. There are three files that might be need to be modified: config.json, cluster.json, and run.sh.  

#### config.json Options  
The config.json file contains the parameters that the snakemake file works with. These parameters are necessary for the snakefile to work under any condition. The following is a list of all parameters in the config.json file:
* **datadir**: The absolute directory that contains your fasta file and all of your fastq files.
* **snakedir**: The absolute directory where the whole pipeline respository is stored.
* **fasta**: The name of the fasta file (no directory attached).
* **groups**: The number of chromosomes to separate into.
* **enzyme_sites**: Sequence to detect (HindIII, MboI, Arime).
* **split**: How many times to split the trimmed fastq files (helps with efficiency).
* **threads**: The number of threads to use for each job.  

#### cluster.json Options  
If you plan on running snakemake through a cluster, you need to specify the required parameters needed to submit the job. The cluster.json file contains arrays of objects for each rule in the snakemake file, which allows for enhanced resource management. Each array is in the following form,  
```
"rule name":{  
    "var1": "value",  
    "var2": value,  
    . . .  
}
```
The rule name section should be left unmodified, along with the left-hand object variable names. The values of the objects, on the right-hand of the object variable names, can be modified to your liking. The following is a list of all the variable names in the cluster.json and their description,  
* **account**: The cluster account name that the program will run under (string).  
* **memory**: The amount of memory (e.g. RAM) to allocate to your jobs (string, e.g. "10g").  
* **name**: The name to submit your jobs under. A default one has been provided, but it can be changed (string).  
* **ncpus**: The number of cpus to allow for your jobs (integer).  
* **nodes**: The number of nodes to allocate for your jobs (integer).  
* **partition**: The subset cluster to use in the account (string).  
* **time**: The max time to allow your submissions to run (string).  
* **output**: The name of the file to output the stdout of each individual submitted job. A default one has been provided, but it can be changed.  
* **error**: The name of the file to output the stderr of each individual submitted job. A default one has been provided, but it can be changed.  

For additional cluster options, refer to the run.sh section.  

#### run.sh Options  
The run.sh file is what calls the necessary command(s) to submit the pipeline into the cluster. If you plan on using a cluster for the pipeline, there are a few parameters to take note of:  
* **-s**: The path of the SnakeFile
* **--cluster-config**: The path of the cluster.json file.
* **--cluster**: The parameters to submit to the cluster.
* **--max-jobs-per-second**: The max number of jobs to perform on at one second.
* **--max-status-checks-per-second**: The max number of status checks to perform at one second.
* **--jobs**: The number of jobs to have submitted at a time.
* **--latency-wait**: Time allowed for snakemake to wait to check if the appropriate output file(s) was/were created.
* **--cluster-status**: The directory to a program that checks the status of the cluster (optional, can be removed).  

If you require any additional cluster parameters to be submitted, you will have to enter them manually within the "--cluster" parameter string.  
