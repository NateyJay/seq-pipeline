# seq-pipeline
A toolset for storage-safe RNA-seq analyses


### Why?

This is meant to help with analyses to save hard-drive space.

### What does it do?

Basically, this is just a set of wrappers for common RNA-seq tools. Here, I standardize some pipelines and give an appropriate download/analysis location.

Finally, the real key of this tool is the '''cleanup''' function. This will remove the large files downloaded from the NCBI, saving valuable disk space while retaining the analyses.

All of the analyses are set to run in single-core mode, with the idea that you would clone this process across many libraries in parallel.



### Requirements and installation

This tool is just a script, which should be run with python3. This can be run in the git folder or made executable and put in your path:

```
## from git folder
python3 ./seq-pipeline.py

## making executable
chmod +x ./seq-pipeline.py
./seq-pipeline.py
```

As a wrapper, this just uses other tools. It is not sophisticated enough to confirm they function, so you must download, install, and have command-line executable from your PATH the following tools:

```
## sra-tools from NCBI
prefetch
fasterq-dump

## fastqc (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
fastqc

## cutadapt (https://cutadapt.readthedocs.io/en/stable/)
cutadapt 

## kallisto (https://pachterlab.github.io/kallisto/)
kallisto quant
kallisto index
```

The tool uses several python packages you may not have, which can be installed through the python package manager pip, i.e.
```
pip install click
```


### How to use it?

Several functions are generic, with the endpoint functions (kallisto, chip) being interchangeable

#### -> Quantification
This downloads, trims, performs QC, and quantifies single- or paired-end RNA seq data based on SRR code. The cleanup step removes the NCBI files, including trimmed and untrimmed libraries.

```
./seq-pipeline.py download SRR6369351
./seq-pipeline.py trim SRR6369351
./seq-pipeline.py kallisto SRR6369351 transcriptome.fa
./seq-pipeline.py cleanup SRR6369351
```

#### -> ChIP peak finding

This requires ChIP libraries, which may have replicates and treatment/control designations. This process requires treatments, while control (input) libraries are optional for macs2.

```
## treatment library
seq-pipeline.py download SRR799817
seq-pipeline.py trim SRR799817

## control library
seq-pipeline.py download SRR799819
seq-pipeline.py trim SRR799819

## chip analysis
seq-pipeline.py chip -t SRR799817 -c SRR799819 ./Solyc.GCF_000188115.4_SL3.0.fa -o chip_output

## cleanup
seq-pipeline.py cleanup SRR799817
seq-pipeline.py cleanup SRR799819
```


### Problems

1) **Not very smart**. You may have problems or errors in steps. If so, the script will likely produce a bad file which the next step will have an error with. You will need to look at the output of the run to identify where these broken files are and remedy them. Each step just looks for its expected input, and does little to confirm it is what it actually should be.
2) **Non NCBI accessions**. If you need to use an accession that is not in the NCBI (e.g. chinese accessions CNNxxxxx), you must download them manually. The software should function for trim and subsequent steps if you nest save them in a folder named by their own accession, in this format: ```{accession}/{accession}.fastq```.
3) **Large memory buffer**. While the cleanup steps ensure you won't keep large files around after the analysis, the software is not very respectful of hard drive space *while* analyzing. It produces multiple somewhat redundant bams, fastqs... You should be aware that you will need significant free space to run these analyses.



