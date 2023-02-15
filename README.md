# seq-pipeline
A toolset for storage-safe RNA-seq analyses


### Why?

This is meant to help with analyses to save hard-drive space.

### What does it do?

Basically, this is just a set of wrappers for common RNA-seq tools. Here, I standardize some pipelines and give an appropriate download/analysis location.

Finally, the real key of this tool is the '''cleanup''' function. This will remove the large files downloaded from the NCBI, saving valuable disk space while retaining the analyses.



### Requirements and installation

This tool is just a script, which should be run with python3. This can be run in the git folder or made executable and put in your path:

```
## from git folder
python3 ./seq-pipeline.py

## making executable
chmod +x ./seq-pipeline.py
./seq-pipeline.py
```


### How to use it?

Currently, the tool only supports one pipeline. This downloads, trims, performs QC, and quantifies single- or paired-end RNA seq data based on SRR code.


```
./seq-pipeline.py download SRR6369351
./seq-pipeline.py trim SRR6369351
./seq-pipeline.py kallisto SRR6369351 transcriptome.fa
./seq-pipeline.py cleanup SRR6369351
```