#!/usr/bin/env python3

import os
import click
from subprocess import PIPE, Popen, call
from pathlib import Path
from os.path import isfile, isdir
from pprint import pprint
from time import sleep
import sys


###########################
### Tools
###########################


class Logger(object):
	def __init__(self, file_name):
		self.terminal = sys.stdout
		self.file_name = file_name

		if isfile(file_name):
			self.log = open(file_name, "a")
		else:
			self.log = open(file_name, "w")

		# with open(file_name, "w") as outf:
		# 	outf.write("")

	def clear_ansi(self, message):
		return(message.replace("\033[1m", "").replace("\033[0m",""))

	def write(self, message):
		self.terminal.write(message)
		# with open(self.file_name, 'a') as outf:
		# 	outf.write(message)  
		self.log.write(self.clear_ansi(message))

	def flush(self):
		self.terminal.flush()
		self.log.flush()


def ispaired(srr_accession):
	if isfile(f"{srr_accession}/{srr_accession}.fastq"):
		return(False)
	elif isfile(f"{srr_accession}/{srr_accession}_1.fastq") and isfile(f"{srr_accession}/{srr_accession}_2.fastq"):
		return(True)
	else:
		sys.exit("unknown if file is paired or single end!")



###########################
### Functions
###########################

def prefetch(srr_accession):
	print()
	print('prefetch...')

	call = ['prefetch', srr_accession, '-v']
	print("command:", " ".join(map(str,call)))
	print()
	p = Popen(call, stderr=PIPE, encoding='utf-8')

	for line in p.stderr:
		print("[prefetch]", line.strip())
	p.wait()


def fasterq_dump(srr_accession):
	print()
	print('fasterq-dump...')

	fastq_file = f"{srr_accession}/{srr_accession}.fastq"
	if isfile(fastq_file):
		print("[fasterq-dump]", f"fastq file found ({srr_accession}/{fastq_file}), skipping...")
		return


	call = ['fasterq-dump', srr_accession, '-v', '-v', '-O', srr_accession]
	print("command:", " ".join(map(str,call)))
	print()
	p = Popen(call, stderr=PIPE, encoding='utf-8')

	for line in p.stderr:
		print("[fasterq-dump]", line.strip())
	p.wait()


def fastqc(srr_accession, suffix = ".fastq"):

	files = []
	if ispaired(srr_accession):
		files += [f'{srr_accession}/{srr_accession}_1{suffix}']
		files += [f'{srr_accession}/{srr_accession}_2{suffix}']
	else:
		files += [f'{srr_accession}/{srr_accession}{suffix}']

	for file in files:

		print()
		print('fastqc...')

		call = ['fastqc', file, '--outdir', srr_accession]
		print("command:", " ".join(map(str,call)))
		print()
		p = Popen(call, stderr=PIPE, stdout=PIPE, encoding='utf-8')

		for line in p.stderr:
			print("[fastqc]", line.strip())

		for line in p.stdout:
			print("[fastqc]", line.strip())
		p.wait()





def cutadapt(srr_accession, cores):
	print()
	print('cutadapt...')

	adapter = 'TGGAATTC'

	adapters = [
		'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
		'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA',
		'TACACTCTTTCCCTACACGACGCTCTTCCGATCT',
		'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT',
		'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT',
		'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT',
		'AATGATACGGCGACCACCGAGATCTACACTCTTTCCTACACGACGCTCTTCCGATCT',
		'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTT'
	]


	call = ['cutadapt', '-q', '25,25', '-m', '30', '-j', str(cores)] ## multicore doesn't work on my system....
	if not ispaired(srr_accession):
		for adapter in adapters:
			call += ['-b', adapter]
		call += ['-o', f'{srr_accession}/{srr_accession}.trimmed.fastq', f'{srr_accession}/{srr_accession}.fastq']

	else:
		for adapter in adapters:
			call += ['-b', adapter]
		for adapter in adapters:
			call += ['-B', adapter]
		call += [
			'-o', f'{srr_accession}/{srr_accession}_1.trimmed.fastq', 
			'-p', f'{srr_accession}/{srr_accession}_2.trimmed.fastq', 
			f'{srr_accession}/{srr_accession}_1.fastq', 
			f'{srr_accession}/{srr_accession}_2.fastq'
		]

	print("command:", " ".join(map(str,call)))
	print()


	p = Popen(call, stdout=PIPE, encoding='utf-8')


	for line in p.stdout:
		print("[cutadapt]", line.strip())
	p.wait()





def kallisto_index(cdna_file):
	print()
	print('kallisto index...')

	index_file = f'{cdna_file}.kaidx'
	if isfile(index_file):
		print("[kallisto_index]", f"index file found ({index_file}), skipping...")
		return

	call = ['kallisto', 'index','-i', index_file, cdna_file]
	print("command:", " ".join(map(str,call)))
	print()
	p = Popen(call, stderr=PIPE, encoding='utf-8')

	for line in p.stderr:
		print("[kallisto_index]", line.strip())
	p.wait()




def kallisto_quant(srr_accession, cdna_file, bootstraps, length, sd):


	print()
	print('kallisto quant...')

	index_file = f'{cdna_file}.kaidx'
	output_file = f'{srr_accession}/{srr_accession}_kallisto'


	call = ['kallisto', 'quant', 
		'-i', index_file, 
		'-o', output_file, 
		'-b', str(bootstraps)]

	if ispaired(srr_accession):
		call += [
			f"{srr_accession}/{srr_accession}_1.trimmed.fastq",
			f"{srr_accession}/{srr_accession}_2.trimmed.fastq"
		]
	else:
		call += [
			'--single',
			'-l', str(length),
			'-s', str(sd),
			f"{srr_accession}/{srr_accession}.trimmed.fastq"
		]
	print("command:"," ".join(map(str,call)))
	print()
	p = Popen(call, stderr=PIPE, encoding='utf-8')

	for line in p.stderr:
		print("[kallisto_quant]", line.strip())
	p.wait()


def bowtie2_build(genome_file):
	print()
	print('bowtie2-build...')

	call = ['bowtie2-build', genome_file, genome_file]
	print("command:", " ".join(map(str,call)))
	print()
	p = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')

	for line in p.stdout:
		print("[bowtie2-build]", line.strip())
	for line in p.stderr:
		print("[bowtie2-build]", line.strip())

	p.wait()


def bowtie2_align(srr_accession, genome_file, cores):
	expected_index = genome_file + ".1.bt2"

	if not isfile(expected_index):
		print(f"Index not found ({expected_index}), performing bowtie2-build...")
		bowtie2_build(genome_file)


	print()
	print('bowtie2-align...')

	call = ["bowtie2", '-p', str(cores), '-x', genome_file]

	if not ispaired(srr_accession):
		call += ["-U", f"{srr_accession}/{srr_accession}.trimmed.fastq"]
	else:
		call += ["-1", f"{srr_accession}/{srr_accession}_1.trimmed.fastq"]
		call += ["-2", f"{srr_accession}/{srr_accession}_2.trimmed.fastq"]


	bam_file = f"{srr_accession}/{srr_accession}.bam"

	if isfile(bam_file):
		print(f'bam file ({bam_file}) found, skipping alignment')
		return

	print("command:", " ".join(map(str,call)), f"| samtools view -b -h > {bam_file}")
	print()

	with open(bam_file, 'wb') as outf:
		p1 = Popen(call, stdout=PIPE, stderr=PIPE, encoding='utf-8')
		p2 = Popen(['samtools','view','-b', '-h'], stdin=p1.stdout, stdout=outf)
		p2.wait()







def samtools_sort(srr_accession):

	print()
	print('samtools sort...')

	bam_file = f"{srr_accession}/{srr_accession}.bam"
	sorted_bam_file = f"{srr_accession}/{srr_accession}.sorted.bam"

	if isfile(sorted_bam_file):
		print(f'sorted bam file ({sorted_bam_file}) found, skipping sort')
		return


	call = ['samtools','sort',bam_file]
	print("command:", " ".join(map(str,call)), f"> {sorted_bam_file}")
	print()

	with open(sorted_bam_file, 'wb') as outf:
		p = Popen(call, stdout=outf, stderr=PIPE, encoding='utf-8')
		for line in p.stderr:
			print("[samtools_sort]", line.strip())
		p.wait()


def picard_deduplicate(srr_accession):

	bam_file = f"{srr_accession}/{srr_accession}.sorted.bam"
	dedup_file = f"{srr_accession}/{srr_accession}.dedup.bam"
	metrics_file = f"{srr_accession}/{srr_accession}.dedup_metrics.txt"

	print()
	print('picard deduplicate...')


	if isfile(dedup_file):
		print(f'deduplicated bam file ({dedup_file}) found, skipping deduplication...')
		return

	call = ['picard', 'MarkDuplicates', '-I', bam_file, '-O', dedup_file, '-M', metrics_file]
	print("command:", " ".join(map(str,call)))
	print()

	p = Popen(call, stderr=PIPE, stdout=PIPE, encoding='utf-8')

	for line in p.stderr:
		print("[picard_MarkDuplicates]", line.strip())
	for line in p.stdout:
		print("[picard_MarkDuplicates]", line.strip())
	p.wait()

def quality_filter(srr_accession):

	print()
	print('quality filter...')

	dedup_file = f"{srr_accession}/{srr_accession}.dedup.bam"
	filtered_file = f"{srr_accession}/{srr_accession}.filtered.bam"

	if isfile(filtered_file):
		print(f'quality filtered bam file ({filtered_file}) found, skipping filter...')
		return

	call = ['samtools', 'view', '-b', '-q', '10', dedup_file]
	print("command:", " ".join(map(str,call)))
	print()


	with open(filtered_file, 'wb') as outf:
		p = Popen(call, stdout=outf)
		p.wait()

def get_genome_size(output_directory, genome_file):


	print()
	print('measuring genome size...')

	genome_size_file = f"{output_directory}/genome_size.txt"

	if isfile(genome_size_file):
		print(f'genome size file ({genome_size_file}) found, skipping measurement...')
		with open(genome_size_file, 'r') as f:
			return(f.readline().strip().split()[0])


	call = ['wc', '-c', genome_file]
	print("command:", " ".join(map(str,call)))
	print()

	with open(genome_size_file, 'w') as outf:

		p = Popen(call, stdout=PIPE, encoding='utf-8')

		for line in p.stdout:
			print("[measure_genome]", line.strip())
			print(line.strip(), file=outf)
			return(line.strip().split()[0])


def macs2(treatments, controls, genome_file, genome_size, output_directory):

	print()
	print('getting peaks with macs2...')
	call = ["macs2", 'callpeak']

	call.append('-t')
	for t in treatments:
		call.append(f"{t}/{t}.dedup.bam")

	if len(controls) > 0:
		call.append('-c')
		for c in controls:
			call.append(f"{c}/{c}.dedup.bam")


	if ispaired(treatments[0]):
		call += ['-f', 'BAMPE']

	call += ['-g', genome_size]
	call += ['-n', 'macs_analysis']
	call += ['--outdir', output_directory]

	print("command:", " ".join(map(str,call)))
	print()

	p = Popen(call, stderr=PIPE, encoding='utf-8')
	for line in p.stderr:
			print("[macs2_callpeak]", line.strip())

	p.wait()





	# 	  "+file+" -c ./5-qual/qFiltered_SID01.SRR16905154_dedup.bam -n "+id2+".input --outdir 6-macs"]

	# command = "macs2 callpeak -t "+file+" -c ./5-qual/qFiltered_SID01.SRR16905154_dedup.bam -f BAMPE -n "+id2+".input.bampe --outdir 6-macs"

	# command = "macs2 callpeak -t "+file+" -f BAMPE -n "+id2+".bampe --outdir 6-macs"


###########################
### CLI commands
###########################




@click.group()
def cli():
	'''Tool chain for performing RNA-seq quantification with no leftovers.'''
	pass




@cli.command()
@click.argument("srr_accession")

def download(srr_accession):
	'''Download and QC of libraries.'''



	Path(f"{srr_accession}").mkdir(parents=True, exist_ok=True)

	log_file = f"{srr_accession}/{srr_accession}.log.txt"
	with open(log_file, 'w') as outf:
		outf.write("")
	sys.stdout = Logger(log_file)


	prefetch(srr_accession)
	fasterq_dump(srr_accession)
	fastqc(srr_accession)








@cli.command()
@click.argument("srr_accession")
@click.option('-c','--cores',
	default=1,
	help='Number of cores used for cutadapt.')

def trim(srr_accession, cores):
	'''Trimming and QC of libraries.'''


	log_file = f"{srr_accession}/{srr_accession}.log.txt"
	sys.stdout = Logger(log_file)

	cutadapt(srr_accession, cores=cores)
	fastqc(srr_accession, ".trimmed.fastq")







@cli.command()
@click.argument("srr_accession")
@click.argument("cdna_file")
@click.option('-b', '--bootstraps',
	default=100,
	help='Bootstrap count for kallisto quant.')

@click.option('--length',
	default=400,
	help='Average expected fragment length (used only in single-end mode).')

@click.option('--sd',
	default=100,
	help='Fragment length standard deviation (used only in single-end mode).')


def kallisto(srr_accession, cdna_file, bootstraps, length, sd):
	'''Quantification of RNAseq reads using kallisto.'''

	log_file = f"{srr_accession}/{srr_accession}.log.txt"
	sys.stdout = Logger(log_file)

	kallisto_index(cdna_file)
	kallisto_quant(srr_accession, cdna_file, bootstraps, length, sd)








@cli.command()
@click.argument("genome_file")
@click.option('--cores',
	default=1,
	help='Number of cores used for bowtie2-align.')

@click.option('-o','--output_directory',
	required=True,
	help='Location to save chip experiment data')

@click.option('-t','--treatments',
	required=True,
	multiple=True,
	help='Treatment library SRR accessions, add multiple by using multiple flags')

@click.option('-c','--controls',
	required=False,
	multiple=True,
	help='Control library SRR accessions, add multiple by using multiple flags')

def chip(genome_file, cores, treatments, controls, output_directory):
	'''Analysis of chip data using macs2'''

	for srr_accession in treatments + controls:

		print()
		print(f"PROCESSING {srr_accession}")

		log_file = f"{srr_accession}/{srr_accession}.log.txt"
		sys.stdout = Logger(log_file)


		bowtie2_align(srr_accession, genome_file, cores)
		samtools_sort(srr_accession)
		picard_deduplicate(srr_accession)
		quality_filter(srr_accession)

	Path(f"{output_directory}").mkdir(parents=True, exist_ok=True)

	genome_size = get_genome_size(output_directory, genome_file)


	macs2(treatments, controls, genome_file, genome_size, output_directory)





@cli.command()
@click.argument("srr_accession")

def cleanup(srr_accession):
	'''Cleaning up heavy files'''

	log_file = f"{srr_accession}/{srr_accession}.log.txt"
	sys.stdout = Logger(log_file)

	files = [
		f'{srr_accession}/{srr_accession}.fastq',
		f'{srr_accession}/{srr_accession}_1.fastq',
		f'{srr_accession}/{srr_accession}_2.fastq',
		f'{srr_accession}/{srr_accession}.sra',
		f'{srr_accession}/{srr_accession}.bam',
		f'{srr_accession}/{srr_accession}.sorted.bam',
		f'{srr_accession}/{srr_accession}.filtered.bam',
		f'{srr_accession}/{srr_accession}.dedup.bam'
	]

	print(files)

	files = [f for f in files if isfile(f)]

	print("cleanup...")

	if len(files) == 0:
		print(" -> no files to remove!")
		return

	print()
	print("Removing:")
	for file in files:
		print(" ", file)

	counter_line = "[     ]"
	print(counter_line, flush=False, end='\r')

	for wait in range(6):
		counter_line = "[" + wait*"." + (5-wait)*" " + "]"
		print(counter_line, end='\r')


		# print(".", end='', flush=True)
		sleep(1)

	print()

	for file in files:
		os.remove(file)


	print('done!')





if __name__ == '__main__':
    cli()






