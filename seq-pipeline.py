#!/usr/bin/env python3

import os
import click
from subprocess import PIPE, Popen, call
from pathlib import Path
from os.path import isfile, isdir
from pprint import pprint
from time import sleep



def ispaired(srr_accession):
	if isfile(f"{srr_accession}/{srr_accession}.fastq"):
		return(False)
	elif isfile(f"{srr_accession}/{srr_accession}_1.fastq") and isfile(f"{srr_accession}/{srr_accession}_2.fastq"):
		return(True)
	else:
		sys.exit("unknown if file is paired or single end!")

def prefetch(srr_accession):
	print()
	print('prefetch...')

	call = ['prefetch', srr_accession, '-v']
	print("command:", " ".join(map(str,call)))
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
		p = Popen(call, stderr=PIPE, stdout=PIPE, encoding='utf-8')

		for line in p.stderr:
			print("[fastqc]", line.strip())

		for line in p.stdout:
			print("[fastqc]", line.strip())
		p.wait()





def cutadapt(srr_accession):
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


	if not ispaired(srr_accession):
		call = ['cutadapt', '-q', '25,25', '-m', '30']
		for adapter in adapters:
			call += ['-b', adapter]
		call += ['-o', f'{srr_accession}/{srr_accession}.trimmed.fastq', f'{srr_accession}/{srr_accession}.fastq']

	else:
		call = ['cutadapt', '-q', '25,25', '-m', '30']
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
	p = Popen(call, stderr=PIPE, encoding='utf-8')

	for line in p.stderr:
		print("[kallisto_quant]", line.strip())
	p.wait()


@click.group()
def cli():
	'''Tool chain for performing RNA-seq quantification with no leftovers.'''
	pass




@cli.command()
@click.argument("srr_accession")

def download(srr_accession):
	'''Download and QC of libraries.'''

	prefetch(srr_accession)
	fasterq_dump(srr_accession)
	fastqc(srr_accession)








@cli.command()
@click.argument("srr_accession")

def trim(srr_accession):
	'''Trimming and QC of libraries.'''
	cutadapt(srr_accession)
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
	kallisto_index(cdna_file)
	kallisto_quant(srr_accession, cdna_file, bootstraps, length, sd)








@cli.command()
@click.argument("srr_accession")
@click.argument("cdna_file")

def chip(srr_accession, genome_file):
	'''Analysis of chip data using macs2'''

	print('not yet implemented')







@cli.command()
@click.argument("srr_accession")

def cleanup(srr_accession):
	'''Cleaning up heavy files'''


	files = [
		f'{srr_accession}/{srr_accession}.fastq',
		f'{srr_accession}/{srr_accession}_1.fastq',
		f'{srr_accession}/{srr_accession}_2.fastq',
		f'{srr_accession}/{srr_accession}.sra'
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






