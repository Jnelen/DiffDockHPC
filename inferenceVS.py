# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 10:49:33 2023

@author: Jochem Nelen (jnelen@ucam.edu)
"""

import datetime
import os
import shutil
import subprocess
import sys
import time

import itertools
import glob

import argparse
from argparse import ArgumentParser

parser = ArgumentParser()
  
parser.add_argument('--protein_path', '-r', '-p', required=True, type=str, default='', help='Path to the protein/receptor .pdb file')
parser.add_argument('--ligand', '-l', required=True, type=str, default='', help='The path to the directory of (separate) mol2/sdf ligand files')
parser.add_argument('--out_dir', '-out', '-o', required=True,type=str, default='', help='Directory where the output structures will be saved to')
parser.add_argument('--jobs', '-j', required=True, type=int, default=1, help='Number of jobs to use')
parser.add_argument('--time', '-t', '-tj', required=False, default="", help='Amount of time each job can run')
parser.add_argument('--queue', '-qu', type=str, default="", help='On which node to launch the jobs. The default value is the default queue for the user. Might need to be specified if there is no default queue configured')
parser.add_argument('--mem', '-m', type=str, default="4G", help='How much memory to use for each job. The default value is `4GB')
parser.add_argument('--gpu', '-gpu', '-GPU', '--GPU', action="store_true", default=False, help='Use GPU resources. This will accelerate docking calculations if a compatible GPU is available.')
parser.add_argument('--cores', '-c', type=int, default=None, help='How many cores to use for each job. The default value is 1 when used with the GPU option enabled, otherwise it defaults to 4 cores')
parser.add_argument('--num_outputs', '-n', type=int, default=1, help='How many structures to output per compound. The default value is 1')
parser.add_argument('--remove_hs', action='store_true', default=False, help='Remove the hydrogens in the final output structures')
parser.add_argument('--keep_local_structures', action='store_true', default=False, help='Keeps the local structure when specifying an input with 3D coordinates instead of generating them with RDKit')
parser.add_argument('--keep_cache', action='store_true', default=False, help='Keep the Cache directories after finishing the calculations (Not recommended)')
parser.add_argument('--no_slurm', '-ns', action='store_true', default=False, help='Don\'t use slurm to handle the resources. This will run all samples on 1 GPU. Other Slurm arguments such as the amount memory, time limit, ... will also be ignored')


args = parser.parse_args()

## Check if Singularity image is present and ask to download it
if not os.path.exists("singularity/DiffDockHPC.sif"):
	print("The Singularity image doesn't seem to be present. Please follow the installation instructions and download it from there, or build it manually using the def file")
	answer = input("Would you like to download it automatically? (y/n) ").lower()
	while answer not in ("y", "n", "yes", "no"):
		print("Invalid input. Please enter y(es) or n(o).")
		answer = input("Would you like to download the Singularity image automatically? (y/n) ").lower()
	if answer == "y" or answer == "yes":
		subprocess.run('wget --no-check-certificate -r "https://drive.usercontent.google.com/download?id=1eo_--K6qZoiaphsTK5G4dA8kikZLHG4y&confirm=t" -O singularity/DiffDockHPC.sif ', shell=True)			
	else:
		sys.exit("Please download or build the Singularity image manually and try again")

## Determine the amount of cores if not defined by the user
if args.cores is None:
	if args.gpu:
		args.cores = 1
	else:
		args.cores = 4
		
## If --no_slurm is set, always only use 1 job
if args.no_slurm:
	args.jobs = 1
	
if args.time == "":
	timeArg = ""
else:
	timeArg = f" --time {args.time} "
	
outputPath, outputDirName = os.path.split(args.out_dir)

currentDateNow = datetime.datetime.now()

if not outputPath == "":
	outputPath += "/"
outputDir = outputPath + "_".join(["VS_DD", outputDirName, str(currentDateNow.year), str(currentDateNow.month), str(currentDateNow.day)])

## Check if the output directory already exists, and asks the user what to do if it does
if os.path.isdir(outputDir):
	print(f"The directory {outputDir} already exists. To continue you must delete this directory or choose another outputname.")
	answer = input("Do you want to remove it? (y/n) ").lower()
	while answer not in ("y", "n", "yes", "no"):
		print("Invalid input. Please enter y(es) or n(o).")
		answer = input("Do you want to remove or overwrite it? (y/n) ").lower()
	if answer == "y" or answer == "yes":
		shutil.rmtree(outputDir, ignore_errors=False)			
	else:
		sys.exit()
			
os.mkdir(outputDir)

os.mkdir(outputDir + "/molecules")
os.mkdir(outputDir + "/csvs")
os.mkdir(outputDir + "/jobs_out")
os.mkdir(outputDir + "/jobs")

## Copy protein file to the directory
shutil.copy(args.protein_path, outputDir)

## Check if protein embeddings exist already, if not generate them
proteinName = os.path.basename(args.protein_path)
if not os.path.isdir("data/esm2_output/"):
	os.mkdir("data/esm2_output/")
embeddingPathList = glob.glob(f"data/esm2_output/{proteinName}*.pt")

if len(embeddingPathList) == 0:
	print("No protein embeddings present, generating now..")
	nvArgument = ""
	if args.gpu == True:
		nvArgument = "--nv"
	subprocess.run(f"singularity run {nvArgument} --bind $PWD singularity/DiffDockHPC.sif python proteinEmbedding.py {args.protein_path}", shell=True)
print("Launching jobs now..")	

ligandPaths = glob.glob(f"{args.ligand}/*.sdf") + glob.glob(f"{args.ligand}/*.mol2")

## Code to distribute the query ligands among the amount of jobs 
def split(a, n):
    if n > len(a):
        print("more jobs than files, launching 1 job per file")
        return [a[i:i+1] for i in range(len(a))]
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

ligandPathsSplit = list(split(ligandPaths, args.jobs))

queueArgument = ""
if not args.queue == "":
	queueArgument = " -p " + args.queue

keep_original_struct = ""
if args.keep_local_structures == True:
	 keep_original_struct = "--keep_local_structures"

remove_hs = ""
if args.remove_hs == True:
	 remove_hs = "--remove_hs"
	 
keep_cache = ""
if args.keep_cache == False:
	 keep_cache = "--delete_cache"
	 
for i, jobLigands in enumerate(ligandPathsSplit):
	csvFilePath = f"{outputDir}/csvs/job_csv_{str(i+1)}.csv"
	with open(csvFilePath, 'w') as jobCSV:
		jobCSV.write("protein_path\tligand_description\n")
		for jobLigand in jobLigands:
			jobCSV.write(f"{args.protein_path}\t{jobLigand}\n")

	jobCSV.close()

	if not args.no_slurm:
		## Execute command using singularity and sbatch wrap giving the csv as an input, and passing the input variables as well
		if args.gpu == True:
			jobCMD = f'sbatch --wrap="singularity run --nv --bind $PWD singularity/DiffDockHPC.sif python3 inference.py --protein_ligand_csv {csvFilePath} --samples_per_complex {args.num_outputs} {remove_hs} --out_dir {outputDir}/molecules/ {keep_original_struct} {keep_cache} -c {str(args.cores)}" --mem {args.mem} --output={outputDir}/jobs_out/job_{str(i+1)}_%j.out --gres=gpu:1 --job-name=DiffDockHPC -c {str(args.cores)} {timeArg} {queueArgument}'
		else:
			jobCMD = f'sbatch --wrap="singularity run --bind $PWD singularity/DiffDockHPC.sif python3 inference.py --protein_ligand_csv {csvFilePath} --samples_per_complex {args.num_outputs} {remove_hs} --out_dir {outputDir}/molecules/ {keep_original_struct} {keep_cache} -c {str(args.cores)}" --mem {args.mem} --output={outputDir}/jobs_out/job_{str(i+1)}_%j.out --job-name=DiffDockHPC -c {str(args.cores)} {timeArg} {queueArgument}'
	else:
		if args.gpu == True:
			jobCMD = f'singularity run --nv --bind $PWD singularity/DiffDockHPC.sif python3 inference.py --protein_ligand_csv {csvFilePath} --samples_per_complex {args.num_outputs} {remove_hs} --out_dir {outputDir}/molecules/ {keep_original_struct} {keep_cache} -c {str(args.cores)}'
		else:
			jobCMD = f'singularity run --bind $PWD singularity/DiffDockHPC.sif python3 inference.py --protein_ligand_csv {csvFilePath} --samples_per_complex {args.num_outputs} {remove_hs} --out_dir {outputDir}/molecules/ {keep_original_struct} {keep_cache} -c {str(args.cores)}'
		
	with open(f"{outputDir}/jobs/job_{str(i+1)}.sh", "w") as jobfile:
		jobfile.write("#!/usr/bin/env bash\n")
		jobfile.write(jobCMD)
	subprocess.run(jobCMD, shell=True)
