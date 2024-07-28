#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 13:46:19 2024

@author: Jochem Nelen (jnelen@ucam.edu)
"""

# This script is used to rerun molecules that failed to dock during a DiffDockHPC job. It checks the output directory for missing sdf files and prompts the user to relaunch the failed molecules.
# To use this script, simply run it from the command line and provide the main output directory as an argument:
# python relaunchFailedCompounds.py VS_DD_.../


import glob
import os
import re
import shutil
import subprocess
import sys

if len(sys.argv) < 2:
	sys.exit("You have to put in a DiffDockHPC run as an argument")
	
inputPath = sys.argv[1]

if not os.path.isdir(inputPath):
	sys.exit("The input path doesn't seem to be a valid directory")

successCounter = 0
pathDict = {}
finishedList = []

## Code to distribute the query ligands among the amount of jobs 
def split(a, n):
    if n > len(a):
        print("more jobs than files, launching 1 job per file")
        return [a[i:i+1] for i in range(len(a))]
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

## Check the csv files and store the information to a dict
for path in glob.glob(f"{inputPath}/csvs/*.csv"):
	with open(path) as inputFile:
		inputLines = inputFile.readlines()
		for line in inputLines[1:]:
			lineSplit = line.strip().split(';')
			molName = lineSplit[0]
			molPath = lineSplit[-1]		
			pathDict[molName] = molPath
			
finishedPaths = glob.glob(f"{inputPath}/molecules/*.sdf")

## Identify which compounds finished successfully
for finishedPath in finishedPaths:
	finishedName = finishedPath.split("VS_DD_")[-1].split("_rank")[0]
	finishedList.append(finishedName)

## Loop through the finished compounds and remove the values from the pathDict
for name in finishedList:
	del(pathDict[name])
	successCounter += 1

print(f"Failed to process {len(pathDict)} files, {successCounter} did process successfully")

## Ask with how many jobs it should be rerun (interactively)
answer = input("With how many jobs would you like to relaunch the failed compounds? ").lower()
if not answer.isdigit():
	sys.exit("Not launching any jobs..")
else:
	jobNumber = int(answer)
	print(f"launching {jobNumber} jobs..")

ligandPathsSplit = list(split(list(pathDict.values()), jobNumber))

## Check if the output directory already exists, and asks the user what to do if it does
redoDir = f"{inputPath}/redo/"
if os.path.isdir(redoDir):
	print(f"The directory {redoDir} already exists. To continue you must delete this directory or choose another outputname.")
	answer = input("Do you want to remove it? (y/n) ").lower()
	while answer not in ("y", "n", "yes", "no"):
		print("Invalid input. Please enter y(es) or n(o).")
		answer = input("Do you want to remove or overwrite it? (y/n) ").lower()
	if answer == "y" or answer == "yes":
		shutil.rmtree(redoDir, ignore_errors=False)			
	else:
		sys.exit()
		
## write a "redo" or directory in the output directory, write csvs using this information and relaunch it as jobs using the original settings
os.mkdir(redoDir)
os.mkdir(f"{redoDir}/csvs/")
os.mkdir(f"{redoDir}/jobs_out/")
os.mkdir(f"{redoDir}/jobs/")

## Get a job file to copy the settings automatically
jobPaths = glob.glob(f"{inputPath}/jobs/job_*.sh")
with open(jobPaths[0]) as jobFile:
	jobLine = jobFile.readlines()[1]

## Get the protein path
protein_path = glob.glob(f"{inputPath}/*.pdb")[0]

## Write and launch jobs
for i, jobLigands in enumerate(ligandPathsSplit):
	csvFilePath = f"{redoDir}/csvs/job_csv_{str(i+1)}.csv"
	with open(csvFilePath, 'w') as jobCSV:
		jobCSV.write("complex_name;protein_path;ligand_description\n")
		for jobLigand in jobLigands:
			complexName = os.path.basename(jobLigand).split('.')[0]
			jobCSV.write(f"{complexName};{protein_path};{jobLigand}\n")
	jobCSV.close()

	## Modify the command to use the correct csv and job output path
	jobCMD = re.sub(r"/jobs_out/job_.*\.out", f"/redo/jobs_out/redo_job_{i}_%j.out", jobLine)
	jobCMD = re.sub(r"(--protein_ligand_csv\s+)[^ ]+(\s+--samples_per_complex)", r"\1" + csvFilePath + r"\2", jobCMD)

	## Write the DiffDockHPC job file	
	with open(f"{redoDir}/jobs/redo_job_{str(i+1)}.sh", "w") as jobfile:
		jobfile.write("#!/usr/bin/env bash\n")
		jobfile.write(jobCMD)

	## Run the DiffDockHPC job file
	subprocess.run(jobCMD, shell=True)
