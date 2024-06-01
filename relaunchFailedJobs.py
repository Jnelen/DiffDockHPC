#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 15:58:13 2024

@author: Jochem Nelen (jnelen@ucam.edu)
"""

### Sometimes jobs fail to launch or finish. This can be especially annoying when launching many jobs, and it can be a hassle to run them again.
### This script checks which jobs didn't finish correctly and tries to relaunch them.
### Note: Only use this script after all jobs associated with this run have been finished. Jobs that are still running might be initiated again and cause problems.
### The only argument you should provide is the main output directory as produced by DiffDockHPC. It should contain the jobs and jobs_out directories.

import glob
import os
import subprocess
import sys

if len(sys.argv) < 2:
	sys.exit("You have to put in a DiffDockHPC run as an argument")
	
inputPath = sys.argv[1]

if not os.path.isdir(inputPath):
	sys.exit("The input path doesn't seem to be a valid directory")
	
finishedList = []
jobPaths = []

## Check the .out files to see if the job finished successfully
for path in glob.glob(f"{inputPath}/jobs_out/*.out"):
	with open(path) as inputFile:
		inputLines = inputFile.readlines()
		for line in inputLines:
			if "Calculations finished after " in line:
				## Add this jobnumber to joblist so it can be rerun
				finishedList.append(os.path.basename(path).split(".")[0].split("job_")[-1].split("_")[0])
				break

finishedList = set(finishedList)

## Get all the paths and sort them by jobnumber
jobPaths = sorted(glob.glob(f"{inputPath}/jobs/*.sh"), key=lambda x: int(os.path.basename(x).split("job_")[-1].split("_")[0].split(".")[0]))

if not len(jobPaths) > len(finishedList):
	sys.exit(f"All the {len(finishedList)} jobs have successfully finished")
	
## Check the jobs and relaunch them if they haven't finished
for path in jobPaths:
	jobNumber = str(os.path.basename(path).split("job_")[-1].split("_")[0].split(".")[0])
	
	if not jobNumber in finishedList:
		print(f"Relaunching job_{jobNumber}")
		subprocess.run(f"sh {path}", shell=True)
		