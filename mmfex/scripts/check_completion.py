#!/usr/bin/env python

import os
import shutil
import subprocess
import glob

m = <NUMBER OF THE FIRST MONOMER> #first monomer
n = <NUMBER OF THE LAST MONOMER> #last monomer

# all the files are into the current work directory
source = os.getcwd()

# list of failed couplings
fails = []

# the cycle starts here
for i in range(m, n+1):
  alfa = str(i)
  for j in range(i+1, n+1):
    beta = str(j)
    
#   define the destination directory for the files
    dest = source + '/coup_' + alfa + '_' + beta
    os.chdir(dest)
    
#   checks for "Happy Landing" 
    if (subprocess.Popen('grep Happy *out', shell=True, stdout=subprocess.PIPE,stderr=subprocess.DEVNULL, text=True).communicate()[0] == ''):
      fails.append(dest)
    
    os.chdir(source)
    
if not fails:
  print('All done, good job!\n')
else:
  print('Some work still to be done. I will submit the corresponding calculations, if required\n')
  i = -1
  #queu = ['batch-l', 'long-l', 'pcomp-h', 'pcomp-l']
  #n = len(queu)
  runningjobs = subprocess.Popen('squeue -h -o "%20j" -u devicol', shell=True, stdout=subprocess.PIPE,stderr=subprocess.DEVNULL, text=True).communicate()[0].split()
  for dirs in fails:
    tocheck = 'coupling_' + dirs.split('/')[-1][-7:]
    if (tocheck in runningjobs): # the job is still running
      print(tocheck + ' is still running or in the queue')
    else:
      i = (i + 1) % n
      os.chdir(dirs)
      inputfile = glob.glob('*.input')[0]
      print('I will re-submit ' + tocheck)
      os.popen(<FILL HERE WITH YOUR USUAL COMMAND LINE RUN>)
      #os.popen('submit_OMP_2024.sh ' + inputfile + ' 1 50 ' + queu[i])

