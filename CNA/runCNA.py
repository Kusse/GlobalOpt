 # 1 #!/bin/sh
 # 2 
#inFile="AlFcc.r"
inFile="Au55_chiral_min.r"
 # 5 ######################################
#outFile="AlFcc.CN"
outFile="Au55_chiral_min.CN"
rcut="3.6"
 # 9 cnA.exc -v -R $rcut -I $inFile -F baskes -O $outFile

import subprocess
import os
#1 Compile CNA
subprocess.call(["make"])

#2 Make it executable
cmd = "chmod +x cnA.exc"
subprocess.call(cmd, shell=True)

#3 Export path
cwd = os.getcwd()
# export PATH=/Users/KSB/git/GlobalOptimization/GlobalOpt/CNA:$PATH
cmd = "export PATH={}/:$PATH".format(cwd)
print(cmd)
subprocess.call(cmd, shell=True)
#4 Run CNA

#cmd = "./RunCNA.scr"
cmd = "./cnA.exc -v -R {} -I {} -F baskes -O {}".format(rcut, inFile, outFile)
subprocess.call(cmd, shell=True)

#5 Run CNA_post_processing
cmd = "python cnbond.py"
subprocess.call(cmd, shell=True)
#6 Done

