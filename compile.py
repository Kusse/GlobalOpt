import subprocess
import os
cwd = os.getcwd()
cmd = "cmake -D CMAKE_Fortran_COMPILER='/global/apps/gcc/4.8.4/bin/gfortran' .."
#subprocess.call("rm -rf build; mkdir build; cd build", shell=True)
#subprocess.call("rm -rf build", shell=True)
#subprocess.call("mkdir build", shell=True) 
#subprocess.call("cd build", shell=True)
subprocess.call(cmd, shell=True)

#cmake -D CMAKE_Fortran_COMPILER="/global/apps/gcc/4.8.4/bin/gfortran" ..
