import subprocess
import os
cwd = os.getcwd()
import platform
platform = platform.platform().split("-")[0]
print("platform=", platform)
if(platform == "Darwin" or platform == "darwin"):
	cmd = "cmake .."
elif(platform=="Linux" or platform == "linux"):
	cmd = "cmake -D CMAKE_Fortran_COMPILER='/global/apps/gcc/4.8.4/bin/gfortran .."
subprocess.call("rm -rf build", shell=True)
subprocess.call("rm -rf lib", shell=True)
subprocess.call("rm -rf bin", shell=True)
subprocess.call("mkdir build", shell=True)
path = "{}/build".format(cwd)
os.chdir(path)
#subprocess.call("make", shell=True)

#print("cwd = ", os.getcwd())
#cmd = "cmake -D CMAKE_Fortran_COMPILER='/global/apps/gcc/4.8.4/bin/gfortran' .."
#subprocess.call("rm -rf build; mkdir build; cd build", shell=True)
#subprocess.call("rm -rf build", shell=True)
#subprocess.call("mkdir build", shell=True) 
#subprocess.call("cd build", shell=True)
subprocess.call(cmd, shell=True)
subprocess.call("make", shell=True)
os.chdir(cwd)
#cmake -D CMAKE_Fortran_COMPILER="/global/apps/gcc/4.8.4/bin/gfortran" ..
