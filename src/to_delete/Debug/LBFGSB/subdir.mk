################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../LBFGSB/mylbfgsb.f90 

F_SRCS += \
../LBFGSB/blas.f \
../LBFGSB/lbfgsb.f \
../LBFGSB/linpack.f \
../LBFGSB/timer.f 

OBJS += \
./LBFGSB/blas.o \
./LBFGSB/lbfgsb.o \
./LBFGSB/linpack.o \
./LBFGSB/mylbfgsb.o \
./LBFGSB/timer.o 


# Each subdirectory must supply rules for building sources it contributes
LBFGSB/%.o: ../LBFGSB/%.f
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

LBFGSB/blas.o: ../LBFGSB/blas.f

LBFGSB/lbfgsb.o: ../LBFGSB/lbfgsb.f

LBFGSB/linpack.o: ../LBFGSB/linpack.f

LBFGSB/%.o: ../LBFGSB/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

LBFGSB/mylbfgsb.o: ../LBFGSB/mylbfgsb.f90 commons.o ga_population.o objective.o

LBFGSB/timer.o: ../LBFGSB/timer.f


