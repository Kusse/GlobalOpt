################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../commons.f90 \
../comparator.f90 \
../ga_cluster.f90 \
../ga_main.f90 \
../ga_params.f90 \
../ga_population.f90 \
../ga_select.f90 \
../main.f90 \
../modhess.f90 \
../mycpu_time.f90 \
../objective.f90 \
../qmodule.f90 

F_SRCS += \
../Gupta.f \
../dprand.f 

OBJS += \
./Gupta.o \
./commons.o \
./comparator.o \
./dprand.o \
./ga_cluster.o \
./ga_main.o \
./ga_params.o \
./ga_population.o \
./ga_select.o \
./main.o \
./modhess.o \
./mycpu_time.o \
./objective.o \
./qmodule.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

Gupta.o: ../Gupta.f

%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

commons.o: ../commons.f90

comparator.o: ../comparator.f90

dprand.o: ../dprand.f

ga_cluster.o: ../ga_cluster.f90 commons.o ga_params.o ga_population.o

ga_main.o: ../ga_main.f90 LBFGSB/mylbfgsb.o commons.o ga_params.o ga_population.o qmodule.o

ga_params.o: ../ga_params.f90

ga_population.o: ../ga_population.f90 ga_params.o

ga_select.o: ../ga_select.f90 ga_params.o ga_population.o

main.o: ../main.f90

modhess.o: ../modhess.f90 commons.o

mycpu_time.o: ../mycpu_time.f90 commons.o

objective.o: ../objective.f90

qmodule.o: ../qmodule.f90


