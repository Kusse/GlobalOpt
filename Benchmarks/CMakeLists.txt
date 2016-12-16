ENABLE_LANGUAGE(Fortran)
SET(OBJSRC 	${SRC_OBJ}/objective.f90
		${SRC_OBJ}/Potentials/lj.f90
		${SRC_OBJ}/Potentials/morse.f90
		${SRC_OBJ}/Potentials/Gupta.f
		${SRC_OBJ}/Potentials/modhess.f90
		${SRC_OBJ}/ACSTEM/acstem.f90
		${SRC_OBJ}/ACSTEM/filter.f90
		${SRC_OBJ}/ACSTEM/merge_rank.f90
		${SRC_OBJ}/ACSTEM/acstem_parameters.f90
		${SRC_OBJ}/objective_parameters.f90
)
ADD_LIBRARY(${OBJLIB} STATIC ${OBJSRC})
