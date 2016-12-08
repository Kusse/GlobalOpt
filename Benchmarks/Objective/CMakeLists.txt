# ACSTEM LIBRARY
ENABLE_LANGUAGE(Fortran)
SET(OBJ_SRC     ${SRC_OBJ}/Gupta.f
                ${SRC_OBJ}/lj.f90
		${SRC_OBJ}/modhess.f90
                ${SRC_OBJ}/morse.f90
                ${SRC_OBJ}/objective.f90
)
ADD_LIBRARY(${OBJLIB} STATIC ${OBJ_SRC})
