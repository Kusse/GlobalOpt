SET(GA_SRC ${SRC_GA}/commons.f90
		  ${SRC_GA}/comparator.f90
		  ${SRC_GA}/dprand.f
		  ${SRC_GA}/ga_cluster.f90
		  ${SRC_GA}/ga_main.f90
		  ${SRC_GA}/ga_params.f90
		  ${SRC_GA}/ga_population.f90
		  ${SRC_GA}/ga_selection.f90
		  ${SRC_GA}/modhess.f90
		  ${SRC_GA}/mycpu_time.f90
		  ${SRC_GA}/qmodule.f90
)
# Define the executable in terms of the source files
ADD_EXECUTABLE(${GARUN} ${GA_SRC})
#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# This links foo to the bar library
TARGET_LINK_LIBRARIES(${SEARUN} ${ACSTEMLIB})
TARGET_LINK_LIBRARIES(${SEARUN} ${LBFGSBLIB})

# Uncomment if you need to link to BLAS and LAPACK
#TARGET_LINK_LIBRARIES(${FOOEXE} ${BLAS_LIBRARIES}
#                                ${LAPACK_LIBRARIES}
#                                ${CMAKE_THREAD_LIBS_INIT})

# Uncomment if you have parallization
#IF(USE_OPENMP)
#    SET_TARGET_PROPERTIES(${FOOEXE} PROPERTIES
#                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
#                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
#ELSEIF(USE_MPI)
#    SET_TARGET_PROPERTIES(${FOOEXE} PROPERTIES
#                          COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
#                          LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
#    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
#    TARGET_LINK_LIBRARIES(${FOOEXE} ${MPI_Fortran_LIBRARIES})
#ENDIF(USE_OPENMP)

#####################################
# Tell how to install this executable
#####################################

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
INSTALL(TARGETS ${SEARUN} RUNTIME DESTINATION bin)