cmake_minimum_required (VERSION 2.8)
PROJECT(Test Fortran)
ENABLE_LANGUAGE(Fortran)
add_executable(Test main.f90)
TARGET_LINK_LIBRARIES(Test ${ACSTEMLIB})
#ADD_SUBDIRECTORY(${SRC_OBJ} ${LIB})
#ADD_SUBDIRECTORY(${SRC_ACSTEM})
