SET(LBFGSB_SRC ${SRC_LBFGSB}/blas.f
		  ${SRC_LBFGSB}/lbfgsb.f
		  ${SRC_LBFGSB}/linpack.f
		  ${SRC_LBFGSB}/timer.f
)
ADD_LIBRARY(${LBFGSBLIB} STATIC ${LBFGSB_SRC})
