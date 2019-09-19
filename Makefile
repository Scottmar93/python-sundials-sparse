SHELL = sh

prefix       = /home/scott/Projects/ida_test/sundials/instdir
exec_prefix  = /home/scott/Projects/ida_test/sundials/instdir
includedir   = /home/scott/Projects/ida_test/sundials/instdir/include
libdir       = /home/scott/Projects/ida_test/sundials/instdir/lib

CPP      = /usr/bin/c++
CPPFLAGS = -O3 -DNDEBUG
CC       = /usr/bin/c++
CFLAGS   = -O3 -DNDEBUG
LDFLAGS  = 
LIBS     =  -lm /usr/lib/x86_64-linux-gnu/librt.so -lblas

LINKFLAGS = -Wl,-rpath,/home/scott/Projects/ida_test/sundials/instdir/lib

# -----------------------------------------------------------------------------------------

LIBRARIES_LAPACK = -lsundials_sunlinsollapackdense -lsundials_sunlinsollapackband  
LINKFLAGS_LAPACK = ${LINKFLAGS}::

# INCLUDES_KLU  = 
# LIBRARIES_KLU = -lsundials_sunlinsolklu 
# LINKFLAGS_KLU = ${LINKFLAGS}:

# INCLUDES_SLUMT  = 
# LIBRARIES_SLUMT = -lsundials_sunlinsolsuperlumt  
# LINKFLAGS_SLUMT = ${LINKFLAGS}::

TMP_INCS  = ${includedir} 
INCLUDES  = $(addprefix -I, ${TMP_INCS})
# INCLUDES  = $(${TMP_INCS})
LIBRARIES = -lsundials_idas -lsundials_nvecserial ${LIBS}

# -----------------------------------------------------------------------------------------

EXAMPLES = my_simple_example
EXAMPLES_DEPENDENCIES = residual jacobian events

OBJECTS = ${EXAMPLES:=.o}
OBJECTS_DEPENDENCIES = ${EXAMPLES_DEPENDENCIES:=.o}

# -----------------------------------------------------------------------------------------

.SUFFIXES : .o .c

.c.o :
	${CC} ${CFLAGS} ${INCLUDES} -c $<

# -----------------------------------------------------------------------------------------
all: sparse

my_simple_example: ${OBJECTS}
	@for i in ${EXAMPLES} ; do \
	  echo "${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS}" ; \
	  ${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS} ; \
	done

python_shared_example: ${OBJECTS}
	@for i in ${EXAMPLES} ; do \
	  echo "${CC} -FPIC -shared -o $${i}.so $${i}.o ${OBJECTS_DEPENDENCIES} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS}" ; \
	  ${CC} -FPIC -shared -o $${i}.so $${i}.o ${OBJECTS_DEPENDENCIES} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS} ; \
	done

sparse:
	g++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` -I/home/scott/Projects/sparse-test/sundials/instdir/include -I/home/scott/Projects/sparse-test/SuiteSparse/KLU/Include -I/home/scott/Projects/sparse-test/SuiteSparse/AMD/Include -I/home/scott/Projects/sparse-test/SuiteSparse/COLAMD/Include -I/home/scott/Projects/sparse-test/SuiteSparse/BTF/Include -I/home/scott/Projects/sparse-test/SuiteSparse/SuiteSparse_config -o sundials`python3-config --extension-suffix` sundials_sparse.cpp -L/home/scott/Projects/sparse-test/sundials/instdir/lib -lsundials_sunmatrixsparse -lsundials_ida -lsundials_nvecserial -lm /usr/lib/x86_64-linux-gnu/librt.so -lblas -Wl,-rpath,/home/scott/Projects/sparse-test/sundials/instdir/lib -lsundials_sunlinsolklu /home/scott/Projects/sparse-test/SuiteSparse/KLU/Lib/libklu.a /home/scott/Projects/sparse-test/SuiteSparse/AMD/Lib/libamd.a /home/scott/Projects/sparse-test/SuiteSparse/COLAMD/Lib/libcolamd.a /home/scott/Projects/sparse-test/SuiteSparse/BTF/Lib/libbtf.a /home/scott/Projects/sparse-test/SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a;

ida_rob_klu: 
	g++ -O3 -Wall -std=c++11 -fPIC -I/home/scott/Projects/sparse-test/sundials/instdir/include -I/home/scott/Projects/sparse-test/SuiteSparse/KLU/Include -I/home/scott/Projects/sparse-test/SuiteSparse/AMD/Include -I/home/scott/Projects/sparse-test/SuiteSparse/COLAMD/Include -I/home/scott/Projects/sparse-test/SuiteSparse/BTF/Include -I/home/scott/Projects/sparse-test/SuiteSparse/SuiteSparse_config -o ida_robs idaRoberts_klu.c -L/home/scott/Projects/sparse-test/sundials/instdir/lib -lsundials_sunmatrixsparse -lsundials_ida -lsundials_nvecserial -lm /usr/lib/x86_64-linux-gnu/librt.so -lblas -Wl,-rpath,/home/scott/Projects/sparse-test/sundials/instdir/lib -lsundials_sunlinsolklu /home/scott/Projects/sparse-test/SuiteSparse/KLU/Lib/libklu.a /home/scott/Projects/sparse-test/SuiteSparse/AMD/Lib/libamd.a /home/scott/Projects/sparse-test/SuiteSparse/COLAMD/Lib/libcolamd.a /home/scott/Projects/sparse-test/SuiteSparse/BTF/Lib/libbtf.a /home/scott/Projects/sparse-test/SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a;

${OBJECTS}: ${OBJECTS_DEPENDENCIES}

clean:
	rm -f ${OBJECTS_DEPENDENCIES} 
	rm -f ${OBJECTS}
	rm -f ${EXAMPLES}.so



# c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` example.cpp -o example`python3-config --extension-suffix`