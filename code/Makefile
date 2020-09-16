#$Id: Makefile $

.SUFFIXES: .f90 .o

# Compilers

# HPC Strasbourg:
FC = ifort
FFLAGS = -fstack-check -g -O2
LIBS = -L/b/home/quant/bsenjean/LAPACK/lapack-3.7.1 -llapacke -lrefblas

# Laptop:
#FC = gfortran
#FFLAGS = -fp-stack-check -g -O2 -traceback -check all
#LIBS = -L/usr/local/lib -llapack -lrefblas

MAIN1 = beta_and_derivatives.o
OBJS1 = secant.o quadpack_double.o
MODS1 = module.o
EXEC1 = beta_and_derivatives

MAIN2 = schmidt_decomposition.o
OBJS2 = diasym.o SVD.o inverse_matrice.o
EXEC2 = schmidt_decomposition

##################### End of Configurable Options ###############
all: beta_and_derivatives schmidt_decomposition

beta_and_derivatives: ${MODS1} ${OBJS1} ${MAIN1} Makefile 
	${FC} ${FFLAGS} -o ${EXEC1} ${MAIN1} ${OBJS1} ${MODS1} ${LIBS}

schmidt_decomposition: ${OBJS2} ${MAIN2} Makefile 
	${FC} ${FFLAGS} -o ${EXEC2} ${MAIN2} ${OBJS2} ${LIBS}

.f90.o: ${MODS1} Makefile
	${FC} ${FFLAGS} -c $<
