.SUFFIXES: .F90 .f90 .F .f .o .mod

#
# COSA executables
#

ifeq ($(HOST),cirrus)
   FC = mpif90
   PFC = mpif90
   FFLAGS= -I/work/z04/shared/adrianj/cgns/include -cpp -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-132 -fcray-pointer -O3 -ffast-math  -fallow-argument-mismatch
   L_SEQ =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -L/work/z04/shared/adrianj/cgns/lib -lcgns -L/work/z04/shared/adrianj/tecplot -lteciompi -lstdc++ -ldl
else ifeq ($(HOST),oca)
   FC = mpif90 #gfortran
   PFC = mpif90
   FFLAGS= -I/usr/local/cgns-gnu/include -cpp -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-132 -fcray-pointer -O3 -ffast-math
   L_SEQ =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -L/usr/local/cgns-gnu/lib -lcgns_gnu
else ifeq ($(HOST),ngio-intel)
   FC = ifort
   PFC = mpiifort
   FFLAGS=  -g -I/home/nx01/nx01/adrianj/cgns/develop/include -132 -fpp -m64 -O3 -r8 #-traceback -check bounds #-O3
   L_SEQ =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -L/home/nx01/nx01/adrianj/cgns/develop/lib -lcgns
else ifeq ($(HOST),ngio-gfortran)
   FC = mpif90
   PFC = mpif90
   FFLAGS= -g -I/home/nx01/nx01/adrianj/cgns/develop-gnu/include -m64 -fallow-argument-mismatch -fdefault-real-8 -fdefault-double-8 -ffree-line-length-132 -fcray-pointer -cpp -O3  #-fbounds-check -fbacktrace  #-ffast-math -fbounds-check -fbacktrace
   L_SEQ =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -L/home/nx01/nx01/adrianj/cgns/develop-gnu/lib -lcgns -L/home/software/tecplot/2023.1.0/lib -lteciompi -lstdc++
else ifeq ($(HOST),acquarius)
    FC = gfortran
    PFC = mpif90
    FFLAGS= -m64 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-132 -fcray-pointer -cpp -O3 -ffast-math -fbounds-check -fbacktrace #-fimplicit-none
    L_SEQ = -llapack -lrefblas -ltmglib -L/home/sergio/local/lapack-3.5.0
else ifeq ($(HOST),tacchino)
   FC = ifort
   PFC = mpiifort
   FFLAGS= -132 -fpp -m64 -r8 -O3 -I$(MKLROOT)/include
   L_SEQ = -L$(MKLROOT)/lib/intel64 -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm
else ifeq ($(HOST),wayland-2020)
   export OMPI_FC=ifort
   FC = ifort
   PFC = mpif90
   FFLAGS= -132 -fpp -m64 -r8 -O3 -ipo
   L_SEQ = -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm
else ifeq ($(HOST),archer)
   FC = ftn
   PFC = ftn
#  Lapack static linking start
   FFLAGS= -132 -fpp -r8 -m64 -O3 -I$(MKLROOT)/include
   L_SEQ = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm
#  Lapack static linking end
else ifeq ($(HOST),joule)
   FC  = bgxlf
   PFC = mpixlf90_r
   FFLAGS= -q64 -qhot -O3 -qstrict=none:library -qautodbl=dbl4 -qfree=132 -qalias=intptr -qextname=etime
   LDFLAGS= ${FFLAGS}
   L_SEQ = -lesslbg -L${ESSL_LIB}
else ifeq ($(HOST),fermi)
   FC  = bgxlf90
   PFC = mpixlf90
   FFLAGS= -q64 -qhot -O3 -qarch=qp -qtune=qp -qautodbl=dbl4 -qfree=132 -qalias=intptr -qextname=etime #-g #-qstrict
   LDFLAGS= ${FFLAGS}
   L_SEQ = -llapack -L${LAPACK_LIB} -lblas -L${BLAS_LIB}
else ifeq ($(HOST),marconi-xlf)
   FC  = xlf
   PFC = mpixlf
   FFLAGS= -q64 -qhot -O3 -qstrict=none:library -qautodbl=dbl4 -qfree=132 -qalias=intptr -qextname=etime
   LDFLAGS= ${FFLAGS}
   L_SEQ = -lessl -L${ESSL_LIB}
else ifeq ($(HOST),marconi-pgi)
   FC  = nvfortran
   PFC = mpif90
   FFLAGS= -fastsse -Mextend -Minform=warn -g  -D MPI -r8 -Mcray=pointer  -Mpreprocess
   LDFLAGS= ${FFLAGS}
   L_SEQ = -lessl -L${ESSL_LIB}
endif

mpi:   COMP = ${PFC}
mpixl: COMP = ${PFC}
mpi:   CPPFLAGS = -DMPI
mpixl: CPPFLAGS = -qsuffix=cpp=f -WF,-DMPI -WF,-DXL
seq:   COMP = ${FC}
seqxl: COMP = ${FC}
seq:   CPPFLAGS =
seqxl: CPPFLAGS = -qsuffix=cpp=f -WF,-DXL
pos: COMP = ${FC}
pos: CPPFLAGS =

COSA   = cosa_precision.o common_variables.o cosa_variables.o parallelutils.o input.o cosa.o flux.o src.o update.o losp.o \
         grid.o util.o bc.o debug.o output.o 
POS    = cosa_precision.o common_variables.o cosa_variables.o parallelutils.o input.o poscosa.o grid.o util.o flux.o \
         debug.o bc.o

#
# transformation rules
#

.f90.o:       ; ${COMP} ${CPPFLAGS} ${FFLAGS} -c $<
%.o: %.mod

#
# program dependencies
#

#
# programs
#

default:
	@echo 'Please specify seq for sequential or mpi for parallel version.'
	@echo 'Warning: do not specify both targets on the same make invocation.'

mpi:	${COSA}
	${PFC} $(FFLAGS) -o cosa.mpi ${COSA} ${L_SEQ}
#	${FC} ${FFLAGS} -o ../bin/cosa_$(ver).par ${COSA}
mpixl:	${COSA}
	${PFC} $(FFLAGS) -o cosa.mpi ${COSA} ${L_SEQ}

seq:	${COSA}
	${FC} $(FFLAGS) -o cosa.seq ${COSA} ${L_SEQ}
seqxl:	${COSA}
	${FC} $(FFLAGS) -o cosa.seq ${COSA} ${L_SEQ}

pos:	${POS}
	${FC} ${FFLAGS} -o poscosa.seq ${POS}
#	${FC} ${FFLAGS} -o ../bin/poscosa_$(ver).seq ${POS}

#
# clean things up
#

clean:	
	/bin/rm -f *.o *.mod 

