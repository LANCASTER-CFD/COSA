.SUFFIXES: .F90 .f90 .F .f .o

#
# COSA executables
#

ifeq ($(HOST),oca)
   FC = mpif90 #gfortran
   PFC = mpif90
   FFLAGS= -I/usr/local/cgns-gnu/include -cpp -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-132 -fcray-pointer -O3 -ffast-math
   L_SEQ =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -L/usr/local/cgns-gnu/lib -lcgns_gnu
   #FC = ifort
   #PFC = mpiifort
   #FFLAGS=  -g -132 -fpp -m64 -O3 -r8 -I /usr/local/cgns-intel/include
   #L_SEQ =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -L/usr/local/cgns-intel/lib -lcgns
else ifeq ($(HOST),ngio-intel)
   FC = ifort
   PFC = mpiifort
   FFLAGS=  -g -I/home/nx01/nx01/adrianj/cgns/develop/include -132 -fpp -m64 -O3 -r8 #-traceback -check bounds #-O3
   L_SEQ =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -L/home/nx01/nx01/adrianj/cgns/develop/lib -lcgns
else ifeq ($(HOST),ngio-gfortran)
   FC = mpif90
   PFC = mpif90
   FFLAGS= -g -I/home/nx01/nx01/adrianj/cgns/develop-gnu/include -m64 -fallow-argument-mismatch -fdefault-real-8 -fdefault-double-8 -ffree-line-length-132 -fcray-pointer -cpp -O3  #-fbounds-check -fbacktrace  #-ffast-math -fbounds-check -fbacktrace
   L_SEQ =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -L/home/nx01/nx01/adrianj/cgns/develop-gnu/lib -lcgns
else ifeq ($(HOST),acquarius)
    FC = gfortran
    PFC = mpif90
#   FFLAGS= -m64 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-132 -fcray-pointer -cpp -O3 -ffast-math
    FFLAGS= -m64 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-132 -fcray-pointer -cpp -O3 -ffast-math -fbounds-check -fbacktrace #-fimplicit-none
    L_SEQ = -llapack -lrefblas -ltmglib -L/home/sergio/local/lapack-3.5.0
else ifeq ($(HOST),tacchino)
#  export OMPI_FC=ifort
   FC = ifort
   PFC = mpiifort
#  FFLAGS= -132 -fpp -m64 -r8 -O3 -fpe0 -traceback -check bounds -I$(MKLROOT)/include
   FFLAGS= -132 -fpp -m64 -r8 -O3 -I$(MKLROOT)/include
   L_SEQ = -L$(MKLROOT)/lib/intel64 -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm
else ifeq ($(HOST),wayland-2020)
   export OMPI_FC=ifort
   FC = ifort
   PFC = mpif90
#  FFLAGS= -132 -fpp -m64 -r8 -O3 -fpe0 -traceback -check bounds #-C
   FFLAGS= -132 -fpp -m64 -r8 -O3 -ipo
   L_SEQ = -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm
else ifeq ($(HOST),archer)

   FC = ftn
   PFC = ftn

#  Lapack static linking start
   FFLAGS= -132 -fpp -r8 -m64 -O3 -I$(MKLROOT)/include
   L_SEQ = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm
#  Lapack static linking end

#  Lapack dynamic linking start
#  FFLAGS= -132 -fpp -r8 -m64 -O3 -mkl=sequential
#  L_SEQ = -lpthread -lm
#  Lapack dynamic linking end

else ifeq ($(HOST),joule)
   FC  = bgxlf
   PFC = mpixlf90_r
   FFLAGS= -q64 -qhot -O3 -qstrict=none:library -qautodbl=dbl4 -qfree=132 -qalias=intptr -qextname=etime
#    FFLAGS= -q64 -qhot -O3 -qstrict -qautodbl=dbl4 -qfree=132 -qalias=intptr -qextname=etime#-g 
##  serial does not work without -qstrict
##  using essl the code is twice as fast as using lapack and blas
   LDFLAGS= ${FFLAGS}
   L_SEQ = -lesslbg -L${ESSL_LIB}
## L_SEQ = -llapack -lblas -L${LAPACK_LIB}

#  FC = gfortran
#  PFC = mpif90
#  FFLAGS= -cpp -fdefault-real-8 -fdefault-double-8 -ffree-line-length-132 -fcray-pointer -O3 -ffast-math  -fno-underscoring -g
###Remark 1: Lapack for MPI code needs to be compiled using mpif90 and NOT gfortran
###Remark 2: Use -fno-underscoring in case of SEQ, alternatively recompile lapack using gfortran
#  L_SEQ = -llapack -lrefblas -ltmglib -L/gpfs/home/HCBG016/mmc01/jxd02-mmc01/local/lapack-3.4.0
else ifeq ($(HOST),fermi)

   FC  = bgxlf90
   PFC = mpixlf90
   FFLAGS= -q64 -qhot -O3 -qarch=qp -qtune=qp -qautodbl=dbl4 -qfree=132 -qalias=intptr -qextname=etime #-g #-qstrict
#  mpi code works fine without -qstrict; serial does not work without -qstrict
   LDFLAGS= ${FFLAGS}
   L_SEQ = -llapack -L${LAPACK_LIB} -lblas -L${BLAS_LIB}

#  FC = gfortran
#  PFC = mpif90
#  FFLAGS= -cpp -fdefault-real-8 -fdefault-double-8 -ffree-line-length-132 -fcray-pointer -O3 -ffast-math #-fno-underscoring #-g
#  -nounderscoring is needed to make serial code work with FERMI LAPACK; mpi code works without -fno-underscoring
#  L_SEQ = -L/cineca/prod/libraries/lapack/3.4.1/bgq-gnu--4.4.6/lib -llapack -L/cineca/prod/libraries/blas/2007/bgq-gnu--4.4.6/lib -lblas

#  FC  = bgxlf90
#  PFC = mpif90
#  FFLAGS= -q64 -qhot -O3 -qstrict -qautodbl=dbl4 -qfree=132 -qalias=intptr -qextname=etime #-g
#  LDFLAGS= ${FFLAGS}
#  following libraries work for both serial and MPI code
#  L_SEQ = -llapack -L${LAPACK_LIB} -lblas -L${BLAS_LIB}
#  following libraries are to fix LSP prob. but work only in MPI
#  L_SEQ = -L${MASS_LIB} -lmass -lmassv -L${ESSL_LIB} -lesslbg

#      FC = gfortran
#      PFC = mpif90
#      FFLAGS= -cpp -fdefault-real-8 -fdefault-double-8 -ffree-line-length-132 -fcray-pointer -O3 #-g
#   #  following libraries are for MPI code
#      L_SEQ = -L/cineca/prod/libraries/lapack/3.4.1/bgq-gnu--4.4.6/lib -llapack -L/cineca/prod/libraries/blas/2007/bgq-gnu--4.4.6/lib -lblas
#   #  following library if to make LSP code work with MPI code
#   #  L_SEQ = -L${MASS_LIB} -lmass -lmassv
#   #  following libraries are for serial code (they are local!)
#   #  L_SEQ = -llapack -lrefblas -ltmglib -L/fermi/home/userexternal/jdrofeln/local/lapack-3.4.2
#   #  at present serial code does not work with LSP code

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

#  FC = gfortran
#  FFLAGS= -g -m64 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fcray-pointer -cpp -O3

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

.F.o .f.o:     ; ${COMP} ${CPPFLAGS} ${FFLAGS} -c $< 

.f90.o:       ; ${COMP} ${CPPFLAGS} ${FFLAGS} -c $<
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

