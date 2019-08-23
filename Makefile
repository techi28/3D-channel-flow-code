#######################################################################
#                                                                     #
#                        MAKEFILE 3D: Sequential                      #
#                      ---  Miguel Uh Zapata ---                      #
#                                 2013                                #
#                                                                     #
#######################################################################
#
# .__________________________________________________________________.
# |                                                                  |
# |                     Declaration of extension                     |
# |__________________________________________________________________|


.SUFFIXES:.F90

# .__________________________________________________________________.
# |                                                                  |
# |                   Compiler and editor of links                   |
# |__________________________________________________________________|
# |                                                                  |
# | Sequential:   FC  = gfortran                                     |
# | # Compiler Fortran77         : xlf                               |
# |                                                                  |
# | Parallel:     FC  = mpif77                                       |
# | # Compiler Fortran77 + OpenMP: xlf_r / OPAR = -qsmp=omp / KeyOMP |
# |__________________________________________________________________|

FC         = mpif77
FLD        = $(FC)
RM         = /bin/rm -f
MV         = mv

#######################################################################
#                            COMPILATION                              #
#######################################################################

# .__________________________________________________________________.
# |                                                                  |
# |        Optimization options                                      |
# |                                                                  |
# |  Joe    : OLEV  = -qhot -O3 -qtune=auto -qarch=auto -qcache=auto |
# |  Power4 : OLEV  = -qhot -O3 -qtune=pwr4 -qarch=pwr4              |
# |  Power5 : OLEV  = -qhot -O3 -qtune=pwr5 -qarch=pwr5              |
# |__________________________________________________________________|

OLEV       =
FOPT       = -qstrict
FOPT       =
FMAX       = -qipa
FMAX       =

# .__________________________________________________________________.
# |                                                                  |
# |        Parallelization options: OpenMP                           |
# |__________________________________________________________________|

OPAR       = -qsmp=omp
OPAR       =

# .__________________________________________________________________.
# |                                                                  |
# |        Debugging options                                         |
# |__________________________________________________________________|

ODEB       = -qnooptimize
FDEB       = -qcheck -qdbg  \
#            -qflttrap=:ov:und:zero:inv:en \
             -qinitauto=FF \
			 -qhalt=l
#            -qextchk
#            -qfullpath
#            -qundef : implicit none

# .__________________________________________________________________.
# |                                                                  |
# |        Portage options                                           |
# |__________________________________________________________________|

OPRT       = -qautodbl=dbl4 -qdpc=e
OPRT       =

# .__________________________________________________________________.
# |                                                                  |
# |        Inputs/outputs options of the compiler                    |
# |__________________________________________________________________|

OES        = -qfixed=132 -qsuffix=f=f
OES        =

# .__________________________________________________________________.
# |                                                                  |
# |         Listing options                                          |
# |                                                                  |
# |  Sequential:  OLST   = -qlistopt -qreport=hotlist -qsource       |
# |  OpenMP    :  OLST   = -qlistopt -qreport=smplist -qsource       |
# |__________________________________________________________________|

OLST       = -qlistopt -qreport=hotlist -qsource
OLST       =

# .__________________________________________________________________.
# |                                                                  |
# |        Profile options                                           |
# |__________________________________________________________________|

PROF       = -pg -qfullpath -qdbg
PROF       =

# .__________________________________________________________________.
# |                                                                  |
# |        COMPILATION:                                              |
# |                                                                  |
# |        optimize:                                                 |
# |                 FFLAGS =  -fdefault-real-8 -std=f95 \            |
# |                           $(OLEV) $(FOPT) $(FMAX) $(OPAR)\       |
# |                           $(OPRT) $(PROF) $(OES) $(OLST)         |
# |        debug:                                                    |
# |                 FFLAGS =  $(ODEB)$(FDEB) $(FMAX) $(OPAR)\        |
# |                           $(OPRT) $(PROF) $(OES) $(OLST)         |
# |__________________________________________________________________|


FFLAGS = $(OLEV) $(FOPT) $(FMAX) $(OPAR) $(OPRT) $(PROF) $(OES) $(OLST)


#######################################################################
#                         EDITION OF LINKS                            #
#######################################################################

# .__________________________________________________________________.
# |                                                                  |
# |        Include Paths:                                            |
# |        CRIHAN                                                    |
# |        INC       = -I/soft/library/netCDF.32/include/            |
# |        GNU                                                       |
# |__________________________________________________________________|

INC       = -I /usr/include/

# .__________________________________________________________________.
# |                                                                  |
# |       Fortran libraries                                          |
# |__________________________________________________________________|

#LIBS      = -lmass -lmassv -lessl -lm /soft/library/netCDF.32/lib/libnetcdf.a
LIBS      = -lm

# .__________________________________________________________________.
# |                                                                  |
# |       Optimization option                                        |
# |__________________________________________________________________|

LDLEV     = $(OLEV)

# .__________________________________________________________________.
# |                                                                  |
# |        Portage options                                           |
# |__________________________________________________________________|

#LDREN     = -brename:.flush,.flush_ -brename:.etime,.etime_ -brename:.lnblnk,.lnblnk_
#LDREN     = -brename:.etime,.etime_ -brename:.lnblnk,.lnblnk_
#LDREN     = -pedantic -std=std

# .__________________________________________________________________.
# |                                                                  |
# |        Memory options                                            |
# |__________________________________________________________________|

LDMEM     = -bmaxdata:0x80000000 -blpdata
LDMEM     = -blpdata
LDMEM     =

# .__________________________________________________________________.
# |                                                                  |
# |        Edition of links options                                  |
# |__________________________________________________________________|

LDFLAGS   = $(LDLEV) $(LDMEM) $(OPAR) $(LDREN) $(PROF)


#######################################################################
#                         Fortran includes                            #
#######################################################################

# .__________________________________________________________________.
# |                                                                  |
# |        Fortran includes                                          |
# |__________________________________________________________________|

FHDRS      = \
             common.mpf

# .__________________________________________________________________.
# |                                                                  |
# |        Directories: object, source & input                       |
# |__________________________________________________________________|


OBJDIR    = bin
SRCDIR    = src
INPUTDIR  = input
VPATH     = ${SRCDIR}:\
            ${SRCDIR}/Main:\
            ${SRCDIR}/AdvectionDiffusion:\
            ${SRCDIR}/BoundaryConditions:\
            ${SRCDIR}/Parallel:\
            ${SRCDIR}/Interpolation:\
            ${SRCDIR}/Gradient:\
            ${SRCDIR}/LoopUpdate:\
            ${SRCDIR}/Save:\
            ${SRCDIR}/Stats:\
            ${SRCDIR}/LES:

# .__________________________________________________________________.
# |                                                                  |
# |        Objects Fortran                                           |
# |__________________________________________________________________|

FOBJS = \
	${OBJDIR}/module_variables.o\
	${OBJDIR}/module_geometry.o\
	${OBJDIR}/module_interfaces.o\
	${OBJDIR}/module_les.o\
	${OBJDIR}/module_stats.o\
	${OBJDIR}/module_parallel3D.o\
	${OBJDIR}/parallel_distribution.o\
	${OBJDIR}/parallel_input_global.o\
	${OBJDIR}/parallel_input_local.o\
	${OBJDIR}/parallel_parameters.o\
	${OBJDIR}/nsmp3D.o\
	${OBJDIR}/hydro.o\
	${OBJDIR}/input.o\
	${OBJDIR}/input_parameters.o\
	${OBJDIR}/calcul_geometry.o\
	${OBJDIR}/initial.o\
	${OBJDIR}/user_initial.o\
	${OBJDIR}/cfl.o\
	${OBJDIR}/stats.o\
	${OBJDIR}/bicgstab.o\
	${OBJDIR}/SavetecVertex.o\
	${OBJDIR}/SavetecStats.o\
	${OBJDIR}/SaveParaviewVertex.o\
	${OBJDIR}/SavetecCenter.o\
	${OBJDIR}/SavetecVertexCenter.o\
	${OBJDIR}/LoopTimeInitial.o\
	${OBJDIR}/LoopTimeUpdate.o\
	${OBJDIR}/LoopRKInitial.o\
	${OBJDIR}/LoopRKUpdate.o\
	${OBJDIR}/LoopRKFinal.o\
	${OBJDIR}/LoopStats.o\
	${OBJDIR}/interpolation2D.o\
	${OBJDIR}/interpolation3D.o\
	${OBJDIR}/interpolationFace.o\
	${OBJDIR}/interpolationU.o\
	${OBJDIR}/gradientLSM.o\
	${OBJDIR}/gradientLSM2D.o\
	${OBJDIR}/gradientFVM.o\
	${OBJDIR}/gradientP.o\
	${OBJDIR}/advection3D.o\
	${OBJDIR}/diffusion3D.o\
	${OBJDIR}/PDE.o\
	${OBJDIR}/source3D.o\
	${OBJDIR}/newu3D.o\
	${OBJDIR}/AdvDiffNS.o\
	${OBJDIR}/AdvDiffSolve.o\
	${OBJDIR}/BC.o\
	${OBJDIR}/BCoflow.o\
	${OBJDIR}/update_bc.o\
	${OBJDIR}/sgm_les.o\
	${OBJDIR}/Functions.o\


# .__________________________________________________________________.
# |                                                                  |
# |        Sources Fortran                                           |
# |__________________________________________________________________|

INC_COMMON = -I ${SRCDIR}/Include

# .__________________________________________________________________.
# |                                                                  |
# |        Compilation Fortran                                       |
# |                                                                  |
# |        Pour faire sortir les fichiers pre-compile avec cpp       |
# |        au lieu de les effacer ajouter -d apres $(FC)             |
# |__________________________________________________________________|


${OBJDIR}/%.o:  %.F90
	$(FC) $(FFLAGS) $(INC) $(INC_COMMON) -o $@ -c $<


#######################################################################
#                         FINAL COMPILATION                           #
#######################################################################

# .__________________________________________________________________.
# |                                                                  |
# |        Link Fortran                                              |
# |__________________________________________________________________|

EXEC  = nsmp3D
EXECE = nsmp3D.exe

# .__________________________________________________________________.
# |                                                                  |
# |       Destination work directory                                 |
# |__________________________________________________________________|

DEST = ./input

# .__________________________________________________________________.
# |                                                                  |
# |        Compilation                                               |
# |__________________________________________________________________|

all:   $(FOBJS)
	   $(FLD) $(LDFLAGS) $(FOBJS) $(LIBS) -o $(EXEC)
	   $(MV) $(EXEC) $(DEST)
	   $(RM) $(FOBJS) $(EXEC) *.mod

# .__________________________________________________________________.
# |                                                                  |
# |                            Cleaning                              |
# |__________________________________________________________________|

clean:
	   $(RM) $(FOBJS) $(EXEC) *.mod

#######################################################################
#                                                                     #
#                               END                                   #
#                                                                     #
#######################################################################

