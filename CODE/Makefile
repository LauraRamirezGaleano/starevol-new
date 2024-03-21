####################################################################################
#
# 	         makefile for STAREVOL      
#    
#	        Lionel Siess 01/02/2007
#
# commands :
#
#  make man              : display man page
#
#  make                  : generate execution file $(NAME).e
#  make clean            : remove .o files in the working directory $(DIR_CODE)
#  make cleanall         : complete cleaning (including library .o files)
#  make code             : generate starevolc.f
#  make tar              : concatenate sources files+include needed for starevol
#  make utils            : concatenate utilities needed by STAREVOL
#  make lsR              : create ls-R STAREVOL database
#  make alltar           : combine make tar and make utils
#  make DD=on            : compile with debugging options on
#  make SFLAGS="option"  : specify flag to link with static libraries (by default it is on)
#  make FFLAGS="option"  : specify flag options
#  make FC=compiler      : specify fortran compiler
#  make PROF=on          : compile with profiling options
#  make DIR_EXEC=dirname : the execution file will be written in the dirname directory
#  make NAME=exec        : specify output filename : exec.e instead of $(NAME).e
#  make EXT=extension    : specify the output filename extention : starevol_EXT.e instead of $(NAME).e
#  make ARCH=off         : compile without cpu optimizations
#  make PAR=on	         : compile with parallelization optimizations
#  make var              : display system information    
#  make tex              : compile all the .tex sources  
#####################################################################################

###################
# default variables
###################

#hostname
HOSTMACHINE=$(shell basename ${HOST} .graal.univ-montp2.fr)

# starevol type
stype = evol

#root filename
ROOT_NAME = starevol

# compiler flags
SFLAGS=${FORTRAN_SFLAGS}
FFLAGS=${FORTRAN_FLAGS}

# fortran compiler
FC = ${FORTRAN}

# linker, is the same as the compiler
LD = $(FC)

# debugging option
DD=off

# parallelizing
PAR=on

# define code version
star_version = $(shell grep code_version starevol.f |cut -c 22-25)
star_step = $(subst .,,$(star_version))
#star_banner= $(shell grep ersion starevol.f | grep 200 | cut -c 33-36)
star_banner= $(shell grep "   Version" starevol.f | cut -c 33-36)

# processor type
PROC=$(shell cat /proc/cpuinfo | grep name | tail -1 | perl -n -e "s/(^model.*:(.*Mobile)?)(\s*)(\w*)(.*)/\4/g; print | cut" )

#architecture
BIT=$(shell uname -m | grep 64)
ifeq ($(BIT),)
	BIT=x32
	FC = ${FORTRAN32}
else
	BIT=x64
endif

# # cpu architecture optimization
# ARCH=off
# ifeq ($(FC),ifort)
#   ifeq ($(PROC),Intel)
#     FARCH=-xT
#   else
#     ifeq ($(BIT),x32)
#       FARCH=
#     else
#       FARCH=-xW
#     endif
#   endif
# else
  FARCH=
# endif

# additional links
LINKS=

# extension filename
SUFF = ${PROC}${BIT}

# cpu architecture optimization
ifeq ($(DD),on)
   ARCH=off
   PAR=off
endif

ifeq ($(PAR),on)
   #PARALLEL = -parallel
   PARALLEL=
   SFLAGS=
   #PARALLEL = -openmp -fpp
  # PARLINKS = -openmp -lpthread
   PARLINKS = -qopenmp 
#   PARLINKS = -fopenmp 
#-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack -lguide -openmp -lpthread
endif

ifeq ($(ARCH),on)
   FARCH = -ipo
   ifeq ($(SUFF),Intelx32)
     ifeq ($(FC),f90)
	 FARCH = --tp4
     else
	 FARCH = -xN -ipo
     endif
   endif
   ifeq ($(SUFF),Intelx64)
     ifeq ($(FC),f90)
	 FARCH = --tp4
     else
	 FARCH =-xT -ipo
     endif
   endif
   ifeq ($(PROC),AMD)
     ifeq ($(FC),f90)
	 FARCH = --tp4
     else
	 FARCH =-xW -ipo
     endif
   endif
   SUFF:=opt_${SUFF}
else
   SUFF=${BIT}
endif

# compiler debugging flags
ifeq ($(DD),on)
   FARCH=
   PARALLEL=
   PARLINKS=
   ARCH=off
   ifeq ($(FC),f90)
     FFLAGS=${F90_DB}
   else
     ifeq ($(FC),g95)
       FFLAGS=${G95_DB}
     else
       FFLAGS=${FORTRAN_DBFLAGS}
     endif
   endif
   SUFF:=db_${SUFF}
endif

#extension filename
#EXT = ${SUFF}_g${star_version}
EXT = V${star_version}
 
# execution filename
NAME = ${ROOT_NAME}_${EXT}
NAME_EXEC = $(shell basename ${NAME} .e)

ifeq ($(FC),ifort)
    LDFLAGS = $(FFLAGS)
else
    ifeq ($(FC),g95)
	LDFLAGS =
    else
	LDFLAGS = ${SFLAGS}
    endif
endif

# Profiling options
PROF = 

# Profiling notice
#  1 -  $(NAME).e
#  2 -  gprof ./$(NAME).e gmon.out > log
ifdef PROF
	FFLAGS = -pg
	LDFLAGS = -pg
	SUFF:=${SUFF}_prof
	SFLAGS =
	FARCH =
endif

# exec directory
DIR_EXEC = $(DIR_EXE)/

#code source directory (under $DIR_EXEC)
#DIR_SRC=`echo ${SUFF} | tr "[[:lower:]]" "[[:upper:]]"`
DIR_BASE=$(shell basename ${PWD})

# source directory repository
DIR_SRC = $(DIR_EXE)

# Installation directories (default directories)
IDIR_CODE = $(DIR_SRC)/CODE
IDIR_INCL = $(DIR_SRC)/INCL
IDIR_NUC = $(DIR_SRC)/LIBNUC
IDIR_NUM = $(DIR_SRC)/LIBNUM
IDIR_PACK = $(DIR_SRC)/PACKAGES
IDIR_OPA = $(DIR_SRC)/LIB
IDIR_FLAME = $(DIR_SRC)/FLAME
IDIR_BATCH = $(DIR_BATCH)
IDIR_SMPROG = $(DIR_SMONGO)/PROG
IDIR_SMOBS = $(HOME)/SMONGO/OBS
IDIR_UTILS = $(DIR_PROG)
IDIR_HEADERS = $(DIR_RESEVOL)/HEADERS
IDIR_MAN = $(DIR_MAN)

ifeq ($(DIR_BASE),tmp_CODEDIR)
	RDIR_CODE = ../tmp_CODEDIR
	RDIR_INCL = ../tmp_INCLDIR
	RDIR_NUC = ../tmp_NUCDIR
	RDIR_NUM = ../tmp_NUMDIR
	RDIR_PACK = ../tmp_PACKDIR
	RDIR_OPA = ../tmp_OPADIR
	RDIR_FLAME = ../tmp_FLAMEDIR
	RDIR_HEADERS = ../tmp_HEADERS
	RDIR_BATCH = ../tmp_BATCH
	RDIR_SMPROG = ../tmp_SMPROG
	RDIR_SMOBS = ../tmp_SMOBS
	RDIR_UTILS = ../tmp_UTILS
	RDIR_MAN = ../tmp_MAN
else
	RDIR_CODE=$(IDIR_CODE)
	RDIR_INCL = $(IDIR_INCL)
	RDIR_NUC = $(IDIR_NUC)
	RDIR_NUM = $(IDIR_NUM)
	RDIR_PACK = $(IDIR_PACK)
	RDIR_OPA = $(IDIR_OPA)
	RDIR_FLAME = $(IDIR_FLAME)
	RDIR_HEADERS = $(IDIR_HEADERS)
	RDIR_BATCH = $(IDIR_BATCH)
	RDIR_SMPROG = $(IDIR_SMPROG)
	RDIR_SMOBS = $(IDIR_SMOBS)
	RDIR_UTILS = $(IDIR_UTILS)
	RDIR_MAN = $(IDIR_MAN)
endif

# path to the include files, example: -I../ -I$(HOME)/include
INCL_PATH = -I$(RDIR_INCL)

# define date
stardate = $(shell date +%d"_"%m"_"%Y)

# Use MKL library
#LIB = -mkl

EVO_OBJ =
EVO_SRC =

# Compile with GYRE if decided
ifeq ($(GYRE),on)

   # add include directories containing GYRE modules
   INCL_PATH = -I$(RDIR_INCL) -I$(GYRE_INC) -module $(IDIR_PACK) 

   # link librairies for GYRE
   # for gfortran with Linux use
#   LIB = -L$(GYRE_INC) -lgyre -L$(MESASDK_ROOT)/bin `mesasdk_lapack95_link` `mesasdk_lapack_link` `mesasdk_blas_link` `mesasdk_hdf5_link` -lz
   # for gfortran with MacOS use
#   LIB = -L$(GYRE_INC) -lgyre -Wl,-rpath=$(MESASDK_ROOT) -L$(MESASDK_ROOT)/bin `mesasdk_lapack95_link` `mesasdk_lapack_link` `mesasdk_blas_link` `mesasdk_hdf5_link` -lz
    # for ifort with Linux use
#   LIB += $(GYRE_INC)/gyre_lib.a -L$(MESASDK_ROOT)/bin `mesasdk_lapack95_link` `mesasdk_lapack_link` `mesasdk_blas_link` -lz
   # for ifort with MacOS use
   LIB += $(GYRE_INC)/gyre_lib.a -Wl,-rpath,$(MESASDK_ROOT) -L$(MESASDK_ROOT)/bin `mesasdk_lapack95_link` `mesasdk_lapack_link` `mesasdk_blas_link` -lz
    # for static compilation with ifort use
#   LIB = -static -L$(GYRE_INC) -lgyre -L$(MESASDK_ROOT)/bin `mesasdk_lapack95_link` `mesasdk_lapack_link` `mesasdk_blas_link` -lz -lgfortran -lquadmath

   # list modules that need to be recompiled
   EVO_OBJ = $(RDIR_PACK)/evolstell_gyre.o
   EVO_SRC = $(RDIR_PACK)/evolstell_gyre.f90
   # add compilation flag for preprocessor
   FFLAGS := $(FORTRAN_FLAGS) -DGYRE 
endif


####### include files #######

###ATTENTION#####
## Si on calcule avec refsolar=2 (Asplund) il faut :
## 1) copier evolpar_155.star dans evolpar.star
## 2) utiliser livopaint_new.f et lthopa_155.f
### Sinon, garder evolpar.star et utiliser livopaint_opaext.f


SRC_INCL = eoscom.cons evolcom.const \
	evolpar.star eoscom.fermi evolcom.conv evolcom.var \
	eoscom.ion evolcom.data evolcom.mass eoscom.math \
	evolcom.diff evolcom.spec eoscom.par \
	evolcom.lum evolcom.mod evolcom.surf \
	eoscom.phys evolcom.eng evolcom.acc evolcom.eos \
	evolcom.nuc evolcom.therm evolcom.therm2 evolcom.grad \
	evolcom.opa evolcom.cons evolcom.ion evolcom.rot \
	evolcom.transp evolcom.teq evolcom.flame \
	evolcom.igw evolcom.num evolcom.atm

EVOLINCL = $(patsubst %,$(RDIR_INCL)/%,$(SRC_INCL))


####### source files #######

SRC_STAR = runevol.f starevol_init.f starevol.f accret.f atmtp.f boreas.f calcevo.f centeq.f \
	chedif.f convzone.f denucl.f dfdr.f diffinit.f \
	diffsolve_OK.f diffturbul.f diffusion.f evodat.f fdiff.f freac.f fscreen.f \
	guess.f init.f inteq.f interp_CMFGEN.f interpmesh.f interpZ.f kappa.f \
        load_CMFGEN.f load_TLUSTY.f mchange_Vink.f mesh.f mix.f mlt.f netdiff.f network.f neutri.f nuceng.f \
	nucsolve.f nucprod.f nwracc.f nwtraf.f opa_co.f oversh.f \
	pconv.f prvar.f quintic.f resulpr.f rinimod.f rininet.f rmodpar.f \
	screen.f select_table.f structure.f surfeq.f thermo.f tstep.f \
	tstepx.f vit.f ydequil.f yequil.f yneutron.f

SRC_NUC = nulosc.f vitsub_ce2_qinu.f vitsub_167_Nacre2.f
#SRC_NUC = nulosc.f vitsub_ce2_qinu.f vitsub_167_Nacre1.f

SRC_PACK = corrwind.f diffusc_ondes.f overshootc.f igw.f eosc.f jtranspc_ondesexcit.f rotation.f diffinitondesexcit.f
#endif

SRC_FLAME = myengen.f myflamespeed.f mycompo.f

ifeq ($(SOLCOMP),AGSS09)
	SRC_OPA = livopaint_opaAS09.f lthopa_155.f lthopa_155AS09.f lthopa_120.f
endif
ifeq ($(SOLCOMP),AY18)
	SRC_OPA =  lthopa_155AY18.f livopaint_opaAY18.f
endif
ifeq ($(SOLCOMP),GN93)
        SRC_OPA = livopaint_new.f lthopa_155.f
endif

SRC_NUM = exv.f fit.f nwrmat.f interlin.f shearnwt.f \
	lubksb.f ludcmp.f savgol.f lsame.f xerbla.f idamax.f ilaenv.f \
        iparmq.f ieeeck.f \
        nagmod_ok.f sort.f splineatm.f splintatm.f splineatmb.f splintatmb.f \
	dgetrf.f dgetrs.f dtrsm.f dlaswp.f dgetf2.f \
	dswap.f dscal.f dger.f dgemm.f  
#	gift_54.f leqs.f 
#	##ma28ad.f ma28bc.f mc20ad.f ma30ad.f ma30dd.f \
	##mc20ad.f mc24ad.f mc23ad.f mc22ad.f mc21a.f mc13d.f mc21b.f \
	##mc13e.f mc28cd.f ma28cd.f ma30cd.f


######### object files #############

FOR_STAR = $(SRC_STAR)
FOR_NUM = $(patsubst %,$(RDIR_NUM)/%,$(SRC_NUM))
FOR_NUC = $(patsubst %,$(RDIR_NUC)/%,$(SRC_NUC))
FOR_OPA = $(patsubst %,$(RDIR_OPA)/%,$(SRC_OPA))
FOR_PACK = $(patsubst %,$(RDIR_PACK)/%,$(SRC_PACK))
FOR_FLAME = $(patsubst %,$(RDIR_FLAME)/%,$(SRC_FLAME))

OBJ_STAR = $(FOR_STAR:.f=.o)
OBJ_NUM = $(FOR_NUM:.f=.o)
OBJ_NUC = $(FOR_NUC:.f=.o)
OBJ_OPA = $(FOR_OPA:.f=.o)
OBJ_PACK = $(FOR_PACK:.f=.o)
OBJ_FLAME = $(FOR_FLAME:.f=.o)

OBJ_LIB = $(OBJ_NUM) $(OBJ_NUC) $(OBJ_PACK) $(OBJ_FLAME)
OBJECTS = $(OBJ_STAR) $(OBJ_LIB)

######### listing files #############

LST_OBJS = $(OBJECTS:.o=_mod.f90) $(OBJECTS:.o=_mod.mod)
LST_STAR = $(OBJ_STAR:.o=_mod.f90) $(OBJ_STAR:.o=_mod.mod)
LST_OPA = $(OBJ_OPA:.o=_mod.f90) $(OBJ_OPA:.o=_mod.mod)


######### batch files #############

SRC_BATCH = runevol save nexte batche batchesun cpevol deline gensm fusevol \
        header restart renseq evolpar concatevol cgpar genkipsm \
        jcl jclb rcreate catpagb starevolnuc_${stype}$(star_step).par \
	tcshrc ecut efuse ecutfinal mkdirevol lse esh \
	bash_before bash_after bash_after_sun bash_evolfile


######### program files #############

SRC_UTILS = Make_binlist_${stype} binlist_${stype}.f Make_starevolman \
        starevolman.f Grevesse.dat transfomu.f transfomu.link \
        evoldeline.f evolcut.f evolstat.f \
	GN93.dat AG89.dat GNS96.dat GS98.dat AGS05.dat AGSS09.dat AY18.dat


######### headers files #############

SRC_HEADERS =  mevol.c1.gz mevol.c4.gz mevol.s2.gz mevol.v11.gz  \
	mevol.v2.gz mevol.v5.gz mevol.v8.gz mevol.c2.gz mevol.hr.gz \
	mevol.s3.gz mevol.v12.gz mevol.v3.gz mevol.v6.gz mevol.c3.gz \
	mevol.s1.gz mevol.s4.gz mevol.v1.gz mevol.v4.gz mevol.v7.gz


######### manual files #############

SRC_MAN = starevolpar.tex strucfile.tex evolfile.tex \
	glossaireos.tex conservationJ.tex run_log.ps \
	jtransp_notes.tex matrice_complete.tex \
	remplissage_jacobiens_JTRANSP.tex opacity.pdf \
	pulse_def_properties.eps


######### supermongo files #############

SRC_SMPROG = evol.sm Make_nucsm_${stype} nucsm_${stype}.f colours.sm \
	default .sm starevol_${stype}.par network_${stype}.dat \
	mygensm.pl sm.el sm.el.hlp sm.elc .emacs twoD.sm

SRC_SMOBS = Armosky94_M92M15.dat Brown92_M22.dat CostaFrogel96.dat \
	Guandalini06.dat Kawaler.dat LMC_Cstar.dat opaCO.dat vanLoon99.dat \
	atom.dat Carretta05.dat Gratton2000.dat HarrisDB.dat Koester.dat \
	NGC6752_MgAl.dat Pasquini05.dat Yong03Tab1.dat Briley04.dat \
	CarrettaNGC2808.dat guandalini06.dat Ivans1999.dat Letarte06.dat \
	NGC6752_ONa.dat Starevolnucl.dat

######### tar files #############

TAR_INCL = $(patsubst %,tmp_INCLDIR/%,$(SRC_INCL))
TAR_STAR = $(patsubst %,tmp_CODEDIR/%,$(SRC_STAR))
TAR_NUM = $(patsubst %,tmp_NUMDIR/%,$(SRC_NUM))
TAR_NUC = $(patsubst %,tmp_NUCDIR/%,$(SRC_NUC))
TAR_OPA = $(patsubst %,tmp_OPADIR/%,$(SRC_OPA))
TAR_PACK = $(patsubst %,tmp_PACKDIR/%,$(SRC_PACK))
TAR_FLAME = $(patsubst %,tmp_FLAMEDIR/%,$(SRC_FLAME))
TAR_UTILS = $(patsubst %,tmp_UTILS/%,$(SRC_UTILS))
TAR_BATCH = $(patsubst %,tmp_BATCH/%,$(SRC_BATCH))
TAR_SMPROG = $(patsubst %,tmp_SMPROG/%,$(SRC_SMPROG))
TAR_SMOBS = $(patsubst %,tmp_SMOBS/%,$(SRC_SMOBS))
TAR_HEADERS = $(patsubst %,tmp_HEADERS/%,$(SRC_HEADERS))
TAR_MAN = $(patsubst %,tmp_MAN/%,$(SRC_MAN))

######### untar files #############

UNTAR_INCL = $(patsubst %,../tmp_INCLDIR/%,$(SRC_INCL))
UNTAR_STAR = $(patsubst %,../tmp_CODEDIR/%,$(SRC_STAR))
UNTAR_NUM = $(patsubst %,../tmp_NUMDIR/%,$(SRC_NUM))
UNTAR_NUC = $(patsubst %,../tmp_NUCDIR/%,$(SRC_NUC))
UNTAR_OPA = $(patsubst %,../tmp_OPADIR/%,$(SRC_OPA))
UNTAR_PACK = $(patsubst %,../tmp_PACKDIR/%,$(SRC_PACK))
UNTAR_FLAME = $(patsubst %,../tmp_FLAMEDIR/%,$(SRC_FLAME))
UNTAR_UTILS = $(patsubst %,../tmp_UTILS/%,$(SRC_UTILS))
UNTAR_BATCH = $(patsubst %,../tmp_BATCH/%,$(SRC_BATCH))
UNTAR_SMPROG = $(patsubst %,../tmp_SMPROG/%,$(SRC_SMPROG))
UNTAR_SMOBS = $(patsubst %,../tmp_SMOBS/%,$(SRC_SMOBS))
UNTAR_HEADERS = $(patsubst %,../tmp_HEADERS/%,$(SRC_HEADERS))
UNTAR_MAN = $(patsubst %,../tmp_MAN/%,$(SRC_MAN))


######### lsR files #############

LSR_STAR = $(patsubst %,"\\\n"$(IDIR_CODE)/%,$(SRC_STAR))
LSR_NUM = $(patsubst %,"\\\n"$(IDIR_NUM)/%,$(SRC_NUM))
LSR_NUC = $(patsubst %,"\\\n"$(IDIR_NUC)/%,$(SRC_NUC))
LSR_OPA = $(patsubst %,"\\\n"$(IDIR_OPA)/%,$(SRC_OPA))
LSR_PACK = $(patsubst %,"\\\n"$(IDIR_PACK)/%,$(SRC_PACK))
LSR_FLAME = $(patsubst %,"\\\n"$(IDIR_FLAME)/%,$(SRC_FLAME))
LSR_INCL = $(patsubst %,"\\\n"$(IDIR_INCL)/%,$(SRC_INCL))
LSR_TAR = `echo $(LSR_BAN1) $(IDIR_CODE) $(LSR_BAN1) $(LSR_STAR) $(LSR_BAN) $(IDIR_INCL) $(LSR_BAN1) $(LSR_INCL) \
	$(LSR_BAN) $(IDIR_NUM) $(LSR_BAN1) $(LSR_NUM) $(LSR_BAN) $(IDIR_NUC) $(LSR_BAN1) $(LSR_NUC) \
	$(LSR_BAN) $(IDIR_OPA) $(LSR_BAN1) $(LSR_OPA) $(LSR_BAN) $(IDIR_FLAME) $(LSR_BAN1) $(LSR_FLAME) \
	$(LSR_BAN) $(IDIR_PACK) $(LSR_BAN1) $(LSR_PACK) | sed 's/\/\//\//g'`
LSR_UTILS = $(patsubst %,"\\\n"$(IDIR_UTILS)/%,$(SRC_UTILS))
LSR_BATCH = $(patsubst %,"\\\n"$(IDIR_BATCH)/%,$(SRC_BATCH))
LSR_SMPROG = $(patsubst %,"\\\n"$(IDIR_SMPROG)/%,$(SRC_SMPROG))
LSR_SMOBS = $(patsubst %,"\\\n"$(IDIR_SMOBS)/%,$(SRC_SMOBS))
LSR_HEADERS = $(patsubst %,"\\\n"$(IDIR_HEADERS)/%,$(SRC_HEADERS))
LSR_MAN = $(patsubst %,"\\\n"$(IDIR_MAN)/%,$(SRC_MAN))
LSR_BAN ="\\\n\\\n=====================================\\\n"
LSR_BAN1 ="\\\n=====================================\\\n"
LSR_ALL = `echo $(LSR_BAN1) $(IDIR_CODE) $(LSR_BAN1) $(LSR_STAR) $(LSR_BAN) $(IDIR_INCL) $(LSR_BAN1) $(LSR_INCL) \
	$(LSR_BAN) $(IDIR_NUM) $(LSR_BAN1) $(LSR_NUM) $(LSR_BAN) $(IDIR_NUC) $(LSR_BAN1) $(LSR_NUC) \
	$(LSR_BAN) $(IDIR_OPA) $(LSR_BAN1) $(LSR_OPA) $(LSR_BAN) $(IDIR_FLAME) $(LSR_BAN1) $(LSR_FLAME) \
	$(LSR_BAN) $(IDIR_PACK) $(LSR_BAN1) $(LSR_PACK) $(LSR_BAN) $(IDIR_UTILS) $(LSR_BAN1) $(LSR_UTILS) \
	$(LSR_BAN) $(IDIR_BATCH) $(LSR_BAN1) $(LSR_BATCH) $(LSR_BAN) $(IDIR_SMPROG) $(IDIR_SMOBS) $(LSR_BAN1) \
	$(LSR_SMPROG) $(LSR_SMOBS) $(LSR_BAN) $(IDIR_HEADERS) $(LSR_BAN1) $(LSR_HEADERS) $(LSR_BAN) $(IDIR_MAN) \
	$(LSR_BAN1) $(LSR_MAN) | sed 's/\/\//\//g'`


#########################
# generate execution file
#########################

.f.o:
	$(FC) $(INCL_PATH) -c $(FFLAGS) $(FARCH) $(PARALLEL) $(SFLAGS) -o $@ $<

$(NAME_EXEC).e: $(EVO_OBJ) $(OBJECTS) $(OBJ_OPA)
	$(LD) $(SFLAGS) $(FARCH) -o $(DIR_EXEC)/$(NAME_EXEC).e $(EVO_OBJ) $(OBJECTS) $(OBJ_OPA) $(LDFLAGS) $(PARLINKS) $(LINKS) $(LIB)

$(EVO_OBJ): $(EVO_SRC)
	$(FC) $(INCL_PATH) -c $(FFLAGS) $(EVO_SRC) $(PARALLEL) $(PARLINKS) -o $(EVO_OBJ)

	@echo ; echo "$(DIR_EXEC)/$(NAME_EXEC).e generated" ; echo
ifneq ($(star_version),$(star_banner))
	@echo ; echo "WARNING : banner version ($(star_banner)) <> code version ($(star_version))"; echo
endif
ifeq ($(HOSTMACHINE),astropc2)
	cp $(DIR_EXEC)/$(NAME_EXEC).e /c1/${USER}/EVOL/CODE
	@echo
endif
	@echo "$(FC) $(INCL_PATH) -c $(FFLAGS) $(FARCH) $(PARALLEL) $(SFLAGS) ${NAME_EXEC} : on ${HOSTMACHINE} ${PROC} ${BIT} ARCH=${ARCH} " > last_compil
	@chmod u+x ./last_compil

$(RDIR_OPA)/livopaint.o:  $(EVOLINCL)
	$(FC) $(INCL_PATH) -c $(FFLAGS) $(SFLAGS) $(FARCH) $(RDIR_OPA)/livopaint.f -o $(RDIR_OPA)/livopaint.o

$(OBJECTS) : $(EVOLINCL)

##########
# cleaning
##########

clean:
	rm -f $(OBJECTS) $(LST_OBJS) *.core fort.*

cleanall:
	rm -f $(OBJECTS) $(OBJ_OPA) $(LST_OBJS) $(LST_OPA) *.core fort.* \
	starevol.e paq_* *.bin starevolmess \
	$(IDIR_PACK)/evolstell_gyre.mod $(IDIR_PACK)/evolstell_gyre.o

cleancode:
	rm -f $(SRC_STAR) $(OBJ_STAR) $(LST_STAR)


####################
# generate starevolc
####################

code:
	cat $(SRC_STAR) > $(ROOT_NAME)c_${EXT}_$(stardate).f
	@echo
	@echo generate compact source code $(ROOT_NAME)c_${EXT}_$(stardate).f
	@echo



###################
# generate archives
###################

tar:
	@ln -s $(RDIR_CODE) tmp_CODEDIR
	@ln -s $(RDIR_NUM) tmp_NUMDIR
	@ln -s $(RDIR_NUC) tmp_NUCDIR
	@ln -s $(RDIR_OPA) tmp_OPADIR
	@ln -s $(RDIR_PACK) tmp_PACKDIR
	@ln -s $(RDIR_FLAME) tmp_FLAMEDIR
	@ln -s $(RDIR_INCL) tmp_INCLDIR
	@echo -e $(LSR_TAR) > tmp_CODEDIR/$(ROOT_NAME)_$(star_version).lsR
	@tar --one-file-system -cPhzf $(NAME_EXEC)_$(stardate).tgz $(TAR_STAR) $(TAR_NUM) \
	$(TAR_NUC) $(TAR_OPA) $(TAR_PACK) $(TAR_INCL) $(TAR_FLAME) tmp_CODEDIR/Makefile \
	tmp_CODEDIR/starevol_$(star_version).par tmp_CODEDIR/$(ROOT_NAME)_$(star_version).lsR
	@rm -f tmp_CODEDIR
	@rm -f tmp_NUMDIR
	@rm -f tmp_NUCDIR
	@rm -f tmp_OPADIR
	@rm -f tmp_PACKDIR
	@rm -f tmp_FLAMEDIR
	@rm -f tmp_INCLDIR
	@echo
	@echo generate archive file $(NAME_EXEC)_$(stardate).tgz
	@echo

rarutils:
	@ln -s $(RDIR_HEADERS) tmp_HEADERS
	@ln -s $(RDIR_BATCH) tmp_BATCH
	@ln -s $(RDIR_SMPROG) tmp_SMPROG
	@ln -s $(RDIR_SMOBS) tmp_SMOBS
	@ln -s $(RDIR_UTILS) tmp_UTILS
	@ln -s $(RDIR_MAN) tmp_MAN
	tar --one-file-system -cPhzf starutils_$(stardate).tgz $(TAR_UTILS) $(TAR_BATCH) \
	$(TAR_SMPROG) $(TAR_SMOBS) $(TAR_HEADERS) $(TAR_MAN)
	@rm -f tmp_BATCH
	@rm -f tmp_SMPROG
	@rm -f tmp_SMOBS
	@rm -f tmp_UTILS
	@rm -f tmp_HEADERS
	@rm -f tmp_MAN
	@echo
	@echo generate utilities archive file starutils_$(stardate).tgz
	@echo

alltar:
	@make lsR
	@ln -s $(RDIR_CODE) tmp_CODEDIR
	@ln -s $(RDIR_NUM) tmp_NUMDIR
	@ln -s $(RDIR_NUC) tmp_NUCDIR
	@ln -s $(RDIR_OPA) tmp_OPADIR
	@ln -s $(RDIR_PACK) tmp_PACKDIR
	@ln -s $(RDIR_FLAME) tmp_FLAMEDIR
	@ln -s $(RDIR_INCL) tmp_INCLDIR
	@ln -s $(RDIR_BATCH) tmp_BATCH
	@ln -s $(RDIR_SMPROG) tmp_SMPROG
	@ln -s $(RDIR_SMOBS) tmp_SMOBS
	@ln -s $(RDIR_UTILS) tmp_UTILS
	@ln -s $(RDIR_HEADERS) tmp_HEADERS
	@ln -s $(RDIR_MAN) tmp_MAN
	@tar --one-file-system -cPhzf $(NAME_EXEC)_all_$(stardate).tgz tmp_CODEDIR/Makefile \
	tmp_CODEDIR/starevol_$(star_version).par tmp_CODEDIR/$(ROOT_NAME)_$(star_version).lsR \
	$(TAR_STAR) $(TAR_NUM) $(TAR_NUC) $(TAR_OPA) $(TAR_PACK) $(TAR_FLAME) \
	$(TAR_INCL) $(TAR_UTILS) $(TAR_BATCH) $(TAR_SMPROG) $(TAR_SMOBS) $(TAR_HEADERS) $(TAR_MAN)
	@rm -f tmp_CODEDIR
	@rm -f tmp_NUMDIR
	@rm -f tmp_NUCDIR
	@rm -f tmp_OPADIR
	@rm -f tmp_PACKDIR
	@rm -f tmp_FLAMEDIR
	@rm -f tmp_INCLDIR
	@rm -f tmp_BATCH
	@rm -f tmp_SMPROG
	@rm -f tmp_SMOBS
	@rm -f tmp_UTILS
	@rm -f tmp_HEADERS
	@rm -f tmp_MAN
	@echo
	@echo generate entire file archive $(NAME_EXEC)_all_$(stardate).tgz
	@echo

lsR:
	@echo generate lsR database
	@echo -e $(LSR_ALL) > $(ROOT_NAME)_$(star_version).lsR

untar:
	mv -f $(UNTAR_NUM) $(IDIR_NUM)
	mv -f $(UNTAR_NUC) $(IDIR_NUC)
	mv -f $(UNTAR_OPA) $(IDIR_OPA)
	mv -f $(UNTAR_PACK) $(IDIR_PACK)
	mv -f $(UNTAR_FLAME) $(IDIR_FLAME)
	mv -f $(UNTAR_INCL) $(IDIR_INCL)
	mv -f $(UNTAR_STAR) ../tmp_CODEDIR/Makefile ../tmp_CODEDIR/$(ROOT_NAME)_$(star_version).lsR \
	../tmp_CODEDIR/starevol_$(star_version).par $(IDIR_CODE)
	rmdir ../tmp_NUMDIR
	rmdir ../tmp_NUCDIR
	rmdir ../tmp_OPADIR
	rmdir ../tmp_PACKDIR
	rmdir ../tmp_FLAMEDIR
	rmdir ../tmp_INCLDIR
	rmdir ../tmp_CODEDIR

unalltar:
	@make untar
	mv -f $(UNTAR_BATCH) $(IDIR_BATCH)
	mv -f $(UNTAR_SMPROG) $(IDIR_SMPROG)
	mv -f $(UNTAR_SMOBS) $(IDIR_SMOBS)
	mv -f $(UNTAR_UTILS) $(IDIR_UTILS)
	mv -f $(UNTAR_HEADERS) $(IDIR_HEADERS)
	mv -f $(UNTAR_MAN) $(IDIR_MAN)
	rmdir ../tmp_SMPROG
	rmdir ../tmp_SMOBS
	rmdir ../tmp_UTILS
	rmdir ../tmp_HEADERS
	rmdir ../tmp_BATCH
	rmdir ../tmp_MAN

makeall:
	rm -f *~
	make unalltar
	cd $(DIR_SRC) ; make all

all:
	cd $(DIR_SRC) ; make -f ./Makefile cleanall
	cd $(DIR_SRC) ; make -f ./Makefile
	make prog
	make tex

prog:
	cd $(IDIR_UTILS) ; make -f ./Make_binlist_${stype}
	cd $(IDIR_UTILS) ; make -f ./Make_starevolman
	cd $(IDIR_UTILS) ; $(FC) evoldeline.f -o $(DIR_BATCH)/evoldeline.e
	cd $(IDIR_UTILS) ; $(FC) evolcut.f -o $(DIR_BATCH)/evolcut.e
	cd $(IDIR_UTILS) ; $(FC) evolstat.f -o $(DIR_BATCH)/evolstat.e
	cd $(IDIR_SMPROG) ; make -f ./Make_nucsm_${stype}

tex:
	cd $(DIR_MAN) ; latex starevolpar.tex ; latex starevolpar.tex ; makeindex starevolpar ; latex starevolpar
	cd $(DIR_MAN) ; latex strucfile.tex ; latex strucfile.tex ; makeindex strucfile ; latex strucfile
	cd $(DIR_MAN) ; latex starevolpar.tex ; latex starevolpar.tex ; makeindex starevolpar ; latex starevolpar
	cd $(DIR_MAN) ; latex evolfile.tex ; latex evolfile.tex ; makeindex evolfile ; latex evolfile

mklink:
	ln -s $(IDIR_CODE) ../tmp_CODEDIR
	ln -s $(IDIR_NUM) ../tmp_NUMDIR
	ln -s $(IDIR_NUC) ../tmp_NUCDIR
	ln -s $(IDIR_OPA) ../tmp_OPADIR
	ln -s $(IDIR_PACK) ../tmp_PACKDIR
	ln -s $(IDIR_FLAME) ../tmp_FLAMEDIR
	ln -s $(IDIR_INCL) ../tmp_INCLDIR

allfor:
	cat $(SRC_STAR) $(FOR_NUC) $(FOR_PACK) $(FOR_FLAME) $(FOR_OPA) $(FOR_NUM) > $(NAME_EXEC)_all_src_$(stardate).f

var:
	@echo
	@echo "processor type  (PROC) : ${PROC}"
	@echo "cpu architecture (BIT) : ${BIT}"
	@echo "code version           : ${star_version}"
	@echo "exectution filename    : ${NAME_EXEC}.e"
	@echo "fortran compiler       : ${FC}"
	@echo "FLAGS :"
	@echo "   compilation      (FFLAGS) : ${FFLAGS}"
	@echo "   static           (SFLAGS) : ${SFLAGS}"
	@echo "   architecture     (FARCH)  : ${FARCH}"
	@echo "   debugging        (DD)     : ${DD}"
	@echo "   parallelization  (PAR)    : ${PARLINKS}"
	@echo "   additional links (LINKS)  : $(LINKS)"
	@echo
	@echo "compilation command : ${FC} -c ${FFLAGS} ${FARCH} ${PARALLEL} ${SFLAGS}"
	@echo "linking     command : ${LD} ${LDFLAGS} ${SFLAGS} ${FARCH}"
	@echo


man:
	@echo 
	@echo "      *************"
	@echo "      *  MAN PAGE *"
	@echo "      *************"
	@echo
	@echo "   *** COMPILE ***"
	@echo "make                  : generate execution file $(NAME_EXEC).e"
	@echo "make tex              : compile tex sources"
	@echo "make all              : compile all sources (fortran & latex)"
	@echo
	@echo "   *** CLEANING ***"
	@echo "make clean            : remove .o files in the working directory $(RDIR_CODE)"
	@echo "make cleanall         : complete cleaning (including library .o files)"
	@echo "make cleancode        : remove starevolc.f (in $(RDIR_CODE))"
	@echo
	@echo "   *** CONCATENATE SOURCES ***"
	@echo "make code             : generate starevolc.f"
	@echo "make allfor           : combine all fortran sources"
	@echo "make lsR              : create ls-R STAREVOL database"
	@echo
	@echo "   *** ARCHIVES ***"
	@echo "make tar              : concatenate sources files+include needed for starevol"
	@echo "make tarutils         : concatenate utilities needed by STAREVOL"
	@echo "make alltar           : combine make tar and make utils"
	@echo
	@echo "   *** INSTALL ***"
	@echo "make untar            : install code (output of make tar)"
	@echo "make unalltar         : install everything (output of make alltar)"
	@echo "make makeall          : install everything and compile everything"
	@echo "make prog             : install the programs (in ${IDIR_UTILS} and ${IDIR_SMPROG})"
	@echo
	@echo "   *** FLAGS ***"
	@echo "make var              : display the current (machine dependant) option"
	@echo "make NAME=exec        : specify output filename : exec.e instead of $(NAME_EXEC).e"
	@echo "make DD=on            : compile with debugging options on"
	@echo "make PAR=on           : compile with parallelization options on"
	@echo "make SFLAGS=option    : specify flag to link with static libraries (by default it is on)"
	@echo "make FFLAGS=\"option\"  : specify flag options"
	@echo "make FC=compiler      : specify fortran compiler"
	@echo "make PROF=on          : compile with profiling options on, may need to do 'make cleanall' before"
	@echo "                         1 >  run $(NAME).e"
	@echo "                         2 >  gprof ./$(NAME).e gmon.out > file.out, where file.out contains all the profiling information"
	@echo "make ARCH=on          : compile with cpu optimizations"
	@echo "make DIR_EXEC=dirname : the execution file will be written in the \"dirname\" directory"
	@echo "make EXT=extension    : specify the output filename extension : starevol_\"EXT\".e instead of ${ROOT_NAME}_\"$(EXT)\".e"
	@echo
	@echo "    example :          make FC=f77 DIR_EXEC=/tmp FFLAGS=\"-O -v\" EXT=gag_demo"
	@echo
	@echo " will generate the execution file /tmp/starevol_gag_demo.e using the f77 compiler with the compilation option -O -v"
	@echo
	@echo "make var              : display system information    "
	@echo "make tex              : compile all the .tex sources  "
	@echo


#########
# dlink : dgefa,dgesl,dgbfa,dgbsl,daxpy,dcopy,dscal
# jtransp : needs lubksb,ludcmp
# daux : xerrwd,xsetf,xsetun
# dgbsv : dgbsv,dgbtrf,dgbtrs,xerbla,dgbtf2,dger,dlaswp,dtrsm,dgemm,dtbsv,
#          dgemv + dswap,dcopy,dscal,idamax which have been removed
# livopaint : contains spline and splint

