# MAKE FILE to compile SPHDIFF

today_date = `date`
compile_date = $(shell date)

CASE = -fpp -DPremixed

DEBUG ?= 0
ifeq ($(DEBUG), 1)
    COMP_FLAGS = $(CASE) -Ddate='"$(compile_date)"' -fpp -DDEBUG -warn unused -debug full -g -check all -r8  -mcmodel=medium -shared-intel -traceback -ftrapuv -c -o  $@
    LINK_FLAGS = -v -g -debug -check all -traceback -r8 -mcmodel=medium -shared-intel -i-dynamic
else
# COMP_FLAGS = -O3 -r8 -c -mcmodel=medium -shared-intel -o $@
# LINK_FLAGS =  -v -r8 -mcmodel=medium -shared-intel -i-dynamic
#     COMP_FLAGS = -O3 -r8 -c -mcmodel=medium -shared-intel -o $@
#     LINK_FLAGS =  -v -r8 -mcmodel=medium -shared-intel -i-dynamic
    COMP_FLAGS = $(CASE) -Ddate='"$(compile_date)"' -O2 -r8 -c -mcmodel=large -shared-intel -o $@
    LINK_FLAGS = -parallel  -v -r8 -mcmodel=large -shared-intel -i-dynamic
endif

COMPILE_F90 = ifort $(COMP_FLAGS)
LINK        = ifort $(LINK_FLAGS)

#COMPILE_F90 = ifort -debug -g -check bounds -prof-gen -traceback -r8 -c -o  $@
#COMPILE_F90 = ifort -ipo -prof-use -r8 -c -o  $@
#LINK = ifort -v -g -debug -check all -traceback -r8
#LINK = ifort -v -debug -g -check bounds -traceback -r8 -prof-gen
#LINK = ifort -v -r8 -prof-use
#========================================================================================

PROG = Sphdiff.x

#========================================================================================

#TARGETS:
all: executable

clean :; @rm -f *.o *.mod $(PROG)

#----------------------------------------------------------------------------------------
# Defining the object files
OBJ_MAIN = tranlib.o    \
	math.o       \
	dmach.o      \
	cklib.o      \
	fldriv.o     \
 	point.o      \
        prsave.o       \
	gridgen.o       \
	rdkey.o          \
	restrt.o         \
	sphdifb.o        \
	start.o          \
	upcase.o         \
	area.o           \
	optic.o          \
	temp.o           \
	driver.o         \
        inikey.o	\
	cvcm.o		\
	public.o	\
	ddaspk.o       \
       	dlinpk.o	\
	daux.o		\
	publicarray.o	\
	rescvcm.o	\
	cvcminit.o	\
	outcvcm.o	\
	JACCVCM.o	\
	PSOLCVCM.o	\
	drbdpre.o

#----------------------------------------------------------------------------------------
#----------------------  COMPILE SUBROUTINES AND FUNCTIONS  -----------------------------
#----------------------------------------------------------------------------------------

# Compiling source files

sphdifb.o: sphdifb.f
		$(COMPILE_F90) -c sphdifb.f
tranlib.o: tranlib.f
		$(COMPILE_F90) -c tranlib.f
math.o: math.f
		$(COMPILE_F90) -c math.f
dmach.o: dmach.f
		$(COMPILE_F90) -c dmach.f
cklib.o: cklib.f ckstrt.h
		$(COMPILE_F90) -c cklib.f
soot.o: soot.f soot.h
		$(COMPILE_F90) -c soot.f
point.o: point.f
		$(COMPILE_F90) -c point.f
prsave.o: prsave.f
		$(COMPILE_F90) -c prsave.f
rdkey.o: rdkey.f 
		$(COMPILE_F90) -c rdkey.f
restrt.o: restrt.f
		$(COMPILE_F90) -c restrt.f
start.o: start.f 
		$(COMPILE_F90) -c start.f
upcase.o: upcase.f
		$(COMPILE_F90) -c upcase.f
area.o: area.f
		$(COMPILE_F90) -c area.f
optic.o: optic.f
		$(COMPILE_F90) -c optic.f
temp.o: temp.f
		$(COMPILE_F90) -c temp.f
public.o: public.f
		$(COMPILE_F90) -c public.f
driver.o: driver.f 
		$(COMPILE_F90) -c driver.f
cvcm.o: cvcm.f
		$(COMPILE_F90) -c cvcm.f
inikey.o: inikey.f
		$(COMPILE_F90) -c inikey.f
ddaspk.o: ddaspk.f
		$(COMPILE_F90) -c ddaspk.f
daux.o: daux.f
		$(COMPILE_F90) -c daux.f
dlinpk.o: dlinpk.f
		$(COMPILE_F90) -c dlinpk.f
fldriv.o: fldriv.f
		$(COMPILE_F90) -c fldriv.f
publicarray.o: publicarray.f
		$(COMPILE_F90) -c publicarray.f
rescvcm.o: rescvcm.f
		$(COMPILE_F90) -c rescvcm.f
cvcminit.o: cvcminit.f
		$(COMPILE_F90) -c cvcminit.f
outcvcm.o: outcvcm.f
		$(COMPILE_F90) -c outcvcm.f	
JACCVCM.o: JACCVCM.f
		$(COMPILE_F90) -c JACCVCM.f
PSOLCVCM.o: PSOLCVCM.f
		$(COMPILE_F90) -c PSOLCVCM.f
drbdpre.o: drbdpre.f
		$(COMPILE_F90) -c drbdpre.f
#----------------------------------------------------------------------------------------
executable: $(OBJ_MAIN)
	$(LINK) $(OBJ_MAIN) -o $(PROG) 

		echo $(PROG): make complete
