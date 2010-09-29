# Makefile for GIPAW

include ../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../Modules \
          $(MOD_FLAG)../EE $(MOD_FLAG)../PW $(MOD_FLAG).

GIPAWOBJS = \
gipaw_module.o \
gipaw_main.o \
paw_gipaw.o \
stop_code.o \
greenfunction.o \
cgsolve_all.o \
h_psiq.o \
cg_psi.o \
symmetrize_field.o \
ch_psi_all.o \
compute_u_kq.o \
suscept_crystal.o \
j_para.o \
biot_savart.o \
init_gipaw_1.o \
init_gipaw_2.o \
init_us_2_no_phase.o \
g_tensor_crystal.o \
write_tensor_field.o \
velocity.o \
nmr_routines.o \
epr_routines.o \
efg.o hyperfine.o core_relax.o \
util.o ascheq.o

QEMODS = ../Modules/libqemod.a
EEOBJS=../EE/libee.a
PWOBJS = ../PW/libpw.a

TLDEPS=bindir pw mods libs libiotk

all : tldeps gipaw.x

gipaw.x : $(GIPAWOBJS) $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ $(GIPAWOBJS) \
		$(PWOBJS) $(EEOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin; ln -fs ../QE-GIPAW/$@ . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f gipaw.x *.o *.F90 *.d *.mod *.i *.L

include make.depend

# DO NOT DELETE

