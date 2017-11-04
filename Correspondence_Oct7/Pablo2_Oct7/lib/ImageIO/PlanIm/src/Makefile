
.SUFFIXES:	.cxx .h .o

.cxx.o :
	$(CC) $(CFLAGS) -c $<

#
#  BINARY is the name of the thing we are trying to build
#
BINARY = libplanio.a
#
#  INSOPTS is the list of options to install
#
INSOPTS = 
#
#  MAKEFILE is the name of this Makefile - I like makefile, others
#  like Makefile.
#
MAKEFILE = Makefile
#
#  DEPEND chooses the type of dependence munging that will be
#  performed when you say 'make depend'.  See notes further
#  downstream.
#
DEPEND = depend43
#
#  RCSDIR is the absolute pathname of the directory where RCS files
#  are to be found.
#
RCSDIR = ./RCS
#
#  INSDIR is the absolute path name of the directory in which the
#  binary is to be installed
#
ARCH = `arch -k`
INSDIR = $(UP)/lib/$(ARCH)
#
#  CFLAGS are options to the c compiler
#
CFLAGS = $(HFLAGS) -I$(UP)/src/include -DUP=\"$(UP)\" -DPHYS_DAT_DIR=\"$(UP)/phys_dat\"
#
#
#
#  LINTFLAGS are options to lint
#
LINTFLAGS = $(CFLAGS) -abhx
#
#  OBJS is a list of .o (and possibly other) files on which the binary
#  depends
#


LOBJS =	get_accessory.o \
	get_filter.o \
	get_lantis_data.o \
	get_phys_dat_dir.o \
	get_plunc_root_dir.o \
	get_slice.o \
	get_time_calc.o \
	get_tray.o \
	get_unit.o \
	plan_xdr_defs.o \
 	read_2d_grid.o \
	read_accessory.o \
	read_anastruct.o \
	read_beam.o \
	read_byte.o \
	read_contour.o \
	read_filter.o \
	read_grid.o \
	read_image_header.o \
	read_implant.o \
	read_many.o \
	read_map.o \
	read_norm.o \
	read_object.o \
	read_objects.o \
	read_plan.o \
	read_points.o \
	read_sar.o \
	read_scan.o \
	read_time_calc.o \
	read_tray.o \
	read_unit.o \
	read_weights.o \
	realloc.o \
	runlen.o \
	sock_utils.o \
	write_implant.o \
	write_object.o \
	write_objects.o \
	xdr_fd.o \
	xdr_ll_planio.o 

#
#  LIBS is a list of libraries to be handed to ld (in -lfoo form)
#
LIBS = 
#
#  SRCS is a list of .c files for use by 'make depend'
#
LSRCS =	get_accessory.cxx \
	get_filter.cxx \
	get_lantis_data.cxx \
	get_plunc_root_dir.cxx \
	get_phys_dat_dir.cxx \
	get_time_calc.cxx \
	get_slice.cxx \
	get_tray.cxx \
	get_unit.cxx \
	plan_xdr_defs.cxx \
	read_2d_grid.cxx \
	read_accessory.cxx \
	read_anastruct.cxx \
	read_beam.cxx \
	read_byte.cxx \
	read_contour.cxx \
	read_filter.cxx \
	read_grid.cxx \
	read_image_header.cxx \
	read_implant.cxx \
	read_many.cxx \
	read_map.cxx \
	read_norm.cxx \
	read_object.cxx \
	read_objects.cxx \
	read_plan.cxx \
	read_points.cxx \
	read_sar.cxx \
	read_scan.cxx \
	read_time_calc.cxx \
	read_tray.cxx \
	read_unit.cxx \
	read_weights.cxx \
	realloc.cxx \
	runlen.cxx \
	sock_utils.cxx \
	write_implant.cxx \
	write_object.cxx \
	write_objects.cxx \
	xdr_fd.cxx \
	xdr_ll_planio.cxx

#
#  Defines for little rpc helper program
#
ROBJS = rpclib_name.o
RSRCS = rpclib_name.c

#
#  Composite lists
#
#OBJS = $(LOBJS) $(ROBJS)
#SRCS = $(LSRCS) $(RSRCS)
OBJS = $(LOBJS)
SRCS = $(LSRCS)

#
#  HDRS is a list of .h (or other) files that need to be checked out
#  to 'make depend'
#
HDRS =

$(BINARY): $(OBJS)
	-rm -f $(BINARY)
#	ld -r -o $(BINARY) $(OBJS)
	ar rvc $(BINARY) $(OBJS)
	-ar ts $(BINARY)

rpclib_name:	$(ROBJS)
	$(CC) $(CFLAGS) -o rpclib_name $(ROBJS)

$(OBJS):

unlink:
	-rm -f ./*,v
	-set -x; for i in `cat \#colist`; do \
		-rm -f $$i; \
	done
	-rm -f \#colist

link:
#
#  Turn this on to protect against the stupidity of 'make'ing in the
#  RCS directory.
#	@echo RCS files in the current directory are being moved to ./save
#	-mkdir ./save
#	-mv -f ./*,v ./save
	-rm -f ./*,v
	-ln $(RCSDIR)/* .
	cp /dev/null \#colist

#
#  There are two scripts for depend mashing provided here.  They are
#  lifted from the 4.2 and 4.3 kernel Makefiles as indicated.  The 4.3
#  version is the preferred choice but it requires the new-with-4.3 -M
#  flag to the C compiler.  Choose the one you want by defining DEPEND
#
depend: $(HDRS) $(SRCS) $(DEPEND)

depend42:
	grep '^#include' ${SRCS} | grep -v '<' | \
	sed -e 's/:[^"]*"\([^"]*\)".*/: \1/' \
	    -e 's/\.c/.o/' \
	    -e 's/\.s:/.o:/' \
	    -e 's,../[a-zA-Z]*/,,' >xtmp
	awk '$$2~/\/data\// {printf "%s ",$$2}' xtmp >xtmp1
	-if [ -s xtmp1 ]; \
	then \
		grep '^#include' `echo ``cat xtmp1``` | grep -v '<' | \
		sed -e 's/:[^"]*"\([^"]*\)".*/: \1/' \
	    	    -e 's/_data\.c/.o/' \
	    	    -e 's,../[a-zA-Z]*/,,' >>xtmp ; \
	else \
		true; \
	fi
	-rm -f xtmp1
	awk ' { if ($$1 != prev) { print rec; rec = $$0; prev = $$1; } \
		else { if (length(rec $$2) > 78) { print rec; rec = $$0; } \
		       else rec = rec " " $$2 } } \
	      END { print rec } ' > makedep <xtmp
	-rm -f xtmp
	echo '/^# DO NOT DELETE THIS LINE/+1,$$d' >eddep
	echo '$$r makedep' >>eddep
	echo 'w' >>eddep
	cp ${MAKEFILE} ${MAKEFILE}.bak
	ed - ${MAKEFILE} < eddep
	rm eddep makedep

depend43: 
	${CC} -M ${CFLAGS} ${SRCS} | \
	awk ' { if ($$1 != prev) { print rec; rec = $$0; prev = $$1; } \
	        else { if (length(rec $$2) > 78) { print rec; rec = $$0; } \
		       else rec = rec " " $$2 } } \
	END { print rec } ' > makedep
	echo '/^# DO NOT DELETE THIS LINE/+1,$$d' >eddep
	echo '$$r makedep' >>eddep
	echo 'w' >>eddep
	cp ${MAKEFILE} ${MAKEFILE}.bak
	ex - ${MAKEFILE} < eddep
	rm eddep makedep

undepend: 
	echo ' ' > makedep
	echo '/^# DO NOT DELETE THIS LINE/+1,$$d' >eddep
	echo '$$r makedep' >>eddep
	echo 'w' >>eddep
	cp ${MAKEFILE} ${MAKEFILE}.bak
	ex - ${MAKEFILE} < eddep
	rm eddep makedep

install: $(BINARY)
	-rm -f $(INSDIR)/$(BINARY)
	mv $(BINARY) $(INSDIR)
	cp libplanio.h $(UP)/src/include/

lint: link $(HDRS) $(SRCS) llint unlink

llint:
	lint $(LINTFLAGS) $(SRCS)

tags: link $(HDRS) $(SRCS) ttags unlink

ttags:
	etags $(HDRS) $(SRCS)

clean:
	-rm -f $(BINARY) *.o *~ core ./*,v \#colist

# DO NOT DELETE THIS LINE - make depend depends on it
 
