#
# Makefile for virtual widgets and example programs
#
# Copyright 1999 by Larry Charlton/Michael Sweet.
#
# Permission is granted to use, copy, modify, and redistribute this file as
# needed.
#

include ../makeinclude

#
# Make all targets...
#

all:	libflvw.a

#
# Remove object and target files...
#

clean:
	@echo Removing temporary files
	@$(RM) -f *.o core *~ ../FL/*~ ../FL/*.bak

distclean: clean
	@$(RM) -f ../lib/libflvw.a makedepend

install: libflvw.a
	@echo "Installing libraries..."
	@-mkdir -p $(libdir)
	@rm -f $(libdir)/$(LIBRARY)
	@-cp ../lib/$(LIBRARY).a $(libdir)
	@-chmod 644 $(libdir)/$(LIBRARY).a
	@echo "Installing include files..."
	@-mkdir -p $(includedir)/FL
	@-cp -r ../FL/* $(includedir)/FL/

#
# Make the flvw widget library.
#

FLVW_OBJS	= Flv_Style.o Flv_List.o Flv_Table.o Flv_CStyle.o \
		  Flve_Input.o Flve_Check_Button.o Flvt_Edit_Cell.o \
		  Flvt_Edit.o Flvt_Editor.o

FLVW_HDRS:	FL/Flv_Style.H FL/Flv_List.H FL/Flv_Table.H FL/Flv_CStyle.H \
		FL/Flve_Input.H FL/Flve_Check_Button.H FL/Flvt_Edit_Cell.H \
		FL/Flvt_Edit.H

libflvw.a:	$(FLVW_OBJS) $(FLVW_HDRS)
	@echo Building library $@...
	@$(RM) -f ../lib/libflvw.a
	@$(AR) $(ARFLAGS) ../lib/libflvw.a $(FLVW_OBJS)
	@$(RANLIB) ../lib/libflvw.a


#
# End of Makefile.
#
