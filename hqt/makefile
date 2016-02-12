##########################
# Makefile for libhqt.so #
# Author: Gherardo Vita  #
# E-Mail: vita@mit.edu   #
##########################

FORTRANCOMPILER= gfortran -fno-automatic -fno-f2c -O3

FORT= interface_lha.f

FORT+= autohiggs.f besselk.f hankel.f integration.f evhiggs.f \
      subhiggs.f higgsqt.f vegas.f

INCLUDE= const.F scales.F
FILES= $(MAKEFILES) $(FORT) $(INCLUDE) $(LIBFORT)
OBJ = $(patsubst %.f,%.o,$(FORT))
OBJ += $(patsubst %.F,%.o,$(wildcard ./lib/*.F))

#AUTODBL is a flag for automatic double precision type

%.o: %.f $(INCLUDE)
	$(FORTRANCOMPILER) $(AUTODBL) -fPIC -c $<

%.o: %.F $(INCLUDE)
	$(FORTRANCOMPILER) $(AUTODBL) -fPIC -c $< -o $@

libhqt.so: $(OBJ)
	$(FORTRANCOMPILER) $(AUTODBL) -shared -fPIC -o $@ $^

clean:
	rm -f *.o *~ lib/*.o
