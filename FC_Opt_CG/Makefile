# no implicit rules
.SUFFIXES:

# definitions
FC=gfortran
COMPILE.f08 = $(FC) $(FCFLAGS) $(TARGET_ARCH) -c

# sources for the program
SOURCES=dummy_driver.f90 geom_parameters.f90 C60_mod.f90

# main and clean targets
dostuff: $(subst .f90,.o,$(SOURCES))
	$(FC) -o $@ $+

.PHONY: clean
clean:
	-rm -f *.o *.mod *.smod main

# compilation rules
%.o %.mod %.smod: %.f90
	$(COMPILE.f08) -o $*.o $<
	@touch $@

## .o -> .mod of the modules it uses (when submodules present in modules)
#main.o: one.mod
#one@proc.o: two.mod
#two@proc.o: one.mod
