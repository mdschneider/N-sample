# ============================================================
# Makefile for N-sample
# ============================================================

# choose platform
include make.CHARON

# Compiler flags
# Choose from FDEBUG or FOPTIMIZED
FCFLAGS = $(FDEBUG) $(OPTIONS)

# Shared libraries
INC   = $(HDF5_INC) $(HPX_INC) $(FFTW_INC) 
FLIBS = $(HDF5_LIB) $(HPX_LIB) $(FFTW_LIB) 
RLIBS = $(HDF5_RLIB) $(HPX_RLIB) $(FFTW_RLIB) 
LIBS = $(RLIBS) $(FLIBS)

ifeq ($(USEFITPACK), "yes")
FITPACK = fitpack.o
else
FITPACK =
endif

READMILLENNIUM = read_millennium/indexxxi8.o read_millennium/file_path.o read_millennium/readfile.o read_gadget.o

PARAMS = parameters_nth8.o parameters_nth16.o

OBJS = \
types_nsample.o get_free_unit.o $(FITPACK) \
$(PARAMS) \
$(READMILLENNIUM) \
cosmology.o parameters.o \
io_nsample.o fft_nsample.o \
gridtools.o \
zeldovich.o resample.o 

#EXE = addPowerSnapshots.x 
#EXE = addPowerTaylor.x 
EXE = nsample.x powerspectrum.x 

# --- Targets ----------------------------------

sources: ${OBJS}

all: ${OBJS} ${EXE}

clean:
	rm -f *.o *.d *.mod *.x.$(HOSTNAME) *~

# --- Suffix Rules -------------------------------
# First clear out the default suffixes, then declare our own 
# and define the rules

.SUFFIXES:
.SUFFIXES: .F90 .f90 .o .x

# --- Compiling ----------------------------------

%.o: %.f
	@ echo 'Building file: $<'
	@ echo '$(FC) $(FCFLAGS) $(INC) -c -o $@ $<'
	@ $(FC) $(FCFLAGS) $(INC) -c -o $@ $<
	@ echo 'Finished building: $<'
	@ echo

%.o: %.f90
	@ echo 'Building file: $<'
	@ echo '$(FC) $(FCFLAGS) $(INC) -c $< -o $@'
	@ $(FC) $(FCFLAGS) $(INC) -c $< -o $@
	@ echo 'Finished building: $<'
	@ echo ' '

%.o: %.F90
	@ echo 'Building file: $<'
	@ echo '$(FC) $(FCFLAGS) $(INC) -c $< -o $@'
	@ $(FC) $(FCFLAGS) $(INC) -c $< -o $@
	@ echo 'Finished building: $<'
	@ echo ' '

%.x : %.o
	@ echo 'Linking: $@'
	@ echo '$(FC) $(FCFLAGS) $(OBJS) -o $@ $< $(LIBS)'
	@ $(FC) $(FCFLAGS) $(OBJS) -o $@ $< $(LIBS)
	@ mv $@ $@.$(HOSTNAME)$(BINNDX)
	@ echo ' '
