# Compiler and Flags Configuration
FC       := gfortran
FCFLAGS  := -ffpe-summary=none

.PHONY: all
all: Nexp_GLSE Nexp_dwelltime_simulator clean

Nexp_GLSE: src/Nexp_GLSE/Nexp_GLSE.f src/topfit/topfit_marquardt.f src/topfit/lincomb_Poisson_weighting.f src/topfit/glob_par_map_storage.f src/math_lib/chi2_full_set.f src/math_lib/lineqsys.f src/tools/timer.f
	$(FC) $(FCFLAGS) -c src/Nexp_GLSE/Nexp_GLSE.f
	$(FC) $(FCFLAGS) -c src/topfit/topfit_marquardt.f
	$(FC) $(FCFLAGS) -c src/topfit/lincomb_Poisson_weighting.f
	$(FC) $(FCFLAGS) -c src/topfit/glob_par_map_storage.f
	$(FC) $(FCFLAGS) -c src/math_lib/chi2_full_set.f
	$(FC) $(FCFLAGS) -c src/math_lib/lineqsys.f
	$(FC) $(FCFLAGS) -c src/tools/timer.f
	$(FC) $(FCFLAGS) Nexp_GLSE.o topfit_marquardt.o lincomb_Poisson_weighting.o glob_par_map_storage.o chi2_full_set.o lineqsys.o timer.o -o Nexp_GLSE

Nexp_dwelltime_simulator: src/Nexp_GLSE/Nexp_dwelltime_simulator.f
	$(FC) $(FCFLAGS) src/Nexp_GLSE/Nexp_dwelltime_simulator.f -o Nexp_dwelltime_simulator

.PHONY: clean
clean:
	rm Nexp_GLSE.o topfit_marquardt.o lincomb_Poisson_weighting.o glob_par_map_storage.o chi2_full_set.o lineqsys.o timer.o
