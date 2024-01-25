# Timing tests using either 10,000 profiles (with NREPEAT=1, suitable for single-thread run) or
# 40,000 profiles multi-threaded and repeated 10 times (NREPEAT=10, suitable for multi-thread run)
# To get the input file run:
# wget https://zenodo.org/record/7852526/files/ecrad_highres_40k.nc?download=1
# and provide the location below
# Then just run the test e.g.: OMP_NUM_THREADS=64 make test_big_ecckd_tc

INPUT_BIG = ecrad_highres_40k.nc
# INPUT_BIG = /perm/papu/ecrad_inputs/ecrad_highres_40k.nc
# INPUT_BIG = /perm/papu/ecrad_inputs/ecrad_highres_10k.nc
NREPEAT_TIMING = 1
# NREPEAT_TIMING = 10

DRIVER = ../../bin/ecrad
CHANGENAM = ../common/change_namelist.sh

# Various configuration files corresponding to cycles of ECMWF's
# Integrated Forecasting System (IFS)
CONFIG = ../ifs/configCY47R3.nam
#CONFIG = ../ifs/configCY49R1.nam
CONFIG_ECCKD = ../ifs/configCY49R1_ecckd.nam
CONFIG_RRTMGP = ../ifs/configRRTMGP.nam

all: test_big

test_big: test_big_default test_big_tripleclouds test_big_spartacus test_big_ecckd_tc \
	test_big_ecckd_spartacus

# TIMING runs
# RRTMG + McICA
test_big_default:
	$(CHANGENAM) $(CONFIG) config_timing_default.nam \
		do_save_spectral_flux="false" block_derived_types="true" nrepeat="$(NREPEAT_TIMING)" 
	$(DRIVER) config_timing_default.nam $(INPUT_BIG) tmpout.nc

# RRTMG + TripleClouds
test_big_tripleclouds:
	$(CHANGENAM) $(CONFIG) config_timing_tc.nam \
		sw_solver_name='"Tripleclouds"' lw_solver_name='"Tripleclouds"' \
		do_save_spectral_flux="false" block_derived_types="true" nrepeat="$(NREPEAT_TIMING)"
	$(DRIVER) config_timing_tc.nam $(INPUT_BIG) tmpout.nc

# RRTMG + SPARTACUS
test_big_spartacus:
	$(CHANGENAM) $(CONFIG) config_timing_spartacus.nam \
		sw_solver_name='"SPARTACUS"' lw_solver_name='"SPARTACUS"' \
		do_3d_effects="true" \
		do_sw_delta_scaling_with_gases=""false"" \
		do_save_spectral_flux="false" block_derived_types="true" nrepeat="$(NREPEAT_TIMING)"
	$(DRIVER) config_timing_spartacus.nam $(INPUT_BIG) tmpout.nc

# ecCKD + TripleClouds
test_big_ecckd_tc:
	$(CHANGENAM) $(CONFIG_ECCKD) config_timing_ecckd_tc.nam \
		do_save_spectral_flux="false" block_derived_types="true" nrepeat="$(NREPEAT_TIMING)" 
	$(DRIVER) config_timing_ecckd_tc.nam $(INPUT_BIG) tmpout.nc

# ecCKD + SPARTACUS
test_big_ecckd_spartacus:
	$(CHANGENAM) $(CONFIG_ECCKD) config_timing_ecckd_spartacus.nam \
		sw_solver_name='"SPARTACUS"' lw_solver_name='"SPARTACUS"' \
		do_3d_effects="true" \
		do_save_spectral_flux="false" block_derived_types="true" nrepeat="$(NREPEAT_TIMING)" 
	$(DRIVER) config_timing_ecckd_spartacus.nam $(INPUT_BIG) tmpout.nc

# Clean data files
clean:
	rm -f *_out.nc tmp*.nc radiative_properties*.nc inputs.nc sw_diag*nc \
	config_*.nam gmon.out
