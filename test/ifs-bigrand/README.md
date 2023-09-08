To produce the timings in the paper ecRad was first compiled using one of the following (in this example, with the Intel compiler):

**Full optimizations ("OPT3")** when using 32-term ecCKD (compile-time NG_LW, NG_SW must must match gas optics model!):

`make PROFILE=aa_intel  SINGLE_PRECISION=1 GPTL_TIMING=1 OPTIM_CODE=2 NG_SW=32 NG_LW=32`

For RRMTG:
`make PROFILE=aa_intel  SINGLE_PRECISION=1 GPTL_TIMING=1 OPTIM_CODE=2 NG_SW=112 NG_LW=140`

No compile-time ng ("OPT2"):
`make PROFILE=aa_intel  SINGLE_PRECISION=1 GPTL_TIMING=1 OPTIM_CODE=2`

Without optimized main kernels ("OPT1"):
`make PROFILE=aa_intel  SINGLE_PRECISION=1 GPTL_TIMING=1 OPTIM_CODE=1`

**Reference non-optimized ecRad**:
`make PROFILE=aa_intel  SINGLE_PRECISION=1 GPTL_TIMING=1`

Timings were then produced with run_5_tests.sh, changing the command line argument to correspond to different ecRad configurations, for instance:
 
`sh run_5_tests.sh test_big_ecckd_tc`

..for ecCKD+TripleClouds. See Makefile for the different timing tests available. 

Each ecRad run produces a GPTL timing file (timing.*). For each configuration the fastest of the 5 runs was saved
