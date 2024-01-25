To produce the timings in the paper ecRad was first compiled with or without optimizattions (in this example, with the Intel compiler):

**Full optimizations ("OPT3")** when using 32-term ecCKD (compile-time NG_LW, NG_SW must must match gas optics model!):

`make PROFILE=aa_intel  SINGLE_PRECISION=1 GPTL_TIMING=1 OPTIM_CODE=2 NG_SW=32 NG_LW=32`

**Reference non-optimized ecRad**:
`make PROFILE=aa_intel  SINGLE_PRECISION=1 GPTL_TIMING=1`

Timings were then produced using run_5_tests.sh. First, [download the larger input file from Zenodo]( https://zenodo.org/record/7852526/files/ecrad_highres_40k.nc?download=1) into this folder or set the directory in Makefile. Then change the command line argument to try different ecRad configurations, for instance :
 
`sh run_5_tests.sh test_big_ecckd_tc`

..for ecCKD + TripleClouds, which is of main interest. See Makefile for the different ecRad tests available. 

Each ecRad run produces a GPTL timing file (timing.*) if you had compiled with GPTL_TIMING=1 or =2. For each configuration the fastest of the 5 runs was saved


**Other compilation options:**

Full optimizations with RRMTG gas optics (compile-time ng):
`make PROFILE=aa_intel  SINGLE_PRECISION=1 GPTL_TIMING=1 OPTIM_CODE=2 NG_SW=112 NG_LW=140`

No compile-time ng ("OPT2"):
`make PROFILE=aa_intel  SINGLE_PRECISION=1 GPTL_TIMING=1 OPTIM_CODE=2`

Without optimized main kernels ("OPT1"):
`make PROFILE=aa_intel  SINGLE_PRECISION=1 GPTL_TIMING=1 OPTIM_CODE=1`
