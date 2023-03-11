# Code for photosynthesis model in Johnson, Field and Berry (2021) Oecologia

Author: Jen Johnson    
Email: jjohnson@carnegiescience.edu    
Last revised: 2021-10-23    
[![DOI](https://zenodo.org/badge/418789869.svg)](https://zenodo.org/badge/latestdoi/418789869)

### Citation:

Johnson, J. E., Field, C.B., Berry, J.A. 2021. The limiting factors and regulatory 
processes that control the environmental responses of C<sub>3</sub>, 
C<sub>3</sub>–C<sub>4</sub> intermediate, and C<sub>4</sub> photosynthesis.
*Oecologia*, DOI: [10.1007/s00442-021-05062-y](https://doi.org/10.1007/s00442-021-05062-y)

### Notes:

1. This code was written in MATLAB R2020b. It is fully compatible with GNU Octave 6.2.0
for the C<sub>3</sub> implementation of the model. However, the symbolic code used for the 
bundle sheath calculations in the C<sub>3</sub>-C<sub>4</sub> and C<sub>4</sub> 
implementations is not yet compatible with Octave. 

2. This directory includes three example simulations that are described below. These call 
functions in the subdirectory `scripts` and write output to the subdirectory `outputs`.

3. The scripts `run_forward_example1.m`, `run_forward_example2.m`, and 
`run_forward_example3.m` simulate the responses of photosynthesis to light, carbon 
dioxide, and temperature, respectively. Other environmental variables are held constant. 

4. Each simulation can be customized by setting the `pathway_opt` variable to either
 `C3`, `Type-I-C3-C4`, or `NADP-ME-C4`. 

5. The subdirectory `inputs` contains two files with the input data used in Figures 6 and 7. 
For the analysis in Figure 6, the forward model in `scripts/model_fun_c3c4.m` needs to be 
fit to the physiological measurements with an optimization routine. For the analysis in 
Figure 7, the environmental measurements can simply be substituted into the example run 
files to recreate the simulations. 

### Directory structure:
```
├── LICENSE
├── README.md
├── inputs
│   ├── Figure-6.csv
│   └── Figure-7.csv
├── outputs
│   ├── Example1-Light-C3
│   │   ├── Example1-Light-C3-figure.png
│   │   ├── Example1-Light-C3-modelinputs.mat
│   │   └── Example1-Light-C3-modeloutputs.mat
│   ├── Example1-Light-NADP-ME-C4
│   │   ├── Example1-Light-NADP-ME-C4-figure.png
│   │   ├── Example1-Light-NADP-ME-C4-modelinputs.mat
│   │   └── Example1-Light-NADP-ME-C4-modeloutputs.mat
│   ├── Example1-Light-Type-I-C3-C4
│   │   ├── Example1-Light-Type-I-C3-C4-figure.png
│   │   ├── Example1-Light-Type-I-C3-C4-modelinputs.mat
│   │   └── Example1-Light-Type-I-C3-C4-modeloutputs.mat
│   ├── Example2-CO2-C3
│   │   ├── Example2-CO2-C3-figure.png
│   │   ├── Example2-CO2-C3-modelinputs.mat
│   │   └── Example2-CO2-C3-modeloutputs.mat
│   ├── Example2-CO2-NADP-ME-C4
│   │   ├── Example2-CO2-NADP-ME-C4-figure.png
│   │   ├── Example2-CO2-NADP-ME-C4-modelinputs.mat
│   │   └── Example2-CO2-NADP-ME-C4-modeloutputs.mat
│   ├── Example2-CO2-Type-I-C3-C4
│   │   ├── Example2-CO2-Type-I-C3-C4-figure.png
│   │   ├── Example2-CO2-Type-I-C3-C4-modelinputs.mat
│   │   └── Example2-CO2-Type-I-C3-C4-modeloutputs.mat
│   ├── Example3-Temperature-C3
│   │   ├── Example3-Temperature-C3-figure.png
│   │   ├── Example3-Temperature-C3-modelinputs.mat
│   │   └── Example3-Temperature-C3-modeloutputs.mat
│   ├── Example3-Temperature-NADP-ME-C4
│   │   ├── Example3-Temperature-NADP-ME-C4-figure.png
│   │   ├── Example3-Temperature-NADP-ME-C4-modelinputs.mat
│   │   └── Example3-Temperature-NADP-ME-C4-modeloutputs.mat
│   └── Example3-Temperature-Type-I-C3-C4
│       ├── Example3-Temperature-Type-I-C3-C4-figure.png
│       ├── Example3-Temperature-Type-I-C3-C4-modelinputs.mat
│       └── Example3-Temperature-Type-I-C3-C4-modeloutputs.mat
├── run_forward_example1.m
├── run_forward_example2.m
├── run_forward_example3.m
└── scripts
    ├── configure_fun.m
    ├── loadvars_fun.m
    ├── model_fun_c3c4.m
    ├── plotter_forward_fun_c3c4.m
    ├── symsolver_c3c4_fun.m
    └── workspace2struct_fun.m

13 directories, 40 files
```