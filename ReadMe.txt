
ADMM: contains .c code and mex-object of ADMM algorithm for constraint-lasso problem.

construction_functions: contains functions for evaluate SH, SN basis on vertex; compute 
			transition matrix C; Response function matrix.

MEALPix: obtain quadrature points.

toolbox_wavelet_meshes: obtain equal angle grid and plot_spherical_functions (modified).

help_functions: help functions

DWI_generated.m: function that generates DWI data.

SH_ridge.m: apply SH-ridge estimation to DWI data.

superCSD.m: apply superCSD estimation to DWI data with SH_ridge results as input.

SN_lasso.m apply SN-lasso estimation to DWI data.

Examples: contains code for simulation example and read data example

NIfTI: package that loads .nii data into matlab.