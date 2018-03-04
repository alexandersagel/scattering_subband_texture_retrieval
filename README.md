# Scattering Subband Texture Retrieval

### 1. Folder Contents

This folder contains all the data and MATLAB code for reproducing the results from the paper [1]  
The following is a description of the folder contents.

UIUC.mat (MATLAB data file):  D2 Patches for the UIUC [2] experiment arranged in a multidimensional array  
VisTex.mat (MATLAB data file):  D1 Patches for the VisTex [3] experiment arranged in a multidimensional array  
bhatt_weibull.m (MATLAB function):  Computes the Bhattacharryya inspired distance measure between two sets of Weibull parameters  
cross_entropy_ggd.m (MATLAB function):  Computes the cross entropy distance measure between two sets of GGD parameters  
fwt2alphabeta.m (MATLAB function):  Performs subband decompositions and calculate the GGD parameters for each subband  
ggd_beta.m (MATLAB function):  Computes the beta parameter of a GGD from a set of Samples  
normalize_scat.m (MATLAB function):   Converts the regular Scattering Transform of a signal to its Normalized Scattering Transform  
run_simulation.m (MATLAB script):  Performs the experiments  
subbandDecomposition.m (MATLAB function):   Decomposes an image into a multiresolution representation via FWT  


### 2. Software Requirements

a) Statistics Toolbox, https://mathworks.com/products/statistics.html  
b) Image Processing Toolbox, https://mathworks.com/products/image.html  
c) ScatNet, http://www.di.ens.fr/data/software/scatnet/  
d) DT-CWT Toolbox, http://www-sigproc.eng.cam.ac.uk/Main/NGK  


### 3. Running the Experiments

To reproduce the tables from the paper, run run_simulation.m. The according parameters can be set in the first section of the script.

### 4. References

[1] A. Sagel, D. Meyer and H. Shen, *Texture Retrieval Using Scattering Coefficients and Probability Product Kernels*, Int. Conference on Latent Variable Analysis and Signal Separation (LVA/ICA), August 2015.  
[2] Svetlana Lazebnik, Cordelia Schmid, and Jean Ponce. *A Sparse Texture Representation Using Local Affine Regions*. IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 27, no. 8, pp. 1265-1278, August 2005.   
[3] http://vismod.media.mit.edu/vismod/imagery/VisionTexture/vistex.html, see also COPYRIGHT_VISTEX.TXT  
