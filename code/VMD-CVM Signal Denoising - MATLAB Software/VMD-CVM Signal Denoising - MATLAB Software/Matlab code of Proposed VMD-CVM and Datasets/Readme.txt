This matlab code implements a non-parametric signal denoising method currently under review for publication in Signal Processing Journal
Following are the details of the code

1. VMD_CVM_Main: This is the main file from where the whole code can be executed to perform denoising
2. Prop_VMD_CVM: Implements the the proposed denoising algorithm.
3. VMD.m: Performs variational mode decomposition of the noisy signal.
4. cvm.m: Implements the Cramer Von Mises statistic of the selected segment of the data
5. threshvspfa: Estimates threshold for a given Pfa using the rejected noise modes. 
6. cdf_calc.m: matlab code to compute EDF
7. mse.m: Compute mean squared error for the denoised signal
8. snr.m: Compute signal to noise ratio for the denoised signal
9. sofar.mat, taichi.mat and ECGsig.mat contain the real world signals, of the same name, used in this study.