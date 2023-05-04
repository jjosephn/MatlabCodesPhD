# MatlabCodesPhD
Matlab codes and functions used during my PhD.

   ------------------------------------------------------------------------------------------------------------------------------------
   Requirements (Need to download and install separately)
   ------------------------------------------------------------------------------------------------------------------------------------

   Matlab
   cvx toolbox

   KSVD_Matlab_ToolBox
   
   dbListUpdated.xls
   ptb.txt
   
   ------------------------------------------------------------------------------------------------------------------------------------
   PTBDatabaseFull contains whole ptb database.
   
   cvxLinux64 and cvxWindows64 contains cvx tool boxes for linux and windows. cvx tool box needs to be installed depending on the os.

1. Chapter 1: LeadSelectiveMultisacleLinearRegression3To12 : code for running lead multiscale linear regression for whole ptb database.

2. Chapter 2: MultipleCoupledDictionary3To12 : code for running multiple coupled dictionary learning. cvx tool box has to be installed 
   before running this code. It might take some time to finish running as the code is for whole ptb database.
   
3. Chapter 3: dnnAnalysis : code for computing distortion measures after the signals are reconstructed using dnn models. Models are 
   learned in python and the data is saved in .mat format.
   
   ------------------------------------------------------------------------------------------------------------------------------------
   Functions
   ------------------------------------------------------------------------------------------------------------------------------------
   
   1. LPFilter : Low pass filter for baseline removal.
   
   2. wltTfm : computes wavelet transform return a structure of approximation and detailed coefficients.
   
   3. invWltTfm : computes inverse wavelet transform
   
   4. optimizationCVX and optimizationCVXNew : uses cvx tool box for optimization.
   
   5. waveletDist : computes Wavele Energy based Diagnostic Distortion (WEDD).
   
   6. weightedDist : computes Weighted Diagnostic Distortion. Depends upon other functions which are:   
   
      a. featureExtraction and featureExtractionNew : used for extracting shape, amplitude and duration features.
         a.1. PFeatureExtraction, PFeatureExtractionNew : Extracts P wave features.
         a.2. QFeatureExtraction, QFeatureExtractionNew : Extracts Q wave features.
         a.3. TFeatureExtraction, TFeatureExtractionNew : Extracts T wave features.
   
   
