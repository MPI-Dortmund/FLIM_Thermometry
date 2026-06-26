This repository contains exemplary code to generate and fit photon arrival time histograms from time correlated single photon counting (TCSPC) measurements for two TCSPC-raw-data formats (SPC and PTU, 1 sample data file each provided) that were used in [1]. Binary data is imported into MATLAB (currently deployed version R2025b) by the "parsePTUFile.m" and "parseSPCFile.m" functions respectively, as demonstrated in the scripts "fitGlobalHistogramFrom*.m", which also perform fitting and plotting of the sample data.

Files in this repository:  
parsePTUFile.m              -- convert binary PTU-data file of photon arrival times to Matlab structure  
parseSPCFile.m              -- convert binary SPC-data file of photon arrival times to Matlab structure  
tcspcIRF.m                  -- adapted from [2]: 5-parameter IRF-model based on a gaussian with asymmetric flanks  
convTCSPCMonoDecay.m        -- convolution of tcspcIRF.m with a mono-exponential decay, also adapted from [2]  
fitGlobalHistogramFromPTU.m -- sample code to generate a histogram from PTU-data and fit it to 2 different models  
fitGlobalHistogramFromSPC.m -- sample code to generate a histogram from SPC-data and fit it to 2 different models  
SampleFile_01.ptu           -- part of the data used in Fig. 1  
SampleFile_02.spc           -- part of the data used in Fig. 3  

[1] "Robust thermometry-imaging at sub-micrometer and millisecond-resolution by fluorescence lifetime microscopy
    allows for additional acquisition of multiple imaging channels", submitted  
[2] "Precise measurement of protein interacting fractions with fluorescence lifetime imaging microscopy",
    Mol. BioSyst. 2011, 7, 322-336; DOI: 10.1039/c0mb00132e
