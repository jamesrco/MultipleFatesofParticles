# MultipleFatesofParticles
Contains code and (some) data for:

[Collins, J. R., B. R. Edwards, K. Thamatrakoln, J. E. Ossolinski, G. R. DiTullio, K. D. Bidle, S. C. Doney, and B. A. S. Van Mooy. 2015. The multiple fates of sinking particles in the North Atlantic Ocean. Global Biogeochem. Cycles 29: 1471-1494, doi:10.1002/2014GB005037.](http://dx.doi.org/10.1002/2014GB005037)

and 

[Laber. C. P., J. E. Hunter, F. Carvalho, E. J. Hunter, B. M. Schieler, E. Boss, K. More, M. Frada, K. Thamatrakoln, C. M. Brown, L. Haramaty, J. E. Ossolinski, H. F. Fredricks, J. I. Nissimov, R. Vandzura, U. Sheyn, Y. Lehahn, R. J. Chant, A. M. Martins, M. J. L. Coolen, A. Vardi, G. R. DiTullio, B. A. S. Van Mooy, and K. D. Bidle. 2018. Coccolithovirus facilitation of carbon export in the North Atlantic. Nature Microbiology, doi:10.1038/s41564-018-0128-4.](https://doi.org/10.1038/s41564-018-0128-4)

The supporting data files for these scripts (in [SinkingParticles_GBC2015/data](https://github.com/jamesrco/SinkingParticles_GBC2015/tree/master/data)) contain data from the two cruises in various, somewhat disorganized formats. Final data from the cruises, with appropriate metadata, are archived to: http://www.bco-dmo.org/deployment/58787 (KN207-1 cruise data) and http://www.bco-dmo.org/project/2136 (KN207-3 cruise data)

The repo [3H_Leu_BactProd](https://github.com/jamesrco/3H_Leu_BactProd) contains all code for analysis of the bacterial production data.

**Caveat:** Unlike the code I currently have in other repos, the scripts here aren't spiffed up for easy portability. User beware!

Brief descriptions of what's here:

1. [respflux2012_v3.m](https://github.com/jamesrco/SinkingParticles_GBC2015/blob/master/respflux2012_v3.m): MATLAB script used for the sensitivity analyses of the particle flux model attenuation described in the paper. Data inputs for each cruise process station are specified separately, then a common block of code is run. Produces the plots in Fig. 5 of the paper.
   * This script requires quite a few of the data files in [SinkingParticles_GBC2015/data](https://github.com/jamesrco/SinkingParticles_GBC2015/tree/master/data), which are described below, in addition to [Shipcast_metadata_KN207-1.xlsx](https://github.com/jamesrco/3H_Leu_BactProd/blob/master/sample_data_metadata/Shipcast_metadata_KN207-1.xlsx) and [Shipcast_metadata_KN207-3.xlsx](https://github.com/jamesrco/3H_Leu_BactProd/blob/master/sample_data_metadata/Shipcast_metadata_KN207-3.xlsx)

2. [Depth_int_WCR_BP.m](https://github.com/jamesrco/SinkingParticles_GBC2015/blob/master/Depth_int_WCR_BP.m): Calculates depth-integrated rates of water column respiration and bacterial production, along with their associated uncertainties, using a Monte Carlo simulation. The simulation is fed by three kinds of inputs: (1) observed rates of BP and WCR and their uncertainites for 0-50 m depth and at 150 m, (2) estimates of rates between 50 and 150 m interpolated using a power law function, and (3) a range of values for several BP conversion factors, reported from a literaure survey in [Giering et al. 2014](http://www.nature.com/articles/nature13123).
   * The necessary observed data are contained in [Depth-T-BP-Resp Correlations?.xlsx](https://github.com/jamesrco/SinkingParticles_GBC2015/blob/master/data/Depth-T-BP-Resp%20Correlations%3F.xlsx)
   * The power law function used for interpolation is contained in the m-file [Powerlaw_func_y0_at_0.m](https://github.com/jamesrco/SinkingParticles_GBC2015/blob/master/Powerlaw_func_y0_at_0.m)
   * The script uses the dependencies [nlleasqr.m](https://github.com/jamesrco/dependencies-useful-scripts/blob/master/nlleasqr.m) and [dfdp.m](https://github.com/jamesrco/dependencies-useful-scripts/blob/master/dfdp.m) for calculation of function parameters

3. [Depth-T-BP-Resp Correlations.R](https://github.com/jamesrco/SinkingParticles_GBC2015/blob/master/Depth-T-BP-Resp%20Correlations.R) and [Depth-T-BP-Resp-Correlations.m](https://github.com/jamesrco/SinkingParticles_GBC2015/blob/master/Depth-T-BP-Resp-Correlations.m): R and MATLAB script for generating the plots and performing the various regressions reported in Fig. 4 of the paper. The MATLAB script requires [linfit.m](https://github.com/jamesrco/dependencies-useful-scripts/blob/master/linfit.m) and [lsqfitma.m](https://github.com/jamesrco/dependencies-useful-scripts/blob/master/lsqfitma.m)

4. In [SinkingParticles_GBC2015/data](https://github.com/jamesrco/SinkingParticles_GBC2015/tree/master/data): Various, somewhat disorganized data files necessary for some or all of the scripts. _Users looking to access final cruise data are strongly advised to use the links above, which will take you to final, "official" (and more nicely organized) cruise data on BCO-DMO._ The files here include:
   * [BLATZ_VICE_BactProd_calcs_20140115.mat](https://github.com/jamesrco/SinkingParticles_GBC2015/blob/master/data/BLATZ_VICE_BactProd_calcs_20140115.mat): Bacterial production rates, calculated using [3H_Leu_BactProd](https://github.com/jamesrco/3H_Leu_BactProd)
   * [Depth-T-BP-Resp Correlations?.xlsx](https://github.com/jamesrco/SinkingParticles_GBC2015/blob/master/data/Depth-T-BP-Resp%20Correlations%3F.xlsx): Excel file containing temperature, BP, and community respiration data by depth
   * [Depth-T-BP-Resp Correlations.csv](https://github.com/jamesrco/SinkingParticles_GBC2015/blob/master/data/Depth-T-BP-Resp%20Correlations.csv): .csv version of preceding Excel file
   * [KN207-1 and KN207-3 PIT samples.xlsx](https://github.com/jamesrco/SinkingParticles_GBC2015/blob/master/data/KN207-1%20and%20KN207-3%20PIT%20samples.xlsx): Excel file containing PIT (sediment) trap data
   * [VICE and BLATZ net trap data.xlsx](https://github.com/jamesrco/SinkingParticles_GBC2015/blob/master/data/VICE%20and%20BLATZ%20net%20trap%20data.xlsx): Excel file containing net trap data (methodology described in Collins et al. 2015 or [Peterson, M. L., S. G. Wakeham, C. Lee, M. A. Askea, and J. C. Miquel (2005), Novel techniques for collection of sinking particles in the ocean and determining their settling rates, Limnol. Oceanogr. Methods, 3, 520–532](http://onlinelibrary.wiley.com/doi/10.4319/lom.2005.3.520/abstract))
   * [VICE and BLATZ shipboard O2 incubations.xlsx](https://github.com/jamesrco/SinkingParticles_GBC2015/blob/master/data/VICE%20and%20BLATZ%20shipboard%20O2%20incubations.xlsx): Excel file containing data from shipboard incubations of particle material in which dissolved oxygen was measured over time (methods described in Collins et al. 2015)
