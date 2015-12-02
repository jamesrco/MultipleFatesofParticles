% Depth-T-BP-Resp-Correlations.m
% 23 Mar 2015, JRC

% Companion file to R script Depth-T-BP-Resp Correlations.R

% load in data necessary for calculation of regressions

[num_pooledR txt_pooledR raw_pooledR] =...
    xlsread('/Users/jrcollins/Dropbox/KN207-3 NA VICE/Data/2012 North Atlantic resp & flux manny/Depth-T-BP-Resp Correlations?.xlsx',...
    'Depth-T-BP-Resp Correlations.cs'); % read in pooled WCR data

% extract variables

Depth_Resp = num_pooledR(1:25,1);
Sta_Resp = txt_pooledR(2:26,1);
Resp_mg_C_m3_d = num_pooledR(1:25,2);
Uncert_resp_mg_C_m3_d = num_pooledR(1:25,3);

Depth_BP = num_pooledR(1:36,6);
Sta_BP = txt_pooledR(2:37,6);
BP_mg_C_m3_d = num_pooledR(1:36,7);
Uncert_BP_mg_C_m3_d = num_pooledR(1:36,8);

Resp_mg_C_m3_d_comp = num_pooledR(1:16,12);
BP_mg_C_m3_d_comp = num_pooledR(1:16,14);

% --------------------------------------------------------
% resp on depth
% --------------------------------------------------------

% pooled dataset by OLS regression

x=Depth_Resp;
y=Resp_mg_C_m3_d;
sy=Uncert_resp_mg_C_m3_d; %
[a_R sa_R cov_R r_R]=linfit(log(x),y,sy);
[a_R sa_R cov_R r_R]

plot(log(Depth_Resp),Resp_mg_C_m3_d,'o')

% --------------------------------------------------------
% BP on depth
% --------------------------------------------------------

% pooled dataset by OLS regression

x=Depth_BP;
y=BP_mg_C_m3_d;
sy=Uncert_BP_mg_C_m3_d; %
[a_BP sa_BP cov_BP r_bp]=linfit(log(x),y,sy);
[a_BP sa_BP cov_BP r_bp]

plot(log(Depth_BP),BP_mg_C_m3_d,'o')

% --------------------------------------------------------
% BP on resp
% --------------------------------------------------------

% pooled dataset by Type II MA regression using lsqfit ma

% *** yields the same results as the MA regression routine
% in the lmodel2 R package, but lsqfitma.m will give us SE's on the
% parameters rather than confidence intervals

x=log(Resp_mg_C_m3_d_comp);
y=log(BP_mg_C_m3_d_comp);

[m,b,r,sm,sb,xbar,ybar]=lsqfitma(x,y);
