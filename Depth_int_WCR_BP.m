% Depth_int_WCR_BP.m
% 24 Mar 2015, JRC

% *** Improves on the method used to calculate depth-integrated BP and resp
% in respflux2012_v3.m ***'

% This script uses a bootstrap Monte Carlo approach to estimate
% depth-integrated rates of water col. respiration and bacterial
% production, using volumetric rate data for each of those parameters
% obtained at each process station during the KN207-1 and KN207-3 cruises.

% The depth-integrated figures are estimated using a range of values for
% both isotope dilution and the C:leu conversion factor

% Because data are sparse between 50 and 150 m (the range over which we
% want to integrate for Collins et al, 2015, GBC), we will create a series
% on = numsims simulated data sets for each station and parameter, assuming
% the values are drawn from normal distributions with sigmas as calculated
% from replication. We will then use those individual data sets to generate
% n = numsims unique power-law curve fits. We will numerically integrate these
% individual curves over the range of interest. For each station, the xbar
% and sigma of each set of integrated estimates will be used as the
% centroid and uncertainty of the rate at that station.

%------------------------------------------------
% load data & extract variables
%------------------------------------------------

% load in data necessary for calculation of regressions

[num_pooledR txt_pooledR raw_pooledR] =...
    xlsread('/Users/jrcollins/Dropbox/Cruises & projects/KN207-3 NA VICE/Data/2012 North Atlantic resp & flux manny/Depth-T-BP-Resp Correlations?.xlsx',...
    'Depth-T-BP-Resp Correlations.cs');

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

%------------------------------------------------
% specification of constants
%------------------------------------------------

% create a vector to specify our stations

Proc_stations = ['QL1';'QL2';'PS1';'PS2';'PS3';'PS4'];
Proc_stations = cellstr(Proc_stations);

% define number of iterations for our monte carlo bootstraps

numsims = 100000;

% create a vector of depths (0.5 m intervals) for which we want to obtain 
% predicted values

Int_depths = [50:0.5:150];

% randomly generate vector of possible ID's, assuming ID ranges from 1 to 2
% following a uniform distribition

ID_simdist = 1 + (2-1).*rand(numsims,1);

% generate vector of C:leu CF's (units of kg C/mol leu), assuming
% normal distribution of CF's for depths > 50 m according to dataset
% assembled from literature by Giering et al 2014, Nature (Extended Data
% Fig 4)

% for this analysis, we'll have to truncate simulated values at ± 1.5 SD since we
% can't have a CF < 0

v_kgC_mol_leu = 0.44;
v_kgC_mol_leu_SD = 0.27;

v_kgC_mol_leu_simdist = v_kgC_mol_leu + v_kgC_mol_leu_SD*randn(numsims,1);

v_kgC_mol_leu_simdist(v_kgC_mol_leu_simdist<v_kgC_mol_leu-1.5*v_kgC_mol_leu_SD)...
    = v_kgC_mol_leu-1.5*v_kgC_mol_leu_SD;

v_kgC_mol_leu_simdist(v_kgC_mol_leu_simdist>v_kgC_mol_leu+1.5*v_kgC_mol_leu_SD)...
    = v_kgC_mol_leu+1.5*v_kgC_mol_leu_SD;

% also need to specify the CF we used to get the current BP numbers in
% terms of mg_C_m3_d in the other script (respflux2012_v3.m), so we can
% back-convert to a point where we can use the range of values given by
% Giering et al

v_kgC_mol_leu_base = 1.5; % 1.5 kg C/mol leu

%------------------------------------------------
% generation of simulated data
%------------------------------------------------

% generate simulated data for values of each parameter (n = numsims),
% assuming each set of replicate observations reflects a normal
% distribution with means Resp_mg_C_m3_d or BP_mg_C_m3_d and SD's
% Uncert_resp_mg_C_m3_d or Uncert_BP_mg_C_m3_d

% in case of BP, also evaluating a range of ID's (1 to 2, assuming continuous uniform
% dist) and a range of C:leu factors, assuming for > 50 m the mean and SD
% specifed by Giering et al 2014
        
% respiration (easier than BP)

clear Resp_simdists; % in case they already exist from a previous run

for i=1:size(Resp_mg_C_m3_d,1)
    Resp_simdists(i,:) = Resp_mg_C_m3_d(i) + Uncert_resp_mg_C_m3_d(i).*randn(numsims,1);
end

% BP (two separate steps)

% first, generate "base" rate data, using only replication uncertainities,
% as for resp

clear BP_simdists; % in case they already exist from a previous run

for i=1:size(BP_mg_C_m3_d,1)
    BP_simdists(i,:) = BP_mg_C_m3_d(i) + Uncert_BP_mg_C_m3_d(i).*randn(numsims,1);
end

% now, go through and scale by ID and the C:leu CF

for j=1:numsims
    BP_simdists(:,j) = BP_simdists(:,j)/v_kgC_mol_leu_base*v_kgC_mol_leu_simdist(j)*ID_simdist(j);
end

%------------------------------------------------
% the actual bootstrap analysis
%------------------------------------------------

% preallocate matrices for our integration simulation data and for our
% results

Int_sims = nan(numsims,1);

Bootout_depthint_50_150_powerfit_mg_C_m2_d = nan(length(Proc_stations),4);

% for each process station:

% 1. generate n = numsims power-law curve fits using the simulated data
% created above 
% 2. numerically integrate each individual curve from 50 to 150 m
% 3. the xbar and sigma of this final distributon of n = numsims integrated
% estimates per station will be our centroid and uncertainty for each station 

for j=1:length(Proc_stations)
    
    %------------------------------------------------
    % first, respiration
    %------------------------------------------------
    
    % check to make sure we have enough respiration data for this station
    % to proceed
    
    if sum(ismember(Sta_Resp,Proc_stations(j))) > 3
        
        % specify some guesses for nlleasqr to use when fitting curves
        
        pin = [192 -0.75 2];
        
        for k=1:numsims % cycle through each simulated data set
            
            % fit a power-law curve to this data set
            
            x = Depth_Resp(find(ismember(Sta_Resp,Proc_stations(j))));
                
            y = Resp_simdists(find(ismember(Sta_Resp,Proc_stations(j))),k);

            [f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2]=nlleasqr(x,y,pin,'Powerlaw_func_y0_at_0',.0001,50);
            
            % using this curve fit, estimate rate of resp at each 0.5 m
            % depth from 50 to 150 m

            yhats = p(1)*Int_depths.^(p(2))+p(3);
            
            % numerically integrate these rates to get an estimate for this
            % simulation, then store

            Int_sims(k) = trapz(Int_depths,yhats);
        end
    
    end
    
    % obtain an average and SD for our simulated depth-integrated estimates, and
    % record to our matrix

    Bootout_depthint_50_150_powerfit_mg_C_m2_d(j,1) = mean(Int_sims);
    Bootout_depthint_50_150_powerfit_mg_C_m2_d(j,2) = std(Int_sims);

    Int_sims = nan(numsims,1); % reset values in our placeholder matrix, as a precaution
    
    %------------------------------------------------
    % second, BP
    %------------------------------------------------
    
    % check to make sure we have enough BP data for this station
    % to proceed
    
    if sum(ismember(Sta_BP,Proc_stations(j))) > 3
        
        % specify some guesses for nlleasqr to use when fitting curves
        
        pin = [5 -0.09 -2];

        for k=1:numsims  % cycle through each simulated data set
            
            % fit a power-law curve to this data set

            x = Depth_BP(find(ismember(Sta_BP,Proc_stations(j))));

            y = BP_simdists(find(ismember(Sta_BP,Proc_stations(j))),k);

            [f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2]=nlleasqr(x,y,pin,'Powerlaw_func_y0_at_0',.0001,50);

            % using this curve fit, estimate rate of resp at each 0.5 m
            % depth from 50 to 150 m

            yhats = p(1)*Int_depths.^(p(2))+p(3);
            
            % numerically integrate these rates to get an estimate for this
            % simulation, then store

            Int_sims(k) = trapz(Int_depths,yhats);

        end
    
    end
    
    % obtain an average and SD for our depth-integrated estimates, and
    % record to our matrix

    Bootout_depthint_50_150_powerfit_mg_C_m2_d(j,3) = mean(Int_sims);
    Bootout_depthint_50_150_powerfit_mg_C_m2_d(j,4) = std(Int_sims);

    Int_sims = nan(numsims,1); % reset values in our placeholder matrix, as a precaution
    
end
