% respflux2012_v3.m
% 7 Dec 2013, JRC
% 30 Apr 2014, JRC
% Verison 3 code m-file created from v2 code on 28 Jan 2015; uses Monte
% Carlo sims to estimate uncertaintes rather than error prop

% data workup for data from KN207-1 and KN207-3, calculating contribution
% of bacterial respiration/remineralization to particle flux attenuation

% BP data as calculated in script riBPdata.m

clear all; close all; clf;

scrsz = get(0,'ScreenSize'); % define a screen size variable so I can make figures that look decent

%% some assumed constants

RQ = 117/170; % respiratory quotient, i.e., mol CO2 produced : mol O2 consumed; Anderson and Sarmiento, GBC, 1994

% more constants for BP calculations in BP section below

%% load in data from different sources

% first, load in BP data from the riBPdata.m output file

load ('/Users/jrcollins/Dropbox/KN207-3 NA VICE/Data/BP workups/BLATZ_VICE_BactProd_calcs_20140115.mat','BP_results_replicate_averaged_by_CTD');

BPdata.CruiseID=BP_results_replicate_averaged_by_CTD(:,1);
BPdata.CTDStationID=BP_results_replicate_averaged_by_CTD(:,2);
BPdata.CTD_depth=BP_results_replicate_averaged_by_CTD(:,3);
BPdata.Trit_leu_uptake_pmol_L_hr=BP_results_replicate_averaged_by_CTD(:,16);
BPdata.Std_dev_trit_leu_uptake_pmol_L_hr=BP_results_replicate_averaged_by_CTD(:,17);
BPdata.BP_signal_to_noise_ratio=BP_results_replicate_averaged_by_CTD(:,18);
BPdata.CTD_sample_incutemp=BP_results_replicate_averaged_by_CTD(:,19);
BPdata.CTD_station_time=BP_results_replicate_averaged_by_CTD(:,20);
BPdata.CTD_station_time=x2mdate(BPdata.CTD_station_time,1);
BPdata.CTD_station_lat=BP_results_replicate_averaged_by_CTD(:,21);
BPdata.CTD_station_long=BP_results_replicate_averaged_by_CTD(:,22);

% load PIT trap C data

[num_PITdata txt_PITdata raw_PITdata] = xlsread('/Users/jrcollins/Dropbox/KN207-3 NA VICE/Data/Export & prod/KN207-1 and KN207-3 PIT samples.xlsx','Data for MATLAB export');

PITdata.CruiseID=cell2mat(raw_PITdata(4:end,1));
PITdata.ProcStationID=cell2mat(raw_PITdata(4:end,3));
PITdata.Deploy_depth=cell2mat(raw_PITdata(4:end,4));
PITdata.PIC_flux_mgC_m2_d=cell2mat(raw_PITdata(4:end,5));
PITdata.TPC_flux_mgC_m2_d=cell2mat(raw_PITdata(4:end,6));
PITdata.TPC_flux_anal_uncert_mgC_m2_d=cell2mat(raw_PITdata(4:end,7));

% load PIT trap deployment data

[num_PITdeployments txt_PITdeployments raw_PITdeployments] = xlsread('/Users/jrcollins/Dropbox/KN207-3 NA VICE/Data/Export & prod/KN207-1 and KN207-3 PIT samples.xlsx','PIT deployment data');

PITdeployments.CruiseID=cell2mat(raw_PITdeployments(2:end,1));
PITdeployments.ProcStationID=cell2mat(raw_PITdeployments(2:end,3));
PITdeployments.ProcStation_Lat=cell2mat(raw_PITdeployments(2:end,4));
PITdeployments.ProcStation_Long=cell2mat(raw_PITdeployments(2:end,5));
PITdeployments.Deploy_start=x2mdate(cell2mat(raw_PITdeployments(2:end,6)));
PITdeployments.Deploy_end=x2mdate(cell2mat(raw_PITdeployments(2:end,7)));
PITdeployments.Deploy_duration=PITdeployments.Deploy_end-PITdeployments.Deploy_start;

% load O2 incubation data

% net trap particle incubations from KN207-3

[num_nettrincdata_2073 txt_nettrincdata_2073 raw_nettrincdata_2073] = xlsread('/Users/jrcollins/Dropbox/KN207-3 NA VICE/Data/Export & prod/VICE and BLATZ shipboard O2 incubations.xlsx','Net trap particle incubations');

Net_trap_incudata.CruiseID=cell2mat(raw_nettrincdata_2073(4:end,1));
Net_trap_incudata.NetTrapID=cell2mat(raw_nettrincdata_2073(4:end,2));
Net_trap_incudata.Incu_temp=cell2mat(raw_nettrincdata_2073(4:end,3));
Net_trap_incudata.Read_time=x2mdate(cell2mat(raw_nettrincdata_2073(4:end,4)));
Net_trap_incudata.O2_conc_uM=cell2mat(raw_nettrincdata_2073(4:end,5:7));

% net trap particle incubations from KN207-1 (controls from Bethanie's PUA
% addition experiments)

[num_nettrincdata_2071 txt_nettrincdata_2071 raw_nettrincdata_2071] = xlsread('/Users/jrcollins/Dropbox/KN207-3 NA VICE/Data/Export & prod/VICE and BLATZ shipboard O2 incubations.xlsx','BLATZ net tr incs from Bethanie');

Net_trap_incudata2071.CruiseID=cell2mat(raw_nettrincdata_2071(5:end,1));
Net_trap_incudata2071.NetTrapID=cell2mat(raw_nettrincdata_2071(5:end,2));
Net_trap_incudata2071.Incu_temp=cell2mat(raw_nettrincdata_2071(5:end,3));
Net_trap_incudata2071.Split_frac_2071=cell2mat(raw_nettrincdata_2071(5:end,4));
Net_trap_incudata2071.O2_cons_rates_uM_hr=cell2mat(raw_nettrincdata_2071(5:end,5:7));

% water column respiration data, both cruises

[num_watercolincdata txt_watercolincdata raw_watercolincdata] = xlsread('/Users/jrcollins/Dropbox/KN207-3 NA VICE/Data/Export & prod/VICE and BLATZ shipboard O2 incubations.xlsx','Water column incubations');

Water_col_incudata.CruiseID=cell2mat(raw_watercolincdata(2:end,1));
Water_col_incudata.CTD_StationID=cell2mat(raw_watercolincdata(2:end,2));
Water_col_incudata.CTD_depth=cell2mat(raw_watercolincdata(2:end,3));
Water_col_incudata.Incu_temp=cell2mat(raw_watercolincdata(2:end,4));
Water_col_incudata.Read_time=x2mdate(cell2mat(raw_watercolincdata(2:end,5)));
Water_col_incudata.O2_conc_uM=cell2mat(raw_watercolincdata(2:end,6:15));

% load POC/PON data and station data for net trap deployments

[num_nettrapdata txt_nettrapdata raw_nettrapdata] = xlsread('/Users/jrcollins/Dropbox/KN207-3 NA VICE/Data/Export & prod/VICE and BLATZ net trap data.xlsx','Net trap deployment data');

Net_trap_data.CruiseID=cell2mat(raw_nettrapdata(6:end,1));
Net_trap_data.ProcStationID=cell2mat(raw_nettrapdata(6:end,2));
Net_trap_data.NetTrapID=cell2mat(raw_nettrapdata(6:end,3));
Net_trap_data.Deploy_depth=cell2mat(raw_nettrapdata(6:end,4));
Net_trap_data.ProcStation_Lat=cell2mat(raw_nettrapdata(6:end,5));
Net_trap_data.ProcStation_Long=cell2mat(raw_nettrapdata(6:end,6));
Net_trap_data.Deploy_start=x2mdate(cell2mat(raw_nettrapdata(6:end,7)));
Net_trap_data.Deploy_end=x2mdate(cell2mat(raw_nettrapdata(6:end,8)));
Net_trap_data.Deploy_duration=Net_trap_data.Deploy_end-Net_trap_data.Deploy_start;
Net_trap_data.POC_mgC_per_trap=cell2mat(raw_nettrapdata(6:end,16)); % blank-corrected values for entire trap
Net_trap_data.Std_dev_POC_mgC_per_trap=cell2mat(raw_nettrapdata(6:end,17));
Net_trap_data.PON_mgN_per_trap=cell2mat(raw_nettrapdata(6:end,18));
Net_trap_data.Std_dev_PON_mgN_per_trap=cell2mat(raw_nettrapdata(6:end,19));

% load in CTD ship data, both cruises

[num_KN2071data txt_KN2071data  raw_KN2071data ] = xlsread('/Users/jrcollins/Dropbox/KN207-3 NA VICE/Data/BP workups/KN207-1_shipcasts.xlsx');
[num_KN2073data txt_KN2073data  raw_KN2073data ] = xlsread('/Users/jrcollins/Dropbox/KN207-3 NA VICE/Data/BP workups/KN207-3_shipcasts.xlsx');

Ship_data_2071=cell2mat(raw_KN2071data(9:end,[1 5 6 7]));
Ship_data_2073=cell2mat(raw_KN2073data(11:end,[1 5 6 7]));

Ship_data = nan(length(Ship_data_2071)+length(Ship_data_2073),5);

Ship_data(1:length(Ship_data_2071),1)=2071;
Ship_data(length(Ship_data_2071)+1:end,1)=2073;

Ship_data(1:length(Ship_data_2071),2:5)=Ship_data_2071;
Ship_data(length(Ship_data_2071)+1:end,2:5)=Ship_data_2073;

% convert excel dates to MATLAB dates

Ship_data(:,3)=x2mdate(Ship_data(:,3),1); % Ship_data timestamps in 1904 

%% perform some calculations and data manipulation, propagating through uncertainties

% ----------------------------------------------------------------------------
% PIT trap fluxes
% ----------------------------------------------------------------------------

% calculate mean PIT fluxes and associated uncertainties by cruise/process station ID/depth,
% then store in structure PIT_fluxes, in units of mmolC/m2/d

PIT_depths=[50 150 300]; % define PIT trap deployment depths

c_i=1;
for i=1:length(PITdeployments.ProcStationID)  
    thiscruise=PITdeployments.CruiseID(i);
    thisstation=PITdeployments.ProcStationID(i);
    for j=1:length(PIT_depths)
        thisdepth=PIT_depths(j);
        ind=find(PITdata.CruiseID==thiscruise & ...
                PITdata.ProcStationID==thisstation & PITdata.Deploy_depth==thisdepth);
            PIC_data=PITdata.PIC_flux_mgC_m2_d(ind);
            TPC_data=PITdata.TPC_flux_mgC_m2_d(ind);
            TPC_anal_uncert=PITdata.TPC_flux_anal_uncert_mgC_m2_d(ind);
            ind_nan=find(isnan(TPC_anal_uncert));
            TPC_anal_uncert(ind_nan)=[];
            POC=TPC_data-PIC_data;
            ind_nan=find(isnan(PIC_data));
            PIC_data(ind_nan)=[];
            Avg_PIC=mean(PIC_data);
            Uncert_PIC=std(PIC_data);
            ind_nan=find(isnan(TPC_data));
            TPC_data(ind_nan)=[];
            Avg_TPC=mean(TPC_data);
            Std_TPC=std(TPC_data);
            Uncert_TPC=sqrt(sum((TPC_anal_uncert).^2)/length(TPC_anal_uncert)+Std_TPC^2);
            % for TPC, uncertainty is combination of std dev of replicates and pooled analytical
            % uncertainties from vacuum line
            ind_nan=find(isnan(POC));
            POC(ind_nan)=[];
            Avg_POC=mean(POC);
            Uncert_POC=sqrt(Uncert_TPC^2+Uncert_PIC^2);
            % now, convert to mmoles and write data where appropriate
            PIT_fluxes.CruiseID(c_i,1)=thiscruise;
            PIT_fluxes.ProcStationID(c_i,1)=thisstation;
            PIT_fluxes.ProcStation_Lat(c_i,1)=PITdeployments.ProcStation_Lat(i);
            PIT_fluxes.ProcStation_Long(c_i,1)=PITdeployments.ProcStation_Long(i);
            PIT_fluxes.Deploy_start(c_i,1)=PITdeployments.Deploy_start(i);
            PIT_fluxes.Deploy_duration(c_i,1)=PITdeployments.Deploy_duration(i);
            PIT_fluxes.Deploy_depth(c_i,1)=thisdepth;
            PIT_fluxes.Avg_PIC_flux_mmolC_m2_d(c_i,1)=Avg_PIC/12.01;
            PIT_fluxes.Uncert_PIC_flux_mmolC_m2_d(c_i,1)=Uncert_PIC/12.01;
            PIT_fluxes.Avg_POC_flux_mmolC_m2_d(c_i,1)=Avg_POC/12.01;
            PIT_fluxes.Uncert_POC_flux_mmolC_m2_d(c_i,1)=Uncert_POC/12.01;
            PIT_fluxes.POC_PIC_mol_mol(c_i,1)=Avg_POC/Avg_PIC;
            clear POC;
        c_i=c_i+1;
    end
end

clear i j c_i;

% ----------------------------------------------------------------------------
% water column respiration rates
% ----------------------------------------------------------------------------

% calculate mean O2 consumption rates for water column incubations, store
% in structure Watercol_resprates, in units of mmolC/m3/d using assumed RQ

% includes cleanup of data from cruise 2073, station 87 by eliminating
% first O2 reading;=something was wrong with incubator for deeper samples
% incubated at 5.5C (17.3-150 m)

% get some information about the dataset
cruises=unique(Water_col_incudata.CruiseID);
incu_stations=unique(Water_col_incudata.CTD_StationID);
incu_depths=unique(Water_col_incudata.CTD_depth);

c_i=1;
for i=1:length(cruises)
    for j=1:length(incu_stations);
        for k=1:length(incu_depths);
            ind=find(Water_col_incudata.CruiseID==cruises(i) & ...
                Water_col_incudata.CTD_StationID==incu_stations(j) & ...
                Water_col_incudata.CTD_depth==incu_depths(k));
            if numel(ind)>1 % a respiration profile was made with water from station & depth
                % first, ascertain the number of replicates
                numreps=sum(sum(isnan(Water_col_incudata.O2_conc_uM(ind,:)))<length(ind));
                % calculate O2 consumption rates for each set of replicate
                % measurements
                O2data_k=Water_col_incudata.O2_conc_uM(ind,:);
                readtimes_k=Water_col_incudata.Read_time(ind,:);
                for l=1:numreps % cycle through replicates, making sure we
                    % don't have any NaNs before calculating the regression
                    O2data_l=O2data_k(:,l);
                    readtimes_l=readtimes_k;
                    ind_nan=find(isnan(O2data_l));
                    O2data_l(ind_nan)=[];
                    readtimes_l(ind_nan)=[];
                    % do the regression
                    % but first check to see if the data is from the bad incubation
                    if cruises(i)==2073 && incu_stations(j)==87 && ismember(incu_depths(k),[17.3 30 51 150])
                        x=readtimes_l(2:end); % eliminate the first reading from the regression calculation
                        y=O2data_l(2:end);
                    else % it's not one of the bad incubation replicates, so use the whole data range
                    x=readtimes_l;
                    y=O2data_l;
                    end
                    sy=1; % unweighted since we don't know what uncertaintites are
                    [a sa cov r]=linfit(x,y,sy);
                    rate(l)=a(2);
                    rate_uncert(l)=sa(2);
                    corrcoeff(l)=r;
                end
                % find mean rate and uncertainty, convert to units of mmol C/m3 and store to file
                O2_cons_rate_uM_O2_d_avg=mean(rate);
                if cruises(i)==2073 && incu_stations(j)==87 && ismember(incu_depths(k),[17.3 30 51 150])
                    O2_cons_rate_uM_O2_d_uncert=std(rate(1:2))/sqrt(2); % in case of station 87, can't
                     % use S.E. of regression, so estimate using SE of mean
                else
                O2_cons_rate_uM_O2_d_uncert=sqrt(sum((rate_uncert).^2))/length(rate_uncert);
                end
                Watercol_resprates.CruiseID(c_i,1)=cruises(i);
                Watercol_resprates.CTD_StationID(c_i,1)=incu_stations(j);
                Watercol_resprates.CTD_depth(c_i,1)=incu_depths(k);
                Watercol_resprates.Avg_mmolC_m3_d(c_i,1)=O2_cons_rate_uM_O2_d_avg*(1/1000)*(1000)*RQ;
                Watercol_resprates.Uncert_mmolC_m3_d(c_i,1)=O2_cons_rate_uM_O2_d_uncert*(1/1000)*(1000)*RQ;
                Watercol_resprates.R2_mean(c_i,1)=mean(corrcoeff);
                 c_i=c_i+1;
                 clear rate rate_uncert corrcoeff;
            end
        end
    end
end

clear i j k l c_i;

% match and load in location & time information for the water column incubations,
% using ship CTD data

for i=1:length(Watercol_resprates.CruiseID)
    ind=find(Ship_data(:,1)==Watercol_resprates.CruiseID(i) & ...
        Ship_data(:,2)==Watercol_resprates.CTD_StationID(i));
    Watercol_resprates.CTD_Lat(i,1)=Ship_data(ind,4);
    Watercol_resprates.CTD_Long(i,1)=Ship_data(ind,5);
    Watercol_resprates.CTD_Deploy_start(i,1)=Ship_data(ind,3);
end

clear i ind;

% associate water column incubations with the relevant PIT trap deployment (if proximate)
% so we can use for analysis below

% prepopulate target matrix
Watercol_resprates.ProcStationID=nan(length(Watercol_resprates.Avg_mmolC_m3_d),1);

% perform a query with proximity critera ± 0.5 deg lat and long, and ± 12
% hrs
for i=1:length(Watercol_resprates.ProcStationID)
    ind=find(PITdeployments.CruiseID==Watercol_resprates.CruiseID(i) & ...
        abs(PITdeployments.ProcStation_Lat-Watercol_resprates.CTD_Lat(i))<=0.5 & ...
        abs(PITdeployments.ProcStation_Long-Watercol_resprates.CTD_Long(i))<=0.5 & ...
        Watercol_resprates.CTD_Deploy_start(i)>=(PITdeployments.Deploy_start-0.5) & ...
        Watercol_resprates.CTD_Deploy_start(i)<=(PITdeployments.Deploy_end+0.5));
    if isfinite(ind) % there's a match!
        Watercol_resprates.ProcStationID(i,1)=PITdeployments.ProcStationID(ind(1));
    end
end
clear i ind;

% ----------------------------------------------------------------------------
% calculate depth-integrated rates of water column community resp from 50
% to 150 m (can only accomplish for PS2, PS3, PS4)
%
% *** note that the script Depth_int_WCR_BP.m contains a much
% better way of doing this (and which doesn't require interpolation) ***
% ----------------------------------------------------------------------------

% first, will need to interpolate or approximate a value at 50 m in cases where a bottle wasn't fired at that depth

cruises=unique(Watercol_resprates.CruiseID);
CTDstations=unique(Watercol_resprates.CTD_StationID);

clear Watercol_resprates_50_150;

for i=1:length(cruises)
    for j=1:length(CTDstations)
        ind_thisCTD=find(Watercol_resprates.CruiseID==cruises(i) & Watercol_resprates.CTD_StationID...
            ==CTDstations(j)); % index of all analyzed depths at this CTD station
        if isfinite(ind_thisCTD)
        if isfinite(Watercol_resprates.Avg_mmolC_m3_d(ind_thisCTD)) & numel(Watercol_resprates.Avg_mmolC_m3_d(ind_thisCTD))>1 % first check to see if there's
            % even data for this station, and that there's more than one data point; if not, no point in doing any of
            % this, so bail
            if exist('Watercol_resprates_50_150') % define matrix insertion point (should be first row after end of existing BP_int_50_150 data)
                c_i=length(Watercol_resprates_50_150.CruiseID)+1;
            else
                c_i=1; % need this here for the first run, when no data has been inserted yet
            end
            if any(ismember(Watercol_resprates.CTD_depth(ind_thisCTD),50)) % we have a rate for 50 m exactly
                ind_50to150=find(Watercol_resprates.CruiseID==cruises(i) & Watercol_resprates.CTD_StationID...
                ==CTDstations(j) & Watercol_resprates.CTD_depth>=50);
                % index of all measurements from 50m to 150m, inclusive
                % now, feed in data
                startrow=c_i;
                endrow=c_i+length(ind_50to150)-1;
                Watercol_resprates_50_150.CruiseID(startrow:endrow,1)=cruises(i);
                Watercol_resprates_50_150.CTDStationID(startrow:endrow,1)=CTDstations(j);
                Watercol_resprates_50_150.CTD_depth(startrow:endrow,1)=Watercol_resprates.CTD_depth(ind_50to150);
                Watercol_resprates_50_150.Avg_mmolC_m3_d(startrow:endrow,1)=Watercol_resprates.Avg_mmolC_m3_d(ind_50to150);
                Watercol_resprates_50_150.Uncert_mmolC_m3_d(startrow:endrow,1)=Watercol_resprates.Uncert_mmolC_m3_d(ind_50to150);
                Watercol_resprates_50_150.CTD_station_long(startrow:endrow,1)=Watercol_resprates.CTD_Long(ind_50to150);
                Watercol_resprates_50_150.CTD_station_lat(startrow:endrow,1)=Watercol_resprates.CTD_Lat(ind_50to150);
                Watercol_resprates_50_150.CTD_station_time(startrow:endrow,1)=Watercol_resprates.CTD_Deploy_start(ind_50to150);
                Watercol_resprates_50_150.ProcStationID(startrow:endrow,1)=Watercol_resprates.ProcStationID(ind_50to150);
            elseif any(ismember(Watercol_resprates.CTD_depth(ind_thisCTD),[45:1:55])) % we have a value within 5 m 
                % of 50 m which we can directly substitute
                ind_prox=find(Watercol_resprates.CruiseID==cruises(i) & Watercol_resprates.CTD_StationID...
                ==CTDstations(j) & Watercol_resprates.CTD_depth>=45 & Watercol_resprates.CTD_depth<=55); % make a list
                % if there's more than 1, select the closest
                dz=abs(Watercol_resprates.CTD_depth(ind_prox)-50);
                ind_best=find(min(dz));
                ind_sub=ind_prox(ind_best);
                ind_50to150=find(Watercol_resprates.CruiseID==cruises(i) & Watercol_resprates.CTD_StationID...
                ==CTDstations(j) & Watercol_resprates.CTD_depth>=50);
                % now, feed in data, including the additional estimated value at 50 m
                startrow=c_i;
                endrow=c_i+length(ind_50to150);
                Watercol_resprates_50_150.CruiseID(startrow:endrow,1)=cruises(i);
                Watercol_resprates_50_150.CTDStationID(startrow:endrow,1)=CTDstations(j);
                Watercol_resprates_50_150.CTD_depth(startrow:endrow,1)=...
                    [50 Watercol_resprates.CTD_depth(ind_50to150)'];
                Watercol_resprates_50_150.Avg_mmolC_m3_d(startrow:endrow,1)=...
                    [Watercol_resprates.Avg_mmolC_m3_d(ind_sub) Watercol_resprates.Avg_mmolC_m3_d(ind_50to150)'];
                Watercol_resprates_50_150.Uncert_mmolC_m3_d(startrow:endrow,1)=...
                    [Watercol_resprates.Uncert_mmolC_m3_d(ind_sub) Watercol_resprates.Uncert_mmolC_m3_d(ind_50to150)'];
                Watercol_resprates_50_150.CTD_station_long(startrow:endrow,1)=...
                    [Watercol_resprates.CTD_Long(ind_sub) Watercol_resprates.CTD_Long(ind_50to150)'];
                Watercol_resprates_50_150.CTD_station_lat(startrow:endrow,1)=...
                    [Watercol_resprates.CTD_Lat(ind_sub) Watercol_resprates.CTD_Lat(ind_50to150)'];
                Watercol_resprates_50_150.CTD_station_time(startrow:endrow,1)=...
                    [Watercol_resprates.CTD_Deploy_start(ind_sub) Watercol_resprates.CTD_Deploy_start(ind_50to150)'];
                Watercol_resprates_50_150.ProcStationID(startrow:endrow,1)=...
                    [Watercol_resprates.ProcStationID(ind_sub) Watercol_resprates.ProcStationID(ind_50to150)'];
            else % there isn't a rate measurement at anywhere near 50 m, will have to interpolate
                ind_50to150=find(Watercol_resprates.CruiseID==cruises(i) & Watercol_resprates.CTD_StationID...
                ==CTDstations(j) & Watercol_resprates.CTD_depth>=50);
                dz=Watercol_resprates.CTD_depth(ind_thisCTD)-50;
                ind_shallower=find(dz<0); % all measurements shallower than 50 m
                ind_shallower=[ind_thisCTD(ind_shallower) dz(ind_shallower)];
                [C,I]=max(ind_shallower(:,2));
                ind_nextshallowest=ind_shallower(I,1);
                ind_deeper=find(dz>0); % all measurements deeper than 50 m
                ind_deeper=[ind_thisCTD(ind_deeper) dz(ind_deeper)];
                [C,I]=min(ind_deeper(:,2));
                ind_nextdeepest=ind_deeper(I,1);
                % now that we have the two values bracketing 50 m, can
                % interpolate linearly
                x=[Watercol_resprates.CTD_depth(ind_nextshallowest) Watercol_resprates.CTD_depth(ind_nextdeepest)]';
                y=[Watercol_resprates.Avg_mmolC_m3_d(ind_nextshallowest) Watercol_resprates.Avg_mmolC_m3_d(ind_nextdeepest)]';
                % factor in uncertainties
                sy=[Watercol_resprates.Uncert_mmolC_m3_d(ind_nextshallowest) Watercol_resprates.Uncert_mmolC_m3_d(ind_nextdeepest)]';
                % use linfit to run the error-weighted regression
                [a sa cov r]=linfit(x,y,sy);
                % now, estimate value for 50 m and the associated uncertainty
                est_rate_50m=a(2)*50+a(1);
                est_uncert_50m=sa(2)*50+sa(1);
                % now, feed in data, including the additional estimated value at 50 m
                startrow=c_i;
                endrow=c_i+length(ind_50to150);
                Watercol_resprates_50_150.CruiseID(startrow:endrow,1)=cruises(i);
                Watercol_resprates_50_150.CTDStationID(startrow:endrow,1)=CTDstations(j);
                Watercol_resprates_50_150.CTD_depth(startrow:endrow,1)=...
                    [50 Watercol_resprates.CTD_depth(ind_50to150)'];
                Watercol_resprates_50_150.Avg_mmolC_m3_d(startrow:endrow,1)=...
                    [est_rate_50m Watercol_resprates.Avg_mmolC_m3_d(ind_50to150)'];
                Watercol_resprates_50_150.Uncert_mmolC_m3_d(startrow:endrow,1)=...
                    [est_uncert_50m Watercol_resprates.Uncert_mmolC_m3_d(ind_50to150)'];
                Watercol_resprates_50_150.CTD_station_long(startrow:endrow,1)=...
                    [Watercol_resprates.CTD_Long(ind_50to150(1)) Watercol_resprates.CTD_Long(ind_50to150)'];
                Watercol_resprates_50_150.CTD_station_lat(startrow:endrow,1)=...
                    [Watercol_resprates.CTD_Lat(ind_50to150(1)) Watercol_resprates.CTD_Lat(ind_50to150)'];
                Watercol_resprates_50_150.CTD_station_time(startrow:endrow,1)=...
                    [Watercol_resprates.CTD_Deploy_start(ind_50to150(1)) Watercol_resprates.CTD_Deploy_start(ind_50to150)'];
                Watercol_resprates_50_150.ProcStationID(startrow:endrow,1)=...
                    [Watercol_resprates.ProcStationID(ind_50to150(1)) Watercol_resprates.ProcStationID(ind_50to150)'];
            end
        end
        end
    end
end

clear i j c_i;

% now, can actually calculate some depth-integrated water column CR estimates

cruises=unique(Watercol_resprates_50_150.CruiseID);
CTDstations=unique(Watercol_resprates_50_150.CTDStationID);

% since we'll be using a monte carlo analysis to determine the error on our
% depth-integrated estimates, we will have to first generate some vectors
% (n=1000) of random estimates of the respiration rates at each depth and
% station in Watercol_resprates_50_150, assuming observations drawn from normal distribitions

% cycle through each depth and station to generate and store our randomly
% generated, normally distributed estimates of each rate for the monte carlo

Watercol_resprate_montevals_50_150 = nan(length(Watercol_resprates_50_150.Avg_mmolC_m3_d),1000); % preallocate matrix for estimates

for q=1:size(Watercol_resprate_montevals_50_150,1) % cycle through each station and depth
    Watercol_resprate_montevals_50_150(q,:) = Watercol_resprates_50_150.Avg_mmolC_m3_d(q) + Watercol_resprates_50_150.Uncert_mmolC_m3_d(q).*randn(1000,1);
end

Watercolresp_int_monteest = nan(1000,1); % preallocate matrix to hold values for monte carlo-derived estimates as we cycle through

% now, perform integration and estimate error for each using a monte
% carlo

c_i=1;
for i=1:length(cruises)
    for j=1:length(CTDstations)
        ind_thisCTD=find(Watercol_resprates_50_150.CruiseID==cruises(i) & Watercol_resprates_50_150.CTDStationID...
            ==CTDstations(j)); % index of 50 to 150 m BP rate data at this station
        if isfinite(ind_thisCTD) % if data exists for this CTD number
            
            Watercolresp_int_50_150.mmolC_m2_d(c_i,1)=... % compute integral for actual values, and convert to a positive number
                trapz(Watercol_resprates_50_150.CTD_depth(ind_thisCTD),-Watercol_resprates_50_150.Avg_mmolC_m3_d(ind_thisCTD));
            
            % monte carlo to calculate error
            
            % perform the integration 1000 times for each station, using the
            % randomly generated normally distributed estimates for each
            % depth at this station
            
            for r=1:length(Watercol_resprate_montevals_50_150)
                Watercolresp_int_monteest(r) = ...
                    trapz(Watercol_resprates_50_150.CTD_depth(ind_thisCTD),-Watercol_resprate_montevals_50_150(ind_thisCTD,r));
            end
            
            % calculate & store uncertainty for this station using the 1000
            % estimates in our monte carlo
            
            Watercolresp_int_50_150.Uncert_mmolC_m2_d(c_i,1)=...
                std(Watercolresp_int_monteest);
            
            % store other data
            
            Watercolresp_int_50_150.CruiseID(c_i,1)=cruises(i);
            Watercolresp_int_50_150.CTDStationID(c_i,1)=CTDstations(j);
            Watercolresp_int_50_150.CTD_station_long(c_i,1)=Watercol_resprates_50_150.CTD_station_long(ind_thisCTD(1));
            Watercolresp_int_50_150.CTD_station_lat(c_i,1)=Watercol_resprates_50_150.CTD_station_lat(ind_thisCTD(1));
            Watercolresp_int_50_150.CTD_station_time(c_i,1)=Watercol_resprates_50_150.CTD_station_time(ind_thisCTD(1));
            Watercolresp_int_50_150.ProcStationID(c_i,1)=Watercol_resprates_50_150.ProcStationID(ind_thisCTD(1));            
            c_i=c_i+1; % increment our counter
            
            % reset our monte carlo estimate placeholder
            Watercolresp_int_monteest = nan(1000,1);
        end
    end
end

clear i j c_i q r;

% ----------------------------------------------------------------------------
% particle-attached respiration rates (net trap respiration experiments)
% in this section we will eventually calculate k_R, the substrate-specific
% resp rate for sinking particle material
% ----------------------------------------------------------------------------

% calculate mean O2 consumption rates for net trap incubations, store
% in structure Nettrap_resprates, in units of mmol C/m3/d

% first, the bulk of the data from KN207-3

% get some information about the dataset

cruises=unique(Net_trap_incudata.CruiseID);
nettrap_stations=unique(Net_trap_incudata.NetTrapID);

% set the fraction of the total material from each net trap used for the incubations, in
% the case of KN207-3, all incubations were performed in 2 or 3 aliquots from a 1/8 split;
% will need this info to calculate substrate-specific resp rates later

split_frac_2073=.125;

% set BOD bottle size used for these incubations (mL)

BOD_size_2073=350;

c_i=1;
for i=1:length(cruises)
    for j=1:length(nettrap_stations);
            ind=find(Net_trap_incudata.CruiseID==cruises(i) & ...
                Net_trap_incudata.NetTrapID==nettrap_stations(j));
            if numel(ind)>1 % a net trap resp profile was made for this 
                % cruise and net trap ID
                % first, ascertain the number of replicates
                numreps=sum(sum(isnan(Net_trap_incudata.O2_conc_uM(ind,:)))<length(ind));
                % calculate O2 consumption rates for each set of replicate
                % measurements
                O2data_j=Net_trap_incudata.O2_conc_uM(ind,:);
                readtimes_j=Net_trap_incudata.Read_time(ind,:);
                for k=1:numreps % cycle through replicates, making sure we
                    % don't have any NaNs before calculating the regression
                    O2data_k=O2data_j(:,k);
                    readtimes_k=readtimes_j;
                    ind_nan=find(isnan(O2data_k));
                    O2data_k(ind_nan)=[];
                    readtimes_k(ind_nan)=[];
                    % do the regression
                    x=readtimes_k;
                    y=O2data_k;
                    sy=1; % unweighted since we don't know what uncertaintites are
                    [a sa cov r]=linfit(x,y,sy);
                    rate(k)=a(2);
                    rate_uncert(k)=sa(2);
                    corrcoeff(k)=r;
                end
                % find mean rate and uncertainty
                O2_cons_rate_uM_O2_d_avg=mean(rate);
                O2_cons_rate_uM_O2_d_uncert=sqrt(sum((rate_uncert).^2))/length(rate_uncert);
                if isinf(O2_cons_rate_uM_O2_d_uncert)
                % use std dev of replicates as uncertainty estimate in cases where the
                    % regressions didn't have more than two data points to
                    % play with
                    O2_cons_rate_uM_O2_d_uncert=std(rate);
                end
                Nettrap_resprates.CruiseID(c_i,1)=cruises(i);
                Nettrap_resprates.NetTrapID(c_i,1)=nettrap_stations(j);
                Nettrap_resprates.Num_incu_reps(c_i,1)=numreps;
                Nettrap_resprates.Trap_material_split_frac(c_i,1)=split_frac_2073;
                Nettrap_resprates.BOD_bot_size_mL(c_i,1)=BOD_size_2073;
                Nettrap_resprates.Avg_mmolC_m3_d(c_i,1)=O2_cons_rate_uM_O2_d_avg*(1/1000)*(1000)*RQ;
                Nettrap_resprates.Uncert_mmolC_m3_d(c_i,1)=O2_cons_rate_uM_O2_d_uncert*(1/1000)*(1000)*RQ;
                Nettrap_resprates.R2_mean(c_i,1)=mean(corrcoeff);
                 c_i=c_i+1;
                 clear rate rate_uncert corrcoeff;
            end
        end
end

clear i j k c_i;

% now, average the already calculated net trap resp rates from KN207-1 (Bethanie's
% experimental controls), convert to same units, and flow into structure

% in the case of Bethanie's KN207-1 controls, the split fraction differed by
% experiment; will pull data from fractions calculated by information
% Bethanie provided in email dtd 12/12/13

% set BOD bottle size used for these incubations (mL)

BOD_size_2071=150;

c_i=10;

for i=1:length(Net_trap_incudata2071.CruiseID)
    numreps=length(Net_trap_incudata2071.O2_cons_rates_uM_hr(i,1:3));
    O2_cons_rate_uM_O2_d_avg=-mean(Net_trap_incudata2071.O2_cons_rates_uM_hr(i,1:3))*24;
    O2_cons_rate_uM_O2_d_uncert=std(Net_trap_incudata2071.O2_cons_rates_uM_hr(i,1:3))*24;
    Nettrap_resprates.CruiseID(c_i,1)=Net_trap_incudata2071.CruiseID(i);
    Nettrap_resprates.NetTrapID(c_i,1)=Net_trap_incudata2071.NetTrapID(i);
    Nettrap_resprates.Num_incu_reps(c_i,1)=numreps;
    Nettrap_resprates.Trap_material_split_frac(c_i,1)=Net_trap_incudata2071.Split_frac_2071(i);
    Nettrap_resprates.BOD_bot_size_mL(c_i,1)=BOD_size_2071;
    Nettrap_resprates.Avg_mmolC_m3_d(c_i,1)=O2_cons_rate_uM_O2_d_avg*(1/1000)*(1000)*RQ;
    Nettrap_resprates.Uncert_mmolC_m3_d(c_i,1)=O2_cons_rate_uM_O2_d_uncert*(1/1000)*(1000)*RQ;
    c_i=c_i+1;
end

clear c_i i;

% insert some relevant trap information (depth, process station etc) into
% Nettrap_resprates

for i=1:length(Nettrap_resprates.NetTrapID)
    ind=find(Net_trap_data.CruiseID==Nettrap_resprates.CruiseID(i) & ...
        Net_trap_data.NetTrapID==Nettrap_resprates.NetTrapID(i));
    Nettrap_resprates.ProcStationID(i,1)=Net_trap_data.ProcStationID(ind);
    Nettrap_resprates.Deploy_depth(i,1)=Net_trap_data.Deploy_depth(ind);
    Nettrap_resprates.ProcStation_Lat(i,1)=Net_trap_data.ProcStation_Lat(ind);
    Nettrap_resprates.ProcStation_Long(i,1)=Net_trap_data.ProcStation_Long(ind);
    Nettrap_resprates.Deploy_start(i,1)=Net_trap_data.Deploy_start(ind);
    Nettrap_resprates.Deploy_end(i,1)=Net_trap_data.Deploy_end(ind);
    Nettrap_resprates.Deploy_duration(i,1)=Net_trap_data.Deploy_duration(ind);
end

clear i ind;

% estimate some uncertainties for net trap POC concentrations that don't
% already have std dev's attached from analytical procedure, assuming that
% std dev scales proportionally with signal; we will use available data
% from KN207-3 samples for this

% first, perform a regression of POC signal on std dev to get a
% relationship we can use to make estimates; exclude the very high value
% from 2073 trap 9 (since none of the traps for which we want estimates had
% POC values anywhere near this high) and exclude an apparent outlier, 2073
% trap 38 (rel std dev in this case was > 0.5)

ind=find(isfinite(Net_trap_data.POC_mgC_per_trap) & ...
    isfinite(Net_trap_data.Std_dev_POC_mgC_per_trap) & Net_trap_data.POC_mgC_per_trap<100 & ...
    Net_trap_data.Std_dev_POC_mgC_per_trap./Net_trap_data.POC_mgC_per_trap<.25);
x=Net_trap_data.POC_mgC_per_trap(ind);
y=Net_trap_data.Std_dev_POC_mgC_per_trap(ind);
[a sa cov r]=linfit(x,y,1);

% now, estimate some uncertainties for those trap POC values for which we don't have any
% anaytical uncertainties

ind_for_est=find(isfinite(Net_trap_data.POC_mgC_per_trap) & ...
    isnan(Net_trap_data.Std_dev_POC_mgC_per_trap));
x=Net_trap_data.POC_mgC_per_trap(ind_for_est);
yf=a(2)*x+a(1); % make estimates
Net_trap_data.Std_dev_POC_mgC_per_trap(ind_for_est)=yf; % insert back into Net_trap_data.Std_dev_POC_mgC_per_trap

% final step before generating substrate-specific respiration rates: find
% the "closest" (geographically and temporally) water column respiration
% rate for the applicable net trap deployment; will need this rate to
% subtract from the volumetric net trap particle resp rate before scaling

% this is essentially the parameter "r_control" in A. McDonnell (2011) Ph.D.
% thesis, Ch. 5, eq. 1

% two scenarios:
% (1 - best case): use an actual measured respiration rate from incubations or the
% PHORCYS
% (2 - a less ideal case, but necessary for all 300 m deployments and for 
% most of KN207-1 deployments): estimate a rate from a model based on the data we do have

% for case 1, can search by depth and then by timestamp since the cruise proceeded along a contiguous
% track

% prepopulate vectors for the data

Nettrap_resprates.Concurrent_watercol_rate_mmol_C_m3_d=nan(length(Nettrap_resprates.NetTrapID),1);
Nettrap_resprates.Uncert_concurrent_watercol_rate_mmol_C_m3_d=nan(length(Nettrap_resprates.NetTrapID),1);

for i=1:length(Nettrap_resprates.NetTrapID)
    ind=find(abs(Nettrap_resprates.Deploy_depth(i)-Watercol_resprates.CTD_depth(:))<10 & ...
        abs(Nettrap_resprates.Deploy_start(i)-Watercol_resprates.CTD_Deploy_start(:))<3);
    % find the water column incubations conducted within 10 m depth and
    % closest in time to this net trap deployment
    if isfinite(ind) % there is at least one water column incubation meeting the criteria (case #1)
        if length(ind)>1 % if there are more than two incubations in close proximity, choose the closest
            ind_mintime=find(min(abs(Nettrap_resprates.Deploy_depth(i)-Watercol_resprates.CTD_depth(ind))));
            ind=ind(ind_mintime);
        end
        Nettrap_resprates.Concurrent_watercol_rate_mmol_C_m3_d(i,1)=Watercol_resprates.Avg_mmolC_m3_d(ind);
        Nettrap_resprates.Uncert_concurrent_watercol_rate_mmol_C_m3_d(i,1)=Watercol_resprates.Uncert_mmolC_m3_d(ind);
    else % there isn't a water column incubation meeting the critera, so we
        % will have to get results from some sort of model (case #2)
    end
end

% now, generate substrate-specific respiration rates for particle material
% in each net trap, in units of mmol C respired/d/mmol C substrate

for i=1:length(Nettrap_resprates.NetTrapID)
    
    % lookup relevant trap total POC measurement and uncertainty from Net_trap_data.POC_mgC_per_trap
    % and Net_trap_data.Std_dev_POC_mgC_per_trap
    
    % where we don't have a POC measurement for a particular trap, will
    % substitute data from a different trap so long as the other deployment was
    % concurrent spatially and roughly concurrent temporally 
    
    ind = find(Net_trap_data.CruiseID==Nettrap_resprates.CruiseID(i) & ...
        Net_trap_data.NetTrapID==Nettrap_resprates.NetTrapID(i));
    
    if isfinite(Net_trap_data.POC_mgC_per_trap(ind)) % net trap POC data exists for this trap
        
        POC_mmolC_per_trap=Net_trap_data.POC_mgC_per_trap(ind)/12.01;
        Std_dev_POC_mmolC_per_trap=Net_trap_data.Std_dev_POC_mgC_per_trap(ind)/12.01;
        
    else % we'll have to use POC data from a proximate trap for which data was obtained
        
        ind_nan = find(isnan(Net_trap_data.POC_mgC_per_trap));
        
        ind_prox = find(Net_trap_data.CruiseID==Nettrap_resprates.CruiseID(i) & ...
        Net_trap_data.Deploy_depth==Nettrap_resprates.Deploy_depth(i) & ...
        Net_trap_data.ProcStationID==Nettrap_resprates.ProcStationID(i) & ...
        abs(Net_trap_data.ProcStation_Lat-Nettrap_resprates.ProcStation_Lat(i))<0.2 & ...
        abs(Net_trap_data.ProcStation_Long-Nettrap_resprates.ProcStation_Long(i))<0.2 & ...
        abs(Net_trap_data.Deploy_end-Nettrap_resprates.Deploy_end(i))<0.2 & ...
        abs(Net_trap_data.Deploy_duration-Nettrap_resprates.Deploy_duration(i))<0.2);
    
        ind_exclude = find(ismember(ind_prox,ind_nan));
        
        ind_sub = ind_prox;
        
        ind_sub(ind_exclude)=[]; % we have a winner!
        
        % but, will have to adjust the substitute POC value for the difference deployment
        % duration
        
        r_t = Net_trap_data.Deploy_duration(ind)/Net_trap_data.Deploy_duration(ind_sub);
        
        POC_mmolC_per_trap=(Net_trap_data.POC_mgC_per_trap(ind_sub)*r_t)/12.01;
        Std_dev_POC_mmolC_per_trap=(Net_trap_data.Std_dev_POC_mgC_per_trap(ind_sub)*r_t)/12.01;
        
    end
    
    % use paired water column rate to obtain a "particle only" volumetric
    % rate, IF there is a paired water column rate available; otherwise
    % don't correct it
    
    if isfinite(Nettrap_resprates.Concurrent_watercol_rate_mmol_C_m3_d(i))
        
        resp_corrected_mmolC_m3_d=Nettrap_resprates.Avg_mmolC_m3_d(i)-...
            Nettrap_resprates.Concurrent_watercol_rate_mmol_C_m3_d(i);
        
        resp_corrected_mmolC_m3_d_uncert=sqrt(Nettrap_resprates.Uncert_mmolC_m3_d(i)^2+...
            Nettrap_resprates.Uncert_concurrent_watercol_rate_mmol_C_m3_d(i)^2);
    else
        
        resp_corrected_mmolC_m3_d=Nettrap_resprates.Avg_mmolC_m3_d(i);
        
        resp_corrected_mmolC_m3_d_uncert=Nettrap_resprates.Uncert_mmolC_m3_d(i);
        
    end 
    
    % calculate per trap respiration rates and uncertaintites (still
    % volumetric)
    
    mmolC_resp_per_trap_d = resp_corrected_mmolC_m3_d*...
        Nettrap_resprates.BOD_bot_size_mL(i)*(1/10^6)*Nettrap_resprates.Num_incu_reps(i)...
        *1/Nettrap_resprates.Trap_material_split_frac(i);
    
    Std_dev_mmolC_resp_per_trap_d = resp_corrected_mmolC_m3_d_uncert*...
        Nettrap_resprates.BOD_bot_size_mL(i)*(1/10^6)*Nettrap_resprates.Num_incu_reps(i)...
        *1/Nettrap_resprates.Trap_material_split_frac(i);
    
    % calculate substrate-specific rate by combining trap POC data and the
    % per trap volumetric resp rate
    
    Sub_resp_mmolCresp_mmolCsub_d=mmolC_resp_per_trap_d/POC_mmolC_per_trap;
    
    % uncertainty here obeys rules for propagation of error by
    % multiplication/division (* constants)
    
    Sub_resp_mmolCresp_mmolCsub_d_uncert=abs(Sub_resp_mmolCresp_mmolCsub_d)*...
        sqrt((Std_dev_mmolC_resp_per_trap_d/mmolC_resp_per_trap_d)^2+...
        (Std_dev_POC_mmolC_per_trap/POC_mmolC_per_trap)^2);
    
    % insert results into Nettrap_resprates
    
    Nettrap_resprates.mmolC_mmolCsub_d(i,1)=Sub_resp_mmolCresp_mmolCsub_d;
    Nettrap_resprates.Uncert_mmolC_mmolCsub_d(i,1)=Sub_resp_mmolCresp_mmolCsub_d_uncert;
end

clear i ind ind_for_est;

% ----------------------------------------------------------------------------
% depth-integrated estimates of water column BP from 50 to 150 m
%
% *** note that the script Depth_int_WCR_BP.m contains a much
% better way of doing this (and which doesn't require interpolation) ***
% ----------------------------------------------------------------------------

% first, will need to interpolate or approximate a value at 50 m in cases where a bottle wasn't fired at that depth

cruises=unique(BPdata.CruiseID);
CTDstations=unique(BPdata.CTDStationID);

clear BP_rates_50_150;

for i=1:length(cruises)
    for j=1:length(CTDstations)
        ind_thisCTD=find(BPdata.CruiseID==cruises(i) & BPdata.CTDStationID...
            ==CTDstations(j)); % index of all analyzed depths at this CTD station
        if isfinite(ind_thisCTD)
            % first check to see if there's
            % data for this station at all, and if there was a measurement made at 45 m or below; if not, no point in doing any of
            % this, so bail
        if isfinite(BPdata.Trit_leu_uptake_pmol_L_hr(ind_thisCTD))...
                & max(BPdata.CTD_depth(ind_thisCTD))>=45
            if exist('BP_rates_50_150') % define matrix insertion point (should be first row after end of existing BP_int_50_150 data)
                c_i=length(BP_rates_50_150.CruiseID)+1;
            else
                c_i=1; % need this here for the first run, when no data has been inserted yet
            end
            if any(ismember(BPdata.CTD_depth(ind_thisCTD),50)) % we have a rate for 50 m exactly
                ind_50to150=find(BPdata.CruiseID==cruises(i) & BPdata.CTDStationID...
                ==CTDstations(j) & BPdata.CTD_depth>=50);
                % index of all measurements from 50m to 150m, inclusive
                % now, feed in data
                startrow=c_i;
                endrow=c_i+length(ind_50to150)-1;
                BP_rates_50_150.CruiseID(startrow:endrow,1)=cruises(i);
                BP_rates_50_150.CTDStationID(startrow:endrow,1)=CTDstations(j);
                BP_rates_50_150.CTD_depth(startrow:endrow,1)=BPdata.CTD_depth(ind_50to150);
                BP_rates_50_150.Trit_leu_uptake_pmol_L_hr(startrow:endrow,1)=BPdata.Trit_leu_uptake_pmol_L_hr(ind_50to150);
                BP_rates_50_150.Std_dev_trit_leu_uptake_pmol_L_hr(startrow:endrow,1)=BPdata.Std_dev_trit_leu_uptake_pmol_L_hr(ind_50to150);
                BP_rates_50_150.CTD_station_long(startrow:endrow,1)=BPdata.CTD_station_long(ind_50to150);
                BP_rates_50_150.CTD_station_lat(startrow:endrow,1)=BPdata.CTD_station_lat(ind_50to150);
                BP_rates_50_150.CTD_station_time(startrow:endrow,1)=BPdata.CTD_station_time(ind_50to150);
            elseif any(ismember(BPdata.CTD_depth(ind_thisCTD),[45:1:55])) % we have a value within 5 m 
                % of 50 m which we can directly substitute
                ind_prox=find(BPdata.CruiseID==cruises(i) & BPdata.CTDStationID...
                ==CTDstations(j) & BPdata.CTD_depth>=45 & BPdata.CTD_depth<=55); % make a list
                % if there's more than 1, select the closest
                dz=abs(BPdata.CTD_depth(ind_prox)-50);
                ind_best=find(min(dz));
                ind_sub=ind_prox(ind_best);
                ind_50to150=find(BPdata.CruiseID==cruises(i) & BPdata.CTDStationID...
                ==CTDstations(j) & BPdata.CTD_depth>=50);
                % now, feed in data, including the additional estimated value at 50 m
                startrow=c_i;
                endrow=c_i+length(ind_50to150);
                BP_rates_50_150.CruiseID(startrow:endrow,1)=cruises(i);
                BP_rates_50_150.CTDStationID(startrow:endrow,1)=CTDstations(j);
                BP_rates_50_150.CTD_depth(startrow:endrow,1)=...
                    [50 BPdata.CTD_depth(ind_50to150)'];
                BP_rates_50_150.Trit_leu_uptake_pmol_L_hr(startrow:endrow,1)=...
                    [BPdata.Trit_leu_uptake_pmol_L_hr(ind_sub) BPdata.Trit_leu_uptake_pmol_L_hr(ind_50to150)'];
                BP_rates_50_150.Std_dev_trit_leu_uptake_pmol_L_hr(startrow:endrow,1)=...
                    [BPdata.Std_dev_trit_leu_uptake_pmol_L_hr(ind_sub) BPdata.Std_dev_trit_leu_uptake_pmol_L_hr(ind_50to150)'];
                BP_rates_50_150.CTD_station_long(startrow:endrow,1)=...
                    [BPdata.CTD_station_long(ind_sub) BPdata.CTD_station_long(ind_50to150)'];
                BP_rates_50_150.CTD_station_lat(startrow:endrow,1)=...
                    [BPdata.CTD_station_lat(ind_sub) BPdata.CTD_station_lat(ind_50to150)'];
                BP_rates_50_150.CTD_station_time(startrow:endrow,1)=...
                    [BPdata.CTD_station_time(ind_sub) BPdata.CTD_station_time(ind_50to150)'];
            else % there isn't a rate measurement at anywhere near 50 m, will have to interpolate
                ind_50to150=find(BPdata.CruiseID==cruises(i) & BPdata.CTDStationID...
                ==CTDstations(j) & BPdata.CTD_depth>=50);
                dz=BPdata.CTD_depth(ind_thisCTD)-50;
                ind_shallower=find(dz<0); % all measurements shallower than 50 m
                ind_shallower=[ind_thisCTD(ind_shallower) dz(ind_shallower)];
                [C,I]=max(ind_shallower(:,2));
                ind_nextshallowest=ind_shallower(I,1);
                ind_deeper=find(dz>0); % all measurements deeper than 50 m
                ind_deeper=[ind_thisCTD(ind_deeper) dz(ind_deeper)];
                [C,I]=min(ind_deeper(:,2));
                ind_nextdeepest=ind_deeper(I,1);
                % now that we have the two values bracketing 50 m, can
                % interpolate linearly
                x=[BPdata.CTD_depth(ind_nextshallowest) BPdata.CTD_depth(ind_nextdeepest)]';
                y=[BPdata.Trit_leu_uptake_pmol_L_hr(ind_nextshallowest) BPdata.Trit_leu_uptake_pmol_L_hr(ind_nextdeepest)]';
                % factor in uncertainties
                sy=[BPdata.Std_dev_trit_leu_uptake_pmol_L_hr(ind_nextshallowest) BPdata.Std_dev_trit_leu_uptake_pmol_L_hr(ind_nextdeepest)]';
                % use linfit to run the error-weighted regression
                [a sa cov r]=linfit(x,y,sy);
                % now, estimate value for 50 m and the associated uncertainty
                est_rate_50m=a(2)*50+a(1);
                est_uncert_50m=sa(2)*50+sa(1);
                % now, feed in data, including the additional estimated value at 50 m
                startrow=c_i;
                endrow=c_i+length(ind_50to150);
                BP_rates_50_150.CruiseID(startrow:endrow,1)=cruises(i);
                BP_rates_50_150.CTDStationID(startrow:endrow,1)=CTDstations(j);
                BP_rates_50_150.CTD_depth(startrow:endrow,1)=...
                    [50 BPdata.CTD_depth(ind_50to150)'];
                BP_rates_50_150.Trit_leu_uptake_pmol_L_hr(startrow:endrow,1)=...
                    [est_rate_50m BPdata.Trit_leu_uptake_pmol_L_hr(ind_50to150)'];
                BP_rates_50_150.Std_dev_trit_leu_uptake_pmol_L_hr(startrow:endrow,1)=...
                    [est_uncert_50m BPdata.Std_dev_trit_leu_uptake_pmol_L_hr(ind_50to150)'];
                BP_rates_50_150.CTD_station_long(startrow:endrow,1)=...
                    [BPdata.CTD_station_long(ind_50to150(1)) BPdata.CTD_station_long(ind_50to150)'];
                BP_rates_50_150.CTD_station_lat(startrow:endrow,1)=...
                    [BPdata.CTD_station_lat(ind_50to150(1)) BPdata.CTD_station_lat(ind_50to150)'];
                BP_rates_50_150.CTD_station_time(startrow:endrow,1)=...
                    [BPdata.CTD_station_time(ind_50to150(1)) BPdata.CTD_station_time(ind_50to150)'];

            end
        end
        end
    end
end

clear i j c_i;

% now, can actually calculate some depth-integrated BP estimates

cruises=unique(BP_rates_50_150.CruiseID);
CTDstations=unique(BP_rates_50_150.CTDStationID);

% since we'll be using a monte carlo analysis to determine the error on our
% depth-integrated estimates, we will have to first generate some vectors
% (n=1000) of random estimates of the BP rates at each depth and
% station in BP_rates_50_150, assuming observations drawn from normal distribitions

% cycle through each depth and station to generate and store our randomly
% generated, normally distributed estimates of each rate for the monte carlo

BP_montevals_50_150 = nan(length(BP_rates_50_150.Trit_leu_uptake_pmol_L_hr),1000); % preallocate matrix for estimates

for q=1:size(BP_montevals_50_150,1) % cycle through each station and depth
    BP_montevals_50_150(q,:) = BP_rates_50_150.Trit_leu_uptake_pmol_L_hr(q) + BP_rates_50_150.Std_dev_trit_leu_uptake_pmol_L_hr(q).*randn(1000,1);
end

BP_int_monteest = nan(1000,1); % preallocate matrix to hold values for monte carlo-derived estimates as we cycle through

% now, perform integration and estimate error for each using a monte
% carlo

c_i=1;
for i=1:length(cruises)
    for j=1:length(CTDstations)
        ind_thisCTD=find(BP_rates_50_150.CruiseID==cruises(i) & BP_rates_50_150.CTDStationID...
            ==CTDstations(j)); % index of 50 to 150 m BP rate data at this station
        if isfinite(ind_thisCTD) % if data exists for this CTD number
            
            BP_int_50_150.Trit_leu_uptake_mmol_m2_d(c_i,1)=... % compute integral for actual values, after converting per L rates to rates per m3
                24*(1/10^9)*trapz(BP_rates_50_150.CTD_depth(ind_thisCTD),BP_rates_50_150.Trit_leu_uptake_pmol_L_hr(ind_thisCTD)*1000);
            
            % monte carlo to calculate error
            
            % perform the integration 1000 times for each station, using the
            % randomly generated normally distributed estimates for each
            % depth at this station
            
            for r=1:length(BP_montevals_50_150)
                BP_int_monteest(r) = ...
                    24*(1/10^9)*trapz(BP_rates_50_150.CTD_depth(ind_thisCTD),BP_montevals_50_150(ind_thisCTD,r)*1000);
            end
            
            % calculate & store uncertainty for this station using the 1000
            % estimates in our monte carlo
            
            BP_int_50_150.Uncert_trit_leu_uptake_mmol_m2_d(c_i,1)=...
                std(BP_int_monteest);
                        
            % store other data
            
            BP_int_50_150.CruiseID(c_i,1)=cruises(i);
            BP_int_50_150.CTDStationID(c_i,1)=CTDstations(j);
            BP_int_50_150.CTD_station_long(c_i,1)=BP_rates_50_150.CTD_station_long(ind_thisCTD(1));
            BP_int_50_150.CTD_station_lat(c_i,1)=BP_rates_50_150.CTD_station_lat(ind_thisCTD(1));
            BP_int_50_150.CTD_station_time(c_i,1)=BP_rates_50_150.CTD_station_time(ind_thisCTD(1));
            c_i=c_i+1; % increment our counter
            
            % reset our monte carlo estimate placeholder
            BP_int_monteest = nan(1000,1);
        end
    end
end

clear i j c_i q r;

% associate integrated BP estimates with the relevant PIT trap deployment (if proximate)
% so we can compare to observed POC fluxes and the calculations (earlier section) based on
% respiration rates

% prepopulate target matrix
BP_int_50_150.ProcStationID=nan(length(BP_int_50_150.Trit_leu_uptake_mmol_m2_d),1);

% perform a query with proximity critera ± 0.5 deg lat and long, and ± 12
% hrs
for i=1:length(BP_int_50_150.ProcStationID)
    ind=find(PITdeployments.CruiseID==BP_int_50_150.CruiseID(i) & ...
        abs(PITdeployments.ProcStation_Lat-BP_int_50_150.CTD_station_lat(i))<=0.5 & ...
        abs(PITdeployments.ProcStation_Long-BP_int_50_150.CTD_station_long(i))<=0.5 & ...
        BP_int_50_150.CTD_station_time(i)>=(PITdeployments.Deploy_start-0.5) & ...
        BP_int_50_150.CTD_station_time(i)<=(PITdeployments.Deploy_end+0.5));
    if isfinite(ind) % there's a match!
        BP_int_50_150.ProcStationID(i,1)=PITdeployments.ProcStationID(ind(1));
    end
end
clear i ind;

% calculate averages by process station

% prepopulate destination vectors

Int_BP_by_station.Avg_trit_leu_uptake_mmol_m2_d=nan(6,1);
Int_BP_by_station.Uncert_trit_leu_uptake_mmol_m2_d=nan(6,1);
Int_BP_by_station.CruiseID=nan(6,1);
Int_BP_by_station.ProcStationID=nan(6,1);

cruises=unique(PITdeployments.CruiseID);
procstas=unique(PITdeployments.ProcStationID);

c_i=1;
for i=1:length(cruises)
    for j=1:length(procstas)
        ind_exists=find(PITdeployments.CruiseID==cruises(i) & ...
            PITdeployments.ProcStationID==procstas(j)); % determine if this combination
        % of CruiseID and ProcStaID even exists
        if isfinite(ind_exists)
        ind=find(BP_int_50_150.CruiseID==cruises(i) & ...
            BP_int_50_150.ProcStationID==procstas(j));
        avg_leu_uptake=mean(BP_int_50_150.Trit_leu_uptake_mmol_m2_d(ind));
        uncert_uptake=sqrt(sum((BP_int_50_150.Uncert_trit_leu_uptake_mmol_m2_d(ind)).^2))/length(BP_int_50_150.Uncert_trit_leu_uptake_mmol_m2_d(ind));
        Int_BP_by_station.Avg_trit_leu_uptake_mmol_m2_d(c_i,1)=avg_leu_uptake;
        Int_BP_by_station.Uncert_trit_leu_uptake_mmol_m2_d(c_i,1)=uncert_uptake;
        Int_BP_by_station.CruiseID(c_i,1)=cruises(i);
        Int_BP_by_station.ProcStationID(c_i,1)=procstas(j);
        c_i=c_i+1;
        end
    end
end
clear i j ind;

% ----------------------------------------------------------------------------
% depth-integrated estimates of water column BCD (BP + WCR) from 50 to 150
% m, at least at PS-2, PS-3, PS-4
% ----------------------------------------------------------------------------

cruises=unique(Watercolresp_int_50_150.CruiseID);
ProcStations=unique(Watercolresp_int_50_150.ProcStationID);

ID=1; % set ID for this exercise

c_i=1; % set a counter for our insertion point in the destination table

for i=1:length(cruises)
    for j=1:length(ProcStations)
        ind_thisProcSta_WCR=find(Watercolresp_int_50_150.CruiseID==cruises(i) & Watercolresp_int_50_150.ProcStationID...
            ==ProcStations(j)); % index of depth-integrated WCR data at this process station
        ind_thisProcSta_BP=find(Int_BP_by_station.CruiseID==cruises(i) & Int_BP_by_station.ProcStationID...
            ==ProcStations(j)); % index of depth-integrated BP data at this process station

        if isfinite(ind_thisProcSta_WCR) % if WCR data exists for this process station
            
            WC_BCD_int_50_150.mg_C_m2_d(c_i,1)=... % calculate WC BCD in mg C/m2/d, converting units as appropriate
                Int_BP_by_station.Avg_trit_leu_uptake_mmol_m2_d(ind_thisProcSta_BP).*MW_leu*(1/leu_BPAA)*C_prot*ID +...
                Watercolresp_int_50_150.mmolC_m2_d(ind_thisProcSta_WCR)*12.01;
            
            WC_BCD_int_50_150.Uncert_mg_C_m2_d(c_i,1)=... % calculate uncertainty, converting units as appropriate
                sqrt((Int_BP_by_station.Uncert_trit_leu_uptake_mmol_m2_d(ind_thisProcSta_BP).*MW_leu*(1/leu_BPAA)*C_prot*ID)^2 +...
                (Watercolresp_int_50_150.Uncert_mmolC_m2_d(ind_thisProcSta_WCR)*12.01)^2);
            
            % store other data
            
            WC_BCD_int_50_150.CruiseID(c_i,1)=cruises(i);
            WC_BCD_int_50_150.ProcStationID(c_i,1)=ProcStations(j);
            c_i=c_i+1; % increment our counter
            
        end
    end
end

clear c_i i j

%% analysis: first, define constants, data locations, and parameter value ranges for each station

% global constants

k_DOCZ_axes = [0 200 10^-5 1]; % axis ranges for plots in analysis 1a

% define ranges for sensitity analyses

% k_DOCZ = logspace(-5,0,100); % rate constant for loss of OC to water column in model 1a
% W_avg=[.01:.1:200]; % sinking rates, in m/d

k_DOCZ = logspace(-2,0,50); % rate constant for loss of OC to water column in model 1a
W_avg=[2.01:.1:52]; % sinking rates, in m/d

% assumed constants for BP calcs

leu_BPAA = 0.073; % mol fraction leucine : total bacterial amino acids, Simon & Azam 1989
leu_BPAA_uncert = 0.0191; % associated uncertainty, Simon & Azam 1989
C_prot = 0.86; % mass ratio cellular carbon : protein for marine bacteria, Simon & Azam 1989
MW_leu = 131.2; % molecular weight of leucine


%% KN207-1, process station 1

% define data locations and plot parameters for this station

PIT_loc=[1 2 3]; % location of PIT trap data
Flux_est_netresp_loc=[1 2 3]; % location of net trap-based flux estimates
Nettrap_resprate_loc=[10 11 10 11 10 11]; % location of net trap-derived respiration rates
% specify the trap material incubations that should be used for calculations
% at this station, by depth in order: 50m, 50m, 150m, 150m, 300m, 300m
% if only one incubation exists for each depth, specify it three times
Serial_ProcStaID=1; % 1-6, in order of process station from KN207-1 to KN207-3
cmap_axis_offset=.005;
Labelstep=0.2; % increment at which to label monte carlo runs in plot 2
Labelstep_2=0.2; % increment at which to label monte carlo runs in plot 4
heatcontours_WavgkDOCZ_150 = [-25 -10 -5 0 5 6]; % contours at which to place labels in W_avg-k_DOCZ sensitivity analysis heatmap for 150 m
heatcontours_WavgkDOCZ_300 = [-25 -5 -2.5 0 0.5 1 2.5]; % contours at which to place labels in W_avg-k_DOCZ sensitivity analysis heatmap for 300 m

%% KN207-1, process station 2

% define data locations and plot parameters for this station

PIT_loc=[4 5 6]; % location of PIT trap data
Flux_est_netresp_loc=[4 5 6]; % location of net trap-based flux estimates
Nettrap_resprate_loc=[12 12 12 12 12 12]; % location of net trap-derived respiration rates
% specify the trap material incubations that should be used for calculations
% at this station, by depth in order: 50m, 50m, 150m, 150m, 300m, 300m
% if only one incubation exists for each depth, specify it three times
Serial_ProcStaID=2; % 1-6, in order of process station from KN207-1 to KN207-3
cmap_axis_offset=.005;
Labelstep=0.3; % increment at which to label monte carlo runs in plot
Labelstep_2=0.5; % increment at which to label monte carlo runs in plot 4
heatcontours_WavgkDOCZ_150 = [-25 -15 0 10 20 22.5]; % contours at which to place labels in W_avg-k_DOCZ sensitivity analysis heatmap for 150 m
heatcontours_WavgkDOCZ_300 = [-12 -5 0 5 7.5]; % contours at which to place labels in W_avg-k_DOCZ sensitivity analysis heatmap for 300 m
                  
%% KN207-3, process station 2

% define data locations and plot parameters for this station

PIT_loc=[10 11 12]; % location of PIT trap data
Flux_est_netresp_loc=[10 11 12]; % location of net trap-based flux estimates
Nettrap_resprate_loc=[2 3 2 3 4 4]; % location of net trap-derived respiration rates
% specify the trap material incubations that should be used for calculations
% at this station, by depth in order: 50m, 50m, 150m, 150m, 300m, 300m
% if only one incubation exists for each depth, specify it three times
Serial_ProcStaID=4; % 1-6, in order of process station from KN207-1 to KN207-3
cmap_axis_offset=.005;
Labelstep=0.5; % increment at which to label monte carlo runs in plot
Labelstep_2=0.5; % increment at which to label monte carlo runs in plot 4
heatcontours_WavgkDOCZ_150 = [-17.5 -15 0 25 35 40]; % contours at which to place labels in W_avg-k_DOCZ sensitivity analysis heatmap for 150 m
heatcontours_WavgkDOCZ_300 = [-9.5 -5 0 4 5]; % contours at which to place labels in W_avg-k_DOCZ sensitivity analysis heatmap for 300 m

%% KN207-3, process station 3

% define data locations and plot parameters for this station

PIT_loc=[13 14 15]; % location of PIT trap data
Flux_est_netresp_loc=[13 14 15]; % location of net trap-based flux estimates
Nettrap_resprate_loc=[5 5 6 6 7 7]; % location of net trap-derived respiration rates
% specify the trap material incubations that should be used for calculations
% at this station, by depth in order: 50m, 50m, 150m, 150m, 300m, 300m
% if only one incubation exists for each depth, specify it three times
Serial_ProcStaID=5; % 1-6, in order of process station from KN207-1 to KN207-3
cmap_axis_offset=.005;
Labelstep=1; % increment at which to label monte carlo runs in plot
Labelstep_2=1; % increment at which to label monte carlo runs in plot 4
heatcontours_WavgkDOCZ_150 = [-175 -50 0 25 35 40]; % contours at which to place labels in W_avg-k_DOCZ sensitivity analysis heatmap for 150 m
heatcontours_WavgkDOCZ_300 = [-125 -50 0 50 65]; % contours at which to place labels in W_avg-k_DOCZ sensitivity analysis heatmap for 300 m

%% KN207-3, process station 4

% define data locations and plot parameters for this station

PIT_loc=[16 17 18]; % location of PIT trap data
Flux_est_netresp_loc=[16 17 18]; % location of net trap-based flux estimates
Nettrap_resprate_loc=[8 9 8 9 9 9]; % location of net trap-derived respiration rates
% specify the trap material incubations that should be used for calculations
% at this station, by depth in order: 50m, 50m, 150m, 150m, 300m, 300m
% if only one incubation exists for each depth, specify it three times
Serial_ProcStaID=6; % 1-6, in order of process station from KN207-1 to KN207-3
cmap_axis_offset=0.005;
Labelstep=1.4; % increment at which to label monte carlo runs in plot
Labelstep_2=0.5; % increment at which to label monte carlo runs in plot 4
heatcontours_WavgkDOCZ_150 = [-125 -50 0 10 15]; % contours at which to place labels in W_avg-k_DOCZ sensitivity analysis heatmap for 150 m
heatcontours_WavgkDOCZ_300 = [-100 -50 0 30 45]; % contours at which to place labels in W_avg-k_DOCZ sensitivity analysis heatmap for 300 m

%% now, perform analyses and generate plots for the given inputs

clf;

% top set of plots:

% modeled particle flux attentuation as a function of measured particle & water column respiration rates

% this approach requires some assumptions about sinking rate

% will perform a sensitivity analysis to assess effect of a range of possible
% sinking rates on modeled flux attentuation.

% the question: by constraining flux attenuation with measured respiration rates,
% what values of W_avg allow us to explain a reasonable fraction
% of the observed flux attenuation?

% bottom plots:

% constraining flux attenuation a different way, using measured rates of BP and BR

% this approach requires us to make different assumptions - in this case,
% regarding BGE, ID, and the % of het bact presumed to be particle-attached

% includes sensitivity analysis by station to assess effect of assumed parameter values on BCD

% first, analysis (plots in separate section below)

% first, clean up any old results from previous model runs

clear Flux_est_netresp Resp_matrix est_matrix Flux_est_netresp_alt Resp_matrix_alt est_matrix_alt;

if exist('BCD_pa_int_50_150')
    
    if isfield(BCD_pa_int_50_150,'mmol_C_m2_d')
        BCD_pa_int_50_150=rmfield(BCD_pa_int_50_150,'mmol_C_m2_d');
    end

    if isfield(BCD_pa_int_50_150,'Uncert_mmol_C_m2_d')
        BCD_pa_int_50_150=rmfield(BCD_pa_int_50_150,'Uncert_mmol_C_m2_d');
    end
    
    if isfield(BCD_pa_int_50_150,'CruiseID')
        BCD_pa_int_50_150=rmfield(BCD_pa_int_50_150,'CruiseID');
    end

    if isfield(BCD_pa_int_50_150,'ProcStationID')
        BCD_pa_int_50_150=rmfield(BCD_pa_int_50_150,'ProcStationID');
    end

    if isfield(BCD_pa_int_50_150,'ID')
        BCD_pa_int_50_150=rmfield(BCD_pa_int_50_150,'ID');
    end

    if isfield(BCD_pa_int_50_150,'BGE')
        BCD_pa_int_50_150=rmfield(BCD_pa_int_50_150,'BGE');
    end

end

if isfield(Int_BP_by_station,'Uncert_mmol_C_m2_d')
    Int_BP_by_station=rmfield(Int_BP_by_station,'Uncert_mmol_C_m2_d');
end

if isfield(Int_BP_by_station,'mmol_C_m2_d')
    Int_BP_by_station=rmfield(Int_BP_by_station,'mmol_C_m2_d');
end

if isfield(Int_BP_by_station,'ID')
    Int_BP_by_station=rmfield(Int_BP_by_station,'ID');
end

if isfield(Int_BP_by_station,'BGE')
    Int_BP_by_station=rmfield(Int_BP_by_station,'BGE');
end

if isfield(Nettrap_resprates,'L_remin')
    Nettrap_resprates=rmfield(Nettrap_resprates,'L_remin');
end

if exist('BCD_pa_int_50_150_sensBGE')
    
    if isfield(BCD_pa_int_50_150_sensBGE,'mmol_C_m2_d')
        BCD_pa_int_50_150_sensBGE=rmfield(BCD_pa_int_50_150_sensBGE,'mmol_C_m2_d');
    end

    if isfield(BCD_pa_int_50_150_sensBGE,'Uncert_mmol_C_m2_d')
        BCD_pa_int_50_150_sensBGE=rmfield(BCD_pa_int_50_150_sensBGE,'Uncert_mmol_C_m2_d');
    end

    if isfield(BCD_pa_int_50_150_sensBGE,'CruiseID')
        BCD_pa_int_50_150_sensBGE=rmfield(BCD_pa_int_50_150_sensBGE,'CruiseID');
    end

    if isfield(BCD_pa_int_50_150_sensBGE,'ProcStationID')
        BCD_pa_int_50_150_sensBGE=rmfield(BCD_pa_int_50_150_sensBGE,'ProcStationID');
    end

    if isfield(BCD_pa_int_50_150_sensBGE,'ID')
        BCD_pa_int_50_150_sensBGE=rmfield(BCD_pa_int_50_150_sensBGE,'ID');
    end

    if isfield(BCD_pa_int_50_150_sensBGE,'BGE')
        BCD_pa_int_50_150_sensBGE=rmfield(BCD_pa_int_50_150_sensBGE,'BGE');
    end

    if isfield(BCD_pa_int_50_150_sensBGE,'Frac_atten')
        BCD_pa_int_50_150_sensBGE=rmfield(BCD_pa_int_50_150_sensBGE,'Frac_atten');
    end
    
end

if exist('BCD_pa_int_50_150_f_P')
    
    if isfield(BCD_pa_int_50_150_f_P,'mmol_C_m2_d')
        BCD_pa_int_50_150_f_P=rmfield(BCD_pa_int_50_150_f_P,'mmol_C_m2_d');
    end

    if isfield(BCD_pa_int_50_150_f_P,'Uncert_mmol_C_m2_d')
        BCD_pa_int_50_150_f_P=rmfield(BCD_pa_int_50_150_f_P,'Uncert_mmol_C_m2_d');
    end

    if isfield(BCD_pa_int_50_150_f_P,'CruiseID')
        BCD_pa_int_50_150_f_P=rmfield(BCD_pa_int_50_150_f_P,'CruiseID');
    end

    if isfield(BCD_pa_int_50_150_f_P,'ProcStationID')
        BCD_pa_int_50_150_f_P=rmfield(BCD_pa_int_50_150_f_P,'ProcStationID');
    end

    if isfield(BCD_pa_int_50_150_f_P,'f_P')
        BCD_pa_int_50_150_f_P=rmfield(BCD_pa_int_50_150_f_P,'f_P');
    end

    if isfield(BCD_pa_int_50_150_f_P,'BGE')
        BCD_pa_int_50_150_f_P=rmfield(BCD_pa_int_50_150_f_P,'BGE');
    end
    
end       

if exist('Int_BP_by_station_f_P')

    if isfield(Int_BP_by_station_f_P,'Uncert_mmol_C_m2_d')
        Int_BP_by_station_f_P=rmfield(Int_BP_by_station_f_P,'Uncert_mmol_C_m2_d');
    end

    if isfield(Int_BP_by_station_f_P,'mmol_C_m2_d')
        Int_BP_by_station_f_P=rmfield(Int_BP_by_station_f_P,'mmol_C_m2_d');
    end

    if isfield(Int_BP_by_station_f_P,'f_p')
        Int_BP_by_station_f_P=rmfield(Int_BP_by_station_f_P,'f_P');
    end

    if isfield(Int_BP_by_station_f_P,'BGE')
        Int_BP_by_station_f_P=rmfield(Int_BP_by_station_f_P,'BGE');
    end
    
end

if exist('Int_BP_by_station_sensBGE')

    if isfield(Int_BP_by_station_sensBGE,'Uncert_mmol_C_m2_d')
        Int_BP_by_station_sensBGE=rmfield(Int_BP_by_station_sensBGE,'Uncert_mmol_C_m2_d');
    end

    if isfield(Int_BP_by_station_sensBGE,'mmol_C_m2_d')
        Int_BP_by_station_sensBGE=rmfield(Int_BP_by_station_sensBGE,'mmol_C_m2_d');
    end

    if isfield(Int_BP_by_station_sensBGE,'ID')
        Int_BP_by_station_sensBGE=rmfield(Int_BP_by_station_sensBGE,'ID');
    end

    if isfield(Int_BP_by_station_sensBGE,'BGE')
        Int_BP_by_station_sensBGE=rmfield(Int_BP_by_station_sensBGE,'BGE');
    end
    
end

if exist('BCD_pa_int_50_150_sensID')

    if isfield(BCD_pa_int_50_150_sensID,'mmol_C_m2_d')
        BCD_pa_int_50_150_sensID=rmfield(BCD_pa_int_50_150_sensID,'mmol_C_m2_d');
    end

    if isfield(BCD_pa_int_50_150_sensID,'Uncert_mmol_C_m2_d')
        BCD_pa_int_50_150_sensID=rmfield(BCD_pa_int_50_150_sensID,'Uncert_mmol_C_m2_d');
    end

    if isfield(BCD_pa_int_50_150_sensID,'CruiseID')
        BCD_pa_int_50_150_sensID=rmfield(BCD_pa_int_50_150_sensID,'CruiseID');
    end

    if isfield(BCD_pa_int_50_150_sensID,'ProcStationID')
        BCD_pa_int_50_150_sensID=rmfield(BCD_pa_int_50_150_sensID,'ProcStationID');
    end

    if isfield(BCD_pa_int_50_150_sensID,'ID')
        BCD_pa_int_50_150_sensID=rmfield(BCD_pa_int_50_150_sensID,'ID');
    end

    if isfield(BCD_pa_int_50_150_sensID,'BGE')
        BCD_pa_int_50_150_sensID=rmfield(BCD_pa_int_50_150_sensID,'BGE');
    end

    if isfield(BCD_pa_int_50_150_sensID,'Frac_atten')
        BCD_pa_int_50_150_sensID=rmfield(BCD_pa_int_50_150_sensID,'Frac_atten');
    end
    
end

if exist('Int_BP_by_station_sensID')

    if isfield(Int_BP_by_station_sensID,'Uncert_mmol_C_m2_d')
        Int_BP_by_station_sensID=rmfield(Int_BP_by_station_sensID,'Uncert_mmol_C_m2_d');
    end

    if isfield(Int_BP_by_station_sensID,'mmol_C_m2_d')
        Int_BP_by_station_sensID=rmfield(Int_BP_by_station_sensID,'mmol_C_m2_d');
    end

    if isfield(Int_BP_by_station_sensID,'ID')
        Int_BP_by_station_sensID=rmfield(Int_BP_by_station_sensID,'ID');
    end

    if isfield(Int_BP_by_station_sensID,'BGE')
        Int_BP_by_station_sensID=rmfield(Int_BP_by_station_sensID,'BGE');
    end
    
end

%%
% ----------------------------------------------------------------------------
% flux attenuation constrained by particle/water col
% incubation measurements, with additional consideration of k_DOCZ,
% accounting for loss of OC to water column through both solubilization by attached bacteria and
% zooplankton activities
% ----------------------------------------------------------------------------

% applying an exponential decay model based on a remineralization length
% scale, this time assuming there are other loss processes in addition to
% particle-attached respiration
% 
% generate flux estimates for each assumed sinking rate and k_DOCZ pairing, and uncertainties
% 

c_i=1;
for i=1:length(W_avg)
    for j=1:length(k_DOCZ)

    % estimate at 150 m
    
    % lookup values of parameters for this run & store
    
    F_0=PIT_fluxes.Avg_POC_flux_mmolC_m2_d(PIT_loc(1));
    F_0_sigma=PIT_fluxes.Uncert_POC_flux_mmolC_m2_d(PIT_loc(1));
    k_R_0_a=-Nettrap_resprates.mmolC_mmolCsub_d(Nettrap_resprate_loc(1));
    k_R_0_a_sigma=Nettrap_resprates.Uncert_mmolC_mmolCsub_d(Nettrap_resprate_loc(1));
    k_R_0_b=-Nettrap_resprates.mmolC_mmolCsub_d(Nettrap_resprate_loc(2));
    k_R_0_b_sigma=Nettrap_resprates.Uncert_mmolC_mmolCsub_d(Nettrap_resprate_loc(2));
    k_R_z_a=-Nettrap_resprates.mmolC_mmolCsub_d(Nettrap_resprate_loc(3));
    k_R_z_a_sigma=Nettrap_resprates.Uncert_mmolC_mmolCsub_d(Nettrap_resprate_loc(3));
    k_R_z_b=-Nettrap_resprates.mmolC_mmolCsub_d(Nettrap_resprate_loc(4));
    k_R_z_b_sigma=Nettrap_resprates.Uncert_mmolC_mmolCsub_d(Nettrap_resprate_loc(4));
    z=150;
    z_0=50;
        
    % calculate result
            
    L_remin=W_avg(i)/((mean([mean([k_R_0_a k_R_0_b]) mean([k_R_z_a k_R_z_b])]))+k_DOCZ(j));
    F_z=F_0*exp(-(z-z_0)/L_remin);
      
    k_R_150 = mean([mean([k_R_0_a k_R_0_b]) mean([k_R_z_a k_R_z_b])]);
   
    % store result

    Flux_est_netresp_alt.POC_F_hat_mmolC_m2_d(Flux_est_netresp_loc(2),c_i)=F_z;
    
    % use a Monte Carlo simulation to get uncertainty, assuming all
    % observations were drawn from normal distributions
    
    % generate some vectors (n=1000) of values of our variables, assuming
    % observations drawn from normal distribitions
    
    F_0_montevals = F_0 + F_0_sigma.*randn(1000,1);
    k_R_0_a_montevals = k_R_0_a + k_R_0_a_sigma.*randn(1000,1);
    k_R_0_b_montevals = k_R_0_b + k_R_0_b_sigma.*randn(1000,1);
    k_R_z_a_montevals = k_R_z_a + k_R_z_a_sigma.*randn(1000,1);
    k_R_z_b_montevals = k_R_z_b + k_R_z_b_sigma.*randn(1000,1);
    
    % generate 1000 possible values of F_z, drawing values of variables from
    % our vectors above
    
    for e=1:1000
    
        L_remin_e=W_avg(i)/((mean([mean([k_R_0_a_montevals(e) k_R_0_b_montevals(e)]) mean([k_R_z_a_montevals(e) k_R_z_b_montevals(e)])]))+k_DOCZ(j));
        F_z_monteout(e)=F_0_montevals(e)*exp(-(z-z_0)/L_remin_e);
        
    end     
       
    % calculate & store uncertainty for this run
    
    Flux_est_netresp_alt.POC_F_hat_mmolC_m2_d_uncert(Flux_est_netresp_loc(2),c_i)=std(F_z_monteout);
       
    % now, estimate at 300 m
    
    % lookup values of parameters for this run & store

    F_0=PIT_fluxes.Avg_POC_flux_mmolC_m2_d(PIT_loc(2));
    F_0_sigma=PIT_fluxes.Uncert_POC_flux_mmolC_m2_d(PIT_loc(2));
    k_R_0_a=-Nettrap_resprates.mmolC_mmolCsub_d(Nettrap_resprate_loc(3));
    k_R_0_a_sigma=Nettrap_resprates.Uncert_mmolC_mmolCsub_d(Nettrap_resprate_loc(3));
    k_R_0_b=-Nettrap_resprates.mmolC_mmolCsub_d(Nettrap_resprate_loc(4));
    k_R_0_b_sigma=Nettrap_resprates.Uncert_mmolC_mmolCsub_d(Nettrap_resprate_loc(4));
    k_R_z_a=-Nettrap_resprates.mmolC_mmolCsub_d(Nettrap_resprate_loc(5));
    k_R_z_a_sigma=Nettrap_resprates.Uncert_mmolC_mmolCsub_d(Nettrap_resprate_loc(5));
    k_R_z_b=-Nettrap_resprates.mmolC_mmolCsub_d(Nettrap_resprate_loc(6));
    k_R_z_b_sigma=Nettrap_resprates.Uncert_mmolC_mmolCsub_d(Nettrap_resprate_loc(6));
    z=300;
    z_0=150;    
    
    % calculate result
        
    L_remin=W_avg(i)/((mean([mean([k_R_0_a k_R_0_b]) mean([k_R_z_a k_R_z_b])]))+k_DOCZ(j));
    F_z=F_0*exp(-(z-z_0)/L_remin);
        
    k_R_300 = mean([mean([k_R_0_a k_R_0_b]) mean([k_R_z_a k_R_z_b])]);

    % store result

    Flux_est_netresp_alt.POC_F_hat_mmolC_m2_d(Flux_est_netresp_loc(3),c_i)=F_z;
    
    % use a Monte Carlo simulation to get uncertainty, assuming all
    % observations were drawn from normal distributions
    
    % generate some vectors (n=1000) of values of our variables, assuming
    % observations drawn from normal distribitions
    
    F_0_montevals = F_0 + F_0_sigma.*randn(1000,1);
    k_R_0_a_montevals = k_R_0_a + k_R_0_a_sigma.*randn(1000,1);
    k_R_0_b_montevals = k_R_0_b + k_R_0_b_sigma.*randn(1000,1);
    k_R_z_a_montevals = k_R_z_a + k_R_z_a_sigma.*randn(1000,1);
    k_R_z_b_montevals = k_R_z_b + k_R_z_b_sigma.*randn(1000,1);
    
    % generate 1000 possible values of F_z, drawing values of variables from
    % our vectors above
    
    for e=1:1000
    
        L_remin_e=W_avg(i)/((mean([mean([k_R_0_a_montevals(e) k_R_0_b_montevals(e)]) mean([k_R_z_a_montevals(e) k_R_z_b_montevals(e)])]))+k_DOCZ(j));
        F_z_monteout(e)=F_0_montevals(e)*exp(-(z-z_0)/L_remin_e);
        
    end     
          
    Flux_est_netresp_alt.POC_F_hat_mmolC_m2_d_uncert(Flux_est_netresp_loc(3),c_i)=std(F_z_monteout);
    
    c_i=c_i+1;
    
    % provide some output to user on the progress of the analysis
    
     if  rem(c_i/(length(W_avg)*length(k_DOCZ))*100,.1) == 0 % display progress at 0.1% intervals
        Percent_progress = c_i/(length(W_avg)*length(k_DOCZ))*100
     end
    
    end
    
end

clear i j c_i;

% generate matrix of estimates in form amenable to plotting

for i=1:length(Flux_est_netresp_alt.POC_F_hat_mmolC_m2_d(Flux_est_netresp_loc(2),:))
    est_matrix_alt(i,1)=PIT_fluxes.Avg_POC_flux_mmolC_m2_d(PIT_loc(1));
    est_matrix_alt(i,2)=Flux_est_netresp_alt.POC_F_hat_mmolC_m2_d(Flux_est_netresp_loc(2),i);
    est_matrix_alt(i,3)=Flux_est_netresp_alt.POC_F_hat_mmolC_m2_d(Flux_est_netresp_loc(3),i);
end

x_alt = est_matrix_alt';
Resp_matrix_alt=-(est_matrix_alt-PIT_fluxes.Avg_POC_flux_mmolC_m2_d(PIT_loc(1)));

%%
% ----------------------------------------------------------------------------
% plots for the W_avg-k_DOCZ sensitivity analysis
% ----------------------------------------------------------------------------

% define some needed variables

W_avg_p=W_avg;
k_DOCZ_p=k_DOCZ;

Estflux_150=Flux_est_netresp_alt.POC_F_hat_mmolC_m2_d(Flux_est_netresp_loc(2),:)*12.01;
Estflux_uncert_150=Flux_est_netresp_alt.POC_F_hat_mmolC_m2_d_uncert(Flux_est_netresp_loc(2),:)*12.01;

Estflux_300=Flux_est_netresp_alt.POC_F_hat_mmolC_m2_d(Flux_est_netresp_loc(3),:)*12.01;
Estflux_uncert_300=Flux_est_netresp_alt.POC_F_hat_mmolC_m2_d_uncert(Flux_est_netresp_loc(3),:)*12.01;

dz_150=Estflux_150-PIT_fluxes.Avg_POC_flux_mmolC_m2_d(PIT_loc(2))*12.01;
ind_high=find(dz_150>=1000); % find values that are absurdly high & put in a lower placeholder so visual representation is better
dz_150(ind_high)=1000;
dzmat_150=reshape(dz_150,length(k_DOCZ),length(W_avg));
dz_300=Estflux_300-PIT_fluxes.Avg_POC_flux_mmolC_m2_d(PIT_loc(3))*12.01;
ind_high=find(dz_300>=1000);  % find values that are absurdly high & put in a lower placeholder so visual representation is better
dz_300(ind_high)=1000;
dzmat_300=reshape(dz_300,length(k_DOCZ),length(W_avg));

err_hi_150=PIT_fluxes.Avg_POC_flux_mmolC_m2_d(PIT_loc(2))*12.01+PIT_fluxes.Uncert_POC_flux_mmolC_m2_d(PIT_loc(2))*12.01;
err_lo_150=PIT_fluxes.Avg_POC_flux_mmolC_m2_d(PIT_loc(2))*12.01-PIT_fluxes.Uncert_POC_flux_mmolC_m2_d(PIT_loc(2))*12.01;
err_hi_300=PIT_fluxes.Avg_POC_flux_mmolC_m2_d(PIT_loc(3))*12.01+PIT_fluxes.Uncert_POC_flux_mmolC_m2_d(PIT_loc(3))*12.01;
err_lo_300=PIT_fluxes.Avg_POC_flux_mmolC_m2_d(PIT_loc(3))*12.01-PIT_fluxes.Uncert_POC_flux_mmolC_m2_d(PIT_loc(3))*12.01;

% heat map plots of the W_avg-k_DOCZ analysis, with contour labels

% figure('Position',[1 scrsz(4)*8/10 scrsz(3)*(10/10) scrsz(4)*5/10])

figure('Units','Centimeters','Position',[0 0 32 8])

set(gcf,'PaperPositionMode','auto')

subplot(1,2,1)

% 150 m

hold on;
[c1, h6]=contourf(W_avg,k_DOCZ,dzmat_150,200);
[c10, h10]=contour(W_avg,k_DOCZ,dzmat_150,heatcontours_WavgkDOCZ_150);
clabel(c10,h10,'Color',[147/255 146/255 149/255],'LabelSpacing',100);
set(h10,'LineStyle','-','LineColor',[147/255 146/255 149/255],...
    'LineWidth',1);
set(h6,'LineStyle','none');
load('MyColormaps','cmap_respflux2012');
set(gcf,'Colormap',cmap_respflux2012);
caxis([-max(abs(dz_150))+cmap_axis_offset*12.01 max(abs(dz_150))-cmap_axis_offset*12.01]);
xlabel('Average particle sinking velocity (m d-1)','FontSize',16);
ylabel('k_{DOC+Z}','FontSize',16);
h7=colorbar;
labeltext=strcat(['Deviation from observed POC flux' char(10) '(mg C m^{-2} d^{-1}']);
ylabel(h7,labeltext,'FontSize',16);
set(gca,'FontSize',16,'YScale','log','YTick',[10^-5 10^-4 10^-3 10^-2 10^-1 10^0])
axis(k_DOCZ_axes)

% first, look for statistical significance - i.e., do the modeled and
% observed errors overlap?

clear sens_sig;

for i=1:length(Estflux_150)
    
if Estflux_150(i)>PIT_fluxes.Avg_POC_flux_mmolC_m2_d(PIT_loc(2))*12.01 % modeled value is greater than observed
    
    if Estflux_150(i)-Estflux_uncert_150(i)>err_hi_150
        sens_sig(i)=0; % model result is significant
    else
        sens_sig(i)=1; % null hypothesis holds, model result not significantly different than observed
    end
    
elseif Estflux_150(i)<PIT_fluxes.Avg_POC_flux_mmolC_m2_d(PIT_loc(2))*12.01 % modeled value is less than observed

    if Estflux_150(i)+Estflux_uncert_150(i)<err_lo_150
        sens_sig(i)=0; % model result is significant
    else
        sens_sig(i)=1; % null hypothesis holds, model result not significantly different than observed
    end
    
end
end

% now, overlay significance as shaded area

% reshape significance indicators (1 or 0) into appropriately dimensioned matrix

sens_sig_150 = reshape(sens_sig,length(k_DOCZ),length(W_avg));

[c8, h8]=contourf(W_avg,k_DOCZ,sens_sig_150,[0 1]);
set(h8,'LineStyle','none')

% lay down a second axis to show ratio of k_DOCZ to k_R

posAx_150 = get(gca,'Position');
altax_150 = axes;
set(altax_150,'yaxislocation','left','xtick',[],'Color','none','Position',posAx_150+[-.04 0 0 0],'FontSize',16)
ylabel(altax_150,'k_{DOC+Z}/k_{R}','FontSize',16)
set(altax_150,'YLim',[k_DOCZ_axes(3)/k_R_150 k_DOCZ_axes(4)/k_R_150],'YScale','log',...
    'YTick',[10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);


hold off;

% 300 m

subplot(1,2,2)

hold on;
[c1, h6]=contourf(W_avg,k_DOCZ,dzmat_300,200);
[c11, h11]=contour(W_avg,k_DOCZ,dzmat_300,heatcontours_WavgkDOCZ_300);
clabel(c11,h11,'Color',[147/255 146/255 149/255],'LabelSpacing',200);
set(h11,'LineStyle','-','LineColor',[147/255 146/255 149/255],...
    'LineWidth',1);
set(h6,'LineStyle','none')
load('MyColormaps','cmap_respflux2012')
set(gcf,'Colormap',cmap_respflux2012)
caxis([-max(abs(dz_300))+cmap_axis_offset*12.01 max(abs(dz_300))-cmap_axis_offset*12.01])
xlabel('Average particle sinking velocity (m d-1)','FontSize',16);
ylabel('k_{DOC+Z}','FontSize',16);
h7=colorbar;
labeltext=strcat(['Deviation from observed POC flux' char(10) '(mg C m^{-2} d^{-1}']);
ylabel(h7,labeltext,'FontSize',16);
set(gca,'FontSize',16,'YScale','log','YTick',[10^-5 10^-4 10^-3 10^-2 10^-1 10^0])
axis(k_DOCZ_axes)

% first, look for statistical significance - i.e., do the modeled and
% observed errors overlap?

clear sens_sig;

for i=1:length(Estflux_300)
    
if Estflux_300(i)>PIT_fluxes.Avg_POC_flux_mmolC_m2_d(PIT_loc(3))*12.01 % modeled value is greater than observed
    
    if Estflux_300(i)-Estflux_uncert_300(i)>err_hi_300
        sens_sig(i)=0; % model result is significant
    else
        sens_sig(i)=1; % null hypothesis holds, model result not significantly different than observed
    end
    
elseif Estflux_300(i)<PIT_fluxes.Avg_POC_flux_mmolC_m2_d(PIT_loc(3))*12.01 % modeled value is less than observed

    if Estflux_300(i)+Estflux_uncert_300(i)<err_lo_300
        sens_sig(i)=0; % model result is significant
    else
        sens_sig(i)=1; % null hypothesis holds, model result not significantly different than observed
    end
    
end
end

% now, overlay significance as shaded area

% reshape significance indicators (1 or 0) into appropriately dimensioned matrix

sens_sig_300 = reshape(sens_sig,length(k_DOCZ),length(W_avg));

[c8, h8]=contourf(W_avg,k_DOCZ,sens_sig_300,[0 1]);
set(h8,'LineStyle','none')
% 
% lay down a second axis to show ratio of k_DOC+Z to k_R

posAx_300 = get(gca,'Position');
altax_300 = axes;
set(altax_300,'yaxislocation','left','xtick',[],'Color','none','Position',posAx_300+[-.04 0 0 0],'FontSize',16)
ylabel(altax_300,'k_{DOC+Z}/k_{R}','FontSize',16)
set(altax_300,'YLim',[k_DOCZ_axes(3)/k_R_300 k_DOCZ_axes(4)/k_R_300],'YScale','log',...
    'YTick',[10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);

hold off;

% save the file

fn=strcat(['KN' num2str(PIT_fluxes.CruiseID(PIT_loc(1))) '-' num2str(PIT_fluxes.ProcStationID(PIT_loc(1))) '_W_avg-k_DOC+Z.eps']);
print('-depsc2',fn);


%% Reporting of diagnosed values

%----------------------------------------------------------------
% best-guess diagnosis of kDOCZ, with upper/lower bounds
%----------------------------------------------------------------

% specify W_avg for which we want to report values of k_DOCZ; in this case, the
% median (50 m d-1) from our literature dataset of N Atl sinking speeds

W_avg_lookup_val = 50;

[c ind_W_avg] = min(abs(W_avg-W_avg_lookup_val)); % return index of W_avg (our
% sensitivity analysis vector) that is closest to W_avg_lookup_val

% 150 m

[c ind_kDOCZ] = min(abs(dzmat_150(:,ind_W_avg))); % use this index to find the
% minimum model-obs deviation for this W_avg; the row number corresponding
% to this minimum will be the index to our diagnosed value for k_DOCZ

best_k_DOCZ_150 = k_DOCZ(ind_kDOCZ) % return the value of k_DOCZ corresponding to this index

% now, use similar approach comparing ind_W_avg to our statistical significance matrix
% to get upper and lower bounds

% "upper" and "lower" are counterintuitive here since the sens_sig matrixes
% will superimposed in transpose on the plot beginning at 0,0

% upper bound

ind_kDOCZ_hi = find(sens_sig_150(:,ind_W_avg)==1,1,'last');
kDOCZ_hi_150 = k_DOCZ(ind_kDOCZ_hi)

% lower bound

ind_kDOCZ_lo = find(sens_sig_150(:,ind_W_avg)==1,1);
kDOCZ_lo_150 = k_DOCZ(ind_kDOCZ_lo)

% 300 m

[c ind_kDOCZ] = min(abs(dzmat_300(:,ind_W_avg))); % use this index to find the
% minimum model-obs deviation for this W_avg; the row number corresponding
% to this minimum will be the index to our diagnosed value for k_DOCZ

best_k_DOCZ_300 = k_DOCZ(ind_kDOCZ) % return the value of k_DOCZ corresponding to this index

% now, use similar approach comparing ind_W_avg to our statistical significance matrix
% to get upper and lower bounds

% upper bound

ind_kDOCZ_hi = find(sens_sig_300(:,ind_W_avg)==1,1,'last');
kDOCZ_hi_300 = k_DOCZ(ind_kDOCZ_hi)

% lower bound

ind_kDOCZ_lo = find(sens_sig_300(:,ind_W_avg)==1,1);
kDOCZ_lo_300 = k_DOCZ(ind_kDOCZ_lo)

%----------------------------------------------------------------
% minimum constraint on W_avg
%----------------------------------------------------------------

% 150 m

[ind_W_avg_min_row ind_W_avg_min_col] = find(sens_sig_150(:,:)==1,1);
W_avg_min_150 = W_avg(ind_W_avg_min_col)

% 300 m

[ind_W_avg_min_row ind_W_avg_min_col] = find(sens_sig_300(:,:)==1,1);
W_avg_min_300 = W_avg(ind_W_avg_min_col)

