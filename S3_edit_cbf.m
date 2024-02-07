% if add/delete obs field, make sure to change here 
% and change loademptyCBF
% and change edit_fluxnetCBFuncs
use_Xu_ABGB = 0; % if 0 use Yan_1km_FLUXNET ABGB
use_05deg_soildepth = 0; % if 0 use load_1km_soildepth_Pelletier
% load worldmesh
[mlat,mlon,marea] = loadworldmesh(0.5);

% load ERA5 met data
[met] = Load_05deg_ERA5_met(); % takes less than 10s

metnames = ["time","T2M_MIN","T2M_MAX","SSRD","STRD","CO2",...
    "DOY","BURNED_AREA","VPD","TOTAL_PREC","SNOWFALL","SKT"];

% load obs SOC from HWSD on CARDAMOM google drive
% -- -- -- load obs CARDAMOM 05deg data, Xu ABGB,EWT,SOM,soildepth -- -- --
[obs] = Load_05deg_obs();
%%
% load obs SCF and LAI from MODIS, prepared by Yan
% read from Yan's high resolution data
obsdir2 = '/Users/shuangma/RESEARCH/DATA/CARDAMOM/MODIS_LAI_SCF_FLUXNET_Yan/';
obs.LAI = readtable([obsdir2 'shuang_LAI.csv']);
obs.LAI.Properties.RowNames = obs.LAI.SITE_ID;
obs.LAI('US-MAC',:) = [];
obs.SCF = readtable([obsdir2 'shuang_snowcover.csv']);
obs.SCF.Properties.RowNames = obs.SCF.SITE_ID;
obs.SCF('US-MAC',:) = []; 
% read from GRACE 05deg data, organized by Shuang
obsdir2 = '/Users/shuangma/RESEARCH/DATA/CARDAMOM/GRACE_05deg_FLUXNET/';
obs.EWT = readtable([obsdir2 'GRACE_fluxnetch4_grids.csv']);
obs.EWT.Properties.RowNames = obs.EWT.SITE_ID;

% load coordinates info for all the fluxnet ch4 sites
coord=readtable('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S2_LATLON2CBF_edited.csv');
% edited csv file move around coords for site 12,31,52,64 to be land pixels
rowv = string(coord.r);
colv = string(coord.c);
coordv=string(coord.cbfname);
siteid=string(coord.SITE_ID);
siteindex=string(coord.nsite);

string1 = strings(length(coord.cbfname),1);
string2 = strings(length(coord.cbfname),1);

for i=1:length(coord.cbfname)
    irow = str2num(rowv(i));
    icol = str2num(colv(i));
    isiteindex = str2num(siteindex(i));
    isiteid = (siteid(i));
    yendremove=0;
% load site obs file, these files were processed to include data before end of 2021 
    b=sprintf('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S1/%s_%s_obs.csv',siteindex(i),siteid(i));
    obsdata=readtable(b);
    if size(obsdata,1)==0
        string1(i) = sprintf('%s%s has no obs data before 2021',siteindex(i),siteid(i));
    elseif sum(isnan(squeeze(met.(metnames(2))(irow,icol,:))))~=0
        string2(i) = sprintf('%s%s is too close to the sea',siteindex(i),siteid(i));
    else
        
%% step load empty CBF with basic info ready, exact nom for each field, lat
    % the timeframe of loaded CBF depend on obsdata
    % cut time length according to year from fluxnet obs
      % need to chage MET, OBS, nodays
        ystart = min(obsdata.year)-5; % run 5 years before fusing data
        yend = max(obsdata.year)+3;   % run 3 more year after fusing data
        if ystart<2001   % earliest FLUXNET is in 2006, never smaller than 2001+3
            ystart=2001;
        end
        if yend > 2021
            if yend == (2019+3)
                yendremove=1;
            elseif yend == (2020+3)
                yendremove=2;
            elseif yend == (2021+3)
                yendremove=3;
            end
            yend=2021;
        end        
% read startyear
    obsstartyear = min(obsdata.year);
% read endyear 
    obsendyear = max(obsdata.year);
% read lat
    latitude = mlon(str2double(rowv(i)),1);
% define the MODELID, 1100 is the default, change it in S4 script not here ~!
    MODELID=1100;
% -- -- -- -- -- -- load empty CBF -- -- -- -- -- --
    [CBF_exp] = loademptyCBF(ystart,yend,latitude,MODELID);
            
%% step edit met for CBF
% the loaded met span 2001-2021
    selectmetrows = ((ystart-2001)*12+1):((yend-2000)*12);
    for imet=[2:6, 9:12]
      CBF_exp.(metnames(imet)).values = squeeze(met.(metnames(imet))(irow,icol,selectmetrows));
      CBF_exp.(metnames(imet)).reference_mean = nanmean(CBF_exp.(metnames(imet)).values);
    end
      CBF_exp.BURNED_AREA.values = zeros([length(selectmetrows) 1]);
      CBF_exp.BURNED_AREA.reference_mean = nanmean(CBF_exp.BURNED_AREA.values);
%% step edit obs values for CBF
    yearmark = repelem((ystart:yend),12);
% at some sites observation start/end in the middle of years
    % check if obs start from Jan
    obstartmonth = obsdata.month(1);
    obsendmonth = obsdata.month(end);
    spaceholder0 = repelem(NaN,(obstartmonth-1)); % if obs month start in the middle of year
    spaceholder00= repelem(NaN,(12-obsendmonth)); % if obs month end in the middle of year
    % add 5yrs of NA record to FLUXNET data, as space holder
    spaceholder1 = repelem(NaN,12*5);
    % add 3yrs to the end if within 2021
    if yendremove==0
        spaceholder2 = repelem(NaN,12*3);
    elseif yendremove==1
        spaceholder2 = repelem(NaN,12*2);
    elseif yendremove==2
        spaceholder2 = repelem(NaN,12*1);
    else
        spaceholder2 = [];
    end
            
%   save no obs CBF to fill in later
    CBF_exp_noobs = CBF_exp; 
% fill in all available values for CBF
    CBF_exp.NBE.values = [spaceholder0';spaceholder1';table2array(obsdata(:,8));spaceholder2';spaceholder00'];
    CBF_exp.NBE.values(CBF_exp.NBE.values==-9999)=nan;
    CBF_exp.CH4.values = [spaceholder0';spaceholder1';table2array(obsdata(:,4));spaceholder2';spaceholder00'];
    CBF_exp.CH4.values(CBF_exp.CH4.values==-9999)=nan;
    CBF_exp.ET.values = [spaceholder0';spaceholder1';table2array(obsdata(:,15));spaceholder2';spaceholder00'];
    CBF_exp.ET.values(CBF_exp.ET.values==-9999)=nan;
    CBF_exp.ET.values(CBF_exp.ET.values<0)=nan;
    CBF_exp.GPP.values = [spaceholder0';spaceholder1';table2array(obsdata(:,11));spaceholder2';spaceholder00'];
    CBF_exp.GPP.values(CBF_exp.GPP.values==-9999)=nan;
    CBF_exp.GPP.values(CBF_exp.GPP.values<0)=nan;   
    CBF_exp.ER.values = [spaceholder0';spaceholder1';table2array(obsdata(:,13));spaceholder2';spaceholder00'];
    CBF_exp.ER.values(CBF_exp.ER.values==-9999)=nan;
    CBF_exp.ER.values(CBF_exp.ER.values<0)=nan;       
% the MODIS SCF and LAI span 2001-2021, select same rows as did for
% met,convert from 100% to fraction
    whole_SCF = (table2array(obs.SCF(isiteid,7:end)))'*0.01;
    CBF_exp.SCF.values = whole_SCF(selectmetrows);
% put in time series data and use filter=1 to match mean of whenever there's obs lai data
    whole_LAI = (table2array(obs.LAI(isiteid,7:end)))';
    CBF_exp.LAI.values = whole_LAI(selectmetrows);
% GRACE EWT 2002-2021
    whole_EWT = (table2array(obs.EWT(isiteid,11:end)))'; 
    whole_EWT = [NaN(12,1); whole_EWT];
    CBF_exp.EWT.values = whole_EWT(selectmetrows);    
% add a test: Mean_LAI
   % CBF_exp.Mean_LAI.values = nanmean(whole_LAI(selectmetrows));
% HWSD SOM    
    CBF_exp.PEQ_iniSOM.values = obs.SOM(irow,icol);
% Soil depth, use soil depth database to add maximum limit to PAW_z 
    if use_05deg_soildepth==1
        CBF_exp.PEQ_PAW_z.max_value = obs.soildepth(irow,icol);
    else
        [soil_depth_value] = load_1km_soildepth_Pelletier(isiteid);
        % the script treated -1m locations as 1m
        CBF_exp.PEQ_PAW_z.max_value = soil_depth_value;
    end
% ABGB    
    if use_Xu_ABGB ==1
    % Xu et al 2021 ABGB, available 2001-2019
        if yend < 2020
            smoyABGB = ((ystart-2001)*12+1):((yend-2001)*12+12);
        else
            smoyABGB = ((ystart-2001)*12+1):((2019-2001)*12+12);
        end
        CBF_exp.ABGB.values(1:(min(yend,2019)-ystart+1)*12) = squeeze(obs.ABGB(irow,icol,smoyABGB));
    else
    % Yan FLUXNET 1km available 2020, 2015
        [Yan_2015,Yan_2020] = load_1km_ABGB_Yan(isiteid);
        if ystart<=2015 && yend>=2015
            CBF_exp.ABGB.values(((2015-ystart)*12+1):((2015-ystart)*12+12)) = Yan_2015;
        end
        if ystart<=2020 && yend>=2020
            CBF_exp.ABGB.values(((2020-ystart)*12+1):((2020-ystart)*12+12)) = Yan_2020;
        end
    end
% -- -- -- -- -- fill in MCMC setups -- -- -- -- --
    [CBF_exp] = edit_fluxnetCBFuncs(CBF_exp);

            %% Step 3. save other obs not used in MDF
            clear otherobs;
            otherobs.GPP_DT = [spaceholder0';spaceholder1';table2array(obsdata(:,11));spaceholder2';spaceholder00'];
            otherobs.ER_DT = [spaceholder0';spaceholder1';table2array(obsdata(:,13));spaceholder2';spaceholder00'];
            otherobs.GPP_NT = [spaceholder0';spaceholder1';table2array(obsdata(:,12));spaceholder2';spaceholder00'];
            otherobs.ER_NT = [spaceholder0';spaceholder1';table2array(obsdata(:,14));spaceholder2';spaceholder00'];
            otherobs.NBE_F = [spaceholder0';spaceholder1';table2array(obsdata(:,9));spaceholder2';spaceholder00'];
            otherobs.CH4_F = [spaceholder0';spaceholder1';table2array(obsdata(:,5));spaceholder2';spaceholder00'];
            otherobs.ET_F = [spaceholder0';spaceholder1';table2array(obsdata(:,16));spaceholder2';spaceholder00'];
            otherobs.NBE_ANN = [spaceholder0';spaceholder1';table2array(obsdata(:,10));spaceholder2';spaceholder00'];
            otherobs.CH4_ANN = [spaceholder0';spaceholder1';table2array(obsdata(:,6));spaceholder2';spaceholder00'];
            otherobs.ET_ANN = [spaceholder0';spaceholder1';table2array(obsdata(:,17));spaceholder2';spaceholder00'];

            otherobs_cell = struct2cell(otherobs);
            clear N;
            N = [CBF_exp.NBE.values,CBF_exp.CH4.values,CBF_exp.ET.values, ...
                CBF_exp.LAI.values,CBF_exp.SCF.values,CBF_exp.EWT.values,...
                CBF_exp.ABGB.values, ...
                otherobs_cell{1:length(otherobs_cell)}];
            TN = array2table(N);
            TN.Properties.VariableNames = {'NBE', 'CH4','ET','LAI','SCF','EWT','ABGB', ...
                'GPP_DT','ER_DT','GPP_NT','ER_NT','NBE_F','CH4_F','ET_F','NBE_ANN','CH4_ANN','ET_ANN'};
            
            pathname = sprintf('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/obs/');
            if not(isfolder(pathname))
                mkdir(pathname)
            end
            filename = [pathname num2str(isiteindex) '_' char(isiteid) '_cbf_obs.csv'];
%             filename = [pathname char(isiteid) '_cbf_obs.csv'];
            writetable(TN,filename);  % 
            
            %% Step 4. save met to csv files
%             metnames = ["time","T2M_MIN","T2M_MAX","SSRD","CO2",...
%     "DOY","BURNED_AREA","VPD","TOTAL_PREC","SNOWFALL"];
            for isave = 1:length(metnames)
                savemet{isave} = CBF_exp.(metnames(isave)).values;
            end
            savetable = [yearmark' savemet{1:length(metnames)}];
            tsavetable = array2table(savetable);
            tsavetable.Properties.VariableNames = ['yearmark', metnames(1:length(metnames))];
            ometdir = 'RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/met/';
            if not(isfolder(ometdir))
                mkdir(ometdir)
            end
            writetable(tsavetable,[ometdir num2str(isiteindex) '_' char(isiteid) '_met.csv'])             
            %% Step 5. save to .cbf.nc
            pathname2=sprintf('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/input/');
            if not(isfolder(pathname2))
                mkdir(pathname2)
            end
            CARDAMOM_WRITE_NC_CBF_FILE(CBF_exp,[pathname2 num2str(isiteindex) '_' char(isiteid) '.cbf.nc']);
            CARDAMOM_WRITE_NC_CBF_FILE(CBF_exp_noobs,[pathname2 num2str(isiteindex) '_' char(isiteid) '_noobs.cbf.nc']);
    end
end

