nop = 'S4_P25';
PROJNAME = ['fluxnet_ch4_' nop];
wd = (['/Users/shuangma/GATTACA_LOCAL_MIRROR/' PROJNAME '/']);
pdfpath = ['/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/' nop '/figures/'];
csvpath = ['/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/' nop '/likelihood/'];
opath = ['/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/' nop '/output/'];
opath_hist = ['/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/' nop '/output_hist/'];
if not(isfolder(pdfpath))
    mkdir(pdfpath)
end
if not(isfolder(csvpath))
    mkdir(csvpath)
end
if not(isfolder(opath))
    mkdir(opath)
end
if not(isfolder(opath_hist))
    mkdir(opath_hist)
end
%%
% mask0 = ncread('/Users/shuangma/CARDAMOM-MAPS-DATASETS/Mask/CARDAMOM-MAPS_GC4x5_LAND_SEA_MASK.nc','data');
% size(mask0)
% % figure(1);plotglobal(mask0);
% GC=GEOSChem_xygrids('4x5'); %$$ for extracting LAT.value
% [mh pixelid]=GC4x5_pull_pixelid_AM();
%% read in independent data GRACE EWT at FLUXNET locations
load '/Users/shuangma/RESEARCH/DATA/CARDAMOM/GRACE_05deg_FLUXNET/GRACE_fluxnetch4_grids.mat'
%% load coordinates info for all the fluxnet ch4 sites
coord=readtable('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S2_LATLON2CBF_edited.csv');
rowv = string(coord.r);
colv = string(coord.c);
coordv=string(coord.cbfname);
siteid=string(coord.SITE_ID);
siteindex=string(coord.nsite);

% nfrow=6; % number of figure rows
% nfcol=4; % number of figure cols
string3 = strings(length(coord.cbfname),1); %nominally define a size for string3
istr=1;
wbmatrix = nan(length(coord.cbfname),3);
%%
for i=1:length(coord.cbfname)
% for i=1:40
    i
    irow = str2num(rowv(i));
    icol = str2num(colv(i));
    isiteindex = str2num(siteindex(i));
    isiteid = (siteid(i));    
%   for iexp=[2 3 4 5 6] %$$
% for iexp=[2 3] %$$
   for iexp=[1 2 6 7]
      cbfdir = [wd 'CBF_FILES/' num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '.cbf.nc'];
      cbrdir = [wd 'CBR_FILES/' num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '.cbf.cbr'];
%       cbrdir = {[wd 'cbr_files/TR.1032.311_' char(pid) '_1.cbr'],[wd 'cbr_files/TR.1032.311_' char(pid) '_2.cbr'],[wd 'cbr_files/TR.1032.311_' char(pid) '_3.cbr'],[wd 'cbr_files/TR.1032.311_' char(pid) '_4.cbr']};
%       try
        CBR=CARDAMOM_RUN_MODEL(cbfdir,cbrdir);
% %         mCBR=CARDAMOM_RUN_MODEL(cbfdir,CBR.PARS(1,:)); % added this line
% %         to test cost function set ups, corrected minmax_value prerequisites;
% %         CBRtest=CARDAMOM_RUN_MDF(CBF); % used this line to test of the
% %         cost function corrections worked, short chain test iexp5
% %         may have a problem here, site 1 iexp=3 has 100 posteriors only
% %         CBF.PEQ_PAW_z.max_value=10; % adjust constraint to debug
%       catch ME
%           fprintf([char(isiteid) num2str(iexp) 'get CBR wo success']);
%           string3(istr) = sprintf([char(isiteid) num2str(iexp) 'get CBR wo success']);istr=istr+1;
%           continue;
%       end
      disp('read CBR')
      xvalue = 1:size(CBR.FLUXES,2);
      CBF=CARDAMOM_READ_NC_CBF_FILE(cbfdir);
      obs=CBF;
      MD=CARDAMOM_MODEL_LIBRARY(CBF.ID.values,[],1);
      disp('read model library')

      [mod,obs] = get_cardamom4plot(CBF,CBR,MD);
      
      parnames  = string(fieldnames(MD.PARAMETER_IDs))';
%       varnames1 = string(fieldnames(mod.mean))';
%       varnames2 = string(fieldnames(mod.residence_time))';
%       varnames2 = strcat('RT_', string(fieldnames(mod.residence_time))');
%% save posterior of mean fluxes, posterior of parameters, to csv files
      tsavetable1 = [struct2table(mod.mean),struct2table(mod.residence_time)];
      tsavetable2 = array2table(mod.ppd);
      tsavetable2.Properties.VariableNames = parnames;

      writetable(tsavetable1,[opath_hist num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '_fluxes.csv'])
      writetable(tsavetable2,[opath_hist num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '_pars.csv'])
   end
end
