% main_202306 is the staged branch

% ./gattaca2_local_mirror/fluxnet_ch4_S4_P45/RSYNC_CARDAMOM_PROJECT.sh

nop = 'S4_P46';
PROJNAME = ['fluxnet_ch4_' nop];
wd = (['/Users/shuangma/GATTACA2_LOCAL_MIRROR/' PROJNAME '/']);
pdfpath = ['/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/' nop '/figures/'];
csvpath = ['/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/' nop '/likelihood/'];
opath = ['/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/' nop '/output/'];
if not(isfolder(pdfpath))
    mkdir(pdfpath)
end
if not(isfolder(csvpath))
    mkdir(csvpath)
end
if not(isfolder(opath))
    mkdir(opath)
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
% coord=readtable('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S2_LATLON2CBF_HL.csv');
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
LIKELIHOOD_colnames = {'ABGB','CH4','DOM','ET','EWT','GPP','SIF','LAI',...
    'NBE','ROFF','SCF','Mean_ABGB','Mean_FIR','Mean_GPP','Mean_LAI',...
    'PEQ_Cefficiency','PEQ_CUE','PEQ_iniSnow','PEQ_iniSOM','PEQ_C3frac','PEQ_Vcmax25','PEQ_LCMA',...
    'PEQ_r_ch4','PEQ_S_fv','PEQ_rhch4_rhco2'};
%%
% find(coord.cbfname==317049)
% for i=1:length(coord.cbfname)
% for i=[1 12 61]
  for i=58
    irow = str2num(rowv(i));
    icol = str2num(colv(i));
    isiteindex = str2num(siteindex(i));
    isiteid = (siteid(i));    
  for iexp=[8] %$$
% for iexp=[2 3] %$$
%    for iexp= [6]
      cbfdir = [wd 'CBF_FILES/' num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '.cbf.nc'];
      cbfdir_allobs = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P45/CBF_FILES/' num2str(isiteindex) '_' char(isiteid) '_exp6.cbf.nc'];
      cbrdir = [wd 'CBR_FILES/' num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '.cbf.cbr'];
%       cbrdir = {[wd 'cbr_files/TR.1032.311_' char(pid) '_1.cbr'],[wd 'cbr_files/TR.1032.311_' char(pid) '_2.cbr'],[wd 'cbr_files/TR.1032.311_' char(pid) '_3.cbr'],[wd 'cbr_files/TR.1032.311_' char(pid) '_4.cbr']};
      try
        CBR=CARDAMOM_RUN_MODEL(cbfdir,cbrdir);
%         CBR.PARS(:,MD.PARAMETER_IDs.scf_scalar) = 15;
%         mCBR=CARDAMOM_RUN_MODEL(cbfdir,CBR.PARS); % added this line
%         
%         SWE=CBR.POOLS(:,:,MD.POOL_IDs.H2O_SWE);
%         mSWE=mCBR.POOLS(:,:,MD.POOL_IDs.H2O_SWE);
%         SCF=CBR.POOLS(:,:,MD.POOL_IDs.D_SCF);
%         mSCF=mCBR.POOLS(:,:,MD.POOL_IDs.D_SCF);
%         
%         LIKELIHOOD_table = array2table(CBR.LIKELIHOODS);
%         LIKELIHOOD_table.Properties.VariableNames = LIKELIHOOD_colnames;
%         figure(1); plotunc(SWE);hold on;plotmultilines(mSWE);
%         figure(2); plotunc(SCF);hold on;plotmultilines(mSCF);
        
%         MD.PARAMETER_IDs.scf_scalar
%         MD.POOL_IDs.H2O_SWE
%         to test cost function set ups, corrected minmax_value prerequisites;
%         CBRtest=CARDAMOM_RUN_MDF(CBF); % used this line to test of the
%         cost function corrections worked, short chain test iexp5
%         may have a problem here, site 1 iexp=3 has 100 posteriors only
%         CBF.PEQ_PAW_z.max_value=10; % adjust constraint to debug
      catch ME
          fprintf([char(isiteid) num2str(iexp) 'get CBR wo success']);
          string3(istr) = sprintf([char(isiteid) num2str(iexp) 'get CBR wo success']);istr=istr+1;
          continue;
      end
      disp('read CBR')
      xvalue = 1:size(CBR.FLUXES,2);
      CBF=CARDAMOM_READ_NC_CBF_FILE(cbfdir);
      CBF_allobs=CARDAMOM_READ_NC_CBF_FILE(cbfdir_allobs);

      MD=CARDAMOM_MODEL_LIBRARY(CBF.ID.values,[],1);
      disp('read model library')

      [mod,obs] = get_cardamom4plot_main(CBF,CBR,MD);
      [mod_nomimal,obs_allobs] = get_cardamom4plot_main(CBF_allobs,CBR,MD);      
%% read in independent data for the site
      rawdata = fluxnetcoords(fluxnetcoords.SITE_ID==isiteid,:);
      fluxnet_startyr = sscanf(obs.time.info,'Time since January of %d');
      GRACE_startmon = ((fluxnet_startyr-2001)*12+1); 
      GRACE_endmon   = GRACE_startmon + size(obs.time.values,1) - 1;
      fluxnet_endyr = fluxnet_startyr + size(obs.time.values,1)/12 -1;
      ind.monthstamp = repmat(1:12,1,size(obs.time.values,1)/12)';
      ind.yearstamp = repelem(fluxnet_startyr:fluxnet_endyr,12)';
      ind.GRACE_EWT = table2array(rawdata(:,11))';   % 2002-2021
      ind.GRACE_EWT = [repelem(nan,12)'; ind.GRACE_EWT]; % extend to 2001-2021
      ind.GRACE_EWT = ind.GRACE_EWT(GRACE_startmon:GRACE_endmon);
%% save independent validation data to csv files
      indnames = string(fieldnames(ind))';
      tsavetable = struct2table(ind);
      tsavetable.Properties.VariableNames = indnames;
      oinddir = 'RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/ind/';
      if not(isfolder(oinddir))
          mkdir(oinddir)
      end
      writetable(tsavetable,[oinddir num2str(isiteindex) '_' char(isiteid) '_ind.csv'])
%%            
      % check Water balance:
      % input: mod, CBF, CBR, MD
      % output: three numbers difference of LHS and RHS for PAW, PUW, SWE
%       [WB] = check_water_balance(mod, CBR, MD);
%       wbmatrix(i,:) = [WB.PAW WB.PUW WB.SWE];
      
      disp('assigned mod')
      outdata = {mod.NBE,mod.NBEano,mod.rh_ch4,mod.ABGB,mod.ET,mod.SCF,mod.swe,...
          mod.rhch4_rhco2,mod.rhco2_rh,mod.mean.rhch4_rhco2,mod.mean.rhco2_rh,...
          mod.qsurf,mod.melt,mod.infil,mod.prec,mod.snowfall,...
          mod.EWT,mod.H2O_LY1,mod.H2O_LY2,mod.SM_LY1,mod.SM_LY2,...
          mod.GPP,mod.rh_co2,mod.Rh,mod.NPP,mod.labC,mod.leafC,mod.woodC,mod.rootC,mod.SOM}; % ,mod.SOM,mod.LAI,mod.leafC,mod.rh_co2,mod.rh_ch4,mod.soil_moist,mod.NPP,mod.Ra
      valdata = {obs.NBE.values,obs.NBEano,obs.CH4.values,obs.ABGB.values,obs.ET.values,obs.SCF.values}; % ,obs.SOM.values(1),obs.Mean_LAI.values
      valdata_allobs = {obs_allobs.NBE.values,obs_allobs.NBEano,obs_allobs.CH4.values,obs_allobs.ABGB.values,obs_allobs.ET.values,obs_allobs.SCF.values}; % ,obs.SOM.values(1),obs.Mean_LAI.values
      varnames= {'NBE','NBEano','CH4','ABGB','ET','SCF','SWE'...
          'rhch4_rhco2','rhco2_rh','mean.rhch4_rhco2','mean.rhco2_rh',...
          'qsurf','melt','infil','prec','snowfall',...
          'EWT','H2O_LY1','H2O_LY2','SM_LY1','SM_LY2',...
          'GPP','rh_co2','Rh','NPP','labC','leafC','woodC','rootC','SOM'};

% save data to output folder      
      for idata=1:length(varnames)
          iymat = outdata{idata};
          v95 = quantile(iymat,0.95,1);
          v75 = quantile(iymat,0.75,1);
          v50 = quantile(iymat,0.50,1);
          v25 = quantile(iymat,0.25,1);
          v05 = quantile(iymat,0.05,1);
          vmean = nanmean(iymat,1);
          savedata = [vmean; v05; v25; v50; v75; v95];
          savetable = array2table(savedata);
          savetable.Properties.RowNames(1:size(savetable,1)) = ...
    {'mean', 'v05','v25','v50','v75','v95'};
          ofilename = [opath 'site' num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '_' char(varnames(idata)) '.csv'];
          writetable(savetable,ofilename,'WriteRowNames', true,'WriteVariableNames', false)
      end
% save parameter distributions to outputfolder
          savetable = array2table(CBR.PARS);
          savetable.Properties.VariableNames(1:size(savetable,2)) = string(fieldnames(MD.PARAMETER_IDs));
          ofilename = [opath 'site' num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '_pdf.csv'];
          writetable(savetable,ofilename,'WriteRowNames', false,'WriteVariableNames', true)
% save plots to figures folder
      f1 = figure(1);%ifig=ifig+1;
      isub=1;
      disp('right before plotting figure 1')
      for o=1:6%(length(varnames)-15)
          nfrow=7; % number of figure rows
          nfcol=4; % number of figure cols
          if o==4 || o==6 % i think this is to compare state at beginning of month?
              yhat = nanmedian(outdata{o}(:,1:(end-1)),1)';
          else
              yhat = nanmedian(outdata{o},1)';
          end
          y = valdata{o};
          y_allobs = valdata_allobs{o};
%           y = y(1:xmax);
          % get stats metric
          RMSE = sqrt(nanmean((y - yhat).^2));  % Root Mean Squared Error
          mdl = fitlm(y,yhat);
          R2=mdl.Rsquared.Ordinary;
%           xymin = min(min(y,[],1),min(yhat,[],1));
%           xymax = max(max(y,[],1),max(yhat,[],1));
%             ylim([xymin xymax]);xlim([xymin xymax]);
          subplot(nfrow,nfcol,isub);isub=isub+1;                  
          % green will be covered by red lines if the obs is actually used
          % in the mcmc
          plotunc(outdata{o});plot(xvalue,y_allobs,'g-','LineWidth',2);plot(xvalue,y,'r-','LineWidth',2);title([varnames{o}]);   

          subplot(nfrow,nfcol,isub);isub=isub+1; 
          scatter(yhat,y);title(['exp' num2str(iexp) ' ' varnames{o}]);
          ylabel('obs');xlabel('mod');
          
          NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
          SW = [min(xlim) min(ylim)]+[diff(xlim) diff(ylim)]*0.05;
          NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
          SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
          txt1 = ['rmse = ', num2str(RMSE)];
          text(NW(1), NW(2),txt1);
          txt2 = ['r2 = ', num2str(R2)];
          text(SW(1), SW(2),txt2);
      end
      
% continue plotting model outputs that doesn't need to calculate RMSE R2
      disp('before SOM plot')
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.SOM);title(['SOM']);
%       plotunc(mod.SOM);title(['exp' num2str(iexp) ' SOM']);
      if isfield(obs, 'PEQ_iniSOM')
          plot(1,obs.PEQ_iniSOM.values,'r.','markersize', 18);
      end
          subplot(nfrow,nfcol,isub);isub=isub+1;    
      if isfield(obs, 'PEQ_iniSOM')
          a = median(mod.SOM);
          scatter(a(1),obs.PEQ_iniSOM.values);title(['SOM']);
      end
      disp('after SOM plot')
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.LAI);
      nonnanIndexes = ~isnan(obs.LAI.values);  
      theMean = mean(mod.LAI(nonnanIndexes));
      theMean_wholetime = mean(mod.LAI,[1 2]);
      if sum(isnan(obs.LAI.values))<length(obs.LAI.values) % if there's obs LAI
        yline(nanmean(obs.LAI.values),'r-','LineWidth',3);
        yline(theMean,'g-','LineWidth',3);
      elseif ~isnan(obs.Mean_LAI.values)
        yline(obs.Mean_LAI.values,'r-','LineWidth',3);
        yline(theMean_wholetime,'g-','LineWidth',3);
      end
      title(['mod.LAI obs.MeanLAI']);
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      scatter(mean(median(mod.LAI)),obs.Mean_LAI.values');
      title(['mod.LAI obs.MeanLAI']);
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.leafC);
      title(['leafC']);      
      
%       subplot(nfrow,nfcol,isub);isub=isub+1;      
%       plotunc(mod.rootC);
%       title(['exp' num2str(iexp) ' rootC']);      
%       
%       subplot(nfrow,nfcol,isub);isub=isub+1;      
%       plotunc(mod.woodC);
%       title(['exp' num2str(iexp) ' woodC']);  
            
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.rh_co2);
      title(['rh co2']);  
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.rh_co2+mod.rh_ch4);
      title(['Rh']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;           
      plotunc(mod.Ra);
      title(['Ra']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;           
      plotunc(mod.GPP);
      title(['GPP']); 
      
      if iexp==2 || iexp==6 % when use in situ data
          subplot(nfrow,nfcol,isub);isub=isub+1;                  
          plotunc(outdata{1});plot(xvalue,valdata{1},'r-','LineWidth',2);title([varnames{1}]);
          ylim([min(valdata{1}) max(valdata{1})]);
          subplot(nfrow,nfcol,isub);isub=isub+1;                  
          plotunc(outdata{2});plot(xvalue,valdata{2},'r-','LineWidth',2);title([varnames{2}]);  
          ylim([min(valdata{2}) max(valdata{2})]);
          subplot(nfrow,nfcol,isub);isub=isub+1;                  
          plotunc(outdata{3});plot(xvalue,valdata{3},'r-','LineWidth',2);title([varnames{3}]);  
          ylim([min(valdata{3}) max(valdata{3})]);
      end
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.SCF);
      title(['SCF']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      hist(CBR.PARS(:,MD.PARAMETER_IDs.r_ch4));
%       xlim([MD.parmin(MD.PARAMETER_IDs.r_ch4) MD.parmax(MD.PARAMETER_IDs.r_ch4)]);
      xlim([0.001 0.9])
      title(['r_ch4 (0.001-0.9)']); 
      ylimits = ylim; ymax = ylimits(2);
      format shortg
      text(10,ymax*3/4,['median=',num2str(round(median(CBR.PARS(:,MD.PARAMETER_IDs.r_ch4)),2))]);
      
%       subplot(nfrow,nfcol,isub);isub=isub+1;      
%       plotunc(mod.fV);
%       title(['volumetric frac of ae soil']); 
      
%       subplot(nfrow,nfcol,isub);isub=isub+1;      
%       plotunc(mod.fW);
%       title(['sm factor in soil_resp']);
%       
%       subplot(nfrow,nfcol,isub);isub=isub+1;      
%       plotunc(mod.fT);
%       title(['temp factor in co2 prod']);      
%       
%       subplot(nfrow,nfcol,isub);isub=isub+1;      
%       plotunc(mod.fTch4);
%       title(['temp factor in ch4 prod']);
%       
%       subplot(nfrow,nfcol,isub);isub=isub+1;      
%       plotunc(mod.fCH4);
%       title(['frac of anaeco2 convt to ch4']);
      
      disp('after soil moist')
      sgtitle(['site' num2str(isiteindex) '_' char(isiteid) ' exp' num2str(iexp)]);
      %% save to pdf
      %% define page and save a pdf file
clear figure_property;
figure_property.units = 'inches';
figure_property.format = 'pdf';
figure_property.Preview= 'none';
figure_property.Width= '18'; % Figure width on canvas
figure_property.Height= '18'; % Figure height on canvas
figure_property.Units= 'inches';
figure_property.Color= 'rgb';
figure_property.Background= 'w';
figure_property.FixedfontSize= '6';
figure_property.ScaledfontSize= 'auto';
figure_property.FontMode= 'scaled';
figure_property.FontSizeMin= '6';
figure_property.FixedLineWidth= '2';
figure_property.ScaledLineWidth= 'auto';
figure_property.LineMode= 'none';
figure_property.LineWidthMin= '0.5';
figure_property.FontName= 'Times New Roman';% Might want to change this to something that is available
figure_property.FontWeight= 'auto';
figure_property.FontAngle= 'auto';
figure_property.FontEncoding= 'latin1';
figure_property.PSLevel= '3';
figure_property.Renderer= 'painters';
figure_property.Resolution= '100';
figure_property.LineStyleMap= 'none';
figure_property.ApplyStyle= '0';
figure_property.Bounds= 'tight';
figure_property.LockAxes= 'off';
figure_property.LockAxesTicks= 'off';
figure_property.ShowUI= 'off';
figure_property.SeparateText= 'off';
chosen_figure=gcf;
set(chosen_figure,'PaperUnits','inches');
set(chosen_figure,'PaperPositionMode','auto');
set(chosen_figure,'PaperSize',[str2num(figure_property.Width) str2num(figure_property.Height)]); % Canvas Size
set(chosen_figure,'Units','inches');
hgexport(gcf,[pdfpath 'validation_1100_5e5_site' num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '.pdf'],figure_property); %Set desired file name
clf(f1);
clf(f1,'reset')


%%
      % save water plots to figures folder
      f2 = figure(2);%ifig=ifig+1;
      nfrow = 5;
      nfcol = 6;
      isub=1;
      disp('right before plotting figure 2')

      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.swe);
      title(['SWE (mm)']);
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.H2O_LY1);
      title(['H2O_LY1 (mm)']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.H2O_LY2);
      title(['H2O_LY2 (mm)']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.H2O_LY3);
      title(['H2O_LY3 (mm)']); 

      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.SM_LY1);
      title(['SM_LY1 (mm)']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.SM_LY2);
      title(['SM_LY2 (mm)']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.SM_LY3);
      title(['SM_LY3 (mm)']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.LF_LY1);
      title(['LF_LY1 (mm)']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.LF_LY2);
      title(['LF_LY2 (mm)']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.LF_LY3);
      title(['LF_LY3 (mm)']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.E_LY1);
      title(['E_LY1 (mm)']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.E_LY2);
      title(['E_LY2 (mm)']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.E_LY3);
      title(['E_LY3 (mm)']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.TEMP_LY1);
      title(['TEMP_LY1 (mm)']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.TEMP_LY2);
      title(['TEMP_LY2 (mm)']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.TEMP_LY3);
      title(['TEMP_LY3 (mm)']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.H2O_LY2);
      title(['SCF (mm)']); 
      
%       subplot(nfrow,nfcol,isub);isub=isub+1;      
%       plotunc(mod.soil_moist);
%       title(['soil moist (0-1)']); % sm_PAW
      
%       subplot(nfrow,nfcol,isub);isub=isub+1;      
%       plotunc(mod.sm_PAW);
%       title(['sm_PAW (0-1)']);
%       
%       subplot(nfrow,nfcol,isub);isub=isub+1;      
%       plotunc(mod.sm_PUW);
%       title(['sm_PUW (0-1)']);
% -----------------------------------------------------
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.evap);
      title(['evap (mm/day)']);  
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.qsurf);
      title(['qsurf (mm/day)']);
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.transp1);
      title(['transp1 (mm/day)']);
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.transp2);
      title(['transp2 (mm/day)']);
      
%       subplot(nfrow,nfcol,isub);isub=isub+1;      
%       plotunc(mod.paw2puw);
%       title(['paw2puw (mm/day)']);
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.melt);
      title(['melt (mm/day)']);
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.infil);
      title(['infil (mm/day)']);
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.prec);
      title(['prec (mm/day)']);
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.snowfall);
      title(['snowfall (mm/day)']);
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      hist(1000*CBR.PARS(:,MD.PARAMETER_IDs.LY1_z));xlim([1000*MD.parmin(MD.PARAMETER_IDs.LY1_z) 1000*MD.parmax(MD.PARAMETER_IDs.LY1_z)]);
      title(['LY1_z (converted to mm)']);
      ylimits = ylim; ymax = ylimits(2);
      format shortg
      text(10,ymax*3/4,['median=',num2str(round(median(CBR.PARS(:,MD.PARAMETER_IDs.LY1_z)),2))]);
      text(10,ymax*2/3,['obs=',num2str(round(obs.PEQ_PAW_z.max_value,2))]);

      subplot(nfrow,nfcol,isub);isub=isub+1;      
      hist(CBR.PARS(:,MD.PARAMETER_IDs.LY1_por));xlim([MD.parmin(MD.PARAMETER_IDs.LY1_por) MD.parmax(MD.PARAMETER_IDs.LY1_por)]);
      title(['LY1 por (0-1)']);
      ylimits = ylim; ymax = ylimits(2);
      format shortg
      text(10,ymax*3/4,['median=',num2str(round(median(CBR.PARS(:,MD.PARAMETER_IDs.LY1_por)),2))]);
      
%       subplot(nfrow,nfcol,isub);isub=isub+1; 
%       A=1000.*CBR.PARS(:,MD.PARAMETER_IDs.LY1_z).*CBR.PARS(:,MD.PARAMETER_IDs.LY1_por);
%       B=repmat(A,1,size(mod.ABGB,2));
%       plotunc(B);
%       title(['LY1 por*LY1 depth in mm']); 
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      hist(CBR.PARS(:,MD.PARAMETER_IDs.S_fv));
%       xlim([MD.parmin(MD.PARAMETER_IDs.LY1_por) MD.parmax(MD.PARAMETER_IDs.LY1_por)]);
      xlim([1 100])
      title(['S_fv (1-100)']); 
      ylimits = ylim; ymax = ylimits(2);
      format shortg
      text(10,ymax*3/4,['median=',num2str(round(median(CBR.PARS(:,MD.PARAMETER_IDs.S_fv)),2))]);
      
      subplot(nfrow,nfcol,isub);isub=isub+1;      
      plotunc(mod.EWT); hold on; plot(ind.GRACE_EWT,'r-','LineWidth',3);title(['EWT (mm)']);
          yhat = median(mod.EWT(:,1:(end-1)),1)';
          y = ind.GRACE_EWT;
%           y = y(1:xmax);
          % get stats metric
          RMSE = sqrt(nanmean((y - yhat).^2));  % Root Mean Squared Error
          mdl = fitlm(y,yhat);
          R2=mdl.Rsquared.Ordinary;

          subplot(nfrow,nfcol,isub);isub=isub+1; 
          scatter(yhat,y);title(['exp' num2str(iexp) ' EWT']);
          ylabel('obs');xlabel('mod');
          
          NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
          SW = [min(xlim) min(ylim)]+[diff(xlim) diff(ylim)]*0.05;
          NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
          SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
          txt1 = ['rmse = ', num2str(RMSE)];
          text(NW(1), NW(2),txt1);
          txt2 = ['r2 = ', num2str(R2)];
          text(SW(1), SW(2),txt2);
          
      sgtitle(['site' num2str(isiteindex) '_' char(isiteid) ' exp' num2str(iexp)]); 

      %%
      clear figure_property;
figure_property.units = 'inches';
figure_property.format = 'pdf';
figure_property.Preview= 'none';
figure_property.Width= '25'; % Figure width on canvas
figure_property.Height= '15'; % Figure height on canvas
figure_property.Units= 'inches';
figure_property.Color= 'rgb';
figure_property.Background= 'w';
figure_property.FixedfontSize= '6';
figure_property.ScaledfontSize= 'auto';
figure_property.FontMode= 'scaled';
figure_property.FontSizeMin= '6';
figure_property.FixedLineWidth= '2';
figure_property.ScaledLineWidth= 'auto';
figure_property.LineMode= 'none';
figure_property.LineWidthMin= '0.5';
figure_property.FontName= 'Times New Roman';% Might want to change this to something that is available
figure_property.FontWeight= 'auto';
figure_property.FontAngle= 'auto';
figure_property.FontEncoding= 'latin1';
figure_property.PSLevel= '3';
figure_property.Renderer= 'painters';
figure_property.Resolution= '100';
figure_property.LineStyleMap= 'none';
figure_property.ApplyStyle= '0';
figure_property.Bounds= 'tight';
figure_property.LockAxes= 'off';
figure_property.LockAxesTicks= 'off';
figure_property.ShowUI= 'off';
figure_property.SeparateText= 'off';
chosen_figure=gcf;
set(chosen_figure,'PaperUnits','inches');
set(chosen_figure,'PaperPositionMode','auto');
set(chosen_figure,'PaperSize',[str2num(figure_property.Width) str2num(figure_property.Height)]); % Canvas Size
set(chosen_figure,'Units','inches');
hgexport(gcf,[pdfpath 'validation_water_site' num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '.pdf'],figure_property); %Set desired file name
clf(f2);
clf(f2,'reset')
%%

% save likelihood
        likelihood = array2table(CBR.LIKELIHOODS);
        likelihood.Properties.VariableNames(1:size(likelihood,2)) = LIKELIHOOD_colnames;
        ofilename = [csvpath num2str(isiteindex) '_' char(isiteid) '.xlsx'];
        writetable(likelihood,ofilename,'Sheet',['exp' num2str(iexp)])
        temp_names = LIKELIHOOD_colnames;
% save chi2
        chi2 = CBR.LIKELIHOODS;
        for icol=1:size(CBR.LIKELIHOODS,2)
            if sum(CBR.LIKELIHOODS(:,icol))~=0 % if mod-data mismatch has values
%             if isfield(CBF,temp_names{icol})
                noo = sum(~isnan(CBF.(temp_names{icol}).values)); % calculate number of observations
                chi2(:,icol) = chi2(:,icol).*2./noo;
            end
        end
        t_chi2 = array2table(chi2);
        t_chi2.Properties.VariableNames(1:size(t_chi2,2)) = LIKELIHOOD_colnames;
%         ofilename = [csvpath 'site' pidchar '.xlsx'];
        writetable(t_chi2,ofilename,'Sheet',['chi2_exp' num2str(iexp)])

% save mean of likelihood, each site has a row
        mean_likehood.(['exp' num2str(iexp)])(i,:) = nanmean(CBR.LIKELIHOODS,1);
% save mean of chi2, each site has a row
        mean_chi2.(['exp' num2str(iexp)])(i,:) = nanmean(chi2,1);
  end
end

%% save water balance check table

% savetable = array2table(savedata);
%           savetable.Properties.RowNames(1:size(savetable,1)) = ...
%     {'mean', 'v05','v25','v50','v75','v95'};
%           ofilename = [opath 'site' num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '_' char(varnames(idata)) '.csv'];
%           writetable(savetable,ofilename,'WriteRowNames', true,'WriteVariableNames', false)

%%


% for iexp = [2 3 4 5 6]
% for iexp = [2 3]
for iexp = [8]
    ofilename = [csvpath 'A1_hist.xlsx'];
    savetable = array2table(mean_likehood.(['exp' num2str(iexp)]));
    savetable.Properties.VariableNames(1:size(savetable,2)) = LIKELIHOOD_colnames;
    writetable(savetable,ofilename,'Sheet',['likelihood_exp' num2str(iexp)])

    savetable = array2table(mean_chi2.(['exp' num2str(iexp)]));
    savetable.Properties.VariableNames(1:size(savetable,2)) = LIKELIHOOD_colnames;
    writetable(savetable,ofilename,'Sheet',['chi2_exp' num2str(iexp)])
end
