nop = 'S4_P61';
PROJNAME = ['fluxnet_ch4_' nop];
wd = (['/Users/shuangma/GATTACA2_LOCAL_MIRROR/' PROJNAME '/']);
odir = ['/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S9_sensitivity_test/'];
MD=CARDAMOM_MODEL_LIBRARY(1100,[],1);
disp('read model library')
fluxnames = fieldnames(MD.FLUX_IDs);
poolnames = fieldnames(MD.POOL_IDs);
parnames = fieldnames(MD.PARAMETER_IDs);
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
iexp=6;
delta=0.001; % a smaller delta value reduce the chance of model to crush due to change of par values 
sen_fluxes = nan(length(coord.cbfname),(MD.nopars),(MD.nofluxes));
sen_pools = nan(length(coord.cbfname),(MD.nopars),(MD.nopools));
%%
% find(coord.cbfname==317049)
for i=1:length(coord.cbfname)
%   for i=58
    irow = str2num(rowv(i));
    icol = str2num(colv(i));
    isiteindex = str2num(siteindex(i));
    isiteid = (siteid(i));    
     %$$
% for iexp=[2 3] %$$
%    for iexp= [6]
      cbfdir = [wd 'CBF_FILES/' num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '.cbf.nc'];
      cbfdir_allobs = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P45/CBF_FILES/' num2str(isiteindex) '_' char(isiteid) '_exp6.cbf.nc'];
      cbrdir = [wd 'CBR_FILES/' num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '.cbf.cbr'];
%       cbrdir = {[wd 'cbr_files/TR.1032.311_' char(pid) '_1.cbr'],[wd 'cbr_files/TR.1032.311_' char(pid) '_2.cbr'],[wd 'cbr_files/TR.1032.311_' char(pid) '_3.cbr'],[wd 'cbr_files/TR.1032.311_' char(pid) '_4.cbr']};
      try
        CBR=CARDAMOM_RUN_MODEL(cbfdir,cbrdir);
      catch ME
          fprintf([char(isiteid) num2str(iexp) 'get CBR wo success']);
          string3(istr) = sprintf([char(isiteid) num2str(iexp) 'get CBR wo success']);istr=istr+1;
          continue;
      end
      disp('read CBR')
      xvalue = 1:size(CBR.FLUXES,2);
      CBF=CARDAMOM_READ_NC_CBF_FILE(cbfdir);
      CBF_allobs=CARDAMOM_READ_NC_CBF_FILE(cbfdir_allobs);

      [mod,obs] = get_cardamom4plot_main(CBF,CBR,MD);
      [mod_nomimal,obs_allobs] = get_cardamom4plot_main(CBF_allobs,CBR,MD);   
      
      %% perturbation starts
      CBR_ori = CARDAMOM_RUN_MODEL(cbfdir,CBR.PARS(end,:));
      
      for ipar = 1:MD.nopars
          CBR_purb.PARS = CBR.PARS;
      % single MCMC perturbation
%           CBR_purb.PARS(end,ipar) = CBR.PARS(end,ipar) * (1+delta); 
%           CBR_purbed=CARDAMOM_RUN_MODEL(cbfdir,CBR_purb.PARS(end,:));
%           sen_fluxes(i,ipar,:) = squeeze(nanmean(squeeze(CBR_purbed.FLUXES-CBR_ori.FLUXES),1))/(CBR.PARS(end,ipar) * delta);
%           sen_pools(i,ipar,:) = squeeze(nanmean(squeeze((CBR_purbed.POOLS-CBR_ori.POOLS)),1))/(CBR.PARS(end,ipar) * delta);

% full MCMC perturbation
          CBR_purb.PARS(:,ipar) = CBR.PARS(:,ipar) * (1+delta); 
          CBR_purbed=CARDAMOM_RUN_MODEL(cbfdir,CBR_purb.PARS);
          a1 = squeeze(nanmean(CBR_purbed.FLUXES-CBR.FLUXES,2)) ./ (CBR.PARS(:,ipar) * delta);
          a2 = squeeze(nanmean(CBR_purbed.POOLS-CBR.POOLS,2))   ./ (CBR.PARS(:,ipar) * delta);
% remove outliers          
          idx1 = a1 < -1e2;
          idx2 = a1 > 1e2;
          a1(idx1)=NaN; a1(idx2)=NaN;
          
          idx1 = a2 < -1e2;
          idx2 = a2 > 1e2;
          a2(idx1)=NaN; a2(idx2)=NaN;
          
          sen_fluxes(i,ipar,:) = nanmean(a1,1);
          sen_pools(i,ipar,:) = nanmean(a2,1);
%% start of making figures
%           f1 = figure(1);
%           subplot(2,1,1);plot(a1');
%           subplot(2,1,2);plot(a2');
%           %           clear figure_property;
%           figure_property.units = 'inches';
%           figure_property.format = 'pdf';
%           figure_property.Preview= 'none';
%           figure_property.Width= '18'; % Figure width on canvas
%           figure_property.Height= '18'; % Figure height on canvas
%           figure_property.Units= 'inches';
%           figure_property.Color= 'rgb';
%           figure_property.Background= 'w';
%           figure_property.FixedfontSize= '6';
%           chosen_figure=gcf;
%           set(chosen_figure,'PaperUnits','inches');
%           set(chosen_figure,'PaperPositionMode','auto');
%           set(chosen_figure,'PaperSize',[str2num(figure_property.Width) str2num(figure_property.Height)]); % Canvas Size
%           set(chosen_figure,'Units','inches');
%           hgexport(gcf,[odir num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '_' char(parnames(ipar)) '.pdf'],figure_property); %Set desired file name
%           clf(f1);
%           clf(f1,'reset')
%% end of making figures
      end
end


% rh_CH4 is 55
% Q10_CH4 is 61


%% save csv

  for ifl = 1:MD.nofluxes
      savetable = array2table(sen_fluxes(:,:,ifl));
      savetable.Properties.VariableNames = string(parnames);
      savetable.Properties.RowNames = string(siteid);
      writetable(savetable,[odir 'site_par_' fluxnames{ifl} '.csv'],'WriteRowNames',true,'WriteVariableNames', true);
  end
  for ipl = 1:MD.nopools
      savetable = array2table(sen_pools(:,:,ipl));
      savetable.Properties.VariableNames = string(fieldnames(MD.PARAMETER_IDs));
      savetable.Properties.RowNames = string(siteid);
      writetable(savetable,[odir 'site_par_' poolnames{ipl} '.csv'],'WriteRowNames',true,'WriteVariableNames', true);
  end
  