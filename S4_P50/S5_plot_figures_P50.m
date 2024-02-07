% P50 uses snow_only_shuang branch to fix SCF problem, EDC is on
% main_202306 is the staged branch
nop = 'S4_P50';
PROJNAME = ['fluxnet_ch4_' nop];
wd = (['/Users/shuangma/GATTACA_LOCAL_MIRROR/' PROJNAME '/']);
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
%% load coordinates info for all the fluxnet ch4 sites
coord=readtable('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S2_LATLON2CBF_edited.csv');
% coord=readtable('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S2_LATLON2CBF_HL.csv');
rowv = string(coord.r);
colv = string(coord.c);
coordv=string(coord.cbfname);
siteid=string(coord.SITE_ID);
siteindex=string(coord.nsite);

PROBmatrix=NaN(1000,80);
%%
MD=CARDAMOM_MODEL_LIBRARY(1100,[],1);

% find(coord.cbfname==317049)
for i=1:length(coord.cbfname)
% for i=[1 12 61]
%     i=58;
        irow = str2num(rowv(i));
        icol = str2num(colv(i));
        isiteindex = str2num(siteindex(i));
        isiteid = (siteid(i));    
        iexp = 6;
        cbfdir = [wd 'CBF_FILES/' num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '.cbf.nc'];
        cbrdir = [wd 'CBR_FILES/' num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '.cbf.cbr'];
        
        CBF=CARDAMOM_READ_NC_CBF_FILE(cbfdir);
        CBR=CARDAMOM_RUN_MODEL(cbfdir,cbrdir);

        PROBmatrix(:,i)=CBR.PROB;
        mSCF = CBR.POOLS(:,1:end,MD.POOL_IDs.D_SCF);
        oSCF = CBF.SCF.values(1:end);
        mSWE = CBR.POOLS(:,1:end,MD.POOL_IDs.H2O_SWE);
        
% save results to csv
        outdata = {mSCF}; % ,mod.SOM,mod.LAI,mod.leafC,mod.rh_co2,mod.rh_ch4,mod.soil_moist,mod.NPP,mod.Ra
        valdata = {oSCF}; % ,obs.SOM.values(1),obs.Mean_LAI.values
        varnames= {'SCF'};
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

% ploting
        f1=figure(1); nfrow=4; nfcol=1;isub=1;
%         subplot(nfrow,nfcol,isub);isub=isub+1;
%         plotunc(mSCF);hold on;plotmultilines(mSWE/200);hold on;plot(oSCF,'r*');
        subplot(nfrow,nfcol,isub);isub=isub+1;
        plotunc(mSCF);hold on;plot(oSCF,'r*');
        subplot(nfrow,nfcol,isub);isub=isub+1;
        plotmultilines(mSWE);
%         subplot(nfrow,nfcol,isub);isub=isub+1;
%         plot(mSWE(100,:),mSCF(100,:),'r*');
      % CBR=CARDAMOM_RUN_MODEL(cbfdir,cbrdir);
        subplot(nfrow,nfcol,isub);isub=isub+1;
        scatter(mSWE(1000,:),mSCF(1000,:));
        subplot(nfrow,nfcol,isub);isub=isub+1;
        scatter(mSWE(100,:),mSCF(100,:));
%         subplot(nfrow,nfcol,isub);isub=isub+1;
%         hist(nanmean(mSCF,2));
      
      sgtitle(['site' num2str(isiteindex) ' ' char(isiteid) ' exp' num2str(iexp)]);
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
hgexport(gcf,[pdfpath 'snow_only_gen3_site' num2str(isiteindex) '_' char(isiteid) '_exp' num2str(iexp) '.pdf'],figure_property); %Set desired file name
clf(f1);
clf(f1,'reset')

end  



% save parameter distributions to outputfolder
savetable = array2table(PROBmatrix);
savetable.Properties.VariableNames = string(siteid);
ofilename = [csvpath 'PROB_stats.csv'];
writetable(savetable,ofilename,'WriteRowNames', false,'WriteVariableNames', true)

