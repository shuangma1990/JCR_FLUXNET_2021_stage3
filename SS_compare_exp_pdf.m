MD=CARDAMOM_MODEL_LIBRARY(1100,[],1);
disp('read model library')
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

% for i=1:length(coord.cbfname)
for i=20%57

%     i=58;
    disp(['site' num2str(i)]);

    irow = str2num(rowv(i));
    icol = str2num(colv(i));
    isiteindex = str2num(siteindex(i));
    isiteid = (siteid(i));

    cbfdir_allobs = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P62/CBF_FILES/' num2str(isiteindex) '_' char(isiteid) 'CBF_exp13.cbf.nc'];
    cbfdirexp6 = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P61/CBF_FILES/' num2str(isiteindex) '_' char(isiteid) '_exp6.cbf.nc'];
    cbfdirexp11 = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P62/CBF_FILES/' num2str(isiteindex) '_' char(isiteid) 'CBF_exp11.cbf.nc'];
    cbfdirexp16 = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P64/CBF_FILES/' num2str(isiteindex) '_' char(isiteid) 'CBF_exp16.cbf.nc'];

    cbrdirexp6 = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P61/CBR_FILES/' num2str(isiteindex) '_' char(isiteid) '_exp6.cbf.cbr'];
    cbrdirexp11 = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P62/CBR_FILES/' num2str(isiteindex) '_' char(isiteid) 'CBF_exp11.cbf.cbr'];
    cbrdirexp16 = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P64/CBR_FILES/' num2str(isiteindex) '_' char(isiteid) 'CBF_exp16.cbf.cbr'];
    CBRexp6=CARDAMOM_RUN_MODEL(cbfdirexp6,cbrdirexp6);
    CBRexp11=CARDAMOM_RUN_MODEL(cbfdirexp11,cbrdirexp11);
    CBRexp16=CARDAMOM_RUN_MODEL(cbfdirexp16,cbrdirexp16);

    CBF_allobs=CARDAMOM_READ_NC_CBF_FILE(cbfdir_allobs);
    CBFexp6=CARDAMOM_READ_NC_CBF_FILE(cbfdirexp6);
    CBFexp11=CARDAMOM_READ_NC_CBF_FILE(cbfdirexp11);
    CBFexp16=CARDAMOM_READ_NC_CBF_FILE(cbfdirexp16);

    [mod_exp6,obs_exp6] = get_cardamom4plot_main(CBFexp6,CBRexp6,MD);
    [mod_exp11,obs_exp11] = get_cardamom4plot_main(CBFexp11,CBRexp11,MD);
    [mod_exp16,obs_exp16] = get_cardamom4plot_main(CBFexp16,CBRexp16,MD);
    
    [mod_nomimal,obs_allobs] = get_cardamom4plot_main(CBF_allobs,CBRexp6,MD);

    nrow=4;ncol=6;count=1;
    figure(1);clf
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp6.PARS(:,MD.PARAMETER_IDs.fwc)));title(['fwc w CH4 obs']);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp11.PARS(:,MD.PARAMETER_IDs.fwc)));title(['fwc wo CH4 obs']);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp6.PARS(:,MD.PARAMETER_IDs.Q10ch4)));title(['Q10ch4 w CH4 obs']);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp11.PARS(:,MD.PARAMETER_IDs.Q10ch4)));title(['Q10ch4 wo CH4 obs']);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp6.PARS(:,MD.PARAMETER_IDs.Q10rhco2)));title(['Q10rhco2 w CH4 obs']);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp11.PARS(:,MD.PARAMETER_IDs.Q10rhco2)));title(['Q10rhco2 wo CH4 obs']);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp6.PARS(:,MD.PARAMETER_IDs.r_ch4)));title(['r_ch4 w CH4 obs']);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp11.PARS(:,MD.PARAMETER_IDs.r_ch4)));title(['r_ch4 wo CH4 obs']);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp6.PARS(:,MD.PARAMETER_IDs.S_fv)));title(['S_fv w CH4 obs']);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp11.PARS(:,MD.PARAMETER_IDs.S_fv)));title(['S_fv wo CH4 obs']);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp6.PARS(:,MD.PARAMETER_IDs.max_infil)));title(['max_infil w CH4 obs']);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp11.PARS(:,MD.PARAMETER_IDs.max_infil)));title(['max_infil wo CH4 obs']);    
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp6.PARS(:,MD.PARAMETER_IDs.LY1_z)));title(['LY1_z w CH4 obs']);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp11.PARS(:,MD.PARAMETER_IDs.LY1_z)));title(['LY1_z wo CH4 obs']); 
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp6.PARS(:,MD.PARAMETER_IDs.LY1_por)));title(['LY1_por w CH4 obs']);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp11.PARS(:,MD.PARAMETER_IDs.LY1_por)));title(['LY1_por wo CH4 obs']); 
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp6.PARS(:,MD.PARAMETER_IDs.thermal_cond_surf)));title(['thermal_cond_surf w CH4 obs']);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp11.PARS(:,MD.PARAMETER_IDs.thermal_cond_surf)));title(['thermal_cond_surf wo CH4 obs']);     
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp6.PARS(:,MD.PARAMETER_IDs.thermal_cond)));title(['thermal_cond w CH4 obs']);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp11.PARS(:,MD.PARAMETER_IDs.thermal_cond)));title(['thermal_cond wo CH4 obs']); %     subplot(nrow,ncol,9);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp6.PARS(:,MD.PARAMETER_IDs.LY1_vhc)));title(['LY1_vhc w CH4 obs']);
    subplot(nrow,ncol,count);count=count+1;
    hist(log(CBRexp11.PARS(:,MD.PARAMETER_IDs.LY1_vhc)));title(['LY1_vhc wo CH4 obs']); %     subplot(nrow,ncol,9);

    
%     fontsizev = 20;
%     set(gca,'fontsize', fontsizev);
x0=10;
y0=10;
width=950;
height=400;
    set(gcf,'position',[x0,y0,width,height])
    export_fig(['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output_pdf/' char(isiteid) '_pdf.pdf']); 
    %     append_pdfs(['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/figure.pdf'],['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/figure.pdf'], ['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/figure'  num2str(iplot) '.pdf'])
end
