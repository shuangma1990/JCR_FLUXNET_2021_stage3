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
for i=57

%     i=58;
    disp(['site' num2str(i)]);

    irow = str2num(rowv(i));
    icol = str2num(colv(i));
    isiteindex = str2num(siteindex(i));
    isiteid = (siteid(i));

    cbfdir_allobs = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P62/CBF_FILES/' num2str(isiteindex) '_' char(isiteid) 'CBF_exp13.cbf.nc'];
    cbfdirexp6 = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P61/CBF_FILES/' num2str(isiteindex) '_' char(isiteid) '_exp6.cbf.nc'];
    cbfdirexp11 = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P62/CBF_FILES/' num2str(isiteindex) '_' char(isiteid) 'CBF_exp11.cbf.nc'];
    cbrdirexp6 = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P61/CBR_FILES/' num2str(isiteindex) '_' char(isiteid) '_exp6.cbf.cbr'];
    cbrdirexp11 = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P62/CBR_FILES/' num2str(isiteindex) '_' char(isiteid) 'CBF_exp11.cbf.cbr'];
    CBRexp6=CARDAMOM_RUN_MODEL(cbfdirexp6,cbrdirexp6);
    CBRexp11=CARDAMOM_RUN_MODEL(cbfdirexp11,cbrdirexp11);
    CBF_allobs=CARDAMOM_READ_NC_CBF_FILE(cbfdir_allobs);
    CBFexp6=CARDAMOM_READ_NC_CBF_FILE(cbfdirexp6);
    CBFexp11=CARDAMOM_READ_NC_CBF_FILE(cbfdirexp11);
    [mod_exp6,obs_exp6] = get_cardamom4plot_main(CBFexp6,CBRexp6,MD);
    [mod_exp11,obs_exp11] = get_cardamom4plot_main(CBFexp11,CBRexp11,MD);
    [mod_nomimal,obs_allobs] = get_cardamom4plot_main(CBF_allobs,CBRexp6,MD);

    fluxnames =  {'NBE','ET','rh_ch4','rh_co2','GPP','Ra','runoff','infil',...
        'melt','qsurf','q_ly1','ly1xly2','ly2xly3','LF_LY1','LF_LY2','LF_LY3',...
        'SM_LY1','SM_LY2','SM_LY3','TEMP_LY1','TEMP_LY2','TEMP_LY3'};
    obsnames = {'NBE','ET','CH4'};
    % fluxname = 'runoff';

    for iplot=1:length(fluxnames) 
        fluxname = fluxnames{iplot};
        fontsizev = 20;
        figure(1);clf
        a=mod_exp6.(fluxname);
        b=mod_exp11.(fluxname);
        subplot(2,2,1);plotunc(a);title([fluxname ' w CH4 obs']);set(gca,'fontsize', fontsizev);
        if iplot==1 || iplot==2 || iplot==3
            hold on;
            plot(1:length(obs_exp6.(char(obsnames(iplot))).values),obs_exp6.(char(obsnames(iplot))).values,'LineWidth',2);
        end
        subplot(2,2,2);plotunc(b);title([fluxname ' wo CH4 obs']);set(gca,'fontsize', fontsizev);
        if iplot==1 || iplot==2 || iplot==3
            hold on;
            plot(1:length(obs_exp11.(char(obsnames(iplot))).values),obs_exp11.(char(obsnames(iplot))).values,'LineWidth',2);
        end 
        if iplot==14 || iplot==15 || iplot==16
            hold on;
            yline(0,'b--');
        end
        subplot(2,2,3); plot(median(a(end-999:end,:)),median(b(end-999:end,:)),'o');title('Comparison');set(gca,'fontsize', fontsizev);
        refline(1,0); ylabel([fluxname ' wo CH4 obs']); xlabel([fluxname ' w CH4 obs'])
        xrange = prctile(mean(a,2),[1 99],"all");
        subplot(4,4,11);hist(mean(a,2));xlim(xrange);title([fluxname ' w CH4 obs']);set(gca,'fontsize', fontsizev);
        xrange = prctile(mean(b,2),[1 99],"all");
        subplot(4,4,12);hist(mean(b,2));xlim(xrange);title([fluxname ' wo CH4 obs']);
        set(gca,'fontsize', fontsizev);
        mkdir(['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/' char(isiteid)]);
        export_fig(['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/' char(isiteid) '/figure_'  fluxname '.pdf']); 
    %     append_pdfs(['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/figure.pdf'],['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/figure.pdf'], ['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/figure'  num2str(iplot) '.pdf'])

    end

    pdfnames = [strcat('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/', char(isiteid), '/figure_',fluxnames,'.pdf')];
    append_pdfs(['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/' char(isiteid) '/figure.pdf'],pdfnames);

end








CBRexp11.PARS(:,MD.PARAMETER_IDs.max_infil)=CBRexp6.PARS(:,MD.PARAMETER_IDs.max_infil);

CBRexp11=CARDAMOM_RUN_MODEL(cbfdirexp11,cbrdirexp11);
CBRexp11_modi=CARDAMOM_RUN_MODEL(cbfdirexp11,CBRexp11.PARS);

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
for i=57

%     i=58;
    disp(['site' num2str(i)]);

    irow = str2num(rowv(i));
    icol = str2num(colv(i));
    isiteindex = str2num(siteindex(i));
    isiteid = (siteid(i));

    cbfdir_allobs = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P62/CBF_FILES/' num2str(isiteindex) '_' char(isiteid) 'CBF_exp13.cbf.nc'];
    cbfdirexp6 = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P61/CBF_FILES/' num2str(isiteindex) '_' char(isiteid) '_exp6.cbf.nc'];
    cbfdirexp11 = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P62/CBF_FILES/' num2str(isiteindex) '_' char(isiteid) 'CBF_exp11.cbf.nc'];
    cbrdirexp6 = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P61/CBR_FILES/' num2str(isiteindex) '_' char(isiteid) '_exp6.cbf.cbr'];
    cbrdirexp11 = ['/Users/shuangma/GATTACA2_LOCAL_MIRROR/fluxnet_ch4_S4_P62/CBR_FILES/' num2str(isiteindex) '_' char(isiteid) 'CBF_exp11.cbf.cbr'];
    CBRexp6=CARDAMOM_RUN_MODEL(cbfdirexp6,cbrdirexp6);
    CBRexp11=CARDAMOM_RUN_MODEL(cbfdirexp11,cbrdirexp11);
    CBF_allobs=CARDAMOM_READ_NC_CBF_FILE(cbfdir_allobs);
    CBFexp6=CARDAMOM_READ_NC_CBF_FILE(cbfdirexp6);
    CBFexp11=CARDAMOM_READ_NC_CBF_FILE(cbfdirexp11);
    [mod_exp6,obs_exp6] = get_cardamom4plot_main(CBFexp6,CBRexp6,MD);
    [mod_exp11,obs_exp11] = get_cardamom4plot_main(CBFexp11,CBRexp11,MD);
    [mod_nomimal,obs_allobs] = get_cardamom4plot_main(CBF_allobs,CBRexp6,MD);

    fluxnames =  {'NBE','ET','rh_ch4','rh_co2','GPP','Ra','runoff','infil',...
        'melt','qsurf','q_ly1','ly1xly2','ly2xly3','LF_LY1','LF_LY2','LF_LY3',...
        'SM_LY1','SM_LY2','SM_LY3','TEMP_LY1','TEMP_LY2','TEMP_LY3'};
    obsnames = {'NBE','ET','CH4'};
    % fluxname = 'runoff';

    for iplot=1:length(fluxnames) 
        fluxname = fluxnames{iplot};
        fontsizev = 20;
        figure(1);clf
        a=mod_exp6.(fluxname);
        b=mod_exp11.(fluxname);
        subplot(2,2,1);plotunc(a);title([fluxname ' w CH4 obs']);set(gca,'fontsize', fontsizev);
        if iplot==1 || iplot==2 || iplot==3
            hold on;
            plot(1:length(obs_exp6.(char(obsnames(iplot))).values),obs_exp6.(char(obsnames(iplot))).values,'LineWidth',2);
        end
        subplot(2,2,2);plotunc(b);title([fluxname ' wo CH4 obs']);set(gca,'fontsize', fontsizev);
        if iplot==1 || iplot==2 || iplot==3
            hold on;
            plot(1:length(obs_exp11.(char(obsnames(iplot))).values),obs_exp11.(char(obsnames(iplot))).values,'LineWidth',2);
        end 
        if iplot==14 || iplot==15 || iplot==16
            hold on;
            yline(0,'b--');
        end
        subplot(2,2,3); plot(median(a(end-999:end,:)),median(b(end-999:end,:)),'o');title('Comparison');set(gca,'fontsize', fontsizev);
        refline(1,0); ylabel([fluxname ' wo CH4 obs']); xlabel([fluxname ' w CH4 obs'])
        xrange = prctile(mean(a,2),[1 99],"all");
        subplot(4,4,11);hist(mean(a,2));xlim(xrange);title([fluxname ' w CH4 obs']);set(gca,'fontsize', fontsizev);
        xrange = prctile(mean(b,2),[1 99],"all");
        subplot(4,4,12);hist(mean(b,2));xlim(xrange);title([fluxname ' wo CH4 obs']);
        set(gca,'fontsize', fontsizev);
        mkdir(['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/' char(isiteid)]);
        export_fig(['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/' char(isiteid) '/figure_'  fluxname '.pdf']); 
    %     append_pdfs(['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/figure.pdf'],['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/figure.pdf'], ['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/figure'  num2str(iplot) '.pdf'])

    end

    pdfnames = [strcat('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/', char(isiteid), '/figure_',fluxnames,'.pdf')];
    append_pdfs(['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/' char(isiteid) '/figure.pdf'],pdfnames);

end

CBRexp11=CARDAMOM_RUN_MODEL(cbfdirexp11,cbrdirexp11);
CBRexp11.PARS(:,MD.PARAMETER_IDs.max_infil)=CBRexp6.PARS(:,MD.PARAMETER_IDs.max_infil);
CBRexp11_modi=CARDAMOM_RUN_MODEL(cbfdirexp11,CBRexp11.PARS);

[mod_exp11_modi,obs_exp11_modi] = get_cardamom4plot_main(CBFexp11,CBRexp11_modi,MD);




% replace pdf
fluxnames =  {'NBE','ET','rh_ch4','rh_co2','GPP','Ra','runoff','infil',...
        'melt','qsurf','q_ly1','ly1xly2','ly2xly3','LF_LY1','LF_LY2','LF_LY3',...
        'SM_LY1','SM_LY2','SM_LY3','TEMP_LY1','TEMP_LY2','TEMP_LY3'};
    obsnames = {'NBE','ET','CH4'};
    % fluxname = 'runoff';

    for iplot=1:length(fluxnames) 
        fluxname = fluxnames{iplot};
        fontsizev = 20;
        figure(1);clf
        a=mod_exp6.(fluxname);
        b=mod_exp11_modi.(fluxname);
        subplot(2,2,1);plotunc(a);title([fluxname ' w CH4 obs']);set(gca,'fontsize', fontsizev);
        if iplot==1 || iplot==2 || iplot==3
            hold on;
            plot(1:length(obs_exp6.(char(obsnames(iplot))).values),obs_exp6.(char(obsnames(iplot))).values,'LineWidth',2);
        end
        subplot(2,2,2);plotunc(b);title([fluxname ' wo CH4 obs']);set(gca,'fontsize', fontsizev);
        if iplot==1 || iplot==2 || iplot==3
            hold on;
            plot(1:length(obs_exp11.(char(obsnames(iplot))).values),obs_exp11.(char(obsnames(iplot))).values,'LineWidth',2);
        end 
        if iplot==14 || iplot==15 || iplot==16
            hold on;
            yline(0,'b--');
        end
        subplot(2,2,3); plot(median(a(end-999:end,:)),median(b(end-999:end,:)),'o');title('Comparison');set(gca,'fontsize', fontsizev);
        refline(1,0); ylabel([fluxname ' wo CH4 obs']); xlabel([fluxname ' w CH4 obs'])
        xrange = prctile(mean(a,2),[1 99],"all");
        subplot(4,4,11);hist(mean(a,2));xlim(xrange);title([fluxname ' w CH4 obs']);set(gca,'fontsize', fontsizev);
        xrange = prctile(mean(b,2),[1 99],"all");
        subplot(4,4,12);hist(mean(b,2));xlim(xrange);title([fluxname ' wo CH4 obs']);
        set(gca,'fontsize', fontsizev);
        mkdir(['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/' char(isiteid)]);
        export_fig(['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/' char(isiteid) '/figure_'  fluxname '.pdf']); 
    %     append_pdfs(['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/figure.pdf'],['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/figure.pdf'], ['RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/SS_output/figure'  num2str(iplot) '.pdf'])

    end
    
    
% scatter plot    
    CBRexp6=CARDAMOM_RUN_MODEL(cbfdirexp6,cbrdirexp6);
    CBRexp11=CARDAMOM_RUN_MODEL(cbfdirexp11,cbrdirexp11);
    CBFexp6=CARDAMOM_READ_NC_CBF_FILE(cbfdirexp6);
    CBFexp11=CARDAMOM_READ_NC_CBF_FILE(cbfdirexp11);
    
    [mod_exp6,obs_exp6] = get_cardamom4plot_main(CBFexp6,CBRexp6,MD);
    [mod_exp11,obs_exp11] = get_cardamom4plot_main(CBFexp11,CBRexp11,MD);
    
    subplot(2,2,1); 
    scatter(CBRexp6.PARS(:,MD.PARAMETER_IDs.max_infil),mean(mod_exp6.infil,2));title('x: max infil    y: infil')
    subplot(2,2,2); 
    scatter(mean(mod_exp6.infil,2),mean(mod_exp6.melt,2));title('x: infil    y: melt');
    
    subplot(2,2,3); 
    scatter(CBRexp11.PARS(:,MD.PARAMETER_IDs.max_infil),mean(mod_exp11.infil,2));title('x: max infil    y: infil')
    subplot(2,2,4); 
    scatter(mean(mod_exp11.infil,2),mean(mod_exp11.melt,2));title('x: infil    y: melt');
    