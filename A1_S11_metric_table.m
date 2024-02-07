%% io paths % P45 had the exp1267, copied exp8 over to P45 for conveniency
obsdir = '/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/obs/';
% moddir = '/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S4_P45/output/';
moddir = '/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S4_P62/output/';

varnames = {'NBE','CH4','ET','SCF'};

%% load coordinates info for all the fluxnet ch4 sites
coord=readtable('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S2_LATLON2CBF_edited.csv');
% coord=readtable('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S2_LATLON2CBF_HL.csv');
% rowv = string(coord.r);
% colv = string(coord.c);
coordv=string(coord.cbfname);
siteid=string(coord.SITE_ID);
siteindex=string(coord.nsite);
%%
% (80,5,4,3*5)  80 sites, iexp, DIV RMSE R2 * 5 quantiles, fluxes
X = NaN(80,5,4,3*5); % P45 had the exp1267, copied exp8 over to P45 for conveniency
% a=squeeze(X(1,1,:,:))
for i=1:length(coord.cbfname)
%     for i=1:1

    disp(['site' num2str(i)]);
%     irow = str2num(rowv(i));
%     icol = str2num(colv(i));
    isiteindex = str2num(siteindex(i));
    isiteid = char(siteid(i));    
  
    countiexp=0;
    for iexp=[9:13] %$$
        countiexp = countiexp+1;
        try
        mod.NBE = csvread([moddir 'site' num2str(i) '_' isiteid '_exp' num2str(iexp) '_NBE.csv'],0,1);
            catch ME
        fprintf([num2str(i) '_' char(isiteid) num2str(iexp) ' does not exist\n']);    
            continue;
        end    
        mod.CH4 = csvread([moddir 'site' num2str(i) '_' isiteid '_exp' num2str(iexp) '_CH4.csv'],0,1);
        mod.ET = csvread([moddir 'site' num2str(i) '_' isiteid '_exp' num2str(iexp) '_ET.csv'],0,1);
        mod.SCF = csvread([moddir 'site' num2str(i) '_' isiteid '_exp' num2str(iexp) '_SCF.csv'],0,1);
        mod.GPP = csvread([moddir 'site' num2str(i) '_' isiteid '_exp' num2str(iexp) '_GPP.csv'],0,1);
        obsdata = readtable([obsdir num2str(i) '_' isiteid '_cbf_obs.csv']);
        obs.NBE = table2array(obsdata(:,1));
        obs.CH4 = table2array(obsdata(:,2));
        obs.ET  = table2array(obsdata(:,3));
        obs.SCF = table2array(obsdata(:,5)); % this is where bug is
        obs.GPP = table2array(obsdata(:,8));
        
        outdata = {mod.NBE,mod.CH4,mod.ET,mod.SCF,mod.GPP};
        valdata = {obs.NBE,obs.CH4,obs.ET,obs.SCF,obs.GPP};

        for o=1:5 % o is fluxes
              if o==4% i think this is to compare state at beginning of month?
                  yhat = outdata{o}(2:6,1:(end-1)); % first row is mean
              else
                  yhat = outdata{o}(2:6,:); % first row is mean
              end
%               y = permute(repmat(valdata{o},1,5),[2 1]);
              y = valdata{o}';
              % get stats metric
              count = 0;
              for iqt=1:5
                  count=count+1;
                  RMSE = sqrt(nanmean((y - yhat(iqt,:)).^2));  % Root Mean Squared Error
                  X(i,countiexp,o,count)=RMSE;
                  
                  count=count+1;
                  mdl = fitlm(y,yhat(iqt,:));
                  R2=mdl.Rsquared.Ordinary;
                  X(i,countiexp,o,count)=R2;
                  
                  count=count+1;
                  temp = yhat(iqt,:)./y;
                  temp(isinf(temp))=NaN;
                  DIV = nanmean(temp);
                  X(i,countiexp,o,count)=DIV; 
              end
        end
    end

end

% (80site,5iexp,4fluxes,3*5metric)
X1 = squeeze(nanmean(X,1));
X64 = squeeze(X(64,:,:,:));
X80 = squeeze(X(80,:,:,:));

metric.NBE = squeeze(X1(:,1,:))';
metric.CH4 = squeeze(X1(:,2,:))';
metric.ET = squeeze(X1(:,3,:))';
metric.SCF = squeeze(X1(:,4,:))';
metric.GPP = squeeze(X1(:,5,:))';

metric.SCF64 = squeeze(X64(:,4,:))';
metric.SCF80 = squeeze(X80(:,4,:))';