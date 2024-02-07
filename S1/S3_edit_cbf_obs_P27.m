% due to the length of CBF drier file, FLUXNET data are cut after 2016
% P22 use 1011 to run at 67 fluxnet sites that has obs before 2016
% P23 add demcmc and et gpp ch4 threshold, based on P22
% P24 run obs LAI CH4 job, the rest is the same with P23
% P25 run larger unc cases (loosest), uses obs CH4(1),LAI CH4(7), and LAI(13) job,the rest is the same with P24
% P26 use ID 1005 and nbe lai 3 13 job, 1 unc to see if that's a problem
% with my model structure
% L1 is local quick test project, using 1011 but unc threshold from Yan
% abandoned P26 and L1
% P27 fixes the MLAI MGPP MFIRE with values and branched from P22
% demcmc
projectpurpose = 'P27'; % 
pp = projectpurpose;
% step 1 load cbfname file
coord=readtable('RESEARCH/WORKFLOW/JCR_FLUXNET_2021/S2_LATLON2CBF.csv');
coordv=string(coord.cbfname);
siteid=string(coord.SITE_ID);
siteindex=string(coord.nsite);
for i=1:length(coord.cbfname)
yendremove=0;
%Step 1. Load ori CBF file 
a=sprintf('RESEARCH/DATA/CARDAMOM_MET_0.5/GL05RUN_AUG19_%s.cbf',coordv(i));
CBF=CARDAMOM_READ_BINARY_FILEFORMAT(a);
% load site obs file
b=sprintf('RESEARCH/WORKFLOW/JCR_FLUXNET_2021/S1/%s_%s_obs.csv',siteindex(i),siteid(i));
obsdata=readtable(b);
    if size(obsdata,1)==0
        fprintf('%s%s has no obs data before 2016',siteindex(i),siteid(i));
    else
        if isempty(CBF.OBS.EWT) ==1
            CBF.OBS.EWT = repelem(-9999,192)';
            fprintf('this cbf file does not come with EWT obs');
        else

            % cut time length according to year from fluxnet obs
            % need to chage MET, OBS, nodays, original cbf span 2001-2016, 12*16=192
            ystart = min(obsdata.year)-3; % run 3 years before fusing data
            yend = max(obsdata.year)+1;   % run 1 more year after fusing data
            if ystart<2001   % earliest FLUXNET is in 2006, never smaller than 2001+3
                ystart=2001;
            end
            if yend >2016
                yend=2016;
                yendremove=1; % it'll make changes to CBF.OBS
            end
            CBF.nodays = (yend-ystart+1)*12;
            selectMET = CBF.MET(((ystart-2001)*12+1):((yend-2000)*12),:);
            CBF.MET = []; CBF.MET = selectMET;
            CBF.MET(:,7)=0; % NO FIRE
            yearmark = repelem((ystart:yend),12);
            %Step 2. overwrite CBF.ID 
            CBF.ID = 1011;
%             CBF.ID = 1005;
            a1 = table2array(obsdata(:,8));
            nonZeroIndexes = a1 ~= -9999;  
            theMean = mean(a1(nonZeroIndexes));
            CBF.OBSUNC.ET.obs_unc_threshold = max(theMean * 0.1,0.001);
            
            a2 = table2array(obsdata(:,11));
            nonZeroIndexes = a2 ~= -9999;  
            theMean = mean(a2(nonZeroIndexes));
            CBF.OBSUNC.GPP.obs_unc_threshold = theMean * 0.1;
            
            a3 = table2array(obsdata(:,4));
            nonZeroIndexes = a3 ~= -9999;  % 
            theMean = mean(a3(nonZeroIndexes));
            CBF.OBSUNC.CH4.obs_unc_threshold = max(theMean * 0.1,1e-4);
%             
%             CBF.OBSUNC.ET.obs_unc_threshold = table2array(obsdata(:,8));
%             CBF.OBSUNC.GPP.obs_unc_threshold = table2array(obsdata(:,11));
%             CBF.OBSUNC.CH4.obs_unc_threshold = table2array(obsdata(:,4));
            %Step 3. overwrite CBF.OBS
                % check if obs start from Jan
                obstartmonth = obsdata.month(1);
                obsendmonth = obsdata.month(end);
                spaceholder0 = repelem(-9999,(obstartmonth-1)); % if obs month start in the middle of year
                spaceholder00= repelem(-9999,(12-obsendmonth)); % if obs month end in the middle of year
                % add 3yrs of NA record to FLUXNET data, as space holder
            spaceholder1 = repelem(-9999,36);
            % add another year to the end if within 2016
            if yendremove==0
                spaceholder2 = repelem(-9999,12);
            else
                spaceholder2 = [];
            end
            CBF.OBS.NBE = []; CBF.OBS.NBE = [spaceholder0';spaceholder1';table2array(obsdata(:,8));spaceholder2';spaceholder00'];
            selectSIF = CBF.OBS.GPP(((ystart-2001)*12+1):((yend-2000)*12),:);% ori cbf obs is SIF
            CBF.OBS.GPP = []; CBF.OBS.GPP = [spaceholder0';spaceholder1';table2array(obsdata(:,11));spaceholder2';spaceholder00']; %using DayTime partition method
            CBF.OBSUNC.GPP.gppabs=1; % set flag to actual GPP, not SIF
            CBF.OBS.CH4 = []; CBF.OBS.CH4 = [spaceholder0';spaceholder1';table2array(obsdata(:,4));spaceholder2';spaceholder00'];
            selectLAI = CBF.OBS.LAI(((ystart-2001)*12+1):((yend-2000)*12),:);
            CBF.OBS.LAI = [];CBF.OBS.LAI = selectLAI;
            selectEWT = CBF.OBS.EWT(((ystart-2001)*12+1):((yend-2000)*12),:);
            CBF.OBS.EWT = [];CBF.OBS.EWT = selectEWT;
            
            % save other obs not used in MDF
            OBS.ER = [spaceholder0';spaceholder1';table2array(obsdata(:,13));spaceholder2';spaceholder00'];
%           save all available obs in cbf NBE GPP CH4 from tower,
%           EWT AND LAI from space
            N = [CBF.OBS.NBE,CBF.OBS.GPP,CBF.OBS.CH4,CBF.OBS.LAI,CBF.OBS.EWT,selectSIF,OBS.ER];
            TN = array2table(N);
            TN.Properties.VariableNames(1:7) = {'NBE', 'GPP','CH4','LAI','EWT','SIF','ER'};
            
            pathname = sprintf('RESEARCH/WORKFLOW/JCR_FLUXNET_2021/S3_%s',pp);
            if not(isfolder(pathname))
                mkdir(pathname)
            end
            filename = sprintf('%s/%s_modified_cbf_obs.csv',pathname,siteid(i));
            writetable(TN,filename);  % SIF, LAI, EWT
            
            %Step 4. save the new CBF stucture into a new binary file .cbf
            pathname2=sprintf('RESEARCH/WORKFLOW/JCR_FLUXNET_2021/S3_%s/input',pp);
            if not(isfolder(pathname2))
                mkdir(pathname2)
            end
            CARDAMOM_WRITE_BINARY_FILEFORMAT(CBF,sprintf('RESEARCH/WORKFLOW/JCR_FLUXNET_2021/S3_%s/input/%s_%s_%s.cbf',pp,siteid(i),siteindex(i),num2str(CBF.ID)));
            % Step 5. save to csv files 
            writematrix([yearmark' CBF.MET],sprintf('RESEARCH/WORKFLOW/JCR_FLUXNET_2021/S3_%s/%s_met.csv',pp,siteid(i))) 
        end
    end
end 
