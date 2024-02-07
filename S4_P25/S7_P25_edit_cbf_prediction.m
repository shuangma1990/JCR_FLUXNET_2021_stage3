
nop = 'S4_P25_prediction_1'; 
PROJNAME = ['fluxnet_ch4_' nop];
pathname2=sprintf('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/input/');
pathname3 = ['/Users/shuangma/GATTACA_LOCAL_MIRROR/fluxnet_ch4_S4_P25/']; % previous MCMC io dir
pathname4 = ['/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/' nop]; % current MCMC dir
ncdir = ['/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/' nop '/input/'];
if not(isfolder(ncdir))
    mkdir(ncdir)
end
if not(isfolder(pathname4))
    mkdir(pathname4)
end          
% load coordinates info for all the fluxnet ch4 sites
% coord=readtable('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S2_LATLON2CBF_edited.csv');
coord=readtable('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/siteid_5yrs_goodfit.csv');

rowv = string(coord.r);
colv = string(coord.c);
coordv=string(coord.cbfname);
siteid=string(coord.SITE_ID);
siteindex=string(coord.nsite);

for i=1:length(coord.cbfname)
% for i=1:5
    irow = str2num(rowv(i));
    icol = str2num(colv(i));
    isiteindex = str2num(siteindex(i));
    isiteid = (siteid(i));
% read preprocessed full obs nc files
    dinfo1 = dir([pathname2 '*_' char(isiteid) '.cbf.nc']);
    CBF_exp = CARDAMOM_READ_NC_CBF_FILE([char({dinfo1.folder}) '/' char(dinfo1.name)]);
    dinfo2 = dir([pathname2 '*_' char(isiteid) '_noobs.cbf.nc']);
    CBF_exp_noobs = CARDAMOM_READ_NC_CBF_FILE([char({dinfo2.folder}) '/' char(dinfo2.name)]);
% change model ID here if needed:
%     CBF_exp.ID.values = 1100;
%     CBF_exp_noobs.ID.values = 1100;
%%  insert info for PEQs from previous MCMC
%   exp6 was constrained by RS and EC, use its par posteriors to inform:
%   exp1 no ons prediction, exp6 RS+EC, exp7 RS only, exp6 is a control
    dinfo3 = dir([pathname3 'CBR_FILES/' '*_' char(isiteid) '_exp6.cbf.cbr']);
    cbrdir = [char({dinfo3.folder}) '/' char(dinfo3.name)];
    dinfo4 = dir([pathname3 'CBF_FILES/' '*_' char(isiteid) '_exp6.cbf.nc']);
    cbfdir = [char({dinfo4.folder}) '/' char(dinfo4.name)];
    pre_CBR=CARDAMOM_RUN_MODEL(cbfdir,cbrdir);
    pause(10);
%    pre_CBR=CARDAMOM_RUN_MDF(cbfdir,cbrdir);
    meanvec = exp(mean(log(pre_CBR.PARS),1)); % exp(mean(log(pars))% log() is ln, loge
    sdvec = exp(std(log(pre_CBR.PARS),1));
% -- -- -- --  get CBF_exps ready and submit to gattaca  -- --  --  -- 
    [CBF_cell,string_CBF_cell] = edit_fluxnetexps_used_in_P23_P25_prediction(CBF_exp,CBF_exp_noobs,meanvec,sdvec); 
    
    nc_cell={};
%       for icell=1:length(CBF_cell)
%         for icell=[2 3 4 5 6]
        for icell=[1 6 7]
          fname = [num2str(isiteindex) '_' char(isiteid) char(string_CBF_cell(icell)) '.cbf.nc'];
          nc_cell{icell} = [ncdir fname];
          CARDAMOM_WRITE_NC_CBF_FILE((CBF_cell{icell}),[ncdir fname]);
%% uncomment this line to submit to server (use foldername instead)
%           pause(2)
%           PARRFUN_CARDAMOM_SUBMIT_TO_GATTACA_lg_edge_shuangma(PROJNAME,[nc_cell{icell}]);
%% uncomment this line to test on laptop
%           CBRtest=CARDAMOM_RUN_MDF([nc_cell{icell}]);
%             CBRtest=CARDAMOM_RUN_MDF(CBF_exp_noobs);
         end
end

%% uncomment this line to submit to server with updated JPL server submission script
 PARRFUN_CARDAMOM_SUBMIT_TO_GATTACA_lg_edge_shuangma(PROJNAME,ncdir);
%   PARRFUN_CARDAMOM_SUBMIT_TO_GATTACA_shuangma(PROJNAME,ncdir);