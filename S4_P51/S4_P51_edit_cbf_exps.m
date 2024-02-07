
% P51 uses snow_only_shuang branch to fix SCF problem, EDC is off
% P50 uses snow_only_shuang branch to fix SCF problem
% P40,P41 is the same as P25 but uses main branch DALEC1100 for debugging and testing
% P41 tested main on Oct 9
% P42 tested main on Jan 31
% P43 tested main on apr 27
nop = 'S4_P51'; % P21 first test of min max_value field in PEQs, after correcting JCR C balance,
% replaced 25 with meantemp, and added SKT as a met driver
                % P22 changed DALEC_OBSERVATION_OPERATORS.c:
                % if  (SOBS.value!=DEFAULT_DOUBLE_VAL){(1-D->M_PEQ_CUE)=D->M_PARS[O->CUE_PARAM];}
                % both P21 and P22 are wrong;
                % P23, and does not change anything in DALEC1100 (P.f_auto is an index)
%                 if  (SOBS.value!=DEFAULT_DOUBLE_VAL){(D->M_PEQ_CUE)=1-D->M_PARS[O->CUE_PARAM];}
%               % P24, soil moisture is still low (a bit better though),
%               somewhere in the model want sm stay dry, let's try if CH4
%               constrain can counteract and bring up the sm, by using 
%               S_fv[3 100] instead of [1,100]
PROJNAME = ['fluxnet_ch4_' nop];
pathname2=sprintf('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/input/');
ncdir = ['/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/' nop '/input/'];
if not(isfolder(ncdir))
    mkdir(ncdir)
end
          
% load coordinates info for all the fluxnet ch4 sites
coord=readtable('RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S2_LATLON2CBF_edited.csv');
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
    CBF_exp = CARDAMOM_READ_NC_CBF_FILE([pathname2 num2str(isiteindex) '_' char(isiteid) '.cbf.nc']);
    CBF_exp_noobs = CARDAMOM_READ_NC_CBF_FILE([pathname2 num2str(isiteindex) '_' char(isiteid) '_noobs.cbf.nc']);
% change model ID here if needed:
%     CBF_exp.ID.values = 1100;
%     CBF_exp_noobs.ID.values = 1100;
% change MCMC setup to shorten run time for snow screening branch only
      CBF_exp_noobs.MCMCID.nITERATIONS=1e5;
      CBF_exp_noobs.MCMCID.nSAMPLES=10;
%       CBF_exp_noobs.MCMCID.nSAMPLES_EDC_SEARCH=2e3;
      CBF_exp_noobs.EDC.values=0;
      CBF_exp.SCF.single_unc=1.2;
      CBF_exp.SCF.min_threshold=0.05;
% -- -- -- --  get CBF_exps ready and submit to gattaca  -- --  --  -- 
    [CBF_cell,string_CBF_cell] = edit_fluxnetexps_used_in_P23_P25(CBF_exp,CBF_exp_noobs); 
    
    nc_cell={};
        for icell=[6]
          fname = [num2str(isiteindex) '_' char(isiteid) char(string_CBF_cell(icell)) '.cbf.nc'];
          nc_cell{icell} = [ncdir fname];
          CARDAMOM_WRITE_NC_CBF_FILE((CBF_cell{icell}),[ncdir fname]);
        end
end

%% uncomment this line to submit to server with updated JPL server submission script
%  PARRFUN_CARDAMOM_SUBMIT_TO_GATTACA_lg_edge_shuangma(PROJNAME,ncdir);
   PARRFUN_CARDAMOM_SUBMIT_TO_GATTACA_shuangma(PROJNAME,ncdir);
%    PARRFUN_CARDAMOM_SUBMIT_TO_GATTACA(PROJNAME,ncdir);

% CBR_IVO = CARDAMOM_RUN_MDF('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S4_P42/input/58_US-Ivo_exp6.cbf.nc');
