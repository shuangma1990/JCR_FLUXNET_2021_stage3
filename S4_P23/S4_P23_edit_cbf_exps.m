% P8 corrected ABGB values (was zeros due to mistake in the load_05deg_obs script)
% changed min_threshold of CH4 from 1% to 10% (1% was likely a typo)
% and runs exp4 and exp5 on edge
nop = 'S4_P23'; % P21 first test of min max_value field in PEQs, after correcting JCR C balance,
% replaced 25 with meantemp, and added SKT as a met driver
                % P22 changed DALEC_OBSERVATION_OPERATORS.c:
                % if  (SOBS.value!=DEFAULT_DOUBLE_VAL){(1-D->M_PEQ_CUE)=D->M_PARS[O->CUE_PARAM];}
                % both P21 and P22 are wrong;
                % P23, and does not change anything in DALEC1100 (P.f_auto is an index)
%                 if  (SOBS.value!=DEFAULT_DOUBLE_VAL){(D->M_PEQ_CUE)=1-D->M_PARS[O->CUE_PARAM];}
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
% -- -- -- --  get CBF_exps ready and submit to gattaca  -- --  --  -- 
    [CBF_cell,string_CBF_cell] = edit_fluxnetexps(CBF_exp,CBF_exp_noobs); 
    
    nc_cell={};
%       for icell=1:length(CBF_cell)
        for icell=[2 3 4 5 6]
          fname = [num2str(isiteindex) '_' char(isiteid) char(string_CBF_cell(icell)) '.cbf.nc'];
          nc_cell{icell} = [ncdir fname];
          CARDAMOM_WRITE_NC_CBF_FILE((CBF_cell{icell}),[ncdir fname]);
%% uncomment this line to submit to server
%           pause(2)
%           PARRFUN_CARDAMOM_SUBMIT_TO_GATTACA_lg_edge_shuangma(PROJNAME,[nc_cell{icell}]);
%% uncomment this line to test on laptop
%           CBRtest=CARDAMOM_RUN_MDF([nc_cell{icell}]);
%             CBF_exp_noobs.EDC.values=1;
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'PEQ_Cefficiency');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'PEQ_CUE');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'PEQ_PAW_z');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'PEQ_iniSOM');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'ABGB');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'ABGBunc');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'NBE');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'NBEunc');  
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'LAI');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'LAIunc');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'CH4');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'CH4unc');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'ET');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'ETunc');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'GPP');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'GPPunc');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'EWT');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'EWTunc');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'Mean_Biomass');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'Mean_Fire');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'Mean_GPP');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'Mean_LAI');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'DOM');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'DOMunc');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'SCF');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'SCFunc');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'STRD');
%             CBF_exp_noobs = rmfield(CBF_exp_noobs,'SKT');
%             CARDAMOM_WRITE_NC_CBF_FILE(CBF_exp_noobs,'/Users/shuangma/RESEARCH/WORKFLOW/shuang_debug.nc')
%             CBRtest=CARDAMOM_RUN_MDF(CBF_exp_noobs);
         end
end

%% uncomment this line to submit to server with updated JPL server submission script
 PARRFUN_CARDAMOM_SUBMIT_TO_GATTACA_lg_edge_shuangma(PROJNAME,ncdir);
%   PARRFUN_CARDAMOM_SUBMIT_TO_GATTACA_shuangma(PROJNAME,ncdir);