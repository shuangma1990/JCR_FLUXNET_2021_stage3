% P40,P41 is the same as P25 but uses main branch DALEC1100 for debugging and testing
% P41 tested main on Oct 9
% P42 tested main on Jan 31
% P43 tested main on apr 27
% P44 tested main on branch main_202306, Sep8 2023
% P45 tested main on branch main_202306_addPEQ4above, exp 1 2 6 7 8, 8 is new, remote + pMCMC Sep19 2023
% P46 tested is same as P45, after correcting a mistake in the script DALEC_OBSERVATION_OPERATORS.c
  %  and I changed sfv unc from 1.5 to 1.2 cause I really want it become wetland
% P51 tested on main pulled on 2023Nov3, after Eren had SCF fixes
% P62 uses the main pulled Nov 2, has exp9-13, P61 has exp12678
% P63 has 14-15
  nop = 'S4_P63'; % P21 first test of min max_value field in PEQs, after correcting JCR C balance,
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
% for i=1:2
    irow = str2num(rowv(i));
    icol = str2num(colv(i));
    isiteindex = str2num(siteindex(i));
    isiteid = (siteid(i));
% read preprocessed full obs nc files
    CBF_exp = CARDAMOM_READ_NC_CBF_FILE([pathname2 num2str(isiteindex) '_' char(isiteid) '.cbf.nc']);
%     CBRtest=CARDAMOM_RUN_MDF(CBF_exp);
    CBF_exp_noobs = CARDAMOM_READ_NC_CBF_FILE([pathname2 num2str(isiteindex) '_' char(isiteid) '_noobs.cbf.nc']);
% change model ID here if needed:
%     CBF_exp.ID.values = 1100;
%     CBF_exp_noobs.ID.values = 1100;
% -- -- -- --  get CBF_exps ready and submit to gattaca  -- --  --  -- 
    [CBF_cell,string_CBF_cell] = edit_fluxnetexps_used_in_P46(CBF_exp,CBF_exp_noobs); 
    
    nc_cell={};
%       for icell=1:length(CBF_cell)
%         for icell=[2 3 4 5 6]
%         for icell=[9 10 11 12 13 14 15]
        for icell=[14 15]
          fname = [num2str(isiteindex) '_' char(isiteid) char(string_CBF_cell(icell)) '.cbf.nc'];
          nc_cell{icell} = [ncdir fname];
          CARDAMOM_WRITE_NC_CBF_FILE((CBF_cell{icell}),[ncdir fname]);
%% uncomment this line to submit to server
%           pause(2)
%           PARRFUN_CARDAMOM_SUBMIT_TO_GATTACA_lg_edge_shuangma(PROJNAME,[nc_cell{icell}]);
%% uncomment this line to test on laptop
% CBR_IVO = CARDAMOM_RUN_MDF('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S4_P42/input/58_US-Ivo_exp6.cbf.nc');
% filename = '/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S4_P42/local_test.mat';
% save(filename,'CBR_IVO');

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
%  PARRFUN_CARDAMOM_SUBMIT_TO_GATTACA_lg_edge_shuangma(PROJNAME,ncdir);
 PARRFUN_CARDAMOM_SUBMIT_TO_GATTACA2(PROJNAME,ncdir,'72:00:00');
%    PARRFUN_CARDAMOM_SUBMIT_TO_LONESTAR6(PROJNAME,ncdir,48:00:00);
%    PARRFUN_CARDAMOM_SUBMIT_TO_GATTACA(PROJNAME,ncdir);

% CBR_IVO = CARDAMOM_RUN_MDF('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S4_P42/input/58_US-Ivo_exp6.cbf.nc');
