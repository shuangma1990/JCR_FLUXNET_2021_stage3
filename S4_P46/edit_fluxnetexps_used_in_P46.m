
function [CBF_cell,string_CBF_cell] = edit_fluxnetexps_used_in_P46(CBF_exp,CBF_exp_noobs)
%     P45 removes PAW_z in line with the new model structure, since PAW_z has been removed 
%     for fluxnet sites, obs include
%     obsnames = ["PEQ_CUE","NBE","CH4","ET","SCF","LAI","SOM","ABGB"];
    
    CBF_exp_noobs.PEQ_CUE = CBF_exp.PEQ_CUE; % now PEQ_CUE is everywhere
        
    % ------------------  pMCMC I added them here istead of S3_edit_cbf,
    % shouln't affect previous projects that call foor edit_fluxnetCBFunc.m
    % because the values are constant, may change to be variable depending
    % on inundation situation
    CBF_exp.PEQ_r_ch4.values = 0.9; % values are from the posteriors of 4yr sites exp2
    CBF_exp.PEQ_r_ch4.opt_unc_type=1; % $$$          
    CBF_exp.PEQ_r_ch4.unc=1.5;%$$
    
    CBF_exp.PEQ_S_fv.values = 3.55;
    CBF_exp.PEQ_S_fv.opt_unc_type=1; % $$$          
    CBF_exp.PEQ_S_fv.unc=1.2;%$$
    
    CBF_exp.PEQ_rhch4_rhco2.values = exp(-5);
    CBF_exp.PEQ_rhch4_rhco2.opt_unc_type=1; % $$$          
    CBF_exp.PEQ_rhch4_rhco2.unc=1.2;%$$
    % ------------------  pMCMC
    
    CBF_exp1 = CBF_exp_noobs;
    CBF_exp2 = CBF_exp_noobs;
    CBF_exp3 = CBF_exp_noobs;
    CBF_exp4 = CBF_exp_noobs;
    CBF_exp5 = CBF_exp_noobs;
    CBF_exp6 = CBF_exp_noobs;
    CBF_exp7 = CBF_exp_noobs;
    CBF_exp8 = CBF_exp_noobs;
    CBF_exp9 = CBF_exp_noobs;
    CBF_exp10 = CBF_exp_noobs;
    CBF_exp11 = CBF_exp_noobs;
    CBF_exp12 = CBF_exp_noobs;
    CBF_exp13 = CBF_exp_noobs;
    CBF_exp14 = CBF_exp_noobs;
    CBF_exp15 = CBF_exp_noobs;
    CBF_exp16 = CBF_exp_noobs;
    CBF_exp17 = CBF_exp_noobs;
    CBF_exp18 = CBF_exp_noobs;
    CBF_exp19 = CBF_exp_noobs;
    CBF_exp20 = CBF_exp_noobs;

%     CBF_exp1.ABGB = CBF_exp.ABGB;
    
    CBF_exp2.NBE = CBF_exp.NBE;
    CBF_exp2.CH4 = CBF_exp.CH4;
    CBF_exp2.ET = CBF_exp.ET;
%     CBF_exp2.SCF = CBF_exp.SCF; % used this in P7, proved that single SCF data stream doesn't have problems
    
%     CBF_exp3.NBE = CBF_exp.NBE;
%     CBF_exp3.CH4 = CBF_exp.CH4;
%     CBF_exp3.ET = CBF_exp.ET;

    CBF_exp3 = CBF_exp2;
%     CBF_exp3.PEQ_PAW_z=CBF_exp.PEQ_PAW_z;   % exp3 is now the same as exp2 since PAW_z is obsolete
%     CBF_exp3.Mean_LAI = CBF_exp.Mean_LAI;
    
    CBF_exp4.NBE = CBF_exp.NBE;
    CBF_exp4.CH4 = CBF_exp.CH4;
    CBF_exp4.ET = CBF_exp.ET;
    CBF_exp4.ABGB = CBF_exp.ABGB;
    CBF_exp4.PEQ_iniSOM=CBF_exp.PEQ_iniSOM;
%     CBF_exp4.PEQ_PAW_z=CBF_exp.PEQ_PAW_z;
    
    CBF_exp5=CBF_exp4;
    CBF_exp5.SCF = CBF_exp.SCF;
% both in situ and RS
    CBF_exp6=CBF_exp5;
%     CBF_exp6.PEQ_iniSOM=CBF_exp.PEQ_iniSOM;
%     CBF_exp6.DOM=CBF_exp.DOM;
    CBF_exp6.LAI=CBF_exp.LAI;
% RS only
    % all remote sensing data (P9 but LAI)
    CBF_exp7.ABGB = CBF_exp.ABGB;
    CBF_exp7.SCF = CBF_exp.SCF;
    CBF_exp7.LAI=CBF_exp.LAI;
    CBF_exp7.PEQ_iniSOM=CBF_exp.PEQ_iniSOM;
%     CBF_exp7.PEQ_PAW_z=CBF_exp.PEQ_PAW_z;
    
    CBF_exp8 = CBF_exp7;
    CBF_exp8.PEQ_r_ch4 = CBF_exp.PEQ_r_ch4;
    CBF_exp8.PEQ_S_fv = CBF_exp.PEQ_S_fv;
    CBF_exp8.PEQ_rhch4_rhco2 = CBF_exp.PEQ_rhch4_rhco2;
%   exp 9 is same as exp8 but added in situ GPP
    CBF_exp9 = CBF_exp8;
    CBF_exp9.GPP = CBF_exp.GPP;
%   exp 10 is same as exp8 but added in situ GPP and NBE
    CBF_exp10 = CBF_exp8;
    CBF_exp10.GPP = CBF_exp.GPP;
    CBF_exp10.NBE = CBF_exp.NBE;
%   exp 11 is same as exp6 but took CH4 out
    CBF_exp11 = CBF_exp6;
    CBF_exp11.CH4=CBF_exp_noobs.CH4;
%   exp 12 is same as exp6 but use GPP instead of NBE
    CBF_exp12 = CBF_exp6;
    CBF_exp12.NBE=CBF_exp_noobs.NBE;    
    CBF_exp12.GPP=CBF_exp.GPP;
%   exp 13 is same as exp6 but use GPP and NBE
    CBF_exp13 = CBF_exp6;
    CBF_exp13.GPP=CBF_exp.GPP;   
%   exp 14 has all RS, PEQ, in situ except CH4
    CBF_exp14 = CBF_exp10;
    CBF_exp14.ET = CBF_exp.ET;
%   exp 15 exp6 without CH4
    CBF_exp15 = CBF_exp13;
    CBF_exp15.CH4 = CBF_exp_noobs.CH4;   
%   exp 16 exp6 but uses mean CH4, w NBE ET    
    CBF_exp16 = CBF_exp6;
    % same data but uses mean where data is available, here are the
    % changes:
    CBF_exp16.CH4.opt_filter=1; % 1: mean where data is available
	CBF_exp16.CH4.opt_normalization=0;% none
	CBF_exp16.CH4.opt_unc_type=1;
    CBF_exp16.CH4.single_mean_unc=1.5; % only used with optfilter=1
    CBF_exp16.CH4.single_unc=-9999;
%     CBF_exp16.CH4.single_unc=1.5; % only used with optfilter=0

%   exp 17 exp6 but uses mean CH4, wo NBE ET    
    CBF_exp17 = CBF_exp16;
    % remove NBE and ET
    fields = {'NBE','ET'};
    CBF_exp17 = rmfield(CBF_exp17,fields);

%   exp 18 exp8 but uses mean CH4, wo NBE ET  
    CBF_exp18 = CBF_exp8;
    % exp8+EWT+meanCH4,without NBE and ET:
    CBF_exp18.EWT = CBF_exp.EWT;
    CBF_exp18.CH4 = CBF_exp16.CH4;
    
%   exp 19 RS+EWT+GPP+NBE+ET+PEQ
    CBF_exp19 = CBF_exp14;
    CBF_exp19.EWT = CBF_exp.EWT;
    
%   ep20 RS+EWT+GPP+ET+PEQ
    CBF_exp20 = CBF_exp19;
    fields = {'NBE'};
    CBF_exp20 = rmfield(CBF_exp20,fields);
%%  save nc input files and submit to gattaca
      CBF_cell={CBF_exp1,CBF_exp2,CBF_exp3,CBF_exp4,CBF_exp5,CBF_exp6,CBF_exp7,CBF_exp8,...
          CBF_exp9,CBF_exp10,CBF_exp11,CBF_exp12,CBF_exp13,CBF_exp14,CBF_exp15,CBF_exp16,CBF_exp17,CBF_exp18,CBF_exp19,CBF_exp20}; %$$
      string_CBF_cell=["_exp1","_exp2","_exp3","_exp4","_exp5","_exp6","_exp7","_exp8",...
          "CBF_exp9","CBF_exp10","CBF_exp11","CBF_exp12","CBF_exp13","CBF_exp14","CBF_exp15","CBF_exp16","CBF_exp17","CBF_exp18","CBF_exp19","CBF_exp20"];
end
      