
function [CBF_cell,string_CBF_cell] = edit_fluxnetexps_used_in_P45(CBF_exp,CBF_exp_noobs)
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
    CBF_exp.PEQ_S_fv.unc=1.5;%$$
    
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

%%  save nc input files and submit to gattaca
      CBF_cell={CBF_exp1,CBF_exp2,CBF_exp3,CBF_exp4,CBF_exp5,CBF_exp6,CBF_exp7,CBF_exp8}; %$$
      string_CBF_cell=["_exp1","_exp2","_exp3","_exp4","_exp5","_exp6","_exp7","_exp8"];
end
      