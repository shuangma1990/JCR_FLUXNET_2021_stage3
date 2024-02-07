
function [CBF_cell,string_CBF_cell] = edit_fluxnetexps_used_in_P23_P25(CBF_exp,CBF_exp_noobs)
        % for fluxnet sites, obs include
%     obsnames = ["PEQ_CUE","NBE","CH4","ET","SCF","LAI","SOM","ABGB"];
    
    CBF_exp_noobs.PEQ_CUE = CBF_exp.PEQ_CUE; % now PEQ is everywhere
    
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
    CBF_exp3.PEQ_PAW_z=CBF_exp.PEQ_PAW_z;
%     CBF_exp3.Mean_LAI = CBF_exp.Mean_LAI;
    
    CBF_exp4.NBE = CBF_exp.NBE;
    CBF_exp4.CH4 = CBF_exp.CH4;
    CBF_exp4.ET = CBF_exp.ET;
    CBF_exp4.ABGB = CBF_exp.ABGB;
    CBF_exp4.PEQ_iniSOM=CBF_exp.PEQ_iniSOM;
    CBF_exp4.PEQ_PAW_z=CBF_exp.PEQ_PAW_z;
    
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
    CBF_exp7.PEQ_PAW_z=CBF_exp.PEQ_PAW_z;
%%  save nc input files and submit to gattaca
      CBF_cell={CBF_exp1,CBF_exp2,CBF_exp3,CBF_exp4,CBF_exp5,CBF_exp6,CBF_exp7}; %$$
      string_CBF_cell=["_exp1","_exp2","_exp3","_exp4","_exp5","_exp6","_exp7"];
end
      