% cbrdir = '/Users/shuangma/GATTACA_LOCAL_MIRROR/fluxnet_ch4_S4_P26/CBR_FILES/1_AT-Neu_exp5.cbf.cbr';
% cbfdir = '/Users/shuangma/GATTACA_LOCAL_MIRROR/fluxnet_ch4_S4_P26/CBF_FILES/1_AT-Neu_exp5.cbf.nc';
cbrdir = '/Users/shuangma/GATTACA_LOCAL_MIRROR/fluxnet_ch4_S4_P29/CBR_FILES/1_AT-Neu_exp5.cbf.cbr';
cbfdir = '/Users/shuangma/GATTACA_LOCAL_MIRROR/fluxnet_ch4_S4_P29/CBF_FILES/1_AT-Neu_exp5.cbf.nc';
CBF = CARDAMOM_READ_NC_CBF_FILE(cbfdir);
CBR = CARDAMOM_RUN_MODEL(cbfdir,cbrdir);
MD=CARDAMOM_MODEL_LIBRARY(CBF.ID.values,[],1);

CBR_end = CARDAMOM_RUN_MODEL(cbfdir,CBR.PARS(end,:));

CBR.PARS(end,MD.PARAMETER_IDs.i_PAW)

CBR.PARS(end,MD.PARAMETER_IDs.PAW_por)

CBR.PARS(end,MD.PARAMETER_IDs.PAW_z)

CBR.FLUXES(end,:,MD.FLUX_IDs.sm_PAW)

% sm = PAW/(1000*PAW_z*PAW_por);

i_sm_PAW = CBR.PARS(end,MD.PARAMETER_IDs.i_PAW)/(1000*CBR.PARS(end,MD.PARAMETER_IDs.PAW_z)*CBR.PARS(end,MD.PARAMETER_IDs.PAW_por))

endfirst_sm_PAW = CBR.POOLS(end,1,MD.POOL_IDs.H2O_PAW)/(1000*CBR.PARS(end,MD.PARAMETER_IDs.PAW_z)*CBR.PARS(end,MD.PARAMETER_IDs.PAW_por))
figure(1);
nfrow=2;
nfcol=4;
isub=1;
endpoint = size(CBR.FLUXES,2);
subplot(nfrow,nfcol,isub);isub=isub+1;  
plot(CBR.FLUXES(end,1:endpoint,MD.FLUX_IDs.sm_PAW)); title('SM PAW');
subplot(nfrow,nfcol,isub);isub=isub+1;  
plot(CBR.POOLS(end,1:endpoint,MD.POOL_IDs.H2O_PAW)); title('PAW');
subplot(nfrow,nfcol,isub);isub=isub+1;  
plot(CBR.FLUXES(end,1:endpoint,MD.FLUX_IDs.infil)); title('INFIL');
subplot(nfrow,nfcol,isub);isub=isub+1;  
plot(CBR.FLUXES(end,1:endpoint,MD.FLUX_IDs.q_paw)); title('DRAINAGE');
subplot(nfrow,nfcol,isub);isub=isub+1;  
plot(CBR.POOLS(end,1:endpoint,MD.POOL_IDs.H2O_PUW)); title('PUW');
subplot(nfrow,nfcol,isub);isub=isub+1;  
plot(CBR.FLUXES(end,1:endpoint,MD.FLUX_IDs.paw2puw)); title('PAW2PUW');
subplot(nfrow,nfcol,isub);isub=isub+1;  
scatter(CBR.FLUXES(end,:,MD.FLUX_IDs.sm_PAW),CBR.FLUXES(end,:,MD.FLUX_IDs.rh_ch4));
title('x:sm-PAW y:ch4');
subplot(nfrow,nfcol,isub);isub=isub+1;  
scatter(CBR.FLUXES(end,:,MD.FLUX_IDs.sm_PAW),CBR.FLUXES(end,:,MD.FLUX_IDs.rh_ch4));
title('x:sm-PAW y:ch4');xlim([0 1.2]);
% CBR = CARDAMOM_RUN_MODEL(cbfdir,cbrdir);
CBR_end = CARDAMOM_RUN_MODEL(cbfdir,CBR.PARS(end,:));
figure(2);
nfrow=2;
nfcol=4;
isub=1;
endpoint = 15;
subplot(nfrow,nfcol,isub);isub=isub+1;  
plot(CBR.FLUXES(end,:,MD.FLUX_IDs.sm_PAW)); title('sm_PAW');
subplot(nfrow,nfcol,isub);isub=isub+1;  
plot(CBR.FLUXES(end,1:endpoint,MD.FLUX_IDs.sm_PAW)); title('SM PAW');
subplot(nfrow,nfcol,isub);isub=isub+1;  
plot(CBR.POOLS(end,1:endpoint,MD.POOL_IDs.H2O_PAW)); title('PAW');
subplot(nfrow,nfcol,isub);isub=isub+1;  
plot(CBR.FLUXES(end,1:endpoint,MD.FLUX_IDs.infil)); title('INFIL');
subplot(nfrow,nfcol,isub);isub=isub+1;  
plot(CBR.FLUXES(end,1:endpoint,MD.FLUX_IDs.q_paw)); title('DRAINAGE');
subplot(nfrow,nfcol,isub);isub=isub+1;  
plot(CBR.POOLS(end,1:endpoint,MD.POOL_IDs.H2O_PUW)); title('PUW');
subplot(nfrow,nfcol,isub);isub=isub+1;  
plot(CBR.FLUXES(end,1:endpoint,MD.FLUX_IDs.paw2puw)); title('PAW2PUW');