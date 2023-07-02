%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1.BOLD Index calculation
%   2.function data to T1 data registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc
cpwd = pwd;
cd(cpwd)

addpath(genpath('/home/extern/yangw/Work/fMRIprogramm/NIfTI_20140122'));
addpath(genpath('/home/extern/yangw/Work/fMRIprogramm/fmristat'));
addpath(genpath('/home/extern/yangw/Work/fMRIprogramm/spm12'));

%% parameter setting
TR=0.390;
te = [16,33.16,50.32];
Nt = 924;
k = 1;
kk = 11;

dir_tedana = 'D:\Freiburg_server_data\MB_Experiment\s11\echo3\data4tedana_nosmooth\e1\pca_mdl_results_e1masked';
dir_data = 'D:\Freiburg_server_data\MB_Experiment\s11\echo3\data4tedana_nosmooth\e1';
s0_name = 'smc_m_nonlinear_s0.nii.gz';
t2s_name = 'smc_m_nonlinear_t2s.nii.gz';
mask_name = 's11_3e1_smc_mask.nii.gz';

%%
%%%%%%%%%%%%%%%%% parameters adjust region %%%%%%%%%%%%%%%%%%%%%%%%%%%

h0 = figure;
% set(h0,'Position',[98  352  1417  463])
set(h0,'Position',[472 316 735 463])
TEs = 21;

% load mask
cd(dir_data);
mask_temp = load_untouch_nii(mask_name);
mask = mask_temp.img;
[is,js,ks] = ind2sub(size(mask),find(mask > 0));

ME_s0 = load_untouch_nii(s0_name);
ME_t2s = load_untouch_nii(t2s_name);
ME_data = ME_s0.img.*exp(-TEs./ME_t2s.img);
ME_data(isnan(ME_data)) = 0;
s_mean = mean(ME_data,4);

% create BOLD index image
cd([dir_tedana,'/BOLD_Index_oc'])
load beta.mat
load Er2.mat
load Es0.mat
load res.mat

[Nx,Ny,Nz,Nc] = size(beta);

z = beta;
z_org = beta;

z3d = zeros(size(z_org));
for mm = 1:1:Nc
    ic_temp = zeros(size(is));
    for nn = 1:1:length(is)
        ic_temp(nn) = z_org(is(nn),js(nn),ks(nn),mm);
    end
    ic_backup = ic_temp;
    [muHat,sigmaHat] = normfit(ic_temp);

    ic_temp (ic_temp > muHat+3*sigmaHat | ic_temp < muHat-3*sigmaHat ) = 0; % get rid of data out of 3*sigma to get the normal gaussian distrition
    iss = find(ic_temp ~= 0);
    ic_temp_new = zeros(size(iss));
    for nnn = 1:1:length(iss)
        ic_temp_new(nnn) = ic_temp(iss(nnn));
    end
    [muHat,sigmaHat] = normfit(ic_temp_new);
    ic_temp = ic_backup;
    ic_temp (ic_temp < muHat+1.96*sigmaHat & ic_temp > muHat-1.96*sigmaHat ) = 0; % z-test, p<0.05

    ic_spare = ic_temp;
    for nn = 1:1:length(is)
        z3d(is(nn),js(nn),ks(nn),mm) = ic_spare(nn);
    end
end

z3d = abs(z3d);
z3d(z3d ~= 0) = 1;

%%% 
s_mean_Nc = repmat(s_mean,1,1,1,Nc);

Er2_temp = reshape(z3d.*Er2.*s_mean_Nc,[Nx*Ny*Nz, Nc]);
Es0_temp = reshape(z3d.*Es0.*s_mean_Nc,[Nx*Ny*Nz, Nc]);
res_temp = reshape(z3d.*res.*s_mean_Nc,[Nx*Ny*Nz, Nc]);

vari_z = reshape(z3d,[64*64*20, Nc]);
Er2_ind = sum(Er2_temp)./sum(vari_z);
Es0_ind = sum(Es0_temp)./sum(vari_z);
res_ind = sum(res_temp)./sum(vari_z);

BOLD_Index_cluster = 1./(exp((Es0_ind-Er2_ind)./res_ind)+1);

%%% create the figure
plot(0:1:Nc-1,BOLD_Index_cluster,'Color',[0.8500 0.3250 0.0980]);
hold on
legend('BOLD Index')
legend boxoff
plot(0:1:Nc-1,BOLD_Index_cluster,'*','Color',[0.8500 0.3250 0.0980]*0.8)
plot([-1,Nc+1],[0.5,0.5],'Color',[1,1,1]*0.8)
title(['S',num2str(kk),'\ BOLD Index '])
box on

max_window = max(BOLD_Index_cluster)+(max(BOLD_Index_cluster)-min(BOLD_Index_cluster))*0.1;
min_window = min(BOLD_Index_cluster)-(max(BOLD_Index_cluster)-min(BOLD_Index_cluster))*0.1;
if max_window < 0.6
    max_window = 0.6;
end
axis([-0.9,Nc,min_window,max_window])

saveas(h0,['s',num2str(kk),'_BOLD_Index_k',num2str(k),'.jpg'])
% close all

save('BOLD_Index_cluster.mat','BOLD_Index_cluster')


