
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOLD Index calculation
% using OC data to map the network cluster, (OC, nonlinear -> t2s&s0)
% set the BOLD signal @17ms
% #####################################
% by W. Yang @ Freiburg, 12, June. 2022

clear
addpath(genpath('/home/extern/yangw/Work/fMRIprogramm/NIfTI_20140122')); % Linux
addpath(genpath('/home/extern/yangw/Work/fMRIprogramm/tsv_read'));

te = [16,33.16,50.32];

%% parameter setting
dir_data = 'D:\Freiburg_server_data\MB_Experiment\s11\echo3\data4tedana_nosmooth\e1';
e1_name = 's11_3e1_smc.nii.gz';
e2_name = 's11_3e1_smc.nii.gz';
e3_name = 's11_3e1_smc.nii.gz';
mask_name = 's11_3e1_smc_mask.nii.gz';

dir_tedana = 'D:\Freiburg_server_data\MB_Experiment\s11\echo3\data4tedana_nosmooth\e1\pca_mdl_results_e1masked';

%%
cd(dir_data)

mask_temp = load_untouch_nii(mask_name);
mask = mask_temp.img;
[is,js,ks] = ind2sub(size(mask),find(mask > 0));

s_3e1 = load_untouch_nii(e1_name);
s_3e2 = load_untouch_nii(e2_name);
s_3e3 = load_untouch_nii(e3_name);

OC_data_temp = load_untouch_nii('smc_masked_oc.nii.gz');
OC_data = OC_data_temp.img;

%%%%%%%%%%%%%%%%%%%%%%% Enter the folder for data %%%%%%%%%%%%%%%%%%%%%%%%%
cd(dir_tedana)

IC_t_temp = tsvread('desc-ICA_mixing.tsv');
IC_t = IC_t_temp(2:end,:);

%%%%==== contribution of cth component to the vth voxel ====%%%%
% if 0
% B. this is the General Linear Model
% S_v = beta1*IC1+beta2*IC2+.....betan*ICn+Beta0;
% spatially normalized signal amplitude that characterizes the relative
% contribution of the cth component to the signal from the vth voxel

if ~exist('BOLD_Index_oc', 'dir')
    mkdir('BOLD_Index_oc')
end
cd('BOLD_Index_oc')

Nc = size(IC_t,2);
beta = zeros([size(mask),Nc]);
p_value = zeros([size(mask),Nc]);
t_value = zeros([size(mask),Nc]);

for nn = 1:1:length(is)
    s_v = squeeze((OC_data(is(nn),js(nn),ks(nn),:)));
    %         [b,dev] = glmfit(IC_t,s_v,'normal');
    [b,dev,stats] = glmfit(IC_t,s_v,'normal');
    beta(is(nn),js(nn),ks(nn),:) = b(2:end);
    p_value(is(nn),js(nn),ks(nn),:) = stats.p(2:end);
    t_value(is(nn),js(nn),ks(nn),:) = stats.t(2:end);
end

save('p_value.mat','p_value') % the p_value for estimate the significance
save('t_value.mat','t_value') % the p_value for estimate the significance
save('beta.mat','beta')

%%%%==== figure out whether signal is BOLD-like or non-BOLD-like ====%%%%
% first: get contributions of 3 echos

beta1 = zeros([size(mask),Nc]);
beta2 = zeros([size(mask),Nc]);
beta3 = zeros([size(mask),Nc]);

OC_data = single(s_3e1.img);
for nn = 1:1:length(is)
    s_v = squeeze((OC_data(is(nn),js(nn),ks(nn),:)));
    [b,dev] = glmfit(IC_t,s_v,'normal');
    beta1(is(nn),js(nn),ks(nn),:) = b(2:end);
end

OC_data = single(s_3e2.img);
for nn = 1:1:length(is)
    s_v = squeeze((OC_data(is(nn),js(nn),ks(nn),:)));
    [b,dev] = glmfit(IC_t,s_v,'normal');
    beta2(is(nn),js(nn),ks(nn),:) = b(2:end);
end

OC_data = single(s_3e3.img);
for nn = 1:1:length(is)
    s_v = squeeze((OC_data(is(nn),js(nn),ks(nn),:)));
    [b,dev] = glmfit(IC_t,s_v,'normal');
    beta3(is(nn),js(nn),ks(nn),:) = b(2:end);
end

save('beta1.mat','beta1')
save('beta2.mat','beta2')
save('beta3.mat','beta3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% F-test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('linear function fitting')
load beta.mat
load beta1.mat
load beta2.mat
load beta3.mat

Er2 = zeros([size(mask),Nc]);
Es0 = zeros([size(mask),Nc]);
res = zeros([size(mask),Nc]);

p1 = zeros([size(mask),Nc]);
p2 = zeros([size(mask),Nc]);

for ncc = 1:1:Nc

    delt1 = double(beta1(:,:,:,ncc)./mean(s_3e1.img,4));
    delt2 = double(beta2(:,:,:,ncc)./mean(s_3e2.img,4));
    delt3 = double(beta3(:,:,:,ncc)./mean(s_3e3.img,4));

    for nn = 1:1:length(is)

        delts = [delt1(is(nn),js(nn),ks(nn)),delt2(is(nn),js(nn),ks(nn)),delt3(is(nn),js(nn),ks(nn))];

        p =polyfit(te,delts,1);
        p1(is(nn),js(nn),ks(nn),ncc) = p(1);
        p2(is(nn),js(nn),ks(nn),ncc) = p(2);

        %%%%% Energy of BOLD and non_BOLD estimation %%%%%
        p1_c = p(1);
        p2_c = p(2);
        Er2(is(nn),js(nn),ks(nn),ncc) = abs(p1_c)*21;
        Es0(is(nn),js(nn),ks(nn),ncc) = abs(p2_c);
        res(is(nn),js(nn),ks(nn),ncc) = std(p1_c*te + p2_c - delts);

    end

end

save('p1.mat','p1')
save('p2.mat','p2')

save('Er2.mat','Er2')
save('Es0.mat','Es0')
save('res.mat','res')



