%
% compare the s0,t2s calculation from tedana and nonlinearfitting
% conclusion: the s0 and t2s in tedana are from linear fitting
%

clear
addpath(genpath('/home/extern/yangw/Work/fMRIprogramm/NIfTI_20140122')); % Linux

%% parameter setting part
%%% basic info echo3 resting state %%%
TR = 390;
te = [16.00, 33.16, 50.32]';

dir_data = 'D:\Freiburg_server_data\MB_Experiment\s11\echo3\data4tedana_nosmooth\e1';
e1_name = 's11_3e1_smc.nii.gz';
e2_name = 's11_3e1_smc.nii.gz';
e3_name = 's11_3e1_smc.nii.gz';
mask_name = 's11_3e1_smc_mask.nii.gz';

%%
cd(dir_data)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_temp = load_untouch_nii(mask_name);
mask = mask_temp.img;

s_3e1 = load_untouch_nii(e1_name);
s_3e2 = load_untouch_nii(e2_name);
s_3e3 = load_untouch_nii(e3_name);
data_all_temp = zeros([size(s_3e1.img),3]);
data_all_temp(:,:,:,:,1) = s_3e1.img;
data_all_temp(:,:,:,:,2) = s_3e2.img;
data_all_temp(:,:,:,:,3) = s_3e3.img;

data_all = mean(data_all_temp,4);

s0 = zeros(64,64,20);
t2s = zeros(64,64,20);
[is,js,ks] = ind2sub(size(mask),find(mask>0));

%% Here nonlinear fitting

s0_0 = 1000;
t2s_0 = 40;

for nn = 1:1:length(is)
    tseries = squeeze(data_all(is(nn),js(nn),ks(nn),:));

    beta0 = [s0_0,-1/t2s_0];
    beta = nlinfit(te,tseries,@(b,te)(b(1).*exp(b(2).*te)),beta0);

    s0(is(nn),js(nn),ks(nn)) = beta(1);
    t2s(is(nn),js(nn),ks(nn)) = -1/beta(2);
end

mask_temp.img = s0;
file=strcat('smc_masked_nonlinear_s0','.nii.gz');
save_untouch_nii(mask_temp,file);

mask_temp.img = t2s;
file=strcat('smc_masked_nonlinear_t2s','.nii.gz');
save_untouch_nii(mask_temp,file);

%% optimal combination
w1 = te(1)*exp(-te(1)./t2s)./(te(1)*exp(-te(1)./t2s)+te(2)*exp(-te(2)./t2s)+te(3)*exp(-te(3)./t2s));
w2 = te(2)*exp(-te(2)./t2s)./(te(1)*exp(-te(1)./t2s)+te(2)*exp(-te(2)./t2s)+te(3)*exp(-te(3)./t2s));
w3 = te(3)*exp(-te(3)./t2s)./(te(1)*exp(-te(1)./t2s)+te(2)*exp(-te(2)./t2s)+te(3)*exp(-te(3)./t2s));

w1(isnan(w1)) = 0;    w2(isnan(w2)) = 0;    w3(isnan(w3)) = 0;

data_temp = zeros(size(s_3e1.img));
Nt = length(s_3e1.img);
for mm = 1:1:Nt
    data_temp(:,:,:,mm) = w1.*s_3e1.img(:,:,:,mm) + w2.*s_3e2.img(:,:,:,mm) + w3.*s_3e3.img(:,:,:,mm);
end
s_3e1.img = data_temp;

save_untouch_nii(s_3e1,['smc_masked_oc','.nii.gz']);

