
clear
close all
clc
cpwd = pwd;

addpath(genpath('/home/extern/yangw/Work/fMRIprogramm/NIfTI_20140122'));
addpath(genpath('/home/extern/yangw/Work/fMRIprogramm/fmristat'));
addpath(genpath('/home/extern/yangw/Work/fMRIprogramm/spm12'));
cd('/raid/home/extern/yangw/Work/fMRIprogramm/colormap_fsl_render3')
load rend3.mat % fsl rend3 colorbar

%%%%%%%%%%%%%%%%%%%%%% parameters adjust region %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameter setting
te = [14.0, 40.36, 66.72];
TR=0.390;
k = 4; % for BOLD Index display reason

dir_tedana = 'D:\Freiburg_server_data\MB_Experiment\s11\echo3\data4tedana_nosmooth\e1\pca_mdl_results_e1masked';
dir_data = 'D:\Freiburg_server_data\MB_Experiment\s11\echo3\data4tedana_nosmooth\e1';
e1_name = 's11_3e1_smc.nii.gz';
mask_name = 's11_3e1_smc_mask.nii.gz';

%%

h0 = figure;
% get(h0,'Position')
set(h0,'Position',[371  350  668  461])

kk = 11;

kk
% load mask
cd(dir_data);
mask_temp = load_untouch_nii(mask_name);
mask = mask_temp.img;
[is,js,ks] = ind2sub(size(mask),find(mask > 0));

% load background
background_temp = load_untouch_nii(e1_name);
background = background_temp.img(:,:,:,100);
background = single(background);

% load activation map
cd([dir_tedana,'\BOLD_Index_oc'])
load beta.mat
load BOLD_Index_cluster.mat
Nc = size(beta,4);

z = beta;
z_org = beta;

z3d = zeros(size(z_org));
mu_sigma = zeros(Nc,2);
for mm = 1:1:Nc
    ic_temp = zeros(size(is));
    for nn = 1:1:length(is)
        ic_temp(nn) = z_org(is(nn),js(nn),ks(nn),mm);
    end
    ic_backup = ic_temp;
    [muHat,sigmaHat] = normfit(ic_temp);

    ic_temp (ic_temp > muHat+5*sigmaHat | ic_temp < muHat-5*sigmaHat ) = 0; % select this !
    iss = find(ic_temp ~= 0);
    ic_temp_new = zeros(size(iss));
    for nnn = 1:1:length(iss)
        ic_temp_new(nnn) = ic_temp(iss(nnn));
    end
    [muHat,sigmaHat] = normfit(ic_temp_new);
    ic_temp = ic_backup;
    ic_temp (ic_temp < muHat+1.96*sigmaHat & ic_temp > muHat-1.96*sigmaHat ) = 0; % select this !

    mu_sigma(mm,:) = [muHat,sigmaHat]';

    ic_spare = ic_temp;
    for nn = 1:1:length(is)
        z3d(is(nn),js(nn),ks(nn),mm) = ic_spare(nn);
    end
end

z3d = abs(z3d);
z3d(z3d ~= 0) = 1;
mask_act = z3d;

% calculate the voxel based BOLD_Index
load Er2.mat
load Es0.mat
load res.mat

Engergy_var = (Es0 - Er2)./(res*k);
BOLD_Index = 1./(exp(Engergy_var)+1);
BOLD_Index(isnan(BOLD_Index)) = 0;

cd(dir_tedana)
if ~exist('BOLD_Index_oc_fig', 'dir')
    mkdir('BOLD_Index_oc_fig')
end
cd('BOLD_Index_oc_fig')

for ind_Nc = 1:1:Nc

    % ind_Nc = 8;
    yindex3d = BOLD_Index(:,:,:,ind_Nc);

    % get the structure from 3D to 2D
    str3d = background;
    str2d = zeros([64*4,64*5]);
    for ii = 1:1:4
        for jj = 1:1:5
            A = zeros(4,5);
            A(ii,jj) = 1;
            B = rot90(str3d(:,:,21-((ii-1)*5+jj)));
            str2d = str2d + kron(A,B);
        end
    end

    beta3d = beta(:,:,:,ind_Nc).*mask_act(:,:,:,ind_Nc);
    beta2d = zeros([64*4,64*5]);
    for ii = 1:1:4
        for jj = 1:1:5
            A = zeros(4,5);
            A(ii,jj) = 1;
            B = rot90(beta3d(:,:,21-((ii-1)*5+jj)));
            beta2d = beta2d + kron(A,B);
        end
    end

    % get the index from 3D to 2D
    yindex3dm = yindex3d.*mask_act(:,:,:,ind_Nc);
    yindex3dm(yindex3dm==0) = 0.5;
    yindex2d = zeros([64*4,64*5]);
    for ii = 1:1:4
        for jj = 1:1:5
            A = zeros(4,5);
            A(ii,jj) = 1;
            B = rot90(yindex3dm(:,:,21-((ii-1)*5+jj)));
            yindex2d = yindex2d + kron(A,B);
        end
    end

    img = zeros([size(str2d),3]);
    c_str = gray(256);
    mind = 0 - 1e-4;
    maxd = 800;
    str2d(str2d>maxd) = maxd;
    maxd = max(str2d(:)) + 1e-4;

    [M,N] = size(str2d);
    for ii = 1:1:M
        for jj = 1:1:N
            img(ii,jj,:) = c_str( ceil((str2d(ii,jj)-mind)/(maxd-mind)*size(c_str,1)),:);
        end
    end

    img_beta = img;

    act1 = yindex2d;
    coloru = autumn(256);
    colord = winter(256);
    colord = flipud(colord);
    c_act1 = [colord;0,0,0;coloru];
    mind = 0;
    maxd = 1;
    act1(act1 > maxd) = maxd;
    act1(act1 < mind) = mind;
    mind = mind - 1e-4;
    maxd = maxd + 1e-4;
    [is,js] = ind2sub(size(str2d),find( yindex2d<0.5 | yindex2d>0.5 ));
    for nn = 1:1:length(is)
        img(is(nn),js(nn),:) = c_act1( ceil((act1(is(nn),js(nn))-mind)/(maxd-mind)*size(c_act1,1)),:);
    end

    c_beta = rend3;
    act2 = beta2d;
    mind = mu_sigma(ind_Nc,1)-10*mu_sigma(ind_Nc,2);
    maxd = mu_sigma(ind_Nc,1)+10*mu_sigma(ind_Nc,2);
    act2(act2 > maxd) = maxd;
    act2(act2 < mind) = mind;
    mind = mind - 1e-4;
    maxd = maxd + 1e-4;
    [is,js] = ind2sub(size(str2d),find( beta2d < mu_sigma(ind_Nc,1)-mu_sigma(ind_Nc,2) |  beta2d >mu_sigma(ind_Nc,1)+mu_sigma(ind_Nc,2) ));
    for nn = 1:1:length(is)
        img_beta(is(nn),js(nn),:) = c_beta( ceil((act2(is(nn),js(nn))-mind)/(maxd-mind)*size(c_beta,1)),:);
    end

    %%%%%%%%%%%%%%%%%% figure %%%%%%%%%%%%%%%%%%%
    image(img)
    axis off

    colorbar
    colormap(c_act1)
    caxis([0, 1]);

    title(['s',num2str(kk),' tedana comp',num2str(ind_Nc-1),', BOLD Index'])
    text(100,270, ['BOLD Index mapping (',num2str(BOLD_Index_cluster(ind_Nc)),')'])
    saveas(h0,['s',num2str(kk),'_tedana_comp',num2str(ind_Nc-1),'_BOLD_Index_k',num2str(k),'.png'])
    pause(0.5)
    clf

    %%%%%%%%%%%%%%%%%%%%% image beta map  %%%%%%%%%%%%%%%%%%%%%
    image(img_beta)
    axis off

    colorbar
    colormap(c_beta)
    caxis([mind, maxd]);

    title(['s',num2str(kk),' tedana comp',num2str(ind_Nc-1)])
    text(100,270, [num2str(length(is)),' voxels involved in the component'])
    saveas(h0,['s',num2str(kk),'_tedana_comp',num2str(ind_Nc-1),'.png'])
    pause(0.5)
    clf

end



