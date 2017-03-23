%% Example for real data
%% First you need to add the pathes that contains all the needed functions

clear all
path = '/Users/hao/Dropbox/DMRI_code/';

addpath(path);
addpath(strcat(path,'real_data'));
addpath(strcat(path,'ADMM'));
addpath(strcat(path,'construction_functions'));
addpath(strcat(path,'toolbox_wavelet_meshes'));
addpath(strcat(path,'toolbox_wavelet_meshes/toolbox/'));
addpath(strcat(path,'MEALPix/'));
addpath(strcat(path,'help_functions/'));
addpath(strcat(path,'NIfTI/'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load .nii data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load nii data (raw)
bvec = load('027_S_2245_2011-06-06_16_21_14.0_S110933.rotated.bvec'); %027_S_2245_2011-06-06_16_21_14.0_S110933.rotated.bvec
bval = load('027-S-2245_060611_S110933.bval');
nii_eddy = load_nii('027-S-2245_060611_S110933_data.nii.gz');

%% load fsl processed data
nii_FA = load_nii('dti_FA.nii.gz');
nii_MO = load_nii('dti_MO.nii.gz');
nii_MD = load_nii('dti_MD.nii.gz');
nii_S0 = load_nii('dti_S0.nii.gz');

nii_V1 = load_nii('dti_V1.nii.gz');
nii_V2 = load_nii('dti_V2.nii.gz');
nii_V3 = load_nii('dti_V3.nii.gz');

nii_L1 = load_nii('dti_L1.nii.gz');
nii_L2 = load_nii('dti_L2.nii.gz');
nii_L3 = load_nii('dti_L3.nii.gz');


%% transfer to matrix
img_all = nii_eddy.img;
size(img_all)

img_FA_fsl_all = nii_FA.img;
img_MD_fsl_all = nii_MD.img;
img_S0_fsl_all = nii_S0.img;

img_V1_fsl_all = nii_V1.img;
img_V2_fsl_all = nii_V2.img;
img_V3_fsl_all = nii_V3.img;

img_L1_fsl_all = nii_L1.img;
img_L2_fsl_all = nii_L2.img;
img_L3_fsl_all = nii_L3.img;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% real data analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% real data range
x_range = 108:123; %120:124;  %100:155    101:160
y_range = 124:139; %137:141;  %87:148     91:180;   
z_range = 37:42; %38:40;    %21:47      21:50;

img_temp = img_all(x_range,y_range,z_range,:);

img_S0_fsl_temp = img_S0_fsl_all(x_range,y_range,z_range);
img_MD_fsl_temp = img_MD_fsl_all(x_range,y_range,z_range);
img_FA_fsl_temp = img_FA_fsl_all(x_range,y_range,z_range);
img_V1_fsl_temp = img_V1_fsl_all(x_range,y_range,z_range,:);

[k1, k2, k3, k4] = size(img_temp);


%% MLE estimation of S0 and sigma in each voxel

Est_sigma_all = zeros(k1,k2,k3);
Est_S0_all = zeros(k1,k2,k3);
% Est_sigma_all = [];
% Est_S0_all = [];
for i=1:k1
    for j=1:k2
        for k=1:k3
            
            y = squeeze(img_temp(i,j,k,1:5));
            if(min(y)>0&&max(abs(y))~=Inf)
                options = optimoptions('fminunc','Algorithm','quasi-newton'); %,'TolFun',1e-20);
                options.Display = 'off';
                x0 = double([var(y),mean(y)]);
%                 x00 = double([3600 4000]);
                f = @(x) parameterfun(x,y);
                [x, fval, exitflag,grad] = fminunc(f,x0,options);
    %             [x, fval, exitflag, output] = fminunc(f,x0,options);
%                 Est_sigma_all(i,j,k) = sqrt(x(1));
%                 Est_S0_all(i,j,k) = x(2);
                Est_sigma_all(i,j,k) = sqrt(x(1));
                Est_S0_all(i,j,k) = (x(2));
            end
        end
    end
   	display(i);
end

Est_sigma_var_all = zeros(k1,k2,k3);
Est_S0_mean_all = zeros(k1,k2,k3);
for i=1:k1
    for j=1:k2
        for k=1:k3
            y = squeeze(img_temp(i,j,k,1:5));
            Est_sigma_var_all(i,j,k) = sqrt(var(y));
            Est_S0_mean_all(i,j,k) = mean(y);
        end
    end
end


figure
subplot(4,2,1)  % histogram of MLE est. sigma
hist(reshape(Est_sigma_all,1,k1*k2*k3))
title('Sigma')
subplot(4,2,2)  % histogram of MLE est. S0
hist(reshape(Est_S0_all,1,k1*k2*k3))
title('S0')
subplot(4,2,3)  % MLE sigma vs. Variance est. sigma
scatter(reshape(Est_sigma_all,1,k1*k2*k3),reshape(Est_sigma_var_all,1,k1*k2*k3))
title('var vs. Sigma')
subplot(4,2,4)  % MLE S0 vs. Mean est. S0
scatter(reshape(Est_S0_all,1,k1*k2*k3),reshape(Est_S0_mean_all,1,k1*k2*k3))
title('mean vs S0')
subplot(4,2,5)  % MLE S0 vs. MLE sigma
scatter(reshape(Est_S0_all,1,k1*k2*k3), reshape(Est_sigma_all,1,k1*k2*k3))
title('Sigma vs. S0')
subplot(4,2,6)  % fsl MD vs MLE S0
scatter( reshape(img_MD_fsl_temp,1,k1*k2*k3),reshape(Est_S0_all,1,k1*k2*k3))
title('S0 vs. MD')
subplot(4,2,7)  % fsl MD vs MLE sigma
scatter( reshape(img_MD_fsl_temp,1,k1*k2*k3),reshape(Est_sigma_all,1,k1*k2*k3))
title('Sigma vs. MD')
subplot(4,2,8) % SNR as ratio of MLE S0 over MLE sigma
hist(reshape(Est_S0_all./Est_sigma_all,1,k1*k2*k3),20)
title('SNR')

%% SNR statistics
% median(reshape(Est_S0_all,1,k1*k2*k3))/median(reshape(Est_sigma_all,1,k1*k2*k3))
Est_S0_all_vec = reshape(Est_S0_all,1,k1*k2*k3);
Est_sigma_all_vec = reshape(Est_sigma_all,1,k1*k2*k3);

SNR = Est_S0_all./Est_sigma_all;
SNR_vec = reshape(SNR,1,k1*k2*k3);

% SNR_vec((isnan(SNR_vec))) = 0; 

min(SNR_vec(SNR_vec>0))
max(SNR_vec(SNR_vec>0))
mean(SNR_vec(SNR_vec>0))
median(SNR_vec(SNR_vec>0))

quantile(SNR_vec,0.05)
quantile(SNR_vec,0.95)
figure
hist(SNR_vec,50)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Single tensor model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = k1*k2*k3;
n = 41;
S_min = min(min(img_temp(img_temp>0)));
S_max = max(max(img_temp(img_temp<Inf)));

S_all_temp = zeros(41,N);
for i = 1:41
%     S_all_temp(i,:) = reshape(max (S_min,img_temp(:,:,:,i+5)),1,N);
%     S_all_temp(i,:) = reshape(min(S_max,img_temp(:,:,:,i+5)),1,N);
    S_all_temp(i,:) = reshape(img_temp(:,:,:,i+5),1,N);
    S_all_temp(i,:) = S_all_temp(i,:)./reshape(Est_S0_all,1,N);  %% use S0_fsl
end

index_negative = [];
index_infinity = [];
for i=1:N
    if(min(S_all_temp(:,i))<0)
        index_negative = [index_negative i];
    end
    if(sum(abs(S_all_temp(:,i)==Inf))>0)
        index_infinity = [index_infinity i];
    end
end
index_temp = setdiff(1:N,union(index_negative,index_infinity));




X = zeros(n,6);
X(:,1) = bvec(1,6:46).^2;
X(:,2) = bvec(2,6:46).^2;
X(:,3) = bvec(3,6:46).^2;
X(:,4) = 2*bvec(1,6:46).*bvec(2,6:46);
X(:,5) = 2*bvec(1,6:46).*bvec(3,6:46);
X(:,6) = 2*bvec(2,6:46).*bvec(3,6:46);

% N = length(index_temp);

FA_temp = zeros(1,N);
eval_temp = zeros(3,N);
evec_temp = zeros(3,3,N);

%% single tensor model
for i=1:N   %length(index_temp)
    %i = index_temp(j);
    l_S = log(S_all_temp(:,i));
    D_est_temp = -inv((X'*X))*X'*l_S;
    [D_nl, iter, DWI_est] = LM_dti(X,S_all_temp(:,i),D_est_temp,1e-10);
    D_matrix_nl = [D_nl(1) D_nl(4) D_nl(5);D_nl(4) D_nl(2) D_nl(6);D_nl(5) D_nl(6) D_nl(3)];
    if(sum(isnan(D_matrix_nl))==0)
      [egvec_nl, egval_nl] = eig(D_matrix_nl);
    end
    eval_temp(:,i) = diag(egval_nl);
    evec_temp(:,:,i) = egvec_nl;
    FA_temp(i) = sqrt(1/2)*sqrt(((eval_temp(1,i)-eval_temp(2,i))^2+(eval_temp(1,i)-eval_temp(3,i))^2+(eval_temp(3,i)-eval_temp(2,i))^2))/sqrt((eval_temp(1,i)^2+eval_temp(2,i)^2+eval_temp(3,i)^2)); 
end


MD_temp = sum(eval_temp,1)./3;
eval_ratio23_temp = eval_temp(2,:)./eval_temp(1,:);
eval_ratio_temp = eval_temp(3,:)*2./(eval_temp(1,:)+eval_temp(2,:));
index_ttemp = find(FA_temp<1&MD_temp>=0);
data_FA = reshape(FA_temp,k1,k2,k3);
data_MD = reshape(MD_temp,k1,k2,k3);

FA_show = FA_temp(index_ttemp);
MD_show = MD_temp(index_ttemp);
eval_ratio_show = eval_ratio_temp(index_ttemp);
eval3_show = eval_temp(3,index_ttemp);
eval_show = eval_temp(:,index_ttemp);
eval_ratio23_show = eval_ratio23_temp(index_ttemp);
SNR_show = reshape(Est_S0_all./Est_sigma_all,1,k1*k2*k3);
median(SNR_show(~isnan(SNR_show)))
mean(SNR_show(~isnan(SNR_show)))

idx_response = find(FA_show<1&FA_show>0.8&eval_show(1,:)>0&eval_show(2,:)>0&eval_show(3,:)>0&eval_ratio23_show<1.5);
size(idx_response)

temp = sort(eval_ratio_show(idx_response),'descend');
ratio_response = temp(floor(size(temp,2)/2));

id_temp = find(eval_ratio_show==ratio_response);
b_factor = median(eval_show(3,id_temp));


figure % compare single tensor model estimation vs fsl estimation
subplot(2,1,1)
scatter(FA_show, reshape(img_FA_fsl_temp,1,numel(img_FA_fsl_temp)))
subplot(2,1,2)
scatter(MD_show, reshape(img_MD_fsl_temp,1,numel(img_MD_fsl_temp)))


figure 
subplot(2,3,1)
hist(reshape(Est_S0_all,1,numel(Est_S0_all)),20)
set(gca, 'FontSize', 15)
title('S0','FontSize',25)
subplot(2,3,2)
hist(reshape(log(Est_sigma_all),1,numel(Est_sigma_all)),20)
set(gca, 'FontSize', 15)
title('log Sigma','FontSize',25)
subplot(2,3,3)
hist(reshape(Est_S0_all./Est_sigma_all,1,k1*k2*k3),20)
set(gca, 'FontSize', 15)
title('SNR','FontSize',25)
subplot(2,3,4)
hist(FA_show,25)
set(gca, 'FontSize', 15)
title('FA','FontSize',25)
subplot(2,3,5)
hist(MD_show,20)
set(gca, 'FontSize', 15)
title('MD','FontSize',25)
subplot(2,3,6)
scatter(FA_show,MD_show)
set(gca, 'FontSize', 15)
title('MD vs. FA','FontSize',25)
% subplot(2,4,6)
% hist(eval3_show,20)
% title('Eval. 3','FontSize',15)
% subplot(2,4,7)
% Ctrs = [0:2:20];
% Xtrs = hist(eval_ratio_show,Ctrs);
% bar(Ctrs, Xtrs)
% title('Eval. Ratio','FontSize',15)

% eval_ratio23_fsl = img_L2_temp./img_L3_temp;
% eval_ratio_fsl = 2*img_L1_temp./(img_L3_temp+img_L2_temp);

%% set uniform color range for heat map later
FA_temp_restrict = min(1,FA_temp);
MD_temp_restrict = max(0,MD_temp);
FA_top = max(FA_temp_restrict);
FA_bottom = min(FA_temp_restrict);
MD_top = max(MD_temp_restrict);
MD_bottom = min(MD_temp_restrict);

% for colormap
Eig1_FAcorrected = zeros(k1,k2,k3,3);
Eig1_fsl_FAcorrected = zeros(k1,k2,k3,3);

for i = 1:N
    [i1,i2,i3] = ind2sub([k1,k2,k3],i);
    Eig1_FAcorrected(i1,i2,i3,:) = abs(evec_temp(:,3,i)).*FA_temp(i);
    Eig1_fsl_FAcorrected(i1,i2,i3,:) = abs(img_V1_fsl_temp(i1,i2,i3,:)).*img_FA_fsl_temp(i1,i2,i3);
end

%%%%%%%%%%%%%%%%%%% FA, MD and color maps
%{
    %% 
    for k = 1:k3  %% x-perspective; for z-perspetive, change k1, n1, to k3, n3 and indexing to the third one
        figure;
        h    = [];
        h(1) = subplot(2,3,1);
        h(2) = subplot(2,3,2);
        h(3) = subplot(2,3,3);
        h(4) = subplot(2,3,4);
        h(5) = subplot(2,3,5);
        h(6) = subplot(2,3,6);

        %% FA map
        FA_temp_map = dti_data_reindex(data_FA(:,:,k));
        imagesc(FA_temp_map,'Parent',h(1));
        set(gca,'YDir','normal')
        axis off
        colormap('gray'); %% gray scale 
        caxis manual  %% 
        caxis([FA_bottom FA_top]);  %% fix scale such that it is comparable across slices
%         colorbar;

        %% MD map
        MD_temp_map = dti_data_reindex(data_MD(:,:,k));
        imagesc(MD_temp_map,'Parent',h(2));
        colormap('gray');
        caxis manual
        caxis([MD_bottom MD_top]);
%         colorbar;

        %% Color leading eigenvector map
        Eig1_FAcorrected_temp = dti_data_reindex(squeeze(Eig1_FAcorrected(:,:,k,:)));
%         Eig1_fig = figure;
        axis_FA_fig = axes;
        imagesc(Eig1_FAcorrected_temp,'Parent',h(3));
        axis(axis_FA_fig,'off')
        title('FA, MD, Colormap; Top est.; Bottom fsl')

%%%%%%%%%%%%%%%%% fsl results for comparison
        %% FA map
        FA_temp_map_fsl = dti_data_reindex(img_FA_fsl_temp(:,:,k));
        imagesc(FA_temp_map_fsl,'Parent',h(4));
        colormap('gray');
        caxis manual
        caxis([FA_bottom FA_top]);
%         colorbar;

        %% MD map
        MD_temp_map_fsl = dti_data_reindex(img_MD_fsl_temp(:,:,k));
        imagesc(MD_temp_map_fsl,'Parent',h(5));
        colormap('gray');
        caxis manual
        caxis([MD_bottom MD_top]);
%         colorbar;

        %% Color leading eigenvector map
        Eig1_FAcorrected_temp_fsl = dti_data_reindex(squeeze(Eig1_fsl_FAcorrected(:,:,k,:)));
%         Eig1_fig = figure;
        axis_FA_fig_fsl = axes;
        imagesc(Eig1_FAcorrected_temp_fsl,'Parent',h(6));
        axis(axis_FA_fig_fsl,'off')

        display(k);
    end


    options.use_axis = 0;
%}




clear options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 1;

[v_p,f_p] = compute_semiregular_sphere(5,options);
pos_p = v_p{end};

%%%plotting options 
options.spherical = 1;
% options for the display
options.use_color = 1;
options.color = 'wavelets';
options.use_elevation = 2;
options.rho = 0.5;
options.scaling = 1.5;
% for draw_fiber
plot_rho = options.rho;
plot_scale = options.scaling;


img_min = min(min(min(min(img_temp(img_temp>0))))); %%set all negative signal to the smallest positive signal
img_st = zeros(k1,k2,k3,n);
for i=1:n
    img_st(:,:,:,i) = img_temp(:,:,:,i+5)./Est_S0_all;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SH-ridge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SH-ridge penalty parameter
lambda_seq_SH = 10.^(linspace(log10(1e-2), log10(1e-6), 100)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmax8 = 8;   % SH order
SH_J5_lmax8 = SH_vertex(5, lmax8, 0);     % SH on dense grid for plotting
% generate SH matrix on sample grids, load pre-stored ones will be much faster
SH_matrix8 = sharmonic(bvec(:,6:46),lmax8);  
% generate response matrix, load pre-stored ones will be much much faster
R_matrix8 = Response_Rmatrix_construction(b_factor*bval(6)/1000,ratio_response,5,lmax8);

%% Input for SH_ridge function
penalty_matrix8 = penalty_matrix_ridge(lmax8);    % penalty matrix
design_SH8 = SH_matrix8*R_matrix8;     % design matrix
% 100 log-equally spaced lambda from 10^-2 to 10^-6


%%%%% Store SH-ridge results
f_est_all8 = zeros(k1,k2,k3,size(design_SH8,2));
dwi_SH_all8 =zeros(k1,k2,k3,size(design_SH8,1));
fod_SH_all8 =zeros(k1,k2,k3,size(pos_p,2));

df_SH_all8 = zeros(k1,k2,k3,length(lambda_seq_SH));
RSS_SH_all8 = zeros(k1,k2,k3,size(lambda_seq_SH,2));
BIC_SH_all8 = zeros(k1,k2,k3,size(lambda_seq_SH,2));
index_sele_SH_all8 = zeros(k1,k2,k3,1);


%% SH-ridge estimation
tic;
for i = 1:N
    
    [i1,i2,i3] = ind2sub([k1,k2,k3],i);
    
    DWI = reshape(img_st(i1,i2,i3,:),numel(img_st(i1,i2,i3,:)),1);
    [f_est_all_SH, dwi_est_all_SH, df_all_SH, RSS_all_SH, BIC_all_SH, index_selected_SH] ...
        = SH_ridge(DWI,design_SH8,penalty_matrix8,lambda_seq_SH);
    
    f_est_all8(i1,i2,i3,:) = f_est_all_SH(:,index_selected_SH);
    dwi_SH_all8(i1,i2,i3,:) = dwi_est_all_SH(:,index_selected_SH);
	fod_SH_all8(i1,i2,i3,:) = SH_J5_lmax8*f_est_all_SH(:,index_selected_SH);

	df_SH_all8(i1,i2,i3,:) = df_all_SH(index_selected_SH);
	RSS_SH_all8(i1,i2,i3,:) = RSS_all_SH(index_selected_SH);
	BIC_SH_all8(i1,i2,i3,:) = BIC_all_SH(index_selected_SH);
	index_sele_SH_all8(i1,i2,i3) = index_selected_SH;		    
end
SHridge8_time = toc;


%{
figure
plot_spherical_function(v_p,f_p,squeeze(fod_SH_all8(i1,i2,i3,:)),options);
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmax12 = 12;   % SH order
SH_J5_lmax12 = SH_vertex(5, lmax12, 0);     % SH on dense grid for plotting
% generate SH matrix on sample grids, load pre-stored ones will be much faster
SH_matrix12 = sharmonic(bvec(:,6:46),lmax12);  
% generate response matrix, load pre-stored ones will be much much faster
R_matrix12 = Response_Rmatrix_construction(b_factor*bval(6)/1000,ratio_response,5,lmax12);

%% Input for SH_ridge function
penalty_matrix12 = penalty_matrix_ridge(lmax12);    % penalty matrix
design_SH12 = SH_matrix12*R_matrix12;     % design matrix
% 100 log-equally spaced lambda from 10^-2 to 10^-6


%%%%% Store SH-ridge results
f_est_all12 = zeros(k1,k2,k3,size(design_SH12,2));
dwi_SH_all12 =zeros(k1,k2,k3,size(design_SH12,1));
fod_SH_all12 =zeros(k1,k2,k3,size(pos_p,2));

df_SH_all12 = zeros(k1,k2,k3,length(lambda_seq_SH));
RSS_SH_all12 = zeros(k1,k2,k3,size(lambda_seq_SH,2));
BIC_SH_all12 = zeros(k1,k2,k3,size(lambda_seq_SH,2));
index_sele_SH_all12 = zeros(k1,k2,k3,1);


%% SH-ridge estimation
tic;
for i = 1:N
    
    [i1,i2,i3] = ind2sub([k1,k2,k3],i);
    
    DWI = reshape(img_st(i1,i2,i3,:),numel(img_st(i1,i2,i3,:)),1);
    [f_est_all_SH, dwi_est_all_SH, df_all_SH, RSS_all_SH, BIC_all_SH, index_selected_SH] ...
        = SH_ridge(DWI,design_SH12,penalty_matrix12,lambda_seq_SH);
    
    f_est_all12(i1,i2,i3,:) = f_est_all_SH(:,index_selected_SH);
    dwi_SH_all12(i1,i2,i3,:) = dwi_est_all_SH(:,index_selected_SH);
	fod_SH_all12(i1,i2,i3,:) = SH_J5_lmax12*f_est_all_SH(:,index_selected_SH);

	df_SH_all12(i1,i2,i3,:) = df_all_SH(index_selected_SH);
	RSS_SH_all12(i1,i2,i3,:) = RSS_all_SH(index_selected_SH);
	BIC_SH_all12(i1,i2,i3,:) = BIC_all_SH(index_selected_SH);
	index_sele_SH_all12(i1,i2,i3) = index_selected_SH;		    
end
SHridge12_time = toc;

%{
figure
plot_spherical_function(v_p,f_p,squeeze(fod_SH_all12(i1,i2,i3,:)),options);
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% sCSD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_SCSD_lmax8_all = zeros(k1,k2,k3,size(design_SH8,2));
fod_SCSD_lmax8_all = zeros(k1,k2,k3,size(pos_p,2));

f_SCSD_lmax12_all = zeros(k1,k2,k3,size(design_SH12,2));
fod_SCSD_lmax12_all = zeros(k1,k2,k3,size(pos_p,2));

lmax_truncation = 4;
L_trucation = (lmax_truncation+1)*(lmax_truncation+2)/2;
lambda_scsd = 1;
thresh_scsd = 1e-3;
tau = 0.1;
maxit_scsd = 20;

tic;
for i = 1:N
    
    [i1,i2,i3] = ind2sub([k1,k2,k3],i);
    
    DWI = reshape(img_st(i1,i2,i3,:),numel(img_st(i1,i2,i3,:)),1);
    temp = reshape(f_est_all8(i1,i2,i3,:),numel(f_est_all8(i1,i2,i3,:)),1);
	fmatrix_initial_csd = temp;
    fmatrix_initial_csd((L_trucation+1):length(temp)) = 0;
    [f_est_lmax8_SCSD,thresh_lmax8_SCSD,iter_lmax8_SCSD] = superCSD(DWI, design_SH8, SH_J5_lmax8, fmatrix_initial_csd, lambda_scsd,  thresh_scsd, tau, maxit_scsd,false,0);
	%%
	dwi_SCSD_lmax8=design_SH8*f_est_lmax8_SCSD;
	FOD_SCSD_lmax8=SH_J5_lmax8*f_est_lmax8_SCSD;

	f_SCSD_lmax8_all(i1,i2,i3,:) = f_est_lmax8_SCSD;
    fod_SCSD_lmax8_all(i1,i2,i3,:) = FOD_SCSD_lmax8;
end
sCSD8_time = toc;

tic;
for i = 1:N
    
    [i1,i2,i3] = ind2sub([k1,k2,k3],i);
    
    DWI = reshape(img_st(i1,i2,i3,:),numel(img_st(i1,i2,i3,:)),1);
    temp = reshape(f_est_all12(i1,i2,i3,:),numel(f_est_all12(i1,i2,i3,:)),1);
	fmatrix_initial_csd = temp;
    fmatrix_initial_csd((L_trucation+1):length(temp)) = 0;
    [f_est_lmax12_SCSD,thresh_lmax12_SCSD,iter_lmax12_SCSD] = superCSD(DWI, design_SH12, SH_J5_lmax12, fmatrix_initial_csd, lambda_scsd,  thresh_scsd, tau, maxit_scsd,false,0);
	%%
	dwi_SCSD_lmax12=design_SH12*f_est_lmax12_SCSD;
	FOD_SCSD_lmax12=SH_J5_lmax12*f_est_lmax12_SCSD;

	f_SCSD_lmax12_all(i1,i2,i3,:) = f_est_lmax12_SCSD;
    fod_SCSD_lmax12_all(i1,i2,i3,:) = FOD_SCSD_lmax12;
end
sCSD12_time = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SN-lasso
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jmax = 3;
BW = 2;
half = 0;
Constraint = SN_vertex_symm(5,jmax,BW,half);  %% constraint matrix: Constraint*beta>=0;
C_trans_symm = C_trans_symm_construction(lmax8);  %% SH = SN*C_trans_symm
C_trans=(C_trans_symm*C_trans_symm')\C_trans_symm; %% f = C_trans*beta
design_SN = design_SH8*C_trans;   
SN_matrix_plot = SN_vertex_symm(5, jmax, BW, 0);

lambda_min_la = 1e-5;
lambda_max_la = 1e-2;
lambda_length_la = 50;
lambda_seq_la = 10.^(linspace(log10(lambda_max_la), log10(lambda_min_la), lambda_length_la));
window_percent = 0.05;
% relativeRSS_thresh = 2e-4;
% relativeRSS_thresh = 1e-4;
relativeRSS_thresh = 1e-4;
% relativeRSS_thresh = 1e-2;

ep_a = 10^-4;
ep_r = 10^-2;
maxit = 5000;
print = 1;

beta_SN_all = zeros(k1,k2,k3,size(design_SN,2), size(lambda_seq_la,2));
dwi_SN_all = zeros(k1,k2,k3,size(design_SH8,1), size(lambda_seq_la,2));
fod_SN_all = zeros(k1,k2,k3,size(pos_p,2), size(lambda_seq_la,2));
df_SN_all = zeros(k1,k2,k3,size(lambda_seq_la,2));
df_rank_SN_all = zeros(k1,k2,k3,size(lambda_seq_la,2));
RSS_SN_all = zeros(k1,k2,k3,size(lambda_seq_la,2));
% BICr_classo_all = zeros(k1,k2,k3,size(lambda_seq_la,2));
% AICr_classo_all = zeros(k1,k2,k3,size(lambda_seq_la,2));

% index_sele_SN_BICr = zeros(k1,k2,k3,1);
% df_sele_SN_BICr = zeros(k1,k2,k3,1);
% fod_sele_SN_BICr = zeros(k1,k2,k3,size(pos_p,2));
% BICr_sele_classo = zeros(k1,k2,k3);
BICr_ISO_classo = zeros(k1,k2,k3);

% index_sele_SN_AICr = zeros(k1,k2,k3,1);
% df_sele_SN_AICr = zeros(k1,k2,k3,1);
% fod_sele_SN_AICr = zeros(k1,k2,k3,size(pos_p,2));

index_sele_SN_RSSdiff = zeros(k1,k2,k3,1);
df_sele_SN_RSSdiff = zeros(k1,k2,k3,1);
fod_sele_SN_RSSdiff = zeros(k1,k2,k3,size(pos_p,2));

tic;
for i = 1:N
    
    [i1,i2,i3] = ind2sub([k1,k2,k3],i);
    DWI = reshape(img_st(i1,i2,i3,:),numel(img_st(i1,i2,i3,:)),1);
    [beta_est_all, z_all, w_all, dwi_est_all_SN, df_all_SN, df_rank_all_SN, RSS_all_SN, stop_criterion, index_selected_SN] ...
        = SN_lasso(DWI,design_SN,Constraint,lambda_seq_la,window_percent,relativeRSS_thresh,ep_a,ep_r,maxit,print);
    
%     idx_admm_BIC_rank_C = find(BIC_admm_rank_C==min(BIC_admm_rank_C(1:SN_stop_index)));
% 	idx_admm_AIC_rank_C = find(AIC_admm_rank_C==min(AIC_admm_rank_C(1:SN_stop_index)));
		    
% 	%% FOD estimation under selected models   
% 	FOD_admm_BIC_rank_C=FOD_admm_all_C(:,idx_admm_BIC_rank_C);
%     % FOD_admm_BIC_rank_C_st = fod_stand(FOD_admm_BIC_rank_C);
% 	dwi_admm_BIC_rank_C = dwi_admm_all_C(:,idx_admm_BIC_rank_C);
% 		    
% 	FOD_admm_AIC_rank_C=FOD_admm_all_C(:,idx_admm_AIC_rank_C);
% 	% FOD_admm_AIC_rank_C_st = fod_stand(FOD_admm_AIC_rank_C);
% 	dwi_admm_AIC_rank_C = dwi_admm_all_C(:,idx_admm_AIC_rank_C);
		    
	FOD_admm_RSSdiff_C=SN_matrix_plot*z_all(:,index_selected_SN);
	% FOD_admm_RSSdiff_C_st = fod_stand(FOD_admm_RSSdiff_C);
    dwi_admm_RSSdiff_C = dwi_est_all_SN(:,index_selected_SN);

	beta_SN_all(i1,i2,i3,:,:) = z_all;
    dwi_SN_all(i1,i2,i3,:,:) = dwi_est_all_SN;
% 	fod_SN_all(i1,i2,i3,:,:) = FOD_admm_all_C;
	df_SN_all(i1,i2,i3,:) = df_all_SN;
	df_rank_SN_all(i1,i2,i3,:) = df_rank_all_SN;
	RSS_SN_all(i1,i2,i3,:) = RSS_all_SN;
% 	BICr_classo_all(i1,i2,i3,:) = BIC_admm_rank_C;
% 	AICr_classo_all(i1,i2,i3,:) = AIC_admm_rank_C;

% 	index_sele_SN_BICr(i1,i2,i3) = idx_admm_BIC_rank_C;
% 	df_sele_SN_BICr(i1,i2,i3) = df_admm_rank_C(idx_admm_BIC_rank_C);
% 	fod_sele_SN_BICr(i1,i2,i3,:) = FOD_admm_BIC_rank_C;
%     BICr_sele_classo(i1,i2,i3) = BIC_admm_rank_C(idx_admm_BIC_rank_C);
%             
% 	index_sele_SN_AICr(i1,i2,i3) = idx_admm_AIC_rank_C;
% 	df_sele_SN_AICr(i1,i2,i3) = df_admm_rank_C(idx_admm_AIC_rank_C);
% 	fod_sele_SN_AICr(i1,i2,i3,:) = FOD_admm_AIC_rank_C;

	index_sele_SN_RSSdiff(i1,i2,i3) = index_selected_SN;
	df_sele_SN_RSSdiff(i1,i2,i3) = df_rank_all_SN(index_selected_SN);
	fod_sele_SN_RSSdiff(i1,i2,i3,:) = FOD_admm_RSSdiff_C;

            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Isotropic model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dwi_iso = ones(n,1)*mean(DWI);
    RSS_iso = sum((dwi_iso-DWI).^2);
    BICr_ISO_classo(i1,i2,i3) = n.*log(RSS_iso)+1.*log(n);
    
end
SNlasso_time = toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% z-perspective
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i3 = 1:1
    figure
    for i2 = 1:k2
        for i1 = 1:k1
            subplot(k2,k1,(k2-i2)*k1+i1)
            plot_spherical_function(v_p,f_p,reshape(fod_SH_all8(i1,i2,i3,:),1,size(pos_p,2)),options);
            view([0 0 1])
        end
    end
%     title(strcat('SH_slice',num2str(k3)));
%     savefig(strcat(save_path,'SH_est_z_','slice', num2str(k3),'.fig'));
end

for i3 = 1:1
    figure
    for i2 = 1:k2
        for i1 = 1:k1
            subplot(k2,k1,(k2-i2)*k1+i1)
            plot_spherical_function(v_p,f_p,reshape(fod_SCSD_lmax8_all(i1,i2,i3,:),1,size(pos_p,2)),options);
            view([0 0 1])
        end
    end
%     title(strcat('SH_slice',num2str(k3)));
%     savefig(strcat(save_path,'SH_est_z_','slice', num2str(k3),'.fig'));
end

for i3 = 1:1
    figure
    for i2 = 1:k2
        for i1 = 1:k1
            subplot(k2,k1,(k2-i2)*k1+i1)
            plot_spherical_function(v_p,f_p,reshape(fod_SCSD_lmax12_all(i1,i2,i3,:),1,size(pos_p,2)),options);
            view([0 0 1])
        end
    end
%     title(strcat('SH_slice',num2str(k3)));
%     savefig(strcat(save_path,'SH_est_z_','slice', num2str(k3),'.fig'));
end

for i3 = 1:1
    figure
    for i2 = 1:k2
        for i1 = 1:k1
            subplot(k2,k1,(k2-i2)*k1+i1)
            plot_spherical_function(v_p,f_p,reshape(fod_sele_SN_RSSdiff(i1,i2,i3,:),size(pos_p,2),1),options);
            view([0 0 1])
        end
    end
%     title(strcat('SH_slice',num2str(k3)));
%     savefig(strcat(save_path,'SH_est_z_','slice', num2str(k3),'.fig'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot estimated FOD (takes around 10 min for 15 by 15 slice), 
%% 4 methods (5 if plot relaxed stopping SN result), n3 (or n1) slices
%% usually SN and SN relaxed stopping resutls are similar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oax_left = 0;
oax_bottom = 0;
oax_width = 0.8;
oax_height = 1;

ax_left = 0.1;
ax_bottom = 0.05;
ax_width = 0.6;
ax_height = 0.9;

width_space = ax_width/k1;
height_space = ax_height/k2;     
            
%% plot est. FOD on top of FA color map
 %{  

        for i3 = 5:5
            
            FA_temp_map = dti_data_reindex(squeeze(data_FA(:,:,i3)));
            %% plot est. FOD on FA background inversely scaled with MD value
            FA_fig = figure('units','normalized','position',[0 0 0.8 1]);
            axis_FA_fig = axes;
            colormap(axis_FA_fig, gray);
            imagesc(FA_temp_map);
%             colorbar;
            caxis(axis_FA_fig,[FA_bottom FA_top]);
            ax = gca;  
            ax.Units = 'normalized';
            ax.OuterPosition = [oax_left, oax_bottom, oax_width, oax_height];
            ax.Position = [ax_left, ax_bottom, ax_width, ax_height];
            axis(axis_FA_fig,'off')

            
            hold on;
            for i1 = 1:k1
                for i2 = 1:k2
                    if(data_MD(i1,i2,i3)<1.5)
                        fig_factor = 0.2;
                    else
                        fig_factor = 0.2+(data_MD(i1,i2,i3)-1.5)/MD_top;
                    end

                    ax_temp_coordinate = [ax_left+(i1-1)*width_space+width_space*(fig_factor)/2, ax_bottom+(i2-1)*height_space+height_space*(fig_factor)/2, width_space*(1-fig_factor), height_space*(1-fig_factor)];
                    ax_temp = axes('Position', ax_temp_coordinate);
                    axis(ax_temp,'off');
                    colormap(ax_temp,jet);
                    options.use_axis = ax_temp_coordinate;
                    plot_spherical_function(v_p,f_p,squeeze(fod_sele_SN_RSSdiff(i1,i2,i3,:)),options);
%                     alpha(0.25)
                    view([0 0 1])
                end
            end
            hold off;
%             savefig(strcat(save_path,'SN_FA_est_z_','slice', num2str(k3),'.fig'));          
            display(i3);
        end
%}
