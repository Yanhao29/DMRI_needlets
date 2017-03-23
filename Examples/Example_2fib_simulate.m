%% Example script for FOD with 2 directions with separation angle 90
%% First you need to add the pathes that contains all the needed functions

clear all
path = '/Users/hao/Dropbox/DMRI_code/';

addpath(path);
addpath(strcat(path,'ADMM'));
addpath(strcat(path,'construction_functions'));
addpath(strcat(path,'toolbox_wavelet_meshes'));
addpath(strcat(path,'toolbox_wavelet_meshes/toolbox/'));
addpath(strcat(path,'MEALPix/'));
addpath(strcat(path,'help_functions/'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% plotting options for toolbox_wavelet_meshes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%vertex construction options 
options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 1;

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
plot_sacle = options.scaling;

% denser grid for interpolation and plots 
[v_p,f_p] = compute_semiregular_sphere(5,options);
pos_p = v_p{end};
phi_p = atan2(pos_p(2,:),pos_p(1,:))/(2*pi);   %%phi: azimuthal  angle, [0,2)
phi_p = phi_p+(phi_p<0);
theta_p = acos(pos_p(3,:))/(pi);             %% theta: polar angle, [0,1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Generate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = 3;  % n_sample = 81
b = [3 3];  % bvalue
ratio = [10 10];    % eigenvalue3*2/(eigenvalue1+eigenvalue2)
weight = [0.5, 0.5];    % weight of each fiber
theta0 = [0, pi/2];     % theta
phi0 = [0, 0];      % phi
sigma = 0.05;       % richian noise 
half = 1;       % generate data on half shpere
seed = 31;      % seed for random rotation of the FOD

tic;
% theta_r, phi_r: angles of fiber(s) after random rotation according to seed
[DWI, theta_r, phi_r]= DWI_generate(J, b, ratio, weight, theta0, phi0, sigma, half, seed);
generate_time = toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SH-ridge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmax = 8;   % SH order
% generate SH matrix on sample grids, load pre-stored ones will be much faster
SH_matrix = SH_vertex(J, lmax, half);  
% generate response matrix, load pre-stored ones will be much much faster
R_matrix = Response_Rmatrix_construction(b(1),ratio(1),5,lmax);

%% Input for SH_ridge function
penalty_matrix = penalty_matrix_ridge(lmax);    % penalty matrix
design_SH = SH_matrix*R_matrix;     % design matrix
% 100 log-equally spaced lambda from 10^-2 to 10^-6
lambda_seq_SH = 10.^(linspace(log10(1e-2), log10(1e-6), 100)); 

%% SH-ridge estimation
tic;
[f_est_all_SH, dwi_est_all_SH, df_all_SH, RSS_all_SH, BIC_all_SH, index_selected_SH] ...
    = SH_ridge(DWI,design_SH,penalty_matrix,lambda_seq_SH);
SHridge8_time = toc;

%% BIC selected SH coefficients
f_est_selected_SH = f_est_all_SH(:,index_selected_SH);

SH_matrix_plot = SH_vertex(5, lmax, 0);     % SH on dense grid for plotting

%% SH-ridge estimated FOD using BIC selected model
FOD_SH = SH_matrix_plot*f_est_selected_SH;
FOD_SH_st = fod_stand(FOD_SH);

%% plot selected FOD estimator
figure
plot_spherical_function(v_p,f_p,FOD_SH,options);
hold on;
draw_fiber(theta_r,phi_r,1.5,0.5*max(FOD_SH));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% superCSD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmax_truncation = 4;
L_trucation = (lmax_truncation+1)*(lmax_truncation+2)/2;
lambda_scsd = 1;
thresh_scsd = 1e-5;
tau = 0.1;
maxit_scsd = 500;
f_ini = reshape(f_est_selected_SH,numel(f_est_selected_SH),1);
f_ini((L_trucation+1):length(f_est_selected_SH)) = 0;
tic;
[f_est_sCSD,diff_sCSD,iter_sCSD] = superCSD(DWI, design_SH, SH_matrix_plot, f_ini, lambda_scsd,  thresh_scsd, tau, maxit_scsd,false,0);
SuperCSD8_time = toc;
%%
FOD_sCSD=SH_matrix_plot*f_est_sCSD;
FOD_sCSD_st = fod_stand(FOD_sCSD);

%% plot sCSD FOD estimator
figure
plot_spherical_function(v_p,f_p,FOD_sCSD,options);
hold on;
draw_fiber(theta_r,phi_r,1.5,0.5*max(FOD_sCSD));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SN-lasso
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create mex object that calls ADMM_classo.c (only need to do it once)
% mex -v -largeArrayDims /Users/hao/Dropbox/DMRI_code/ADMM/ADMM_classo.c -lmwlapack -lmwblas
jmax = 3;
BW = 2;
half = 0;
Constraint = SN_vertex_symm(5,jmax,BW,half);  %% constraint matrix: Constraint*beta>=0;
C_trans_symm = C_trans_symm_construction(lmax);  %% SH = SN*C_trans_symm
C_trans=(C_trans_symm*C_trans_symm')\C_trans_symm; %% f = C_trans*beta
design_SN = design_SH*C_trans;   

lambda_min_la = 1e-5;
lambda_max_la = 1e-2;
lambda_length_la = 50;
lambda_seq_la = 10.^(linspace(log10(lambda_max_la), log10(lambda_min_la), lambda_length_la));
window_percent = 0.05;
relativeRSS_thresh = 2e-4;
% relativeRSS_thresh = 1e-4;
% relativeRSS_thresh = 5e-3;
% relativeRSS_thresh = 1e-2;

ep_a = 10^-4;
ep_r = 10^-2;
maxit = 5000;
print = 1;
tic;
[beta_est_all, z_all, w_all, dwi_est_all_SN, df_all_SN, df_rank_all_SN, RSS_all_SN, stop_criterion, index_selected_SN] ...
    = SN_lasso(DWI,design_SN,Constraint,lambda_seq_la,window_percent,relativeRSS_thresh,ep_a,ep_r,maxit,print);
SNlasso_time = toc;

SN_matrix_plot = SN_vertex_symm(5, jmax, BW, 0);     % SN on dense grid for plotting

FOD_SN=SN_matrix_plot*z_all(:,index_selected_SN);
FOD_SN_st = fod_stand(FOD_SN);

%% plot sCSD FOD estimator
figure
plot_spherical_function(v_p,f_p,FOD_SN,options);
hold on;
draw_fiber(theta_r,phi_r,1.5,0.5*max(FOD_SN));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Hellinger distance
% SH representation
coe_sh1 = Dirac_SH_coe(lmax,theta_r(1),phi_r(1)); %% SH coefficients
coe_sh2 = Dirac_SH_coe(lmax,theta_r(2),phi_r(2));
dirac_sh1 = SH_matrix_plot*coe_sh1;
dirac_sh2 = SH_matrix_plot*coe_sh2;
dirac_sh = weight(1)*dirac_sh1+weight(2)*dirac_sh2;  %%SH representation
dirac_sh_st =fod_stand(dirac_sh); %%standardized to be nonnegative and sum =1
% Hellinger distance
HD_SH = hellinger_dis(dirac_sh_st, FOD_SH_st);
HD_sCSD = hellinger_dis(dirac_sh_st, FOD_sCSD_st);
HD_SN = hellinger_dis(dirac_sh_st, FOD_SN_st);

%%%%%%%%%% peak finding, angle error, separation angle
kmin = 40;
Dis = squareform(pdist(pos_p','cosine'));
cut_thresh = 0.25;

[~, ~, ~, ~, ~, peak_pos_SH_final] = FOD_peak(FOD_SH, Dis, kmin, cut_thresh, pos_p, theta_p, phi_p);
[~, ~, ~, ~, ~, peak_pos_sCSD_final] = FOD_peak(FOD_sCSD, Dis, kmin, cut_thresh, pos_p, theta_p, phi_p);
[~, ~, ~, ~, ~, peak_pos_SN_final] = FOD_peak(FOD_SN, Dis, kmin, cut_thresh, pos_p, theta_p, phi_p);
   
fod1_r = [sin(theta_r(1))*cos(phi_r(1)); sin(theta_r(1))*sin(phi_r(1)); cos(theta_r(1))];
fod2_r = [sin(theta_r(2))*cos(phi_r(2)); sin(theta_r(2))*sin(phi_r(2)); cos(theta_r(2))];

SH_err1 = min(separation_angle(fod1_r, peak_pos_SH_final(:,2)),separation_angle(fod1_r, peak_pos_SH_final(:,1)))*180/pi;
sCSD_err1 = min(separation_angle(fod1_r, peak_pos_sCSD_final(:,2)),separation_angle(fod1_r, peak_pos_sCSD_final(:,1)))*180/pi;
SN_err1 = min(separation_angle(fod1_r, peak_pos_SN_final(:,2)),separation_angle(fod1_r, peak_pos_SN_final(:,1)))*180/pi;

SH_err2 = min(separation_angle(fod2_r, peak_pos_SH_final(:,2)),separation_angle(fod2_r, peak_pos_SH_final(:,1)))*180/pi;
sCSD_err2 = min(separation_angle(fod2_r, peak_pos_sCSD_final(:,2)),separation_angle(fod2_r, peak_pos_sCSD_final(:,1)))*180/pi;
SN_err2 = min(separation_angle(fod2_r, peak_pos_SN_final(:,2)),separation_angle(fod2_r, peak_pos_SN_final(:,1)))*180/pi;

SH_sep = separation_angle(peak_pos_SH_final(:,1), peak_pos_SH_final(:,2))*180/pi;
sCSD_sep = separation_angle(peak_pos_sCSD_final(:,1), peak_pos_sCSD_final(:,2))*180/pi;
SN_sep = separation_angle(peak_pos_SN_final(:,1), peak_pos_SN_final(:,2))*180/pi;


