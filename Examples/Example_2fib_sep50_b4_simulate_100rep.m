%% Example script for FOD with 2 fiber (mimic Raymond's 2 fiber simulation)
%% First you need to add the pathes that contains all the needed functions
% Simulation with 100 replicates

clear all
path = '/Users/hao/Dropbox/DMRI_code/';

addpath(path);
addpath(strcat(path,'ADMM'));
addpath(strcat(path,'construction_functions'));
addpath(strcat(path,'toolbox_wavelet_meshes'));
addpath(strcat(path,'toolbox_wavelet_meshes/toolbox/'));
addpath(strcat(path,'MEALPix/'));
addpath(strcat(path,'help_functions/'));
addpath('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/matrix');
path_save='/Users/hao/Dropbox/stats_project/FOD_codes_simulation/';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define parameters for the simulation and load basis matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = 3; % vertex level (decides number of observation of this voxel) (2.5-->41 data points)
J_use = 2.5;
if(J_use==2)
    n_sample = 21;
elseif(J_use==2.5)
    n_sample = 41;
elseif(J_use==3)
    n_sample = 81;
elseif(J_use==4)
    nsample = 321;
end
b = [3, 3]; % back ground magnetic field strength (1 is same as b=1000)
ratio = [10, 10]; % ratio of the leading eigenvalue to the second eigenvalue in the signal simulation
weight = [0.5 0.5];
half = 1; % generate data on half shpere
lmax = 16;  % SH levels
jmax = 4; % SN levels corresponding to lmax

J_r = 5; % vertex level used for graph and representation purpose (dense)
b_response = b(1);  % b value for the response matrix Rmatrix that we use
ratio_response = ratio(1); % shape of the response function, could be misspecified 
sigma = 0.05;  %noixe level %middle  noise: note S0=1, so SNR=20 


%% FOD 
%% separation angle between the two fibers 
fod1_s = [0.7071 0  0.7071]; %%0,0,1; z-[0 0 1] [ sqrt(3)/2 0  1/2] [0.7071 0  0.7071] [1/2 0  sqrt(3)/2] [0.2588 0   0.9659] [sqrt(3)/2 0 1/2] [0.7660 0  0.6428]
fod2_s = [0 0 1]; %%1,0,0; x-axis
sep=acos(fod1_s*fod2_s');

%%% theta, phi of the fibers
phi0 = atan2([fod1_s(2) fod2_s(2)],[fod1_s(1) fod2_s(1)]);   %%phi: azimuthal  angle, [0,2pi)
phi0 = phi0+(phi0<0)*2*pi;
theta0 = acos([fod1_s(3) fod2_s(3)]);             %% theta: polar angle, [0,pi)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% using equal-angle grid with level J (3->81 direction on half sphere, 4->321; 2->21)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 1;

[vertex,~] = compute_semiregular_sphere(J,options); %%vertex and face of the grid 
pos = vertex{end};  %% x-y-z coordinates of the vertex 

% spherical coordinates of the vertex
phi = atan2(pos(2,:),pos(1,:))/(2*pi);   %%phi: azimuthal  angle, [0,1)
phi = phi+(phi<0);
theta = acos(pos(3,:))/(pi);             %% theta: polar angle, [0,1)

pos_corr = pos'*pos;
pos_pair = zeros(size(pos_corr,1),1);

% find the paired cabature points
for i = 1:size(pos_pair,1)
    pos_pair(i) = find(pos_corr(i,:)<-0.9999999);
end

sampling_index = zeros(size(pos_corr,1),1); % only uses half of the symmetrized needlets
for i = 1:size(pos_pair,1)
    if(pos(3,i)>-1e-15&&pos(3,i)<1e-15)
        if(pos(2,i)>-1e-15&&pos(2,i)<1e-15)
            if(pos(1,i)>0)
                sampling_index(i) = i;
            else
                sampling_index(i) = pos_pair(i);
            end
        elseif(pos(2,i)>1e-15)
            sampling_index(i) = i;
        elseif(pos(2,i)<-1e-15)
            sampling_index(i) = pos_pair(i);
        end
    elseif(pos(3,i)>1e-15)
        sampling_index(i) = i;
    elseif(pos(3,i)<-1e-15)
        sampling_index(i) = pos_pair(i);
    end
end
sampling_index = unique(sampling_index);

if(n_sample == 41)
    pos_sampling_h = pos(:,sampling_index); %% position of the half-sphere grid points 
    phi_h=phi(:,sampling_index)*180; 
    theta_h=theta(:,sampling_index)*180;

    %%% take 40 out of these 81 directions: at each level of theta, take about
    %%% half phi 
    index_1=find(theta_h<10); %%only 1
    n_1=size(index_1,2);

    index_t=find(theta_h>10&theta_h<20);
    n_2=size(index_t,2); %%6 
    [~, I]=sort(phi_h(index_t));
    index_2=index_t(1, I([1 3 5]));

    index_t=find(theta_h>20&theta_h<40);
    n_3=size(index_t,2); %%12
    [~, I]=sort(phi_h(index_t));
    index_3=index_t(1, I([1 3 5 7 9 11]));

    index_t=find(theta_h>40&theta_h<50);
    n_4=size(index_t,2); %%12
    [~, I]=sort(phi_h(index_t));
    index_4=index_t(1, I([1 3 5 7 9 11]));

    index_t=find(theta_h>50&theta_h<70);
    n_5=size(index_t,2); %%20
    [~, I]=sort(phi_h(index_t));
    index_5=index_t(1, I([1 3 5 7 9 11 13 15 17 19]));

    index_t=find(theta_h>70&theta_h<85);
    n_6=size(index_t,2); %%22
    [~, I]=sort(phi_h(index_t));
    index_6=index_t(1, I([1 3 5 7 9 11 13 15 17 19 21]));

    index_t=find(theta_h>85);
    n_7=size(index_t,2); %%8
    [~, I]=sort(phi_h(index_t));
    index_7=index_t(1, I([1 3 5 7]));

    index_s=unique([index_1 index_2 index_3 index_4 index_5 index_6 index_7]);
    sampling_grid_index=sampling_index(index_s); 
else
    sampling_grid_index=sampling_index;
end
clearvars index_1 index_2 index_3 index_4 index_5 index_6 index_7 index_t I i n_1 n_2 n_3 n_4 n_5 n_6 n_7 pos_corr pos_pair sampling_index phi_h theta_h;

pos_sampling = pos(:,sampling_grid_index); %% The x , y , z coordinates of the sampling grid.
phi_sampling = phi(:,sampling_grid_index); %% The sampled phi.
theta_sampling = theta(:,sampling_grid_index); %% The sampled theta.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load matrices
lmax8 = 8;
lmax12 = 12;
lmax16 = 16;

load(strcat('SH_vertex_J', num2str(J_r), '_lmax', num2str(lmax8), '.mat'));  % spherical harmonic basis evaluated on J=5(2562 grid points)
SH_J5_all_lmax8 = SH_vertex;

load(strcat('SH_vertex_J', num2str(J_r), '_lmax', num2str(lmax12), '.mat'));  % spherical harmonic basis evaluated on J=5(2562 grid points)
SH_J5_all_lmax12 = SH_vertex;

load(strcat('SH_vertex_J', num2str(J_r), '_lmax', num2str(lmax16), '.mat'));  % spherical harmonic basis evaluated on J=5(2562 grid points)
SH_J5_all_lmax16 = SH_vertex;

load(strcat('SH_vertex_J', num2str(J), '_lmax', num2str(lmax8), '_h.mat'));  % spherical harmonic basis for observations
SH_vertex_lmax8 = SH_vertex;

load(strcat('SH_vertex_J', num2str(J), '_lmax', num2str(lmax12), '_h.mat'));  % spherical harmonic basis for observations
SH_vertex_lmax12 = SH_vertex;

load(strcat('SH_vertex_J', num2str(J), '_lmax', num2str(lmax16), '_h.mat'));  % spherical harmonic basis for observations
SH_vertex_lmax16 = SH_vertex;

load(strcat('C_symm_lmax', num2str(lmax), '.mat'));  %% SH coefficients of the symmetrized needlets basis
load(strcat('SN_vertex_symm_J5_jmax', num2str(jmax), '.mat'));% symmetrized spherical needlets basis J=5(2562 grid points), jmax=4
load(strcat('Rmatrix_J5','_lmax',num2str(lmax8),'_b',num2str(b_response),'_ratio',num2str(ratio_response),'.mat' ));  %file name of the file that stores the corresponding R matrix
Rmatrix_lmax8 = Rmatrix;
load(strcat('Rmatrix_J5','_lmax',num2str(lmax12),'_b',num2str(b_response),'_ratio',num2str(ratio_response),'.mat' ));  %file name of the file that stores the corresponding R matrix
Rmatrix_lmax12 = Rmatrix;
load(strcat('Rmatrix_J5','_lmax',num2str(lmax16),'_b',num2str(b_response),'_ratio',num2str(ratio_response),'.mat' ));  %file name of the file that stores the corresponding R matrix
Rmatrix_lmax16 = Rmatrix;
load(strcat('Rmatrix_J5','_lmax',num2str(lmax),'_b',num2str(b_response),'_ratio',num2str(ratio_response),'.mat' ));  %file name of the file that stores the corresponding R matrix
clearvars SH_vertex Rmatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% saving path and folder name
save_path = strcat(path_save,'simulation_review/','2fib_sep',num2str(round(sep*180/pi,0)),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'/');
mkdir(save_path);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up design matrices and penalty parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SH+ridge
if(n_sample==41)
    design_SH_lmax8 = SH_vertex_lmax8(index_s,:)*Rmatrix_lmax8;
    design_SH_lmax12 = SH_vertex_lmax12(index_s,:)*Rmatrix_lmax12;
    design_SH_lmax16 = SH_vertex_lmax16(index_s,:)*Rmatrix_lmax16;
else
    design_SH_lmax8 = SH_vertex_lmax8*Rmatrix_lmax8;
    design_SH_lmax12 = SH_vertex_lmax12*Rmatrix_lmax12;
    design_SH_lmax16 = SH_vertex_lmax16*Rmatrix_lmax16;
end

%% Penalty matrix: Laplace Beltrame operator
Penalty_matrix_lmax8 = zeros( (lmax8+1)*(lmax8+2)/2 , 1);
for l = 0:2:lmax8
    for m = (-l):l
        Penalty_matrix_lmax8((l+1)*(l+2)/2 - (l-m)) = l^2*(l+1)^2;  %% more penalty on higher order SH basis
    end
end
Penalty_matrix_lmax8 = diag(Penalty_matrix_lmax8);

Penalty_matrix_lmax12 = zeros( (lmax12+1)*(lmax12+2)/2 , 1);
for l = 0:2:lmax12
    for m = (-l):l
        Penalty_matrix_lmax12((l+1)*(l+2)/2 - (l-m)) = l^2*(l+1)^2;  %% more penalty on higher order SH basis
    end
end
Penalty_matrix_lmax12 = diag(Penalty_matrix_lmax12);

Penalty_matrix_lmax16 = zeros( (lmax16+1)*(lmax16+2)/2 , 1);
for l = 0:2:lmax16
    for m = (-l):l
        Penalty_matrix_lmax16((l+1)*(l+2)/2 - (l-m)) = l^2*(l+1)^2;  %% more penalty on higher order SH basis
    end
end
Penalty_matrix_lmax16 = diag(Penalty_matrix_lmax16);

%% symmetric needlets design matrix
Constraint = SN_vertex_symm;  %% constraint matrix: Constraint*beta>=0;
C_trans=(C_trans_symm*C_trans_symm')\C_trans_symm; % SH*f = SN*C_trans_symm' 
design_SN = design_SH_lmax16*C_trans;   

%% sequences of penalty parameters 
%%% for SH+ridge
lambda_min = 1e-6;  %%smallest lamabda
lambda_max = 1e-2;  %%largest lambda 
lambda_seq=10.^(linspace(log10(lambda_min), log10(lambda_max), 100));
 
%%% for classo
lambda_min_la = 1e-5;
lambda_max_la = 1e-2;
lambda_length_la = 50;
lambda_seq_la=10.^(linspace(log10(lambda_max_la), log10(lambda_min_la), lambda_length_la));
window_percent = 0.05;
relativeRSS_thresh = 2e-4;
%%% admm stoping criterion
ep_r=1e-2;
ep_a=1e-4;
maxit=5000;
print=1;

%%% for super_CSD
lmax_truncation = 4;
L_trucation = (lmax_truncation+1)*(lmax_truncation+2)/2;
lambda_scsd = 1;
thresh_scsd = 1e-5;
tau = 0.1;
maxit_scsd = 500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic;
n_rep = 100; 

for rep = 1:n_rep
    
    display(strcat('Start rep ', num2str(rep)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate data and pre analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% generate dwi signals on the equal-angle grid/gradient-direction grid 
    seed = rep*10+21;
    [DWI, theta_r, phi_r]= DWI_generate(J_use, b, ratio, weight, theta0, phi0, sigma, half, seed);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% lmax = 8
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Method I: SH+ridge regression
    [f_est_all_SH_lmax8, dwi_est_all_SH_lmax8, df_all_SH_lmax8, RSS_all_SH_lmax8, BIC_all_SH_lmax8, index_sel_SH_lmax8] ...
    = SH_ridge(DWI,design_SH_lmax8,Penalty_matrix_lmax8,lambda_seq);
    
    %% estimated FOD using BIC selected lambda
    dwi_est_SH_lmax8 = dwi_est_all_SH_lmax8(:,index_sel_SH_lmax8);
    f_est_lmax8= f_est_all_SH_lmax8(:,index_sel_SH_lmax8);
    FOD_SH_est_lmax8=SH_J5_all_lmax8*f_est_lmax8;
    FOD_SH_est_lmax8_st=fod_stand(FOD_SH_est_lmax8);
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Method II: SuperCsd on SH+ridge regression
    f_ini = reshape(f_est_lmax8,numel(f_est_lmax8),1);
    f_ini((L_trucation+1):length(f_est_lmax8)) = 0;
    [f_est_sCSD_lmax8,diff_sCSD_lmax8,iter_sCSD_lmax8] = superCSD(DWI, design_SH_lmax8, SH_J5_all_lmax8, f_ini, lambda_scsd,  thresh_scsd, tau, maxit_scsd,false,0);
    
    dwi_est_SCSD_lmax8=design_SH_lmax8*f_est_sCSD_lmax8;
    FOD_sCSD_lmax8=SH_J5_all_lmax8*f_est_sCSD_lmax8;
    FOD_sCSD_lmax8_st = fod_stand(FOD_sCSD_lmax8);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% lmax = 12
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Method I: SH+ridge regression
    [f_est_all_SH_lmax12, dwi_est_all_SH_lmax12, df_all_SH_lmax12, RSS_all_SH_lmax12, BIC_all_SH_lmax12, index_sel_SH_lmax12] ...
    = SH_ridge(DWI,design_SH_lmax12,Penalty_matrix_lmax12,lambda_seq);
    
    %% estimated FOD using BIC selected lambda
    dwi_est_SH_lmax12 = dwi_est_all_SH_lmax12(:,index_sel_SH_lmax12);
    f_est_lmax12= f_est_all_SH_lmax12(:,index_sel_SH_lmax12);
    FOD_SH_est_lmax12=SH_J5_all_lmax12*f_est_lmax12;
    FOD_SH_est_lmax12_st=fod_stand(FOD_SH_est_lmax12);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Method II: SuperCsd on SH+ridge regression
    f_ini = reshape(f_est_lmax12,numel(f_est_lmax12),1);
    f_ini((L_trucation+1):length(f_est_lmax12)) = 0;
    [f_est_sCSD_lmax12,diff_sCSD_lmax12,iter_sCSD_lmax12] = superCSD(DWI, design_SH_lmax12, SH_J5_all_lmax12, f_ini, lambda_scsd,  thresh_scsd, tau, maxit_scsd,false,0);
    
    dwi_est_SCSD_lmax12=design_SH_lmax12*f_est_sCSD_lmax12;
    FOD_sCSD_lmax12=SH_J5_all_lmax12*f_est_sCSD_lmax12;
    FOD_sCSD_lmax12_st = fod_stand(FOD_sCSD_lmax12);
    
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% lmax = 16
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Method I: SH+ridge regression
    [f_est_all_SH_lmax16, dwi_est_all_SH_lmax16, df_all_SH_lmax16, RSS_all_SH_lmax16, BIC_all_SH_lmax16, index_sel_SH_lmax16] ...
    = SH_ridge(DWI,design_SH_lmax16,Penalty_matrix_lmax16,lambda_seq);
    
    %% estimated FOD using BIC selected lambda
    dwi_est_SH_lmax16 = dwi_est_all_SH_lmax16(:,index_sel_SH_lmax16);
    f_est_lmax16= f_est_all_SH_lmax16(:,index_sel_SH_lmax16);
    FOD_SH_est_lmax16=SH_J5_all_lmax16*f_est_lmax16;
    FOD_SH_est_lmax16_st=fod_stand(FOD_SH_est_lmax16);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Method II: SuperCsd on SH+ridge regression
    f_ini = reshape(f_est_lmax16,numel(f_est_lmax16),1);
    f_ini((L_trucation+1):length(f_est_lmax16)) = 0;
    [f_est_sCSD_lmax16,diff_sCSD_lmax16,iter_sCSD_lmax16] = superCSD(DWI, design_SH_lmax16, SH_J5_all_lmax16, f_ini, lambda_scsd,  thresh_scsd, tau, maxit_scsd,false,0);
    
    dwi_est_SCSD_lmax16=design_SH_lmax16*f_est_sCSD_lmax16;
    FOD_sCSD_lmax16=SH_J5_all_lmax16*f_est_sCSD_lmax16;
    FOD_sCSD_lmax16_st = fod_stand(FOD_sCSD_lmax16);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% useful quantities for analysis (cabature points coordinate, needlets' norm)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    
    pix1 = cell2mat(pix2vec(1,'nest',false));
    pix2 = cell2mat(pix2vec(2^1,'nest',false));
    pix3 = cell2mat(pix2vec(2^2,'nest',false));
    pix4 = cell2mat(pix2vec(2^3,'nest',false));
    pix5 = cell2mat(pix2vec(2^4,'nest',false));

    pix_1 = cell2mat(pix2ang(1,'nest',false));
    pix_2 = cell2mat(pix2ang(2,'nest',false));
    pix_3 = cell2mat(pix2ang(4,'nest',false));
    pix_4 = cell2mat(pix2ang(8,'nest',false));
    pix_5 = cell2mat(pix2ang(16,'nest',false));

    pix_all = [pix1 pix2 pix3 pix4 pix5];
    pixall = [pix_1 pix_2 pix_3 pix_4 pix_5];

    angle11_p1 = acos(pix1'*fod1');
    angle12_p1 = acos(pix1'*(-fod1)');
    angle1_p1 = min(angle11_p1,angle12_p1);
    angle21_p1 = acos(pix1'*fod2');
    angle22_p1 = acos(pix1'*(-fod2)');
    angle2_p1 = min(angle21_p1,angle22_p1);

    angle11_p2 = acos(pix2'*fod1');
    angle12_p2 = acos(pix2'*(-fod1)');
    angle1_p2 = min(angle11_p2,angle12_p2);
    angle21_p2 = acos(pix2'*fod2');
    angle22_p2 = acos(pix2'*(-fod2)');
    angle2_p2 = min(angle21_p2,angle22_p2);

    angle11_p3 = acos(pix3'*fod1');
    angle12_p3 = acos(pix3'*(-fod1)');
    angle1_p3 = min(angle11_p3,angle12_p3);
    angle21_p3 = acos(pix3'*fod2');
    angle22_p3 = acos(pix3'*(-fod2)');
    angle2_p3 = min(angle21_p3,angle22_p3);

    angle11_p4 = acos(pix4'*fod1');
    angle12_p4 = acos(pix4'*(-fod1)');
    angle1_p4 = min(angle11_p4,angle12_p4);
    angle21_p4 = acos(pix4'*fod2');
    angle22_p4 = acos(pix4'*(-fod2)');
    angle2_p4 = min(angle21_p4,angle22_p4);

    angle11_p5 = acos(pix5'*fod1');
    angle12_p5 = acos(pix5'*(-fod1)');%
    angle1_p5 = min(angle11_p5,angle12_p5);
    angle21_p5 = acos(pix5'*fod2');
    angle22_p5 = acos(pix5'*(-fod2)');
    angle2_p5 = min(angle21_p5,angle22_p5);

    cabature_corr = pix_all'*pix_all;
    cabature_pair = zeros(size(cabature_corr,1),1);

    % find the paired cabature points
    for i = 1:size(cabature_pair,1)
        cabature_pair(i) = find(cabature_corr(i,:)<-0.9999999);
    end

    cabature_use = zeros(size(cabature_corr,1)/2,1); % only uses half of the symmetrized needlets
    count = 1;
    for i = 1:size(cabature_pair,1)
        if cabature_pair(i)>i
            cabature_use(count) = i;
            count=count+1;
        end
    end

    cabature_vertex = abs(pix_all'*pos_plot);
    cabature_vertex = cabature_vertex(cabature_use,:);
    [caba_ver, caba_ver_idx] = sort(cabature_vertex,2,'descend');


    angle1_sn={angle1_p1 angle1_p2 angle1_p3 angle1_p4 angle1_p5};
    angle1_sn=cat(1, angle1_sn{:});
    angle1_sn_use = angle1_sn(cabature_use);

    angle2_sn={angle2_p1 angle2_p2 angle2_p3 angle2_p4 angle2_p5};
    angle2_sn=cat(1, angle2_sn{:});
    angle2_sn_use = angle2_sn(cabature_use);

    % figure;
    % plot(pixall(1,:),pixall(2,:),'.');
    % hold on;
    % plot(theta0,phi0,'.','Color','red','MarkerSize', 20);
    % hold on;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% get each level of the needlets and design_SN matrix
    SN_J5_symm_0 = SN_vertex_symm(:,1);
    SN_J5_symm_1 = SN_vertex_symm(:,2:7);
    SN_J5_symm_2 = SN_vertex_symm(:,8:31);
    SN_J5_symm_3 = SN_vertex_symm(:,32:127);
    SN_J5_symm_4 = SN_vertex_symm(:,128:511);
    SN_J5_symm_5 = SN_vertex_symm(:,512:2047);
    %SN_J5_symm_6 = SN_J5_symm_use(:,2048:8191);

    design_SN_0 = design_SN(:,1);
    design_SN_1 = design_SN(:,2:7);
    design_SN_2 = design_SN(:,8:31);
    design_SN_3 = design_SN(:,32:127);
    design_SN_4 = design_SN(:,128:511);
    design_SN_5 = design_SN(:,512:2047);
    %design_SN_6 = design_SN(:,2048:8191);

    l2_1 = sqrt(sum(SN_vertex_symm(:,1+1).^2)*(4*pi)/length(SN_vertex_symm(:,1)));
    l2_2 = sqrt(sum(SN_vertex_symm(:,1+6+1).^2)*(4*pi)/length(SN_vertex_symm(:,1)));
    l2_3 = sqrt(sum(SN_vertex_symm(:,1+6+24+1).^2)*(4*pi)/length(SN_vertex_symm(:,1)));
    l2_4 = sqrt(sum(SN_vertex_symm(:,1+6+24+96+1).^2)*(4*pi)/length(SN_vertex_symm(:,1)));
    l2_5 = sqrt(sum(SN_vertex_symm(:,1+6+24+96+384+1).^2)*(4*pi)/length(SN_vertex_symm(:,1)));
    % l2_6 = sqrt(sum(SN_vertex_symm(:,1+6+24+96+384+1536+1).^2)*(4*pi)/length(SN_vertex_symm(:,1)));
    
    %}

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Method III: classo+ADMM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [beta_est_all, z_all, w_all, dwi_est_all_SN, df_all_SN, df_rank_all_SN, RSS_all_SN, stop_criterion, index_selected_SN] ...
    = SN_lasso(DWI,design_SN,Constraint,lambda_seq_la,window_percent,relativeRSS_thresh,ep_a,ep_r,maxit,print);

    FOD_SN=SN_vertex_symm*z_all(:,index_selected_SN);
    FOD_SN_st = fod_stand(FOD_SN);
    save_name = strcat('2fib_sep',num2str(round(sep*180/pi,0)),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'_rep',num2str(rep),'.mat');
    save(strcat(save_path,save_name));
    display(strcat('Done with rep ', num2str(rep)));
end;
total_time = toc;
total_time

%{
%% plot sCSD FOD estimator
figure
plot_spherical_function(v_p,f_p,FOD_SN,options);  % FOD_sCSD_lmax12
%FOD_SH_est_lmax16 FOD_sCSD_lmax16
hold on;
draw_fiber(theta_r,phi_r,1.5,0.5*max(FOD_SN));
%}


% 285.2760 min
%{


figure;
subplot(1,7,1)
hold all;
plot_spherical_function(v_plot,f_plot,FOD_SH_est_s_lmax8,options);
draw_fiber(theta0,phi0,1.5,0.5*max(FOD_SH_est_s_lmax8));

subplot(1,7,2)
hold all;
plot_spherical_function(v_plot,f_plot,FOD_SCSD_lmax8,options);
draw_fiber(theta0,phi0,1.5,0.5*max(FOD_SCSD_lmax8));

subplot(1,7,3)
hold all;
plot_spherical_function(v_plot,f_plot,FOD_SCSD_lmax12,options);
draw_fiber(theta0,phi0,1.5,0.5*max(FOD_SCSD_lmax12));

subplot(1,7,4)
hold all;
plot_spherical_function(v_plot,f_plot,FOD_SCSD_lmax16,options);
draw_fiber(theta0,phi0,1.5,0.5*max(FOD_SCSD_lmax16));

subplot(1,7,5)
hold all;
plot_spherical_function(v_plot,f_plot,FOD_admm_BIC_rank_C,options);
draw_fiber(theta0,phi0,1.5,0.5*max(FOD_admm_BIC_rank_C));

subplot(1,7,6)
hold all;
plot_spherical_function(v_plot,f_plot,FOD_admm_AIC_rank_C,options);
draw_fiber(theta0,phi0,1.5,0.5*max(FOD_admm_AIC_rank_C));

subplot(1,7,7)
hold all;
plot_spherical_function(v_plot,f_plot,FOD_admm_RSSdiff_C,options);
draw_fiber(theta0,phi0,1.5,0.5*max(FOD_admm_RSSdiff_C));

diff_logRSS = diff(log10(RSS_admm_C),1);
diff_logRSS_adj = diff_logRSS./stop_spacing;
diff_rela = diff_logRSS_adj./log10(RSS_admm_C(1:499));
diff_logRSS_adj_rela = diff(log10(RSS_admm_C),1)./log10(RSS_admm_C(1:499));

%}