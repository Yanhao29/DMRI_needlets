% Analysis of 2fiber simulation with 100 replicates
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
addpath('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/matrix');
path_save='/Users/hao/Dropbox/stats_project/FOD_codes_simulation/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
plot_scale = options.scaling;

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
b = [4,4]; % back ground magnetic field strength (1 is same as b=1000)
ratio = [10,10]; % ratio of the leading eigenvalue to the second eigenvalue in the signal simulation
weight = [0.7,0.3];
half = 1; % generate data on half shpere
lmax = 8;  % SH levels
jmax = 3; % SN levels corresponding to lmax

J_r = 5; % vertex level used for graph and representation purpose (dense)
b_response = b(1);  % b value for the response matrix Rmatrix that we use
ratio_response = ratio(1); % shape of the response function, could be misspecified 
sigma = 0.05;  %noixe level %middle  noise: note S0=1, so SNR=20 


%% FOD 
%% separation angle between the two fibers 
fod1_s = [1 0 0]; %%0,0,1; z-[0 0 1] [ sqrt(3)/2 0  1/2] [0.7071 0  0.7071] [1/2 0  sqrt(3)/2] [0.2588 0   0.9659] [sqrt(3)/2 0 1/2]
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
design_SN = design_SH_lmax8*C_trans;   

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
lambda_sCSD = 1;
thresh_sCSD = 1e-5;
tau = 0.1;
maxit_sCSD = 500;

N_rep = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SH rep of true FOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coe_sh1_lmax8 = Dirac_SH_coe(lmax8,theta0(1),phi0(1)); %% SH coefficients
coe_sh2_lmax8 = Dirac_SH_coe(lmax8,theta0(2),phi0(2));
dirac_sh1_lmax8 = SH_J5_all_lmax8*coe_sh1_lmax8;
dirac_sh2_lmax8 = SH_J5_all_lmax8*coe_sh2_lmax8;
dirac_sh_lmax8 = weight(1)*dirac_sh1_lmax8+weight(2)*dirac_sh2_lmax8;  %%SH representation
dirac_sh_st_lmax8 =fod_stand(dirac_sh_lmax8); %%standardized to be nonnegative and sum =1

coe_sh1_lmax12 = Dirac_SH_coe(lmax12,theta0(1),phi0(1)); %% SH coefficients
coe_sh2_lmax12 = Dirac_SH_coe(lmax12,theta0(2),phi0(2));
dirac_sh1_lmax12 = SH_J5_all_lmax12*coe_sh1_lmax12;
dirac_sh2_lmax12 = SH_J5_all_lmax12*coe_sh2_lmax12;
dirac_sh_lmax12 = weight(1)*dirac_sh1_lmax12+weight(2)*dirac_sh2_lmax12;  %%SH representation
dirac_sh_st_lmax12 =fod_stand(dirac_sh_lmax12); %%standardized to be nonnegative and sum =1

coe_sh1_lmax16 = Dirac_SH_coe(lmax16,theta0(1),phi0(1)); %% SH coefficients
coe_sh2_lmax16 = Dirac_SH_coe(lmax16,theta0(2),phi0(2));
dirac_sh1_lmax16 = SH_J5_all_lmax16*coe_sh1_lmax16;
dirac_sh2_lmax16 = SH_J5_all_lmax16*coe_sh2_lmax16;
dirac_sh_lmax16 = weight(1)*dirac_sh1_lmax16+weight(2)*dirac_sh2_lmax16;  %%SH representation
dirac_sh_st_lmax16=fod_stand(dirac_sh_lmax16); %%standardized to be nonnegative and sum =1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% analyse replicates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df_SH_all_lmax8 = zeros(N_rep,length(lambda_seq));
RSS_SH_all_lmax8 = zeros(N_rep,length(lambda_seq));
BIC_SH_all_lmax8 = zeros(N_rep,length(lambda_seq));

df_SN_rank_all = zeros(N_rep,length(lambda_seq_la));
RSS_SN_all = zeros(N_rep,length(lambda_seq_la));

df_SN = zeros(N_rep,1);


HD_SH_lmax8 = zeros(N_rep,1);
HD_sCSD_lmax8 = zeros(N_rep,1);
HD_sCSD_lmax12 = zeros(N_rep,1);
HD_sCSD_lmax16 = zeros(N_rep,1);
HD_SN = zeros(N_rep,1);

FOD_SH_lmax8 = zeros(N_rep,2562);
FOD_sCSD_lmax8_all = zeros(N_rep,2562);
FOD_sCSD_lmax12_all = zeros(N_rep,2562);
FOD_sCSD_lmax16_all = zeros(N_rep,2562);
FOD_SN_all = zeros(N_rep,2562);

nfib_SH_lmax8 = zeros(N_rep,1);
nfib_sCSD_lmax8 = zeros(N_rep,1);
nfib_sCSD_lmax12 = zeros(N_rep,1);
nfib_sCSD_lmax16 = zeros(N_rep,1);
nfib_SN = zeros(N_rep,1);

angle1_SH_lmax8 = zeros(N_rep,1);
angle1_sCSD_lmax8 = zeros(N_rep,1);
angle1_sCSD_lmax12 = zeros(N_rep,1);
angle1_sCSD_lmax16 = zeros(N_rep,1);
angle1_SN = zeros(N_rep,1);

angle2_SH_lmax8 = zeros(N_rep,1);
angle2_sCSD_lmax8 = zeros(N_rep,1);
angle2_sCSD_lmax12 = zeros(N_rep,1);
angle2_sCSD_lmax16 = zeros(N_rep,1);
angle2_SN = zeros(N_rep,1);

theta1_SH8 = repmat(-99,N_rep,1);
phi1_SH8 = repmat(-99,N_rep,1);

theta1_sCSD8 = repmat(-99,N_rep,1);
phi1_sCSD8 = repmat(-99,N_rep,1);

theta1_sCSD12 = repmat(-99,N_rep,1);
phi1_sCSD12 = repmat(-99,N_rep,1);

theta1_sCSD16 = repmat(-99,N_rep,1);
phi1_sCSD16 = repmat(-99,N_rep,1);

theta1_SN = repmat(-99,N_rep,1);
phi1_SN = repmat(-99,N_rep,1);

theta2_SH8 = repmat(-99,N_rep,1);
phi2_SH8 = repmat(-99,N_rep,1);

theta2_sCSD8 = repmat(-99,N_rep,1);
phi2_sCSD8 = repmat(-99,N_rep,1);

theta2_sCSD12 = repmat(-99,N_rep,1);
phi2_sCSD12 = repmat(-99,N_rep,1);

theta2_sCSD16 = repmat(-99,N_rep,1);
phi2_sCSD16 = repmat(-99,N_rep,1);

theta2_SN = repmat(-99,N_rep,1);
phi2_SN = repmat(-99,N_rep,1);

sep_SH_lmax8 = zeros(N_rep,1);
sep_sCSD_lmax8 = zeros(N_rep,1);
sep_sCSD_lmax12 = zeros(N_rep,1);
sep_sCSD_lmax16 = zeros(N_rep,1);
sep_SN = zeros(N_rep,1);

SN_stop_index = zeros(N_rep,1);

Dis = squareform(pdist(pos_p','cosine'));
% fod1_s_rota = [fod1_s(1) fod1_s(3) 0];
% rota_M = rotation_matrix(fod1_s, fod1_s_rota);

%%%%%%%%%%%%%%%%%%%%%%%%%%
kmin = 40;
peak_thresh = 0.25;

for rep = 1:N_rep
    
    display(rep);
    
    folder_path = strcat(path_save,'simulation_review/','2fib_sep',num2str(round(sep*180/pi,0)),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'/');
    temp_name = strcat('2fib_sep',num2str(round(sep*180/pi,0)),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'_rep',num2str(rep),'.mat');
    filename = strcat(folder_path, temp_name);
    load(filename);
    
    %%%%% record old stop index
    SN_stop_index(rep) = index_selected_SN;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    df_SH_all_lmax8(rep,:) = df_all_SH_lmax8;
    RSS_SH_all_lmax8(rep,:) = RSS_all_SH_lmax8;
    BIC_SH_all_lmax8(rep,:) = BIC_all_SH_lmax8;
    
    df_SN_rank_all(rep,:) = df_rank_all_SN;
    RSS_SN_all(rep,:) = RSS_all_SN;
    
    df_SN(rep) = df_SN_rank_all(rep,index_selected_SN);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% retrieve the fiber directions and the random rotation
    s = RandStream('mcg16807','Seed',rep*10+21);
    %RandStream.setDefaultStream(s); %%for PC: set random number generator seed
    RandStream.setGlobalStream(s); %%for server: set random number generator seed
    temp=RandOrthMat(3); %%% uniformly generate a 3 by 3 orthogonal matrix
    %%rotate the FODs such that they are aligned with x- and z- directions
    for i = 1: size(pos_p,2)
        pos_c=pos_p(:,i);
        pos_tr=pos_c'*temp;
        inner_c=pos_tr*pos_p;
        index_c=find(inner_c==max(inner_c));
        index_c=index_c(1);
        FOD_SH_lmax8(rep,i) = FOD_SH_est_lmax8(index_c);
        FOD_sCSD_lmax8_all(rep,i) = FOD_sCSD_lmax8(index_c);
        FOD_sCSD_lmax12_all(rep,i) = FOD_sCSD_lmax12(index_c);
        FOD_sCSD_lmax16_all(rep,i) = FOD_sCSD_lmax16(index_c);
        FOD_SN_all(rep,i) = FOD_SN(index_c);
    end
    
    FOD_SH_temp_st = fod_stand(FOD_SH_lmax8(rep,:));
    FOD_lmax8_sCSD_temp_st = fod_stand(FOD_sCSD_lmax8_all(rep,:));
    FOD_lmax12_sCSD_temp_st = fod_stand(FOD_sCSD_lmax12_all(rep,:));
    FOD_lmax16_sCSD_temp_st = fod_stand(FOD_sCSD_lmax16_all(rep,:));
    FOD_SN_temp_st = fod_stand(FOD_SN_all(rep,:));
    
    % Hellinger distance
    HD_SH_lmax8(rep) = hellinger_dis(dirac_sh_st_lmax8, FOD_SH_temp_st);
    HD_sCSD_lmax8(rep) = hellinger_dis(dirac_sh_st_lmax8, FOD_lmax8_sCSD_temp_st);
    HD_sCSD_lmax12(rep) = hellinger_dis(dirac_sh_st_lmax8, FOD_lmax12_sCSD_temp_st);
    HD_sCSD_lmax16(rep) = hellinger_dis(dirac_sh_st_lmax8, FOD_lmax16_sCSD_temp_st);
    HD_SN(rep) = hellinger_dis(dirac_sh_st_lmax8, FOD_SN_temp_st);

    % peak detection
    [~, ~, ~, ~, ~, peak_pos_SH_final_lmax8] = FOD_peak(FOD_SH_lmax8(rep,:), Dis, kmin, peak_thresh, pos_p, theta_p, phi_p);
    [~, ~, ~, ~, ~, peak_pos_sCSD_final_lmax8] = FOD_peak(FOD_sCSD_lmax8_all(rep,:), Dis, kmin, peak_thresh, pos_p, theta_p, phi_p);
    [~, ~, ~, ~, ~, peak_pos_sCSD_final_lmax12] = FOD_peak(FOD_sCSD_lmax12_all(rep,:), Dis, kmin, peak_thresh, pos_p, theta_p, phi_p);
    [~, ~, ~, ~, ~, peak_pos_sCSD_final_lmax16] = FOD_peak(FOD_sCSD_lmax16_all(rep,:), Dis, kmin, peak_thresh, pos_p, theta_p, phi_p);
    [~, ~, ~, ~, ~, peak_pos_SN_final] = FOD_peak(FOD_SN_all(rep,:), Dis, kmin, peak_thresh, pos_p, theta_p, phi_p);

    nfib_SH_lmax8(rep) = size(peak_pos_SH_final_lmax8,2);
    nfib_sCSD_lmax8(rep) = size(peak_pos_sCSD_final_lmax8,2);
    nfib_sCSD_lmax12(rep) = size(peak_pos_sCSD_final_lmax12,2);
    nfib_sCSD_lmax16(rep) = size(peak_pos_sCSD_final_lmax16,2);
    nfib_SN(rep) = size(peak_pos_SN_final,2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% SH lmax8                                                    %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(nfib_SH_lmax8(rep)==2)
        fod_temp = [fod1_s; fod2_s];
        
        ang11 = atan2(norm(cross(fod1_s,peak_pos_SH_final_lmax8(:,1))), dot(fod1_s,peak_pos_SH_final_lmax8(:,1)));
        ang12 = abs(pi - atan2(norm(cross(fod1_s,peak_pos_SH_final_lmax8(:,1))), dot(fod1_s,peak_pos_SH_final_lmax8(:,1))));
        ang111 = min(ang11,ang12);
        ang21 = atan2(norm(cross(fod2_s,peak_pos_SH_final_lmax8(:,1))), dot(fod2_s,peak_pos_SH_final_lmax8(:,1)));
        ang22 = abs(pi - atan2(norm(cross(fod2_s,peak_pos_SH_final_lmax8(:,1))), dot(fod2_s,peak_pos_SH_final_lmax8(:,1))));
        ang222 = min(ang21,ang22);
        [angle1_SH_lmax8(rep,1), index_temp] = min([ang111, ang222]);
        %%%%%
        if(dot(fod_temp(index_temp,:)',peak_pos_SH_final_lmax8(:,1))<0)
            peak_pos_SH_final_lmax8(:,1) = -peak_pos_SH_final_lmax8(:,1);
        end
        phi_temp = atan2(peak_pos_SH_final_lmax8(2,1), peak_pos_SH_final_lmax8(1,1));   %%phi: azimuthal  angle, [0,2pi)
        phi_temp = phi_temp+(phi_temp<0)*2*pi;
        theta_temp = acos(peak_pos_SH_final_lmax8(3,1));
        if(index_temp == 1)
            phi1_SH8(rep) = phi_temp;
            theta1_SH8(rep) = theta_temp;
        else
            phi2_SH8(rep) = phi_temp;
            theta2_SH8(rep) = theta_temp;
        end
        %%%%%
        fod_temp(index_temp,:) = [];
        ang11 = atan2(norm(cross(fod_temp(1,:),peak_pos_SH_final_lmax8(:,2))), dot(fod_temp(1,:),peak_pos_SH_final_lmax8(:,2)));
        ang12 = abs(pi - atan2(norm(cross(fod_temp(1,:),peak_pos_SH_final_lmax8(:,2))), dot(fod_temp(1,:),peak_pos_SH_final_lmax8(:,2))));
        ang111 = min(ang11,ang12);
        angle2_SH_lmax8(rep,1) = ang111;
        ang11 = atan2(norm(cross(peak_pos_SH_final_lmax8(:,1),peak_pos_SH_final_lmax8(:,2))), dot(peak_pos_SH_final_lmax8(:,1),peak_pos_SH_final_lmax8(:,2)));
        ang12 = abs(pi - atan2(norm(cross(peak_pos_SH_final_lmax8(:,1),peak_pos_SH_final_lmax8(:,2))), dot(peak_pos_SH_final_lmax8(:,1),peak_pos_SH_final_lmax8(:,2))));
        sep_SH_lmax8(rep,1) = min(ang11, ang12);
        %%%%%
        if(dot(fod_temp(1,:)',peak_pos_SH_final_lmax8(:,2))<0)
            peak_pos_SH_final_lmax8(:,2) = -peak_pos_SH_final_lmax8(:,2);
        end
        phi_temp = atan2(peak_pos_SH_final_lmax8(2,2), peak_pos_SH_final_lmax8(1,2));   %%phi: azimuthal  angle, [0,2pi)
        phi_temp = phi_temp+(phi_temp<0)*2*pi;
        theta_temp = acos(peak_pos_SH_final_lmax8(3,2));
        if(index_temp == 1)
            phi2_SH8(rep) = phi_temp;
            theta2_SH8(rep) = theta_temp;
        else
            phi1_SH8(rep) = phi_temp;
            theta1_SH8(rep) = theta_temp;
        end
        %%%%%
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% sCSD lmax8                                                  %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(nfib_sCSD_lmax8(rep)==2)
        fod_temp = [fod1_s; fod2_s];
        ang11 = atan2(norm(cross(fod1_s,peak_pos_sCSD_final_lmax8(:,1))), dot(fod1_s,peak_pos_sCSD_final_lmax8(:,1)));
        ang12 = abs(pi - atan2(norm(cross(fod1_s,peak_pos_sCSD_final_lmax8(:,1))), dot(fod1_s,peak_pos_sCSD_final_lmax8(:,1))));
        ang111 = min(ang11,ang12);
        ang21 = atan2(norm(cross(fod2_s,peak_pos_sCSD_final_lmax8(:,1))), dot(fod2_s,peak_pos_sCSD_final_lmax8(:,1)));
        ang22 = abs(pi - atan2(norm(cross(fod2_s,peak_pos_sCSD_final_lmax8(:,1))), dot(fod2_s,peak_pos_sCSD_final_lmax8(:,1))));
        ang222 = min(ang21,ang22);
        [angle1_sCSD_lmax8(rep,1), index_temp] = min([ang111, ang222]);
        %%%%%
        if(dot(fod_temp(index_temp,:)',peak_pos_sCSD_final_lmax8(:,1))<0)
            peak_pos_sCSD_final_lmax8(:,1) = -peak_pos_sCSD_final_lmax8(:,1);
        end
        phi_temp = atan2(peak_pos_sCSD_final_lmax8(2,1), peak_pos_sCSD_final_lmax8(1,1));   %%phi: azimuthal  angle, [0,2pi)
        phi_temp = phi_temp+(phi_temp<0)*2*pi;
        theta_temp = acos(peak_pos_sCSD_final_lmax8(3,1));
        if(index_temp == 1)
            phi1_sCSD8(rep) = phi_temp;
            theta1_sCSD8(rep) = theta_temp;
        else
            phi2_sCSD8(rep) = phi_temp;
            theta2_sCSD8(rep) = theta_temp;
        end
        %%%%%
        fod_temp(index_temp,:) = [];
        ang11 = atan2(norm(cross(fod_temp(1,:),peak_pos_sCSD_final_lmax8(:,2))), dot(fod_temp(1,:),peak_pos_sCSD_final_lmax8(:,2)));
        ang12 = abs(pi - atan2(norm(cross(fod_temp(1,:),peak_pos_sCSD_final_lmax8(:,2))), dot(fod_temp(1,:),peak_pos_sCSD_final_lmax8(:,2))));
        ang111 = min(ang11,ang12);
        angle2_sCSD_lmax8(rep,1) = ang111;
        ang11 = atan2(norm(cross(peak_pos_sCSD_final_lmax8(:,1),peak_pos_sCSD_final_lmax8(:,2))), dot(peak_pos_sCSD_final_lmax8(:,1),peak_pos_sCSD_final_lmax8(:,2)));
        ang12 = abs(pi - atan2(norm(cross(peak_pos_sCSD_final_lmax8(:,1),peak_pos_sCSD_final_lmax8(:,2))), dot(peak_pos_sCSD_final_lmax8(:,1),peak_pos_sCSD_final_lmax8(:,2))));
        sep_sCSD_lmax8(rep,1) = min(ang11, ang12);
        %%%%%
        if(dot(fod_temp(1,:)',peak_pos_sCSD_final_lmax8(:,2))<0)
            peak_pos_sCSD_final_lmax8(:,2) = -peak_pos_sCSD_final_lmax8(:,2);
        end
        phi_temp = atan2(peak_pos_sCSD_final_lmax8(2,2), peak_pos_sCSD_final_lmax8(1,2));   %%phi: azimuthal  angle, [0,2pi)
        phi_temp = phi_temp+(phi_temp<0)*2*pi;
        theta_temp = acos(peak_pos_sCSD_final_lmax8(3,2));
        if(index_temp == 1)
            phi2_sCSD8(rep) = phi_temp;
            theta2_sCSD8(rep) = theta_temp;
        else
            phi1_sCSD8(rep) = phi_temp;
            theta1_sCSD8(rep) = theta_temp;
        end
        %%%%%
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% sCSD lmax12                                                 %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(nfib_sCSD_lmax12(rep)==2)
        fod_temp = [fod1_s; fod2_s];
        ang11 = atan2(norm(cross(fod1_s,peak_pos_sCSD_final_lmax12(:,1))), dot(fod1_s,peak_pos_sCSD_final_lmax12(:,1)));
        ang12 = abs(pi - atan2(norm(cross(fod1_s,peak_pos_sCSD_final_lmax12(:,1))), dot(fod1_s,peak_pos_sCSD_final_lmax12(:,1))));
        ang111 = min(ang11,ang12);
        ang21 = atan2(norm(cross(fod2_s,peak_pos_sCSD_final_lmax12(:,1))), dot(fod2_s,peak_pos_sCSD_final_lmax12(:,1)));
        ang22 = abs(pi - atan2(norm(cross(fod2_s,peak_pos_sCSD_final_lmax12(:,1))), dot(fod2_s,peak_pos_sCSD_final_lmax12(:,1))));
        ang222 = min(ang21,ang22);
        [angle1_sCSD_lmax12(rep,1), index_temp] = min([ang111, ang222]);
        %%%%%
        if(dot(fod_temp(index_temp,:)',peak_pos_sCSD_final_lmax12(:,1))<0)
            peak_pos_sCSD_final_lmax12(:,1) = -peak_pos_sCSD_final_lmax12(:,1);
        end
        phi_temp = atan2(peak_pos_sCSD_final_lmax12(2,1), peak_pos_sCSD_final_lmax12(1,1));   %%phi: azimuthal  angle, [0,2pi)
        phi_temp = phi_temp+(phi_temp<0)*2*pi;
        theta_temp = acos(peak_pos_sCSD_final_lmax12(3,1));
        if(index_temp == 1)
            phi1_sCSD12(rep) = phi_temp;
            theta1_sCSD12(rep) = theta_temp;
        else
            phi2_sCSD12(rep) = phi_temp;
            theta2_sCSD12(rep) = theta_temp;
        end
        %%%%%
        fod_temp(index_temp,:) = [];
        ang11 = atan2(norm(cross(fod_temp(1,:),peak_pos_sCSD_final_lmax12(:,2))), dot(fod_temp(1,:),peak_pos_sCSD_final_lmax12(:,2)));
        ang12 = abs(pi - atan2(norm(cross(fod_temp(1,:),peak_pos_sCSD_final_lmax12(:,2))), dot(fod_temp(1,:),peak_pos_sCSD_final_lmax12(:,2))));
        ang111 = min(ang11,ang12);
        angle2_sCSD_lmax12(rep,1) = ang111;
        ang11 = atan2(norm(cross(peak_pos_sCSD_final_lmax12(:,1),peak_pos_sCSD_final_lmax12(:,2))), dot(peak_pos_sCSD_final_lmax12(:,1),peak_pos_sCSD_final_lmax12(:,2)));
        ang12 = abs(pi - atan2(norm(cross(peak_pos_sCSD_final_lmax12(:,1),peak_pos_sCSD_final_lmax12(:,2))), dot(peak_pos_sCSD_final_lmax12(:,1),peak_pos_sCSD_final_lmax12(:,2))));
        sep_sCSD_lmax12(rep,1) = min(ang11, ang12);
        %%%%%
        if(dot(fod_temp(1,:)',peak_pos_sCSD_final_lmax12(:,2))<0)
            peak_pos_sCSD_final_lmax12(:,2) = -peak_pos_sCSD_final_lmax12(:,2);
        end
        phi_temp = atan2(peak_pos_sCSD_final_lmax12(2,2), peak_pos_sCSD_final_lmax12(1,2));   %%phi: azimuthal  angle, [0,2pi)
        phi_temp = phi_temp+(phi_temp<0)*2*pi;
        theta_temp = acos(peak_pos_sCSD_final_lmax12(3,2));
        if(index_temp == 1)
            phi2_sCSD12(rep) = phi_temp;
            theta2_sCSD12(rep) = theta_temp;
        else
            phi1_sCSD12(rep) = phi_temp;
            theta1_sCSD12(rep) = theta_temp;
        end
        %%%%%
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% sCSD lmax16                                                 %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(nfib_sCSD_lmax16(rep)==2)
        fod_temp = [fod1_s; fod2_s];
        ang11 = atan2(norm(cross(fod1_s,peak_pos_sCSD_final_lmax16(:,1))), dot(fod1_s,peak_pos_sCSD_final_lmax16(:,1)));
        ang12 = abs(pi - atan2(norm(cross(fod1_s,peak_pos_sCSD_final_lmax16(:,1))), dot(fod1_s,peak_pos_sCSD_final_lmax16(:,1))));
        ang111 = min(ang11,ang12);
        ang21 = atan2(norm(cross(fod2_s,peak_pos_sCSD_final_lmax16(:,1))), dot(fod2_s,peak_pos_sCSD_final_lmax16(:,1)));
        ang22 = abs(pi - atan2(norm(cross(fod2_s,peak_pos_sCSD_final_lmax16(:,1))), dot(fod2_s,peak_pos_sCSD_final_lmax16(:,1))));
        ang222 = min(ang21,ang22);
        [angle1_sCSD_lmax16(rep,1), index_temp] = min([ang111, ang222]);
        %%%%%
        if(dot(fod_temp(index_temp,:)',peak_pos_sCSD_final_lmax16(:,1))<0)
            peak_pos_sCSD_final_lmax16(:,1) = -peak_pos_sCSD_final_lmax16(:,1);
        end
        phi_temp = atan2(peak_pos_sCSD_final_lmax16(2,1), peak_pos_sCSD_final_lmax16(1,1));   %%phi: azimuthal  angle, [0,2pi)
        phi_temp = phi_temp+(phi_temp<0)*2*pi;
        theta_temp = acos(peak_pos_sCSD_final_lmax16(3,1));
        if(index_temp == 1)
            phi1_sCSD16(rep) = phi_temp;
            theta1_sCSD16(rep) = theta_temp;
        else
            phi2_sCSD16(rep) = phi_temp;
            theta2_sCSD16(rep) = theta_temp;
        end
        %%%%%
        fod_temp(index_temp,:) = [];
        ang11 = atan2(norm(cross(fod_temp(1,:),peak_pos_sCSD_final_lmax16(:,2))), dot(fod_temp(1,:),peak_pos_sCSD_final_lmax16(:,2)));
        ang12 = abs(pi - atan2(norm(cross(fod_temp(1,:),peak_pos_sCSD_final_lmax16(:,2))), dot(fod_temp(1,:),peak_pos_sCSD_final_lmax16(:,2))));
        ang111 = min(ang11,ang12);
        angle2_sCSD_lmax16(rep,1) = ang111;
        ang11 = atan2(norm(cross(peak_pos_sCSD_final_lmax16(:,1),peak_pos_sCSD_final_lmax16(:,2))), dot(peak_pos_sCSD_final_lmax16(:,1),peak_pos_sCSD_final_lmax16(:,2)));
        ang12 = abs(pi - atan2(norm(cross(peak_pos_sCSD_final_lmax16(:,1),peak_pos_sCSD_final_lmax16(:,2))), dot(peak_pos_sCSD_final_lmax16(:,1),peak_pos_sCSD_final_lmax16(:,2))));
        sep_sCSD_lmax16(rep,1) = min(ang11, ang12);
        %%%%%
        if(dot(fod_temp(1,:)',peak_pos_sCSD_final_lmax16(:,2))<0)
            peak_pos_sCSD_final_lmax16(:,2) = -peak_pos_sCSD_final_lmax16(:,2);
        end
        phi_temp = atan2(peak_pos_sCSD_final_lmax16(2,2), peak_pos_sCSD_final_lmax16(1,2));   %%phi: azimuthal  angle, [0,2pi)
        phi_temp = phi_temp+(phi_temp<0)*2*pi;
        theta_temp = acos(peak_pos_sCSD_final_lmax16(3,2));
        if(index_temp == 1)
            phi2_sCSD16(rep) = phi_temp;
            theta2_sCSD16(rep) = theta_temp;
        else
            phi1_sCSD16(rep) = phi_temp;
            theta1_sCSD16(rep) = theta_temp;
        end
        %%%%%
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% SN RSSdiff                                                  %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(nfib_SN(rep)==2)
        fod_temp = [fod1_s;fod2_s];
        ang11 = atan2(norm(cross(fod1_s,peak_pos_SN_final(:,1))), dot(fod1_s,peak_pos_SN_final(:,1)));
        ang12 = abs(pi - atan2(norm(cross(fod1_s,peak_pos_SN_final(:,1))), dot(fod1_s,peak_pos_SN_final(:,1))));
        ang111 = min(ang11,ang12);
        ang21 = atan2(norm(cross(fod2_s,peak_pos_SN_final(:,1))), dot(fod2_s,peak_pos_SN_final(:,1)));
        ang22 = abs(pi - atan2(norm(cross(fod2_s,peak_pos_SN_final(:,1))), dot(fod2_s,peak_pos_SN_final(:,1))));
        ang222 = min(ang21,ang22);
        [angle1_SN(rep,1), index_temp] = min([ang111, ang222]);
        %%%%%
        if(dot(fod_temp(index_temp,:)',peak_pos_SN_final(:,1))<0)
            peak_pos_SN_final(:,1) = -peak_pos_SN_final(:,1);
        end
        phi_temp = atan2(peak_pos_SN_final(2,1), peak_pos_SN_final(1,1));   %%phi: azimuthal  angle, [0,2pi)
        phi_temp = phi_temp+(phi_temp<0)*2*pi;
        theta_temp = acos(peak_pos_SN_final(3,1));
        if(index_temp == 1)
            phi1_SN(rep) = phi_temp;
            theta1_SN(rep) = theta_temp;
        else
            phi2_SN(rep) = phi_temp;
            theta2_SN(rep) = theta_temp;
        end
        %%%%%
        fod_temp(index_temp,:) = [];
        ang11 = atan2(norm(cross(fod_temp(1,:),peak_pos_SN_final(:,2))), dot(fod_temp(1,:),peak_pos_SN_final(:,2)));
        ang12 = abs(pi - atan2(norm(cross(fod_temp(1,:),peak_pos_SN_final(:,2))), dot(fod_temp(1,:),peak_pos_SN_final(:,2))));
        ang111 = min(ang11,ang12);
        angle2_SN(rep,1) = ang111;
        ang11 = atan2(norm(cross(peak_pos_SN_final(:,1),peak_pos_SN_final(:,2))), dot(peak_pos_SN_final(:,1),peak_pos_SN_final(:,2)));
        ang12 = abs(pi - atan2(norm(cross(peak_pos_SN_final(:,1),peak_pos_SN_final(:,2))), dot(peak_pos_SN_final(:,1),peak_pos_SN_final(:,2))));
        sep_SN(rep,1) = min(ang11, ang12);
        %%%%%
        if(dot(fod_temp(1,:)',peak_pos_SN_final(:,2))<0)
            peak_pos_SN_final(:,2) = -peak_pos_SN_final(:,2);
        end
        phi_temp = atan2(peak_pos_SN_final(2,2), peak_pos_SN_final(1,2));   %%phi: azimuthal  angle, [0,2pi)
        phi_temp = phi_temp+(phi_temp<0)*2*pi;
        theta_temp = acos(peak_pos_SN_final(3,2));
        if(index_temp == 1)
            phi2_SN(rep) = phi_temp;
            theta2_SN(rep) = theta_temp;
        else
            phi1_SN(rep) = phi_temp;
            theta1_SN(rep) = theta_temp;
        end
        %%%%%
    end


end


%% saving path and folder name
save_path = strcat(folder_path,'new_aggregate_thresh',num2str(peak_thresh),'/');
mkdir(save_path);


% BIC_SH = mean(BIC_SH_all,1);
mean_BIC_SH_lmax8 = mean(BIC_SH_all_lmax8,1);
mean_RSS_SN = mean(RSS_SN_all,1);
mean_df_SN_rank = mean(df_SN_rank_all,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(log10(lambda_seq),mean_BIC_SH_lmax8);
title('Harmonics: BIC vs lambda');
savefig(strcat(save_path,'BIC_SH.fig'));

figure;
subplot(1,2,1)
plot(log10(lambda_seq_la),mean_df_SN_rank);
title('Needlets: rank df vs lambda');
subplot(1,2,2)
plot(log10(lambda_seq_la),mean_RSS_SN);
title('Needlets: RSS vs lambda');
%savefig(strcat(save_path,'Model_Selection_SN.fig'));

mean_HD_SH_lmax8 = mean(HD_SH_lmax8);
mean_HD_sCSD_lmax8 = mean(HD_sCSD_lmax8);
mean_HD_sCSD_lmax12 = mean(HD_sCSD_lmax12);
mean_HD_sCSD_lmax16 = mean(HD_sCSD_lmax16);
mean_HD_SN  = mean(HD_SN);

display(strcat('Mean HD SH lmax8 = ', num2str(mean_HD_SH_lmax8)));
display(strcat('Mean HD sCSD lmax8 = ', num2str(mean_HD_sCSD_lmax8)));
display(strcat('Mean HD sCSD lmax12 = ', num2str(mean_HD_sCSD_lmax12)));
display(strcat('Mean HD sCSD lmax16 = ', num2str(mean_HD_sCSD_lmax16)));
display(strcat('Mean HD SN RSSdiff = ', num2str(mean_HD_SN)));

median_HD_SH_lmax8 = median(HD_SH_lmax8);
median_HD_sCSD_lmax8 = median(HD_sCSD_lmax8);
median_HD_sCSD_lmax12 = median(HD_sCSD_lmax12);
median_HD_sCSD_lmax16 = median(HD_sCSD_lmax16);
median_HD_SN = median(HD_SN);

display(strcat('Median HD SH lmax8 = ', num2str(median_HD_SH_lmax8)));
display(strcat('Median HD sCSD lmax8 = ', num2str(median_HD_sCSD_lmax8)));
display(strcat('Median HD sCSD lmax12 = ', num2str(median_HD_sCSD_lmax12)));
display(strcat('Median HD sCSD lmax16 = ', num2str(median_HD_sCSD_lmax16)));
display(strcat('Median HD SN RSSdiff = ', num2str(median_HD_SN )));

sd_HD_SH_lmax8 = std(HD_SH_lmax8);
sd_HD_sCSD_lmax8 = std(HD_sCSD_lmax8);
sd_HD_sCSD_lmax12 = std(HD_sCSD_lmax12);
sd_HD_sCSD_lmax16 = std(HD_sCSD_lmax16);
sd_HD_SN  = std(HD_SN );

display(strcat('STD HD SH lmax8 = ', num2str(sd_HD_SH_lmax8)));
display(strcat('STD HD sCSD lmax8 = ', num2str(sd_HD_sCSD_lmax8)));
display(strcat('STD HD sCSD lmax12 = ', num2str(sd_HD_sCSD_lmax12)));
display(strcat('STD HD sCSD lmax16 = ', num2str(sd_HD_sCSD_lmax16)));
display(strcat('STD HD SN RSSdiff = ', num2str(sd_HD_SN )));

%%% get the average/median FOD estimation across 100 replicates.
FOD_SH_med_lmax8=median(FOD_SH_lmax8,1);
FOD_sCSD_med_lmax8=median(FOD_sCSD_lmax8_all,1);
FOD_sCSD_med_lmax12=median(FOD_sCSD_lmax12_all,1);
FOD_sCSD_med_lmax16=median(FOD_sCSD_lmax16_all,1);
FOD_SN_med=median(FOD_SN,1);

FOD_SH_mean_lmax8=mean(FOD_SH_lmax8,1);
FOD_sCSD_mean_lmax8=mean(FOD_sCSD_lmax8_all,1);
FOD_sCSD_mean_lmax12=mean(FOD_sCSD_lmax12_all,1);
FOD_sCSD_mean_lmax16=mean(FOD_sCSD_lmax16_all,1);
FOD_SN_mean=mean(FOD_SN_all,1);

FOD_SH_sd_lmax8=std(FOD_SH_lmax8,1);
FOD_sCSD_sd_lmax8=std(FOD_sCSD_lmax8_all,1);
FOD_sCSD_sd_lmax12=std(FOD_sCSD_lmax12_all,1);
FOD_sCSD_sd_lmax16=std(FOD_sCSD_lmax16_all,1);
FOD_SN_sd=std(FOD_SN_all,1);

std_multi = 2;
FOD_SH_pop_lmax8 = FOD_SH_mean_lmax8+std_multi.*FOD_SH_sd_lmax8;
FOD_sCSD_pop_lmax8 = FOD_sCSD_mean_lmax8+std_multi.*FOD_sCSD_sd_lmax8;
FOD_sCSD_pop_lmax12 = FOD_sCSD_mean_lmax12+std_multi.*FOD_sCSD_sd_lmax12;
FOD_sCSD_pop_lmax16 = FOD_sCSD_mean_lmax16+std_multi.*FOD_sCSD_sd_lmax16;
FOD_SN_pop = FOD_SN_mean+std_multi.*FOD_SN_sd;

%{
figure;
hold all;
plot_spherical_function(v_p,f_p,dirac_sh_lmax8,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(dirac_sh_lmax8));
title(strcat('SH representation lmax',num2str(8)));
savefig(strcat(save_path,'SH_rep.fig'));

figure;
hold all;
plot_spherical_function(v_p,f_p,dirac_sh_lmax12,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(dirac_sh_lmax12));
title(strcat('SH representation lmax',num2str(12)));
savefig(strcat(save_path,'SH12_rep.fig'));

figure;
hold all;
plot_spherical_function(v_p,f_p,dirac_sh_lmax16,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(dirac_sh_lmax16));
title(strcat('SH representation lmax',num2str(16)));
savefig(strcat(save_path,'SH16_rep.fig'));

figure;
hold all;
% subplot(2,2,1)
plot_spherical_function(v_p,f_p,FOD_SH_pop_lmax8,options);
alpha(0.25)
hold on;
plot_spherical_function(v_p,f_p,FOD_SH_mean_lmax8,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_SH_pop_lmax8));
title('Average SH estimation & +2STD');
savefig(strcat(save_path,'SH_est_lmax8.fig'));


figure;
hold all;
% subplot(2,2,2)
plot_spherical_function(v_p,f_p,FOD_sCSD_pop_lmax8,options);
alpha(0.25)
hold on;
plot_spherical_function(v_p,f_p,FOD_sCSD_mean_lmax8,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_sCSD_pop_lmax8));
title('Average SH+sCSD (lmax8) estimation & +2STD');
savefig(strcat(save_path,'sCSD_est_lmax8.fig'));

figure;
hold all;
% subplot(2,2,2)
plot_spherical_function(v_p,f_p,FOD_sCSD_pop_lmax12,options);
alpha(0.25)
hold on;
plot_spherical_function(v_p,f_p,FOD_sCSD_mean_lmax12,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_sCSD_pop_lmax12));
title('Average SH+sCSD (lmax12) estimation & +2STD');
savefig(strcat(save_path,'sCSD_est_lmax12.fig'));

figure;
hold all;
% subplot(2,2,2)
plot_spherical_function(v_p,f_p,FOD_sCSD_pop_lmax16,options);
alpha(0.25)
hold on;
plot_spherical_function(v_p,f_p,FOD_sCSD_mean_lmax16,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_sCSD_pop_lmax16));
title('Average SH+sCSD (lmax16) estimation & +2STD');
savefig(strcat(save_path,'sCSD_est_lmax16.fig'));

figure;
hold all;
% subplot(2,2,3)
plot_spherical_function(v_p,f_p,FOD_SN_pop,options);
alpha(0.25)
hold on;
plot_spherical_function(v_p,f_p,FOD_SN_mean,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_SN_pop));
title('Average SN RSSdiff estimation & +2STD');
savefig(strcat(save_path,'SN_est.fig'));
%}

%{
figure;
plot(theta0,phi0,'.','Color','red','MarkerSize', 20);
hold on;
plot(theta_SH_lmax8, phi_SH_lmax8,'.','Color','blue','MarkerSize', 15);
savefig(strcat(save_path,'SH_peak_angle.fig'));
figure;
plot(theta0,phi0,'.','Color','red','MarkerSize', 20);
hold on;
plot(theta_sCSD_lmax8,phi_sCSD_lmax8,'.','Color','blue','MarkerSize', 15);
savefig(strcat(save_path,'sCSD_lmax8_peak_angle.fig'));
figure;
plot(theta0,phi0,'.','Color','red','MarkerSize', 20);
hold on;
plot(theta_sCSD_lmax12,phi_sCSD_lmax12,'.','Color','blue','MarkerSize', 15);
savefig(strcat(save_path,'sCSD_lmax12_peak_angle.fig'));
figure;
plot(theta0,phi0,'.','Color','red','MarkerSize', 20);
hold on;
plot(theta_sCSD_lmax16,phi_sCSD_lmax16,'.','Color','blue','MarkerSize', 15);
savefig(strcat(save_path,'sCSD_lmax16_peak_angle.fig'));
figure;
plot(theta0,phi0,'.','Color','red','MarkerSize', 20);
hold on;
plot(theta_SN,phi_SN,'.','Color','blue','MarkerSize', 15);
savefig(strcat(save_path,'SN_peak_angle.fig'));
%}

%%%%%
rate_SH8 = sum(nfib_SH_lmax8==2)/N_rep;
under_SH8 = sum(nfib_SH_lmax8<2)/N_rep;
over_SH8 = sum(nfib_SH_lmax8>2)/N_rep;

rate_sCSD8 = sum(nfib_sCSD_lmax8==2)/N_rep;
under_sCSD8 = sum(nfib_sCSD_lmax8<2)/N_rep;
over_sCSD8 = sum(nfib_sCSD_lmax8>2)/N_rep;

rate_sCSD12 = sum(nfib_sCSD_lmax12==2)/N_rep;
under_sCSD12 = sum(nfib_sCSD_lmax12<2)/N_rep;
over_sCSD12 = sum(nfib_sCSD_lmax12>2)/N_rep;

rate_sCSD16 = sum(nfib_sCSD_lmax16==2)/N_rep;
under_sCSD16 = sum(nfib_sCSD_lmax16<2)/N_rep;
over_sCSD16 = sum(nfib_sCSD_lmax16>2)/N_rep;

rate_SN = sum(nfib_SN==2)/N_rep;
under_SN = sum(nfib_SN<2)/N_rep;
over_SN = sum(nfib_SN>2)/N_rep;

display(strcat('rate SH lmax8 = ', num2str(rate_SH8)));
display(strcat('under SH lmax8 = ', num2str(under_SH8)));
display(strcat('over SH lmax8 = ', num2str(over_SH8)));
display(strcat('rate sCSD lmax8 = ', num2str(rate_sCSD8)));
display(strcat('under sCSD lmax8 = ', num2str(under_sCSD8)));
display(strcat('over sCSD lmax8 = ', num2str(over_sCSD8)));
display(strcat('rate sCSD lmax12 = ', num2str(rate_sCSD12)));
display(strcat('under sCSD lmax12 = ', num2str(under_sCSD12)));
display(strcat('over sCSD lmax12 = ', num2str(over_sCSD12)));
display(strcat('rate sCSD lmax16 = ', num2str(rate_sCSD16)));
display(strcat('under sCSD lmax16 = ', num2str(under_sCSD16)));
display(strcat('over sCSD lmax16 = ', num2str(over_sCSD16)));
display(strcat('rate SN RSSdiff = ', num2str(rate_SN)));
display(strcat('under SN RSSdiff = ', num2str(under_SN)));
display(strcat('over SN RSSdiff = ', num2str(over_SN)));
%%%%%

mean_angle1_SH_lmax8 = mean(angle1_SH_lmax8(nfib_SH_lmax8==2))*180/pi;
mean_angle1_sCSD_lmax8 = mean(angle1_sCSD_lmax8(nfib_sCSD_lmax8==2))*180/pi;
mean_angle1_sCSD_lmax12 = mean(angle1_sCSD_lmax12(nfib_sCSD_lmax12==2))*180/pi;
mean_angle1_sCSD_lmax16 = mean(angle1_sCSD_lmax16(nfib_sCSD_lmax16==2))*180/pi;
mean_angle1_SN = mean(angle1_SN(nfib_SN==2))*180/pi;

mean_angle2_SH_lmax8 = mean(angle2_SH_lmax8(nfib_SH_lmax8==2))*180/pi;
mean_angle2_sCSD_lmax8 = mean(angle2_sCSD_lmax8(nfib_sCSD_lmax8==2))*180/pi;
mean_angle2_sCSD_lmax12 = mean(angle2_sCSD_lmax12(nfib_sCSD_lmax12==2))*180/pi;
mean_angle2_sCSD_lmax16 = mean(angle2_sCSD_lmax16(nfib_sCSD_lmax16==2))*180/pi;
mean_angle2_SN = mean(angle2_SN(nfib_SN==2))*180/pi;

%%%%%
angle1_std_SH_lmax8 = std(angle1_SH_lmax8(nfib_SH_lmax8==1));
angle1_std_sCSD_lmax8 = std(angle1_sCSD_lmax8(nfib_sCSD_lmax8==1));
angle1_std_sCSD_lmax12 = std(angle1_sCSD_lmax12(nfib_sCSD_lmax12==1));
angle1_std_sCSD_lmax16 = std(angle1_sCSD_lmax16(nfib_sCSD_lmax16==1));
angle1_std_SN = std(angle1_SN(nfib_SN==1));

angle2_std_SH_lmax8 = std(angle2_SH_lmax8(nfib_SH_lmax8==2));
angle2_std_sCSD_lmax8 = std(angle2_sCSD_lmax8(nfib_sCSD_lmax8==2));
angle2_std_sCSD_lmax12 = std(angle2_sCSD_lmax12(nfib_sCSD_lmax12==2));
angle2_std_sCSD_lmax16 = std(angle2_sCSD_lmax16(nfib_sCSD_lmax16==2));
angle2_std_SN = std(angle2_SN(nfib_SN==2));


%%%%%

%%%%%
mean_sep_SH_lmax8 = mean(sep_SH_lmax8(nfib_SH_lmax8==2))*180/pi;
mean_sep_sCSD_lmax8 = mean(sep_sCSD_lmax8(nfib_sCSD_lmax8==2))*180/pi;
mean_sep_sCSD_lmax12 = mean(sep_sCSD_lmax12(nfib_sCSD_lmax12==2))*180/pi;
mean_sep_sCSD_lmax16 = mean(sep_sCSD_lmax16(nfib_sCSD_lmax16==2))*180/pi;
mean_sep_SN = mean(sep_SN(nfib_SN==2))*180/pi;

sd_sep_SH_lmax8 = std(sep_SH_lmax8(nfib_SH_lmax8==2))*180/pi;
sd_sep_sCSD_lmax8 = std(sep_sCSD_lmax8(nfib_sCSD_lmax8==2))*180/pi;
sd_sep_sCSD_lmax12 = std(sep_sCSD_lmax12(nfib_sCSD_lmax12==2))*180/pi;
sd_sep_sCSD_lmax16 = std(sep_sCSD_lmax16(nfib_sCSD_lmax16==2))*180/pi;
sd_sep_SN = std(sep_SN(nfib_SN==2))*180/pi;

%%%%%


figure;
subplot(1,6,1)
hold all;
plot_spherical_function(v_p,f_p,dirac_sh_lmax8,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(dirac_sh_lmax8));
view([0 1 0]);
% camroll(90);
% title('SH lmax8 rep.');

subplot(1,6,2);
hold all;
plot_spherical_function(v_p,f_p,FOD_SH_pop_lmax8,options);
alpha(0.25)
hold on;
plot_spherical_function(v_p,f_p,FOD_SH_mean_lmax8,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_SH_pop_lmax8));
view([0 1 0]);
% camroll(90);
% title('SH lmax8 est.');


subplot(1,6,3);
hold all;
plot_spherical_function(v_p,f_p,FOD_sCSD_pop_lmax8,options);
alpha(0.25)
hold on;
plot_spherical_function(v_p,f_p,FOD_sCSD_mean_lmax8,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_sCSD_pop_lmax8));
view([0 1 0]);
% camroll(90);
% title('sCSD lmax8');


subplot(1,6,4);
hold all;
plot_spherical_function(v_p,f_p,FOD_sCSD_pop_lmax12,options);
alpha(0.25)
hold on;
plot_spherical_function(v_p,f_p,FOD_sCSD_mean_lmax12,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_sCSD_pop_lmax12));
view([0 1 0]);
% camroll(90);
% title('sCSD lmax12');

subplot(1,6,5);
hold all;
plot_spherical_function(v_p,f_p,FOD_sCSD_pop_lmax16,options);
alpha(0.25)
hold on;
plot_spherical_function(v_p,f_p,FOD_sCSD_mean_lmax16,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_sCSD_pop_lmax16));
view([0 1 0]);
% camroll(90);
% title('sCSD lmax16');

subplot(1,6,6)
hold all;
plot_spherical_function(v_p,f_p,FOD_SN_pop,options);
alpha(0.25);
hold on;
plot_spherical_function(v_p,f_p,FOD_SN_mean,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_SN_pop));
view([0 1 0]);
% camroll(90);
% title('SN RSSdiff rep.');

% savefig(strcat(save_path,'All_est_fib.fig'));
savefig(strcat(save_path,'2fib_sep',num2str(sep*180/pi),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_N',num2str(n_sample),'.fig'));


%{
subplot(1,8,5);
hold all;
% subplot(2,2,2)
plot_spherical_function(v_p,f_p,FOD_sCSD_pop_lmax16,options);
alpha(0.25)
hold on;
plot_spherical_function(v_p,f_p,FOD_sCSD_mean_lmax16,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_sCSD_pop_lmax16));
view([0 1 0]);
% camroll(90);
title('SH lmax16 rep.');


subplot(1,8,6)
hold all;
plot_spherical_function(v_p,f_p,FOD_SN_BIC_rank_pop,options);
alpha(0.25);
hold on;
plot_spherical_function(v_p,f_p,FOD_SN_BIC_rank_mean,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_SN_BIC_rank_pop));
view([0 1 0]);
% camroll(90);
title('SN AICr rep.');

subplot(1,8,7)
hold all;
plot_spherical_function(v_p,f_p,FOD_SN_AIC_rank_pop,options);
alpha(0.25);
hold on;
plot_spherical_function(v_p,f_p,FOD_SN_AIC_rank_mean,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_SN_AIC_rank_pop));
view([0 1 0]);
% camroll(90);
title('SN AICr rep.');
%}



%{
theta1_diff_SH8 = theta1_SH8 - theta0(1);
theta1_diff_SH8(theta1_diff_SH8>(pi/2)) = theta1_diff_SH8(theta1_diff_SH8>(pi/2))-pi;
theta2_diff_SH8 = theta2_SH8 - theta0(2);
theta2_diff_SH8(theta2_diff_SH8>(pi/2)) = theta2_diff_SH8(theta2_diff_SH8>(pi/2))-pi;
phi1_diff_SH8 = phi1_SH8 - phi0(1);
phi1_diff_SH8(phi1_diff_SH8>pi) = phi1_diff_SH8(phi1_diff_SH8>pi)-2*pi;
phi2_diff_SH8 = phi2_SH8 - phi0(2);
phi2_diff_SH8(phi2_diff_SH8>pi) = phi2_diff_SH8(phi2_diff_SH8>pi)-2*pi;
theta1_diff_sCSD8 = theta1_sCSD8 - theta0(1);
theta1_diff_sCSD8(theta1_diff_sCSD8>(pi/2)) = theta1_diff_sCSD8(theta1_diff_sCSD8>(pi/2))-pi;
theta2_diff_sCSD8 = theta2_sCSD8 - theta0(2);
theta2_diff_sCSD8(theta2_diff_sCSD8>(pi/2)) = theta2_diff_sCSD8(theta2_diff_sCSD8>(pi/2))-pi;
phi1_diff_sCSD8 = phi1_sCSD8 - phi0(1);
phi1_diff_sCSD8(phi1_diff_sCSD8>pi) = phi1_diff_sCSD8(phi1_diff_sCSD8>pi)-2*pi;
phi2_diff_sCSD8 = phi2_sCSD8 - phi0(2);
phi2_diff_sCSD8(phi2_diff_sCSD8>pi) = phi2_diff_sCSD8(phi2_diff_sCSD8>pi)-2*pi;
theta1_diff_sCSD12 = theta1_sCSD12 - theta0(1);
theta1_diff_sCSD12(theta1_diff_sCSD12>(pi/2)) = theta1_diff_sCSD12(theta1_diff_sCSD12>(pi/2))-pi;
theta2_diff_sCSD12 = theta2_sCSD12 - theta0(2);
theta2_diff_sCSD12(theta2_diff_sCSD12>(pi/2)) = theta2_diff_sCSD12(theta2_diff_sCSD12>(pi/2))-pi;
phi1_diff_sCSD12 = phi1_sCSD12 - phi0(1);
phi1_diff_sCSD12(phi1_diff_sCSD12>pi) = phi1_diff_sCSD12(phi1_diff_sCSD12>pi)-2*pi;
phi2_diff_sCSD12 = phi2_sCSD12 - phi0(2);
phi2_diff_sCSD12(phi2_diff_sCSD12>pi) = phi2_diff_sCSD12(phi2_diff_sCSD12>pi)-2*pi;
theta1_diff_sCSD16 = theta1_sCSD16 - theta0(1);
theta1_diff_sCSD16(theta1_diff_sCSD16>(pi/2)) = theta1_diff_sCSD16(theta1_diff_sCSD16>(pi/2))-pi;
theta2_diff_sCSD16 = theta2_sCSD16 - theta0(2);
theta2_diff_sCSD16(theta2_diff_sCSD16>(pi/2)) = theta2_diff_sCSD16(theta2_diff_sCSD16>(pi/2))-pi;
phi1_diff_sCSD16 = phi1_sCSD16 - phi0(1);
phi1_diff_sCSD16(phi1_diff_sCSD16>pi) = phi1_diff_sCSD16(phi1_diff_sCSD16>pi)-2*pi;
phi2_diff_sCSD16 = phi2_sCSD16 - phi0(2);
phi2_diff_sCSD16(phi2_diff_sCSD16>pi) = phi2_diff_sCSD16(phi2_diff_sCSD16>pi)-2*pi;
theta1_diff_SN = theta1_SN - theta0(1);
theta1_diff_SN(theta1_diff_SN>(pi/2)) = theta1_diff_SN(theta1_diff_SN>(pi/2))-pi;
theta2_diff_SN = theta2_SN - theta0(2);
theta2_diff_SN(theta2_diff_SN>(pi/2)) = theta2_diff_SN(theta2_diff_SN>(pi/2))-pi;
phi1_diff_SN = phi1_SN - phi0(1);
phi1_diff_SN(phi1_diff_SN>pi) = phi1_diff_SN(phi1_diff_SN>pi)-2*pi;
phi2_diff_SN = phi2_SN - phi0(2);
phi2_diff_SN(phi2_diff_SN>pi) = phi2_diff_SN(phi2_diff_SN>pi)-2*pi;
figure;
subplot(1,2,1);
plot(0,0,'.','Color','red','MarkerSize', 20);
hold on;
plot(theta1_diff_SH8*180/pi, phi1_diff_SH8*180/pi,'+','Color','blue','MarkerSize', 5);
axis([-90 90 -180 180])
subplot(1,2,2);
plot(0,0,'.','Color','red','MarkerSize', 20);
hold on;
plot(theta2_diff_SH8*180/pi, phi2_diff_SH8*180/pi,'+','Color','blue','MarkerSize', 5);
axis([-90 90 -180 180])
savefig(strcat(save_path,'angle_scatter_SH.fig'));
figure;
subplot(1,2,1);
plot(0,0,'.','Color','red','MarkerSize', 20);
hold on;
plot(theta1_diff_sCSD8*180/pi, phi1_diff_sCSD8*180/pi,'+','Color','blue','MarkerSize', 5);
axis([-90 90 -180 180])
subplot(1,2,2);
plot(0,0,'.','Color','red','MarkerSize', 20);
hold on;
plot(theta2_diff_sCSD8*180/pi, phi2_diff_sCSD8*180/pi,'+','Color','blue','MarkerSize', 5);
axis([-90 90 -180 180])
savefig(strcat(save_path,'angle_scatter_sCSD8.fig'));
figure;
subplot(1,2,1);
plot(0,0,'.','Color','red','MarkerSize', 20);
hold on;
plot(theta1_diff_sCSD12*180/pi, phi1_diff_sCSD12*180/pi,'+','Color','blue','MarkerSize', 5);
axis([-90 90 -180 180])
subplot(1,2,2);
plot(0,0,'.','Color','red','MarkerSize', 20);
hold on;
plot(theta2_diff_sCSD12*180/pi, phi2_diff_sCSD12*180/pi,'+','Color','blue','MarkerSize', 5);
axis([-90 90 -180 180])
savefig(strcat(save_path,'angle_scatter_sCSD12.fig'));
figure;
subplot(1,2,1);
plot(0,0,'.','Color','red','MarkerSize', 20);
hold on;
plot(theta1_diff_sCSD16*180/pi, phi1_diff_sCSD16*180/pi,'+','Color','blue','MarkerSize', 5);
axis([-90 90 -180 180])
subplot(1,2,2);
plot(0,0,'.','Color','red','MarkerSize', 20);
hold on;
plot(theta2_diff_sCSD16*180/pi, phi2_diff_sCSD16*180/pi,'+','Color','blue','MarkerSize', 5);
axis([-90 90 -180 180])
savefig(strcat(save_path,'angle_scatter_sCSD16.fig'));
figure;
subplot(1,2,1);
plot(0,0,'.','Color','red','MarkerSize', 20);
hold on;
plot(theta1_diff_SN*180/pi, phi1_diff_SN*180/pi,'+','Color','blue','MarkerSize', 5);
axis([-90 90 -180 180])
subplot(1,2,2);
plot(0,0,'.','Color','red','MarkerSize', 20);
hold on;
plot(theta2_diff_SN*180/pi, phi2_diff_SN*180/pi,'+','Color','blue','MarkerSize', 5);
axis([-90 90 -180 180])
savefig(strcat(save_path,'angle_scatter_SN.fig'));
rate = [rate_SH8 under_SH8 over_SH8; rate_sCSD8 under_sCSD8 over_sCSD8; rate_sCSD12 under_sCSD12 over_sCSD12; rate_sCSD16 under_sCSD16 over_sCSD16;rate_SN under_SN over_SN];
figure;
bar(rate,'stacked');
savefig(strcat(save_path,'rate.fig'));
%}

%{
nfib_all = [nfib_SH_lmax8'; nfib_sCSD_lmax8'; nfib_sCSD_lmax12'; nfib_SN']';
new_nfib_all = nfib_all;
new_nfib_all(new_nfib_all>=6) = 6;

max_nfib = max(max(new_nfib_all));
min_nfib = min(min(new_nfib_all));

ticks = min_nfib:1:max_nfib;
map = zeros(max_nfib-min_nfib+1,3);
for i=1:max_nfib-min_nfib+1
    
    if(i==nfib-min_nfib+1)
        map(i,:) = [0 0 1];
    elseif(i<nfib-min_nfib+1)
        map(i,:) = [0 1-0.2*(nfib-min_nfib+1-i) 0];
    elseif(i>nfib-min_nfib+1)
        map(i,:) = [1-(0.2*(i-nfib+min_nfib-1)) 0 0];
    end
end

figure;
imagesc(new_nfib_all);        % draw image and scale colormap to values range
colormap(map);   % set colormap
set(gca, 'XTick', [1;2;3;4], 'XTickLabel', cellstr(['SH    ';'sCSD8 ';'sCSD12'; 'SN    ']),'FontSize', 20)
colorbar('Ticks',ticks);          % show color scale
savefig(strcat(save_path,'1fib','_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_N',num2str(n_sample),'_heat.fig'));
%}

% openfig('rate.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% latex table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SH_input = round([rate_SH8 under_SH8 over_SH8 mean_angle1_SH_lmax8 mean_angle2_SH_lmax8 mean_sep_SH_lmax8 mean_HD_SH_lmax8 sd_HD_SH_lmax8],4);
sCSD8_input = round([rate_sCSD8 under_sCSD8 over_sCSD8 mean_angle1_sCSD_lmax8 mean_angle2_sCSD_lmax8 mean_sep_sCSD_lmax8 mean_HD_sCSD_lmax8 sd_HD_sCSD_lmax8],4);
sCSD12_input = round([rate_sCSD12 under_sCSD12 over_sCSD12 mean_angle1_sCSD_lmax12 mean_angle2_sCSD_lmax12 mean_sep_sCSD_lmax12 mean_HD_sCSD_lmax12 sd_HD_sCSD_lmax12],4);
sCSD16_input = round([rate_sCSD16 under_sCSD16 over_sCSD16 mean_angle1_sCSD_lmax16 mean_angle2_sCSD_lmax16 mean_sep_sCSD_lmax16 mean_HD_sCSD_lmax16 sd_HD_sCSD_lmax16],4);
SN_input = round([rate_SN under_SN over_SN mean_angle1_SN mean_angle2_SN mean_sep_SN mean_HD_SN sd_HD_SN],4);

input.data = [SH_input;sCSD8_input;sCSD12_input;SN_input];

input.tableColLabels = {'Correct','Under','Over', 'Error-1', 'Error-2', 'Mean Sep.', 'Mean H-dist.', 'SD H-dist.'};
input.tableRowLabels = {'SH-ridge','SCSD8','SCSD12', 'SN-lasso'};

% Switch transposing/pivoting your table:
input.transposeTable = 0;

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used

input.dataFormat = {'%.2f',6,'%.4f',2}; % three digits precision for first two columns, one digit for the last

% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = 'NAN';

% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';

% Switch table borders on/off:
input.tableBorders = 1;

% LaTex table caption:
input.tableCaption = strcat('b=', num2str(b(1)), ' ratio=', num2str(ratio(1)), ' lmax=', num2str(lmax), ' N=',num2str(n_sample));

% LaTex table label:
input.tableLabel = strcat('b', num2str(b(1)), '_ratio', num2str(ratio(1)), '_lmax', num2str(lmax), '_N', num2str(n_sample));

% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;

% call latexTable:
latex = latexTable(input);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
Ctrs_loglambda = -6:0.2:-2;
Ctrs_df = 0:5:45;

figure
subplot(3,2,1)
hist(df_SN_BIC_rank_s,Ctrs_df);
title('df rank BIC');
subplot(3,2,3)
hist(df_SN_AIC_rank_s,Ctrs_df);
title('df rank AIC');
subplot(3,2,5)
hist(df_SN_RSSdiff_s,Ctrs_df);
title('df RSSdiff');

subplot(3,2,2)
hist(log10(lambda_seq_la(ind_SN_BIC_rank_s)),Ctrs_loglambda);
title('log lambda rank BIC');
subplot(3,2,4)
hist(log10(lambda_seq_la(ind_SN_AIC_rank_s)),Ctrs_loglambda);
title('log lambda rank AIC');
subplot(3,2,6)
hist(log10(lambda_seq_la(ind_SN_RSSdiff_s)),Ctrs_loglambda);
title('log lambda RSSdiff');
savefig('df_loglambda_hist');
%}

save(strcat(save_path,'space.mat'));
%{
figure;
% subplot(2,3,5)
plot_spherical_function(v_p,f_p,FOD_BIC_Re_pop,options);
alpha(0.25)
hold on;
plot_spherical_function(v_p,f_p,FOD_BIC_Re_mean,options);
%draw_direction(theta0,phi0,0.001);
title('BIC Re average')
figure;
% subplot(2,3,6)
plot_spherical_function(v_p,f_p,FOD_AIC_Re_pop,options);
alpha(0.25)
hold on;
plot_spherical_function(v_p,f_p,FOD_AIC_Re_mean,options);
%draw_direction(theta0,phi0,0.001);
title('AIC Re average')
%}
%}
%{
scale = max(FOD_admm_rank_C);
figure;
% axis(temp);
% subplot(2,3,4)
hold all;
plot_spherical_function(v_p,f_p,FOD_admm_rank_C*2,options);
alpha(0.25);
% figure;
% hold on;
% axis(temp);
index = 3;
figure;
hold all;
% axis(temp);
% subplot(2,2,4)
plot_spherical_function(v_p,f_p,FOD_SN(index,:),options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_SN(index,:)));
figure;
hold all;
% axis(temp);
% subplot(2,2,4)
plot_spherical_function(v_p,f_p,FOD_SH(index,:),options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_SN(index,:)));
figure;
hold all;
% axis(temp);
% subplot(2,2,4)
plot_spherical_function(v_p,f_p,FOD_SN_rank(index,:),options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_SN(index,:)));
figure;
hold all;
plot_spherical_function(v_p,f_p,FOD_admm_rank_C,options);
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_admm_rank_C));
figure;
subplot(1,3,1)
plot_spherical_function(v_p,f_p,FOD_lmax8_sCSD,options);
subplot(1,3,2)
plot_spherical_function(v_p,f_p,FOD_lmax12_sCSD,options);
subplot(1,3,3)
plot_spherical_function(v_p,f_p,FOD_lmax16_sCSD,options);
index = 1;
figure;
subplot(2,2,1)
plot_spherical_function(v_p,f_p,FOD_sCSD_lmax8(index,:),options);
hold on
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_sCSD_lmax16(index,:)));
view([0 1 0])
subplot(2,2,2)
plot_spherical_function(v_p,f_p,FOD_SN_rank(index,:),options);
hold on
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_SN_rank(index,:)));
view([0 1 0])
index = 96;
subplot(2,2,3)
plot_spherical_function(v_p,f_p,FOD_sCSD_lmax8(index,:),options);
hold on
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_sCSD_lmax16(index,:)));
view([0 1 0])
subplot(2,2,4)
plot_spherical_function(v_p,f_p,FOD_SN_rank(index,:),options);
hold on
draw_fiber(theta0,phi0,plot_scale,plot_rho*max(FOD_SN_rank(index,:)));
view([0 1 0])
%}


