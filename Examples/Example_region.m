%% Simulation for a rigion with crossing fiber

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
b = [1, 1]; % back ground magnetic field strength (1 is same as b=1000)
ratio = [10, 10]; % ratio of the leading eigenvalue to the second eigenvalue in the signal simulation
weight = [0.5 0.5];
half = 1; % generate data on half shpere
lmax = 8;  % SH levels
jmax = 3; % SN levels corresponding to lmax

J_r = 5; % vertex level used for graph and representation purpose (dense)
b_response = b(1);  % b value for the response matrix Rmatrix that we use
ratio_response = ratio(1); % shape of the response function, could be misspecified 
sigma = 0.05;  %noixe level %middle  noise: note S0=1, so SNR=20 

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
save_path = strcat(path_save,'simulation_review/','region/');
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
lambda_scsd = 1;
thresh_scsd = 1e-5;
tau = 0.1;
maxit_scsd = 500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% set simulation region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n1=10;
n2=10;

%% functions in the region
x = 0 : .001 : 1; 
y = 0 : .001 : 1;
[xx, yy] = meshgrid(x,y);
f1 = sqrt(1-(x-1).^2);
f2 = real(sqrt(0.36-(x-1).^2));
f3 = sqrt(1-(x).^2);
f4 = real(sqrt(0.36-(x).^2));

grid_x = n1;
grid_y = n2;
fib_indi = zeros(grid_x,grid_y);
fib1_area = zeros(grid_x,grid_y);
fib2_area = zeros(grid_x,grid_y);
for i = 1:grid_x
    for j =1:grid_y
        lb = yy((j-1)*100+1,1);
        up = yy(j*100+1,1);
        fib1 = ~((lb>=f1(i*100+1))|(up<=(f2((i-1)*100+1))));
        fib2 = 2*(~(((lb>=f3((i-1)*100+1))|(up<=f4((i)*100+1)))));
        fib_indi(i,j) = fib1+fib2;
        f1_u = min(f1(((i-1)*100+1):100*i),up);
        f1_l = max(f2(((i-1)*100+1):100*i),lb);
        f2_u = min(f3(((i-1)*100+1):100*i),up);
        f2_l = max(f4(((i-1)*100+1):100*i),lb);
        
        fib1_area(i,j) = fib1*trapz(max(f1_u-f1_l,0));
        fib2_area(i,j) = 0.5*fib2*trapz(max(f2_u-f2_l,0));
    end
end

r = 0.05 : .1 : 0.95; 
[r_x, r_y] = meshgrid(r, r); 
r_radius1 = (1-r_x).^2+r_y.^2;
r_radius2 = (r_x).^2+r_y.^2;

slope_f1 = -1./((r_radius1-(r_x-1).^2).^(0.5)).*(r_x-1).*(fib_indi==1|fib_indi==3);
l2_norm1 = sqrt((slope_f1.^2+(ones(size(r_x))).^2));
slope_f2 = -1./((r_radius2-(r_x).^2).^(0.5)).*(r_x).*(fib_indi==2|fib_indi==3);
l2_norm2 = sqrt((slope_f2.^2+(ones(size(r_x))).^2));

%% weights in each voxel
w1_matrix = (fib_indi==1)+0.5.*(fib_indi==3);
w2_matrix = (fib_indi==2)+0.5.*(fib_indi==3);


%% angles in each voxel
theta_fib = zeros(n1,n2,2);
phi_fib = zeros(n1,n2,2);

for i=1:n1
    for j=1:n2   
        theta_fib(i,j,1) = atan(1/slope_f1(i,j));
        theta_fib(i,j,2) = -atan(1/slope_f2(i,j))+pi/2;
        phi_fib(i,j,1) = pi/2;
        phi_fib(i,j,2) = pi/2;
    end
end

dirac_sh_all = zeros(n1,n2,2562);
dirac_sh1_all = zeros(n1,n2,2562);
dirac_sh2_all = zeros(n1,n2,2562);

for i=1:n1
    for j=1:n2
        coe_sh1 = Dirac_SH_coe(lmax8,theta_fib(i,j,1),phi_fib(i,j,1)); %% SH coefficients 
        coe_sh2 = Dirac_SH_coe(lmax8,theta_fib(i,j,2),phi_fib(i,j,2));
        coe_sh=w1_matrix(i,j).*coe_sh1+w2_matrix(i,j).*coe_sh2;

        dirac_sh1 = SH_J5_all_lmax8*coe_sh1; 
        dirac_sh2 = SH_J5_all_lmax8*coe_sh2;
        dirac_sh = w1_matrix(i,j).*dirac_sh1+w2_matrix(i,j).*dirac_sh2;  %%SH representation
        dirac_sh_all(i,j,:) = dirac_sh; 
        dirac_sh1_all(i,j,:) = dirac_sh1; 
        dirac_sh2_all(i,j,:) = dirac_sh2; 
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% generate dwi signals on the equal-angle grid/gradient-direction grid 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DWI = zeros(n1,n2,n_sample);  %% no rotation

for i = 1:n1
    for j = 1:n2
        w = [w1_matrix(i,j),w2_matrix(i,j)];
        [DWI(i,j,:),~,~] = DWI_generate_noRotation(J_use,b,ratio,w,theta_fib(i,j,:),phi_fib(i,j,:),sigma,half,0); 
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Method I: SH+ridge lmax8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_est_all_SH_lmax8=zeros(n1,n2,size(design_SH_lmax8,2));
df_all_SH_lmax8=zeros(n1,n2);

RSS_all_SH_lmax8=zeros(n1,n2);
BIC_all_SH_lmax8=zeros(n1,n2);
index_sel_all_SH_lmax8=zeros(n1,n2);

dwi_est_all_SH_lmax8=zeros(n1,n2,n_sample);
FOD_SH_lmax8_all = zeros(n1,n2,size(SH_J5_all_lmax8,1));
FOD_SH_lmax8_all_st = zeros(n1,n2,size(SH_J5_all_lmax8,1));
%%%%%
dwi_est_all_sCSD_lmax8=zeros(n1,n2,n_sample);
FOD_sCSD_lmax8_all = zeros(n1,n2,size(SH_J5_all_lmax8,1));
FOD_sCSD_lmax8_all_st = zeros(n1,n2,size(SH_J5_all_lmax8,1));

%%%%%%%%%%%
f_est_all_SH_lmax12=zeros(n1,n2,size(design_SH_lmax12,2));
df_all_SH_lmax12=zeros(n1,n2);

RSS_all_SH_lmax12=zeros(n1,n2);
BIC_all_SH_lmax12=zeros(n1,n2);
index_sel_all_SH_lmax12=zeros(n1,n2);

dwi_est_all_SH_lmax12=zeros(n1,n2,n_sample);
FOD_SH_lmax12_all = zeros(n1,n2,size(SH_J5_all_lmax12,1));
FOD_SH_lmax12_all_st = zeros(n1,n2,size(SH_J5_all_lmax12,1));
%%%%%
dwi_est_all_sCSD_lmax12=zeros(n1,n2,n_sample);
FOD_sCSD_lmax12_all = zeros(n1,n2,size(SH_J5_all_lmax12,1));
FOD_sCSD_lmax12_all_st = zeros(n1,n2,size(SH_J5_all_lmax12,1));

%%%%%%%%%%%
f_est_all_SH_lmax16=zeros(n1,n2,size(design_SH_lmax16,2));
df_all_SH_lmax16=zeros(n1,n2);

RSS_all_SH_lmax16=zeros(n1,n2);
BIC_all_SH_lmax16=zeros(n1,n2);
index_sel_all_SH_lmax16=zeros(n1,n2);

dwi_est_all_SH_lmax16=zeros(n1,n2,n_sample);
FOD_SH_lmax16_all = zeros(n1,n2,size(SH_J5_all_lmax16,1));
FOD_SH_lmax16_all_st = zeros(n1,n2,size(SH_J5_all_lmax16,1));
%%%%%
dwi_est_all_sCSD_lmax16=zeros(n1,n2,n_sample);
FOD_sCSD_lmax16_all = zeros(n1,n2,size(SH_J5_all_lmax16,1));
FOD_sCSD_lmax16_all_st = zeros(n1,n2,size(SH_J5_all_lmax16,1));

%%%%%%%%%%
beta_SN_all=zeros(n1,n2,size(design_SN,2));
index_sel_all_SN=zeros(n1,n2);
stop_criterion_all_SN=zeros(n1,n2,lambda_length_la);

dwi_est_all_SN=zeros(n1,n2,n_sample);
FOD_SN_all = zeros(n1,n2,size(SH_J5_all_lmax16,1));
FOD_SN_all_st = zeros(n1,n2,size(SH_J5_all_lmax16,1));


for i = 1:n1
    for j = 1:n2
        
        
        DWI_temp = squeeze(DWI(i,j,:));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% lmax = 8
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Method I: SH+ridge regression
        [f_est_temp, dwi_temp, df_temp, ...
            RSS_temp, BIC_temp, index_sel_all_SH_lmax8(i,j)] ...
        = SH_ridge(DWI_temp,design_SH_lmax8,Penalty_matrix_lmax8,lambda_seq);

        %% estimated FOD using BIC selected lambda
        dwi_est_all_SH_lmax8(i,j,:) = dwi_temp(:,index_sel_all_SH_lmax8(i,j));
        f_est_all_SH_lmax8(i,j,:) = f_est_temp(:,index_sel_all_SH_lmax8(i,j));
        df_all_SH_lmax8(i,j) = df_temp(index_sel_all_SH_lmax8(i,j));
        FOD_SH_lmax8_all(i,j,:) = SH_J5_all_lmax8*squeeze(f_est_all_SH_lmax8(i,j,:));
        FOD_SH_lmax8_all_st(i,j,:) = fod_stand(FOD_SH_lmax8_all(i,j,:));
        RSS_all_SH_lmax8(i,j,:) = RSS_temp(index_sel_all_SH_lmax8(i,j));
        BIC_all_SH_lmax8(i,j,:) = BIC_temp(index_sel_all_SH_lmax8(i,j));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Method II: SuperCsd on SH+ridge regression
        f_ini = squeeze(f_est_all_SH_lmax8(i,j,:));
        f_ini((L_trucation+1):end) = 0;
        [f_est_sCSD_lmax8,diff_sCSD_lmax8,iter_sCSD_lmax8] = superCSD(DWI_temp, design_SH_lmax8,...
            SH_J5_all_lmax8, f_ini, lambda_scsd,  thresh_scsd, tau, maxit_scsd,false,0);

        dwi_est_all_sCSD_lmax8(i,j,:) = design_SH_lmax8*f_est_sCSD_lmax8;
        FOD_sCSD_lmax8_all(i,j,:) = SH_J5_all_lmax8*f_est_sCSD_lmax8;
        FOD_sCSD_lmax8_all_st(i,j,:) = fod_stand(FOD_sCSD_lmax8_all(i,j,:));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% lmax = 12
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Method I: SH+ridge regression
        [f_est_temp, dwi_temp, df_temp, ...
            RSS_temp, BIC_temp, index_sel_all_SH_lmax12(i,j)] ...
        = SH_ridge(DWI_temp,design_SH_lmax12,Penalty_matrix_lmax12,lambda_seq);

        %% estimated FOD using BIC selected lambda
        dwi_est_all_SH_lmax12(i,j,:) = dwi_temp(:,index_sel_all_SH_lmax12(i,j));
        f_est_all_SH_lmax12(i,j,:) = f_est_temp(:,index_sel_all_SH_lmax12(i,j));
        df_all_SH_lmax12(i,j) = df_temp(index_sel_all_SH_lmax12(i,j));
        FOD_SH_lmax12_all(i,j,:) = SH_J5_all_lmax12*squeeze(f_est_all_SH_lmax12(i,j,:));
        FOD_SH_lmax12_all_st(i,j,:) = fod_stand(FOD_SH_lmax12_all(i,j,:));
        RSS_all_SH_lmax12(i,j,:) = RSS_temp(index_sel_all_SH_lmax12(i,j));
        BIC_all_SH_lmax12(i,j,:) = BIC_temp(index_sel_all_SH_lmax12(i,j));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Method II: SuperCsd on SH+ridge regression
        f_ini = squeeze(f_est_all_SH_lmax12(i,j,:));
        f_ini((L_trucation+1):end) = 0;
        [f_est_sCSD_lmax12,diff_sCSD_lmax12,iter_sCSD_lmax12] = superCSD(DWI_temp, design_SH_lmax12,...
            SH_J5_all_lmax12, f_ini, lambda_scsd,  thresh_scsd, tau, maxit_scsd,false,0);

        dwi_est_all_sCSD_lmax12(i,j,:) = design_SH_lmax12*f_est_sCSD_lmax12;
        FOD_sCSD_lmax12_all(i,j,:) = SH_J5_all_lmax12*f_est_sCSD_lmax12;
        FOD_sCSD_lmax12_all_st(i,j,:) = fod_stand(FOD_sCSD_lmax12_all(i,j,:));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% lmax = 16
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Method I: SH+ridge regression
        [f_est_temp, dwi_temp, df_temp, ...
            RSS_temp, BIC_temp, index_sel_all_SH_lmax16(i,j)] ...
        = SH_ridge(DWI_temp,design_SH_lmax16,Penalty_matrix_lmax16,lambda_seq);

        %% estimated FOD using BIC selected lambda
        dwi_est_all_SH_lmax16(i,j,:) = dwi_temp(:,index_sel_all_SH_lmax16(i,j));
        f_est_all_SH_lmax16(i,j,:) = f_est_temp(:,index_sel_all_SH_lmax16(i,j));
        df_all_SH_lmax16(i,j) = df_temp(index_sel_all_SH_lmax16(i,j));
        FOD_SH_lmax16_all(i,j,:) = SH_J5_all_lmax16*squeeze(f_est_all_SH_lmax16(i,j,:));
        FOD_SH_lmax16_all_st(i,j,:) = fod_stand(FOD_SH_lmax16_all(i,j,:));
        RSS_all_SH_lmax16(i,j,:) = RSS_temp(index_sel_all_SH_lmax16(i,j));
        BIC_all_SH_lmax16(i,j,:) = BIC_temp(index_sel_all_SH_lmax16(i,j));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Method II: SuperCsd on SH+ridge regression
        f_ini = squeeze(f_est_all_SH_lmax16(i,j,:));
        f_ini((L_trucation+1):end) = 0;
        [f_est_sCSD_lmax16,diff_sCSD_lmax16,iter_sCSD_lmax16] = superCSD(DWI_temp, design_SH_lmax16,...
            SH_J5_all_lmax16, f_ini, lambda_scsd,  thresh_scsd, tau, maxit_scsd,false,0);

        dwi_est_all_sCSD_lmax16(i,j,:) = design_SH_lmax16*f_est_sCSD_lmax16;
        FOD_sCSD_lmax16_all(i,j,:) = SH_J5_all_lmax16*f_est_sCSD_lmax16;
        FOD_sCSD_lmax16_all_st(i,j,:) = fod_stand(FOD_sCSD_lmax16_all(i,j,:));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Method III: classo+ADMM
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [beta_est_all, z_all, w_all, dwi_temp, df_all_SN, df_rank_all_SN, RSS_all_SN, stop_criterion, index_selected_SN] ...
        = SN_lasso(DWI_temp,design_SN,Constraint,lambda_seq_la,window_percent,relativeRSS_thresh,ep_a,ep_r,maxit,print);
        
        beta_SN_all(i,j,:) = z_all(:,index_selected_SN);
        index_sel_all_SN(i,j) = index_selected_SN;
        stop_criterion_all_SN(i,j,:) = stop_criterion;
        dwi_est_all_SN(i,j,:) = dwi_temp(:,index_selected_SN);
        FOD_SN_all(i,j,:) = SN_vertex_symm*z_all(:,index_selected_SN);
        FOD_SN_all_st = fod_stand(FOD_SN_all(i,j,:));
    end
    display(i);
end



%{
figure
for k2=1:n2
    for k1=1:n1
        subplot(n1,n2,(n2-k2)*n1+k1)
        plot_spherical_function(v_p,f_p,squeeze(FOD_SN_all(k1,k2,:)),options)
        %draw_direction(theta_fib(k1,k2,:),phi_fib(k1,k2,:),0.001);
        view([1,0,0])
    end
    display(k2);
end
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'SH estimation (BIC) b=3','HorizontalAlignment','center','VerticalAlignment', 'top');

figure
for k2=1:n2
    for k1=1:n1
        subplot(n1,n2,(n2-k2)*n1+k1)
        plot_spherical_function(v_p,f_p,squeeze(FOD_sCSD_lmax12_all(k1,k2,:)),options)
        %draw_direction(theta_fib(k1,k2,:),phi_fib(k1,k2,:),0.001);
        view([1,0,0])
    end
    display(k2);
end
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'sCSD lmax8','HorizontalAlignment','center','VerticalAlignment', 'top');
%}
    
%{
fod_norm = zeros(n1,n2);
for k1=1:n1
    for k2=1:n2
        fod_norm(k1,k2) = norm(reshape(fod_classo_all_s(k1,k2,:),2562,1));
    end
end

% %% generate 26 FDD directions
[dir_v, dir_v_norm]=FDD_dir();
% %%% use prob sharing scheme
k=2;
thre=0;
grid_prob=Grid_prob(pos_plot,dir_v_norm, k, thre);

true_FDD = zeros(10,10,26);
sh_FDD = zeros(10,10,26);
sn_FDD = zeros(10,10,26);
sh_FDD_bias = zeros(10,10,26);
sn_FDD_bias = zeros(10,10,26);
sh_scsd_FDD = zeros(10,10,26);
sn_scsd_FDD = zeros(10,10,26);
sh_scsd_FDD_bias = zeros(10,10,26);
sn_scsd_FDD_bias = zeros(10,10,26);
for k1 = 1:n1
    for k2 = 1:n2
        true_FDD(k1,k2,:) = reshape(dirac_sh_all(k1,k2,:),1,2562)*grid_prob;
        thresh = 0.5*median(true_FDD(k1,k2,:));
        true_FDD(k1,k2,:) = (true_FDD(k1,k2,:)>thresh).*true_FDD(k1,k2,:);
        true_FDD(k1,k2,:) = (true_FDD(k1,k2,:))./sum(true_FDD(k1,k2,:));
        
        
        sh_FDD(k1,k2,:) = reshape(FOD_est_s_all(k1,k2,:),1,2562)*grid_prob;
        thresh = 0.5*median(sh_FDD(k1,k2,:));
        sh_FDD(k1,k2,:) = (sh_FDD(k1,k2,:)>thresh).*sh_FDD(k1,k2,:);
        sh_FDD(k1,k2,:) = sh_FDD(k1,k2,:)./sum(sh_FDD(k1,k2,:));
        
        sn_FDD(k1,k2,:) = reshape(fod_classo_all_s(k1,k2,:),1,2562)*grid_prob;
        thresh = 0.5*median(sn_FDD(k1,k2,:));
        sn_FDD(k1,k2,:) = (sn_FDD(k1,k2,:)>thresh).*sn_FDD(k1,k2,:);
        sn_FDD(k1,k2,:) = sn_FDD(k1,k2,:)./sum(sn_FDD(k1,k2,:));
        
        sh_scsd_FDD(k1,k2,:) = reshape(fod_scsd_all(k1,k2,:),1,2562)*grid_prob;
        sh_scsd_FDD(k1,k2,:) = sh_scsd_FDD(k1,k2,:)./sum(sh_scsd_FDD(k1,k2,:));
        
        sn_scsd_FDD(k1,k2,:) = reshape(fod_SN_scsd_all(k1,k2,:),1,2562)*grid_prob;
        sn_scsd_FDD(k1,k2,:) = sn_scsd_FDD(k1,k2,:)./sum(sn_scsd_FDD(k1,k2,:));
        
        sh_FDD_bias(k1,k2,:) = abs(sh_FDD(k1,k2,:)-true_FDD(k1,k2,:));
        sn_FDD_bias(k1,k2,:) = abs(sn_FDD(k1,k2,:)-true_FDD(k1,k2,:));
        sh_scsd_FDD_bias(k1,k2,:) = abs(sh_scsd_FDD(k1,k2,:)-true_FDD(k1,k2,:));
        sn_scsd_FDD_bias(k1,k2,:) = abs(sn_scsd_FDD(k1,k2,:)-true_FDD(k1,k2,:));
    end
end
mean(sh_FDD_bias,3)
mean(sn_FDD_bias,3)
mean(sh_scsd_FDD_bias,3)
mean(sn_scsd_FDD_bias,3)
mean(sh_FDD_bias,3)>mean(sn_FDD_bias,3)
mean(sh_scsd_FDD_bias,3)>mean(sn_scsd_FDD_bias,3)
sum(sum(mean(sh_scsd_FDD_bias,3)>mean(sn_scsd_FDD_bias,3)))

(mean(sh_scsd_FDD_bias,3)-mean(sn_scsd_FDD_bias,3))./mean(sh_scsd_FDD_bias,3)

mean(mean(mean(sh_FDD_bias,3)))
mean(mean(mean(sn_FDD_bias,3)))
mean(mean(mean(sh_scsd_FDD_bias,3)))
mean(mean(mean(sn_scsd_FDD_bias,3)))
%}
save_name = strcat('region_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_N',num2str(n_sample),'/');
save_path = strcat('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_review/','region','/',save_name);
mkdir(save_path);

figure
plot(x,f1)
hold on
plot(x,f2)
hold on 
plot(x,f3)
hold on 
plot(x,f4)
hold on
grid on
% quiver(r_y, r_x, w1_matrix.*ones(size(r_x)).*(fib_indi==1|fib_indi==3)./l2_norm1, w1_matrix.*slope_f1./l2_norm1,0.3)
% quiver(r_x, r_y, w2_matrix.*ones(size(r_x)).*(fib_indi==2|fib_indi==3)./l2_norm2, w2_matrix.*slope_f2./l2_norm2,0.3)
quiver(r_y, r_x, ones(size(r_x)).*(fib_indi==1|fib_indi==3)./l2_norm1, slope_f1./l2_norm1,0.3)
quiver(r_x, r_y, ones(size(r_x)).*(fib_indi==2|fib_indi==3)./l2_norm2, slope_f2./l2_norm2,0.3)
title('Fiber flow')
savefig(strcat(save_path,'directions.fig'));

figure
for k2=1:n2
    for k1=1:n1
        subplot(n1,n2,(n2-k2)*n1+k1)
        plot_spherical_function(v_p,f_p,squeeze(dirac_sh_all(k1,k2,:)),options)
        %draw_direction(theta_fib(k1,k2,:),phi_fib(k1,k2,:),0.001);
        view([1,0,0])
    end
    display(k2);
end
text(0.5, 1,'SH representation of true fiber','HorizontalAlignment','center','VerticalAlignment', 'top');
savefig(strcat(save_path,'SH_lmax',num2str(lmax8),'_b',num2str(b(1)),'_rep.fig'));

figure
for k2=1:n2
    for k1=1:n1
        subplot(n1,n2,(n2-k2)*n1+k1)
        plot_spherical_function(v_p,f_p,squeeze(FOD_SH_lmax8_all(k1,k2,:)),options)
        %draw_direction(theta_fib(k1,k2,:),phi_fib(k1,k2,:),0.001);
        view([1,0,0])
    end
    display(k2);
end
text(0.5, 1,'SH estimation (BIC) b=3','HorizontalAlignment','center','VerticalAlignment', 'top');
savefig(strcat(save_path,'SH_lmax',num2str(lmax8),'_b',num2str(b(1)),'_est.fig'));

figure
for k2=1:n2
    for k1=1:n1
        subplot(n1,n2,(n2-k2)*n1+k1)
        plot_spherical_function(v_p,f_p,squeeze(FOD_sCSD_lmax8_all(k1,k2,:)),options)
        %draw_direction(theta_fib(k1,k2,:),phi_fib(k1,k2,:),0.001);
        view([1,0,0])
    end
    display(k2);
end
text(0.5, 1,'SH scsd b=1 BIC','HorizontalAlignment','center','VerticalAlignment', 'top');
savefig(strcat(save_path,'sCSD_lmax',num2str(lmax8),'_b',num2str(b(1)),'_est.fig'));

figure
for k2=1:n2
    for k1=1:n1
        subplot(n1,n2,(n2-k2)*n1+k1)
        plot_spherical_function(v_p,f_p,squeeze(FOD_sCSD_lmax12_all(k1,k2,:)),options)
        %draw_direction(theta_fib(k1,k2,:),phi_fib(k1,k2,:),0.001);
        view([1,0,0])
    end
    display(k2);
end
text(0.5, 1,'SH scsd b=1 BIC','HorizontalAlignment','center','VerticalAlignment', 'top');
savefig(strcat(save_path,'sCSD_lmax',num2str(lmax12),'_b',num2str(b(1)),'_est.fig'));

figure
for k2=1:n2
    for k1=1:n1
        subplot(n1,n2,(n2-k2)*n1+k1)
        plot_spherical_function(v_p,f_p,squeeze(FOD_sCSD_lmax16_all(k1,k2,:)),options)
        %draw_direction(theta_fib(k1,k2,:),phi_fib(k1,k2,:),0.001);
        view([1,0,0])
    end
    display(k2);
end
text(0.5, 1,'SH scsd b=1 BIC','HorizontalAlignment','center','VerticalAlignment', 'top');
savefig(strcat(save_path,'sCSD_lmax',num2str(lmax16),'_b',num2str(b(1)),'_est.fig'));


figure
for k2=1:n2
    for k1=1:n1
        subplot(n1,n2,(n2-k2)*n1+k1)
        plot_spherical_function(v_p,f_p,squeeze(FOD_SN_all(k1,k2,:)),options)
        %draw_direction(theta_fib(k1,k2,:),phi_fib(k1,k2,:),0.001);
        view([1,0,0])
    end
    display(k2);
end
text(0.5, 1,'SN lasso','HorizontalAlignment','center','VerticalAlignment', 'top');
savefig(strcat(save_path,'SN_est.fig'));


nfib_SH_lmax8 = zeros(n1,n2);
nfib_SCSD_lmax8 = zeros(n1,n2);
nfib_SCSD_lmax12 = zeros(n1,n2);
nfib_SCSD_lmax16 = zeros(n1,n2);
nfib_SN_RSSdiff = zeros(n1,n2);

kmin = 40;
peak_thresh = 0.45;
Dis = squareform(pdist(pos_p','cosine'));

peak_vec_SH = [];
peak_vec_sCSD8 = [];
peak_vec_sCSD12 = [];
peak_vec_sCSD16 = [];
peak_vec_SN = [];

vec_map_SH = [];
vec_map_sCSD8 = [];
vec_map_sCSD12 = [];
vec_map_sCSD16 = [];
vec_map_SN = [];


for i=1:n1
    for j=1:n2
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        [~, ~, ~, ~, ~, peak_pos_SH_final_lmax8] = FOD_peak(FOD_SH_lmax8_all(i,j,:), Dis, kmin, peak_thresh, pos_p, theta_p, phi_p);
        [~, ~, ~, ~, ~, peak_pos_SCSD_final_lmax8] = FOD_peak(FOD_sCSD_lmax8_all(i,j,:), Dis, kmin, peak_thresh, pos_p, theta_p, phi_p);
        [~, ~, ~, ~, ~, peak_pos_SCSD_final_lmax12] = FOD_peak(FOD_sCSD_lmax12_all(i,j,:), Dis, kmin, peak_thresh, pos_p, theta_p, phi_p);
        [~, ~, ~, ~, ~, peak_pos_SCSD_final_lmax16] = FOD_peak(FOD_sCSD_lmax16_all(i,j,:), Dis, kmin, peak_thresh, pos_p, theta_p, phi_p);
        [~, ~, ~, ~, ~, peak_pos_SN_final_RSSdiff] = FOD_peak(FOD_SN_all(i,j,:), Dis, kmin, peak_thresh, pos_p, theta_p, phi_p);

        nfib_SH_lmax8(i,j) = size(peak_pos_SH_final_lmax8,2);
        nfib_SCSD_lmax8(i,j) = size(peak_pos_SCSD_final_lmax8,2);
        nfib_SCSD_lmax12(i,j) = size(peak_pos_SCSD_final_lmax12,2);
        nfib_SCSD_lmax16(i,j) = size(peak_pos_SCSD_final_lmax16,2);
        nfib_SN_RSSdiff(i,j) = size(peak_pos_SN_final_RSSdiff,2);
        
        peak_vec_SH = [peak_vec_SH, peak_pos_SH_final_lmax8];
        peak_vec_sCSD8 = [peak_vec_sCSD8, peak_pos_SCSD_final_lmax8];
        peak_vec_sCSD12 = [peak_vec_sCSD12, peak_pos_SCSD_final_lmax12];
        peak_vec_sCSD16 = [peak_vec_sCSD16, peak_pos_SCSD_final_lmax16];
        peak_vec_SN = [peak_vec_SN, peak_pos_SN_final_RSSdiff];
        
        vec_map_SH = [vec_map_SH, repmat([i;j],1,nfib_SH_lmax8(i,j))];
        vec_map_sCSD8 = [vec_map_sCSD8, repmat([i;j],1,nfib_SCSD_lmax8(i,j))];
        vec_map_sCSD12 = [vec_map_sCSD12, repmat([i;j],1,nfib_SCSD_lmax12(i,j))];
        vec_map_sCSD16 = [vec_map_sCSD16, repmat([i;j],1,nfib_SCSD_lmax16(i,j))];
        vec_map_SN = [vec_map_SN, repmat([i;j],1,nfib_SN_RSSdiff(i,j))];
   
    end
    display(i);
end

nfib = (fib_indi>2);
nfib = nfib+1;

n_1fib = sum(sum(nfib==1))
n_2fib = sum(sum(nfib==2))

SH_1fib_over = sum(sum(nfib_SH_lmax8(nfib==1)>1))
SH_1fib_correct = sum(sum(nfib_SH_lmax8(nfib==1)==1))
SH_1fib_under = sum(sum(nfib_SH_lmax8(nfib==1)<1))

SH_2fib_over = sum(sum(nfib_SH_lmax8(nfib==2)>2))
SH_2fib_correct = sum(sum(nfib_SH_lmax8(nfib==2)==2))
SH_2fib_under = sum(sum(nfib_SH_lmax8(nfib==2)<2))

SCSD8_1fib_over = sum(sum(nfib_SCSD_lmax8(nfib==1)>1))
SCSD8_1fib_correct = sum(sum(nfib_SCSD_lmax8(nfib==1)==1))
SCSD8_1fib_under = sum(sum(nfib_SCSD_lmax8(nfib==1)<1))

SCSD8_2fib_over = sum(sum(nfib_SCSD_lmax8(nfib==2)>2))
SCSD8_2fib_correct = sum(sum(nfib_SCSD_lmax8(nfib==2)==2))
SCSD8_2fib_under = sum(sum(nfib_SCSD_lmax8(nfib==2)<2))

SCSD12_1fib_over = sum(sum(nfib_SCSD_lmax12(nfib==1)>1))
SCSD12_1fib_correct = sum(sum(nfib_SCSD_lmax12(nfib==1)==1))
SCSD12_1fib_under = sum(sum(nfib_SCSD_lmax12(nfib==1)<1))

SCSD12_2fib_over = sum(sum(nfib_SCSD_lmax12(nfib==2)>2))
SCSD12_2fib_correct = sum(sum(nfib_SCSD_lmax12(nfib==2)==2))
SCSD12_2fib_under = sum(sum(nfib_SCSD_lmax12(nfib==2)<2))

SCSD16_1fib_over = sum(sum(nfib_SCSD_lmax16(nfib==1)>1))
SCSD16_1fib_correct = sum(sum(nfib_SCSD_lmax16(nfib==1)==1))
SCSD16_1fib_under = sum(sum(nfib_SCSD_lmax16(nfib==1)<1))

SCSD16_2fib_over = sum(sum(nfib_SCSD_lmax16(nfib==2)>2))
SCSD16_2fib_correct = sum(sum(nfib_SCSD_lmax16(nfib==2)==2))
SCSD16_2fib_under = sum(sum(nfib_SCSD_lmax16(nfib==2)<2))

SN_RSSdiff_1fib_over = sum(sum(nfib_SN_RSSdiff(nfib==1)>1))
SN_RSSdiff_1fib_correct = sum(sum(nfib_SN_RSSdiff(nfib==1)==1))
SN_RSSdiff_1fib_under = sum(sum(nfib_SN_RSSdiff(nfib==1)<1))

SN_RSSdiff_2fib_over = sum(sum(nfib_SN_RSSdiff(nfib==2)>2))
SN_RSSdiff_2fib_correct = sum(sum(nfib_SN_RSSdiff(nfib==2)==2))
SN_RSSdiff_2fib_under = sum(sum(nfib_SN_RSSdiff(nfib==2)<2))


SH_1fib_over_rate = sum(sum(nfib_SH_lmax8(nfib==1)>1))/n_1fib
SH_1fib_correct_rate = sum(sum(nfib_SH_lmax8(nfib==1)==1))/n_1fib
SH_1fib_under_rate = sum(sum(nfib_SH_lmax8(nfib==1)<1))/n_1fib

SH_2fib_over_rate = sum(sum(nfib_SH_lmax8(nfib==2)>2))/n_2fib
SH_2fib_correct_rate = sum(sum(nfib_SH_lmax8(nfib==2)==2))/n_2fib
SH_2fib_under_rate = sum(sum(nfib_SH_lmax8(nfib==2)<2))/n_2fib

SCSD8_1fib_over_rate = sum(sum(nfib_SCSD_lmax8(nfib==1)>1))/n_1fib
SCSD8_1fib_correct_rate = sum(sum(nfib_SCSD_lmax8(nfib==1)==1))/n_1fib
SCSD8_1fib_under_rate = sum(sum(nfib_SCSD_lmax8(nfib==1)<1))/n_1fib

SCSD8_2fib_over_rate = sum(sum(nfib_SCSD_lmax8(nfib==2)>2))/n_2fib
SCSD8_2fib_correct_rate = sum(sum(nfib_SCSD_lmax8(nfib==2)==2))/n_2fib
SCSD8_2fib_under_rate = sum(sum(nfib_SCSD_lmax8(nfib==2)<2))/n_2fib

SCSD12_1fib_over_rate = sum(sum(nfib_SCSD_lmax12(nfib==1)>1))/n_1fib
SCSD12_1fib_correct_rate = sum(sum(nfib_SCSD_lmax12(nfib==1)==1))/n_1fib
SCSD12_1fib_under_rate = sum(sum(nfib_SCSD_lmax12(nfib==1)<1))/n_1fib

SCSD12_2fib_over_rate = sum(sum(nfib_SCSD_lmax12(nfib==2)>2))/n_2fib
SCSD12_2fib_correct_rate = sum(sum(nfib_SCSD_lmax12(nfib==2)==2))/n_2fib
SCSD12_2fib_under_rate = sum(sum(nfib_SCSD_lmax12(nfib==2)<2))/n_2fib

SCSD16_1fib_over_rate = sum(sum(nfib_SCSD_lmax16(nfib==1)>1))/n_1fib
SCSD16_1fib_correct_rate = sum(sum(nfib_SCSD_lmax16(nfib==1)==1))/n_1fib
SCSD16_1fib_under_rate = sum(sum(nfib_SCSD_lmax16(nfib==1)<1))/n_1fib

SCSD16_2fib_over_rate = sum(sum(nfib_SCSD_lmax16(nfib==2)>2))/n_2fib
SCSD16_2fib_correct_rate = sum(sum(nfib_SCSD_lmax16(nfib==2)==2))/n_2fib
SCSD16_2fib_under_rate = sum(sum(nfib_SCSD_lmax16(nfib==2)<2))/n_2fib

SN_RSSdiff_1fib_over_rate = sum(sum(nfib_SN_RSSdiff(nfib==1)>1))/n_1fib
SN_RSSdiff_1fib_correct_rate = sum(sum(nfib_SN_RSSdiff(nfib==1)==1))/n_1fib
SN_RSSdiff_1fib_under_rate = sum(sum(nfib_SN_RSSdiff(nfib==1)<1))/n_1fib

SN_RSSdiff_2fib_over_rate = sum(sum(nfib_SN_RSSdiff(nfib==2)>2))/n_2fib
SN_RSSdiff_2fib_correct_rate = sum(sum(nfib_SN_RSSdiff(nfib==2)==2))/n_2fib
SN_RSSdiff_2fib_under_rate = sum(sum(nfib_SN_RSSdiff(nfib==2)<2))/n_2fib


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% latex table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SH_input = [SH_1fib_correct SH_1fib_under SH_1fib_over SH_2fib_correct SH_2fib_under SH_2fib_over];
SCSD8_input = [SCSD8_1fib_correct SCSD8_1fib_under SCSD8_1fib_over SCSD8_2fib_correct SCSD8_2fib_under SCSD8_2fib_over];
SCSD12_input = [SCSD12_1fib_correct SCSD12_1fib_under SCSD12_1fib_over SCSD12_2fib_correct SCSD12_2fib_under SCSD12_2fib_over];
SCSD16_input = [SCSD16_1fib_correct SCSD16_1fib_under SCSD16_1fib_over SCSD16_2fib_correct SCSD16_2fib_under SCSD16_2fib_over];
SN_RSSdiff_input = [SN_RSSdiff_1fib_correct SN_RSSdiff_1fib_under SN_RSSdiff_1fib_over SN_RSSdiff_2fib_correct SN_RSSdiff_2fib_under SN_RSSdiff_2fib_over];

input.data = [SH_input;SCSD8_input;SCSD12_input;SCSD16_input;SN_RSSdiff_input];

input.tableColLabels = {'1fib correct','1fib under','1fib over', '2fib correct','2fib under','2fib over'};
input.tableRowLabels = {'SH','SCSD8', 'SCSD12','SCSD16', 'SN Rssdiff'};

% Switch transposing/pivoting your table:
input.transposeTable = 0;

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used

input.dataFormat = {'%.0f',6}; % three digits precision for first two columns, one digit for the last

% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = 'NAN';

% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';

% Switch table borders on/off:
input.tableBorders = 1;

% LaTex table caption:
input.tableCaption = strcat(' b=', num2str(b(1)), ' ratio=', num2str(ratio(1)), ' lmax=', num2str(lmax), ' N=',num2str(n_sample));

% LaTex table label:
input.tableLabel = strcat('_b', num2str(b(1)), '_ratio', num2str(ratio(1)), '_lmax', num2str(lmax), '_N', num2str(n_sample));

% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;

% call latexTable:
latex1 = latexTable(input);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% latex table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SH_input = [SH_1fib_correct_rate SH_1fib_under_rate SH_1fib_over_rate SH_2fib_correct_rate SH_2fib_under_rate SH_2fib_over_rate];
SCSD8_input = [SCSD8_1fib_correct_rate SCSD8_1fib_under_rate SCSD8_1fib_over_rate SCSD8_2fib_correct_rate SCSD8_2fib_under_rate SCSD8_2fib_over_rate];
SCSD12_input = [SCSD12_1fib_correct_rate SCSD12_1fib_under_rate SCSD12_1fib_over_rate SCSD12_2fib_correct_rate SCSD12_2fib_under_rate SCSD12_2fib_over_rate];
SCSD16_input = [SCSD16_1fib_correct_rate SCSD16_1fib_under_rate SCSD16_1fib_over_rate SCSD16_2fib_correct_rate SCSD16_2fib_under_rate SCSD16_2fib_over_rate];
SN_RSSdiff_input = [SN_RSSdiff_1fib_correct_rate SN_RSSdiff_1fib_under_rate SN_RSSdiff_1fib_over_rate SN_RSSdiff_2fib_correct_rate SN_RSSdiff_2fib_under_rate SN_RSSdiff_2fib_over_rate];

input.data = [SH_input;SCSD8_input;SCSD12_input;SCSD16_input;SN_RSSdiff_input];

input.tableColLabels = {'1fib correct','1fib under','1fib over', '2fib correct','2fib under','2fib over'};
input.tableRowLabels = {'SH','SCSD8', 'SCSD12','SCSD16', 'SN Rssdiff'};

% Switch transposing/pivoting your table:
input.transposeTable = 0;

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used

input.dataFormat = {'%.2f',6}; % three digits precision for first two columns, one digit for the last

% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = 'NAN';

% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';

% Switch table borders on/off:
input.tableBorders = 1;

% LaTex table caption:
input.tableCaption = strcat(' b=', num2str(b(1)), ' ratio=', num2str(ratio(1)), ' lmax=', num2str(lmax), ' N=',num2str(n_sample));

% LaTex table label:
input.tableLabel = strcat('_b', num2str(b(1)), '_ratio', num2str(ratio(1)), '_lmax', num2str(lmax), '_N', num2str(n_sample));

% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;

% call latexTable:
latex2 = latexTable(input);

save(strcat(save_path,'space_peak',num2str(peak_thresh),'.mat'));
