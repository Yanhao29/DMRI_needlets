%% Generate response

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
plot_scale = options.scaling;



b_response = 1;
ratio_response = 10;

J = 5;
jmax = 3;
lmax8 = 8;
lmax12 = 12;
lmax16 = 16;


load(strcat('SH_vertex_J', num2str(J), '_lmax', num2str(lmax8), '.mat'));  % spherical harmonic basis evaluated on J=5(2562 grid points)
SH_J5_all_lmax8 = SH_vertex;

load(strcat('SH_vertex_J', num2str(J), '_lmax', num2str(lmax12), '.mat'));  % spherical harmonic basis evaluated on J=5(2562 grid points)
SH_J5_all_lmax12 = SH_vertex;

load(strcat('SH_vertex_J', num2str(J), '_lmax', num2str(lmax16), '.mat'));  % spherical harmonic basis evaluated on J=5(2562 grid points)
SH_J5_all_lmax16 = SH_vertex;

load(strcat('C_symm_lmax', num2str(lmax8), '.mat'));  %% SH coefficients of the symmetrized needlets basis
load(strcat('SN_vertex_symm_J5_jmax', num2str(jmax), '.mat'));% symmetrized spherical needlets basis J=5(2562 grid points), jmax=4
load(strcat('Rmatrix_J5','_lmax',num2str(lmax8),'_b',num2str(b_response),'_ratio',num2str(ratio_response),'.mat' ));  %file name of the file that stores the corresponding R matrix
Rmatrix_lmax8 = Rmatrix;
load(strcat('Rmatrix_J5','_lmax',num2str(lmax12),'_b',num2str(b_response),'_ratio',num2str(ratio_response),'.mat' ));  %file name of the file that stores the corresponding R matrix
Rmatrix_lmax12 = Rmatrix;
load(strcat('Rmatrix_J5','_lmax',num2str(lmax16),'_b',num2str(b_response),'_ratio',num2str(ratio_response),'.mat' ));  %file name of the file that stores the corresponding R matrix
Rmatrix_lmax16 = Rmatrix;
clearvars SH_vertex Rmatrix;


theta0_response = 0;  %%align with z-axis: this is always the case 
phi0_response = 0;
w_response = 1;
%% get wavelet grid corresponds to J
options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 1;

[vertex_sampling,face_sampling] = compute_semiregular_sphere(J,options);
pos = vertex_sampling{end};
%options.spherical = 1;
%plot_mesh(pos,face_sampling{end}); axis tight;  %plot the generated mesh
%lighting flat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%from perform_spherical_interpolation: (x,y,z) --> spehrical coordinates (r,\theta, \phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = atan2(pos(2,:),pos(1,:))/(2*pi);   %%phi
phi = phi+(phi<0);
theta = acos(pos(3,:))/(pi);                    %% theta

%%% evaluation of the specified R on the wavelet grid J
R = zeros(size(theta,2) , 1);

for at = 1:size(theta,2) %% the vertex

R(at) = myresponse_crossing(b_response,ratio_response, w_response, theta0_response , phi0_response, theta(at)*pi , phi(at)*2*pi  );

end

%%%%%
Response_SH8_coe = zeros(lmax8/2+1 , 1);
for l = 0:2:lmax8
    Response_SH8_coe(l/2+1) = Rmatrix_lmax8((l+1)*(l+2)/2-l,(l+1)*(l+2)/2-l);  %%get m=0 coefficients 
end
Response_SH8_coe_real = zeros((lmax8+1)*(lmax8+2)/2,1);
current = 0;
for i=0:2:lmax8
    Response_SH8_coe_real(current+i+1,1) = Response_SH8_coe(i/2+1)/sqrt(4*pi / (2*i+1));
    current = current+2*i+1;
end
Response_SH8_rep = SH_J5_all_lmax8*Response_SH8_coe_real;

%%%%%
Response_SH12_coe = zeros(lmax12/2+1 , 1);
for l = 0:2:lmax12
    Response_SH12_coe(l/2+1) = Rmatrix_lmax12((l+1)*(l+2)/2-l,(l+1)*(l+2)/2-l);  %%get m=0 coefficients 
end
Response_SH12_coe_real = zeros((lmax12+1)*(lmax12+2)/2,1);
current = 0;
for i=0:2:lmax12
    Response_SH12_coe_real(current+i+1,1) = Response_SH12_coe(i/2+1)/sqrt(4*pi / (2*i+1));
    current = current+2*i+1;
end
Response_SH12_rep = SH_J5_all_lmax12*Response_SH12_coe_real;

%%%%%
Response_SH16_coe = zeros(lmax16/2+1 , 1);
for l = 0:2:lmax16
    Response_SH16_coe(l/2+1) = Rmatrix_lmax16((l+1)*(l+2)/2-l,(l+1)*(l+2)/2-l);  %%get m=0 coefficients 
end
Response_SH16_coe_real = zeros((lmax16+1)*(lmax16+2)/2,1);
current = 0;
for i=0:2:lmax16
    Response_SH16_coe_real(current+i+1,1) = Response_SH16_coe(i/2+1)/sqrt(4*pi / (2*i+1));
    current = current+2*i+1;
end
Response_SH16_rep = SH_J5_all_lmax16*Response_SH16_coe_real;


figure;
plot_spherical_function(vertex_sampling,face_sampling,R,options) %% 
hold on;
draw_fiber(theta0_response,phi0_response,plot_scale,plot_rho*max(R));
view([1,0,0])

figure;
plot_spherical_function(vertex_sampling,face_sampling,Response_SH8_rep,options) %% 
hold on;
draw_fiber(theta0_response,phi0_response,plot_scale,plot_rho*max(Response_SH8_rep));
view([1,0,0])

sum(abs(R'-Response_SH8_rep))
sum(abs(R'-Response_SH12_rep))
sum(abs(R'-Response_SH16_rep))


