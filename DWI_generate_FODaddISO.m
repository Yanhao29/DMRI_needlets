%% function to generate data accordng to S = Phi*R*f, 

function [DWI, theta0_use, phi0_use, rotationMatrix] = DWI_generate_FODaddISO(J, lmax, b, ratio, weight, theta0, phi0, sigma, ISOweight,half, Rotate, seed)
% J defines sample size, J==2.5 corresponds to J==3 n_sample=41
% lmax defines how fine we represent the dirac function using SH
% b defines b-value
% theta0, phi0 defines dirac direction
% sigma defines noise level, if sigma==0 noiseless
% rotate = 1, a random rotation guided by seed; else no rotation

% theta_r, phi_r: fiber direction angular coordinates after random rotation

%{
J = 2.5;
lmax = 16;
b = [3, 3];
ratio = [10, 10];
weight = [0.5, 0.5];
theta0 = [0, pi/2];
phi0 = [0, 0];
sigma = 0.05;
half = 1;
Rotate = 1;
ISOweight = 15;
seed = 31;
%}
k = size(b,2);
fod = zeros(k,3);
for i=1:k
    fod(i,:) = [cos(phi0(i))*sin(theta0(i)) sin(phi0(i))*sin(theta0(i)) cos(theta0(i))];
end

s = RandStream('mcg16807','Seed',seed);%21
if(Rotate == 1)
    %RandStream.setDefaultStream(s); %%for PC: set random number generator seed
    RandStream.setGlobalStream(s); %%for server: set random number generator seed
    temp=RandOrthMat(3); %%% uniformly generate a 3 by 3 orthogonal matrix
else
    temp = eye(3);
end
rotationMatrix = temp;
    
fod_use = fod*temp;
theta0_use = zeros(k,1);
phi0_use = zeros(k,1);
for i = 1:k
    phi0_use(i) = atan2(fod_use(i,2), fod_use(i,1));
    phi0_use(i) = phi0_use(i)+(phi0_use(i)<0)*2*pi;
    theta0_use(i) = acos(fod_use(i,3));
end

%{
options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 1;
%}

if(J==2.5)
    J_use=3;
else
    J_use = J;
end

%{
[vertex,~] = compute_semiregular_sphere(J_use,options); %%vertex and face of the grid 
pos = vertex{end};  %% x-y-z coordinates of the vertex 

% spherical coordinates of the vertex
phi = atan2(pos(2,:),pos(1,:))/(2*pi);   %%phi: azimuthal  angle, [0,1)
phi = phi+(phi<0);
theta = acos(pos(3,:))/(pi);             %% theta: polar angle, [0,1)
    
if(half==1)
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
    
    if(J == 2.5)
    %     pos_sampling_h = pos(:,sampling_index); %% position of the half-sphere grid points 
        phi_h=phi(:,sampling_index)*180; 
        theta_h=theta(:,sampling_index)*180;

        %%% take 40 out of these 81 directions: at each level of theta, take about
        %%% half phi 
        index_1=find(theta_h<10); %%only 1
    %     n_1=size(index_1,2);

        index_t=find(theta_h>10&theta_h<20);
    %     n_2=size(index_t,2); %%6 
        [~, I]=sort(phi_h(index_t));
        index_2=index_t(1, I([1 3 5]));

        index_t=find(theta_h>20&theta_h<40);
    %     n_3=size(index_t,2); %%12
        [~, I]=sort(phi_h(index_t));
        index_3=index_t(1, I([1 3 5 7 9 11]));

        index_t=find(theta_h>40&theta_h<50);
    %     n_4=size(index_t,2); %%12
        [~, I]=sort(phi_h(index_t));
        index_4=index_t(1, I([1 3 5 7 9 11]));

        index_t=find(theta_h>50&theta_h<70);
    %     n_5=size(index_t,2); %%20
        [~, I]=sort(phi_h(index_t));
        index_5=index_t(1, I([1 3 5 7 9 11 13 15 17 19]));

        index_t=find(theta_h>70&theta_h<85);
    %     n_6=size(index_t,2); %%22
        [~, I]=sort(phi_h(index_t));
        index_6=index_t(1, I([1 3 5 7 9 11 13 15 17 19 21]));

        index_t=find(theta_h>85);
    %     n_7=size(index_t,2); %%8
        [~, I]=sort(phi_h(index_t));
        index_7=index_t(1, I([1 3 5 7]));

        index_s=unique([index_1 index_2 index_3 index_4 index_5 index_6 index_7]);
        sampling_grid_index=sampling_index(index_s); 
    else
        sampling_grid_index=sampling_index;
    end
else
    sampling_grid_index = 1:size(pos,2);
end
%}
coe_sh = zeros(k,(lmax+1)*(lmax+2)/2);
for i=1:k
    coe_sh(i,:) = Dirac_SH_coe(lmax,theta0_use(i),phi0_use(i));
end
coe_sh = coe_sh'* weight';
coe_sh = coe_sh/sqrt(sum(coe_sh.^2));  %% make FOD a distribution

coe_1 = [1 zeros(1,length(coe_sh)-1)]';
coe_sh_test = coe_sh + coe_1*ISOweight;
%coe_sh_test = coe_sh_test/sqrt(sum(coe_sh_test.^2));

Phi = SH_vertex(J, lmax, half);
%R = Response_Rmatrix_construction(b(1),ratio(1),J_use,lmax);
load(strcat('Rmatrix_J5','_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'.mat' ));  %file name of the file that stores the corresponding R matrix

DWI_noiseless_all = Phi*Rmatrix*coe_sh_test;

%%% add Rician noise to get observed DWI 
DWI_simulated_all = add_Rician_noise(DWI_noiseless_all, sigma, seed);

%%% get observed DWI on the sampling grid 
DWI=DWI_simulated_all;
    



