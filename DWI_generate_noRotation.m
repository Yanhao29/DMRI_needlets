%% function to generate data

function [DWI, theta0_use, phi0_use] = DWI_generate_noRotation(J, b, ratio, weight, theta0, phi0, sigma, half, seed)
% J defines sample size, J==2.5 corresponds to J==3 n_sample=41
% theta0, phi0 defines dirac direction
% sigma defines noise level, if sigma==0 noiseless
% rotate = 1, a random rotation guided by seed; else no rotation

% theta_r, phi_r: fiber direction angular coordinates after random rotation

k = size(b,2);
fod = zeros(k,3);
for i=1:k
    fod(i,:) = [cos(phi0(i))*sin(theta0(i)) sin(phi0(i))*sin(theta0(i)) cos(theta0(i))];
end

temp = eye(3);
    
fod_use = fod*temp;
theta0_use = zeros(k,1);
phi0_use = zeros(k,1);
for i = 1:k
    phi0_use(i) = atan2(fod_use(i,2), fod_use(i,1));
    phi0_use(i) = phi0_use(i)+(phi0_use(i)<0)*2*pi;
    theta0_use(i) = acos(fod_use(i,3));
end


options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 1;

if(J==2.5)
    J_use=3;
else
    J_use = J;
end

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

% pos_sampling = pos(:,sampling_grid_index); %% The x , y , z coordinates of the sampling grid.
% phi_sampling = phi(:,sampling_grid_index); %% The sampled phi.
% theta_sampling = theta(:,sampling_grid_index); %% The sampled theta.

DWI_noiseless_all = zeros(size(theta,2),1); 
for at = 1:size(theta, 2)
    DWI_noiseless_all(at) = myresponse_crossing(b,ratio,weight,theta0_use,phi0_use,theta(at)*pi,phi(at)*2*pi); 
end

%%% add Rician noise to get observed DWI 
DWI_simulated_all = add_Rician_noise(DWI_noiseless_all, sigma, seed);

%%% get observed DWI on the sampling grid 
DWI=DWI_simulated_all(sampling_grid_index);
    



