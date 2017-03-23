%% Evaluate sphere harmoics on mesh corresponding to J an lmax (on sphere or half sphere)
function Bmatrix = SH_vertex(J, lmax, half)
%% J = 1, 2, 2.5, 3, 4, 5... when J=2.5, J_use = 3, and we sample half of the 81 grid points
%% get wavelet grid corresponds to J
options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 1;

if(J==2.5)
    J_use = 3;
else
    J_use = J;
end

[vertex_sampling,~] = compute_semiregular_sphere(J_use,options);
pos = vertex_sampling{end};

if(half == 0)
    %%%from perform_spherical_interpolation: (x,y,z) --> spehrical coordinates (r,\theta, \phi)
    phi = atan2(pos(2,:),pos(1,:))/(2*pi);   %%phi
    phi = phi+(phi<0);
    theta = acos(pos(3,:))/(pi);                    %% theta
elseif(half == 1)
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
    sampling_grid_index=sampling_index;

    pos_sampling = pos(:,sampling_grid_index); %% The x , y , z coordinates of the sampling grid.

    phi = atan2(pos_sampling(2,:),pos_sampling(1,:))/(2*pi);   %%phi
    phi = phi+(phi<0);
    theta = acos(pos_sampling(3,:))/(pi);                    %% theta
end

%%%
L = (lmax+1)*(lmax+2)/2;  %% number of real SH basis 
Bmatrix = zeros(size(theta,2) , L);

%%% use the function "spharmonic_eval.m" written by Jun
for at = 1:size(theta,2) %% the vertex
    for l = 0:2:lmax     %% even order SH
        for m = (-l):l       %% SH phase      

        Rindex = (l+1)*(l+2)/2 - (l-m); %% column/SH basis  index

        if (m<0 && m>= (-l))  %% m<0: see Descoteaux 2007 equation 3

            Bmatrix(at , Rindex) = ( (-1)^(m)*spharmonic_eval(l,m,theta(at)*pi , phi(at)*2*pi ) + spharmonic_eval(l,-m,theta(at)*pi, phi(at)*2*pi ))/sqrt(2);

        elseif(m == 0)  %% m=0

            Bmatrix(at , Rindex) = spharmonic_eval(l,m,theta(at)*pi , phi(at)*2*pi );

        else  %%m>0

            Bmatrix(at , Rindex) = 1i*( (-1)^(m+1)*spharmonic_eval(l,m,theta(at)*pi , phi(at)*2*pi ) + spharmonic_eval(l,-m,theta(at)*pi , phi(at)*2*pi ))/sqrt(2);

        end

        end
    end
end

if(J==2.5)
    phi_h=phi*180; 
    theta_h=theta*180;

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
    Bmatrix = real(Bmatrix);
    Bmatrix = Bmatrix(index_s,:);
else
    Bmatrix = real(Bmatrix);
end