function psi = SN_vertex(J, jmax, BW, half)

%% get wavelet grid corresponds to J
options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 1;

[vertex,~] = compute_semiregular_sphere(J,options);
pos = vertex{end};

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

phi = phi.*2*pi;
theta = theta.*pi;

psi = cell(jmax+1,1);
    
for j = 1:(jmax+1)
    psi{j} = zeros(12*(4)^(j-1),size(phi,2));
end

for k = 1:size(phi,2)
    current_ND = spneedlet_eval(theta(k),phi(k),BW,jmax);
    for j = 1:(jmax+1)
        psi{j}(:,k) = current_ND{j};
    end
end

