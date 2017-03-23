function SN_vertex_symm = SN_vertex_symm(J, jmax, BW, half)

%% get wavelet grid corresponds to J
options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 1;

[vertex,~] = compute_semiregular_sphere(J,options);
pos = vertex{end};

    %%%from perform_spherical_interpolation: (x,y,z) --> spehrical coordinates (r,\theta, \phi)
    phi = atan2(pos(2,:),pos(1,:))/(2*pi);   %%phi
    phi = phi+(phi<0);
    theta = acos(pos(3,:))/(pi);                    %% theta

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

SN_temp = cell2mat(psi);

    temp = 12*(1-4^(jmax+1))/(1-4);
    pix_all = zeros(3, temp);
    % location of cabature points for each level
    for i = 1:(jmax+1)
        pix_temp = cell2mat(pix2vec(2^(i-1),'nest',false));
        pix_all(:,(12*(1-4^(i-1))/(1-4)+1):(12*(1-4^(i))/(1-4))) = pix_temp;
    end
        
    % compute correlation of each pair of cabature points
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
    
    SN_symm_all = (SN_temp(1:size(cabature_pair,1),:)+SN_temp(cabature_pair,:))/2;
    SN_symm_all_use = SN_symm_all(cabature_use,:);
    
    SH_00 = zeros(1,size(phi,2));
    for k = 1:size(phi,2)
        SH_00(k) = spharmonic_eval(0,0,theta(k),phi(k));
    end

    SN_vertex_symm_all = [SH_00; SN_symm_all_use];
    
if(half==0)
    SN_vertex_symm = transpose(SN_vertex_symm_all);
elseif(half==1)
    SN_vertex_symm = transpose(SN_vertex_symm_all(:,sampling_index));
end
