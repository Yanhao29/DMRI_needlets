%% lmax give s the SH level

function C_trans_symm = C_trans_symm_construction(lmax)

%     path.cur='/Users/hao/Dropbox/stats_project/DTI/';
%     path.save = (strcat(path.cur,'matrix'));
%     addpath (path.save);
%     addpath (strcat(path.cur,'MEALPix'));
% 
%     filename = strcat('C_lmax', num2str(lmax),'.mat');
%     load(filename);
    C_trans = C_trans_construction(lmax);
    B=2;
    jmax = fix( log(lmax)/log(B));
    
    % location of cabature points for each level
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

    C_trans_symm = (C_trans+C_trans(:,[1 (cabature_pair+1)'])); %% first column is the first SH  basis (constant)
    C_trans_symm(:,1)=C_trans(:,1);
    % For C, the first row corresponds to first SH basis (constant), all
    % cabature points used are shifted by one (index shifted)
    cabature_use_C = 1;
    for i = 1:size(cabature_pair,1)
        if cabature_pair(i)>i
            cabature_use_C = [cabature_use_C i+1];
        end
    end

    C_trans_symm = C_trans_symm(:,cabature_use_C);
%     filename = fullfile(path.save, strcat('C_symm_lmax',num2str(lmax),'.mat' ));
%     save(filename,'C_trans_symm');
end