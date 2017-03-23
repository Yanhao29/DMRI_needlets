%% function to apply SN-lasso to the DWI data

function [beta_est_all, z_all, w_all, dwi_est_all, df_all, df_rank_all, RSS_all, stop_criterion, index_selected] ...
    = SN_lasso(DWI, design_matrix, Constraint, lambda_seq, smoothing_window_percent, relativeRSS_thresh, ...
    ep_a, ep_r, maxit, print)

%%%%% Output
%% beta_est_all: SN coefficients across lambda sequence
%% z_all: SN coefficients with sparsity across lambda sequence
%% dwi_est_all: estimated DWI across lambda sequence
%% df_all: degree of freedom across lambda sequence
%% df_rank_all: degree of freedom using rank definition across sequence of lambda
%% RSS_all: RSS across sequence of lambda
%% stop_criterion: stores moving average relative RSS.
%% index_selected: the index of first lambda that selected by relativeRSS_thresh

%%%%% Input
%% DWI: DWI data
%% design_matrix: design_SH*C_trans where design_SH = SH_vertex*R_matrix
%% Constraint: SN evaluated on a dense grid on the sphere
%% lambda_seq: L1 tunning parameter
%% smoothing_widow_percent: Smoothing window size percent of the length of lambda_seq
%% relativeRSS_thresh: the stopping thresh hold. Stops once the relative RSS change is less than this thresh

    X1 = design_matrix'*design_matrix;
    X2 = eye(size(design_matrix,2));
    X3 = (-Constraint)'*(-Constraint);
    Y = design_matrix'*DWI;
    
    beta_est_all=zeros(size(design_matrix,2), size(lambda_seq,2));
    dwi_est_all=zeros(size(design_matrix,1), size(lambda_seq,2));
    df_all=zeros(size(lambda_seq,2),1);
    df_rank_all=zeros(size(lambda_seq,2),1);
    RSS_all=zeros(size(lambda_seq,2),1);
    stop_criterion=zeros(size(lambda_seq,2),1);

    z_all = zeros(size(design_matrix,2),size(lambda_seq,2));
    w_all = zeros(size(Constraint,1),size(lambda_seq,2));
    
    window = floor(smoothing_window_percent*length(lambda_seq));
    stop_spacing = log10(lambda_seq(2))-log10(lambda_seq(1));
    
    for i=1:length(lambda_seq)
        display(i);    
        %%% set current initial 
        if i==1 
            z=z_all(:,1);
            w=w_all(:,1);
        else
            z=z_all(:,i-1);
            w=w_all(:,i-1);
        end
        lambda_c=lambda_seq(1,i); %%current lambda
        X = X1+lambda_c*X2+lambda_c*X3;
        X_inv = inv(X);
        [beta_est_all(:,i), z_all(:,i), w_all(:,i)] = ADMM_classo(Y, X_inv, -Constraint, (-Constraint)', z, w,lambda_c, ep_r,ep_a, lambda_c, maxit, print);
        
        idx_temp = find(abs(z_all(:,i))>0); %% set cutoff manually
        df_rank_all(i,1)=rank(design_matrix(:,idx_temp));
        df_all(i,1)=size(idx_temp,1);
        dwi_est_all(:,i) = design_matrix*z_all(:,i);  %%fitted dwi
        RSS_all(i,1)=sum((dwi_est_all(:,i)-DWI).^2);
%         dwi_temp = dwi_est_all(:,i);
        
        if(i>window)
            rela_diff_temp = diff(log10(RSS_all((i-window):i,1)),1)./stop_spacing;
            stop_criterion(i) = mean(abs(rela_diff_temp));
            if(stop_criterion(i)<relativeRSS_thresh)
                index_selected = i;
                break;
            else
                index_selected = i;
            end
        end
    end
    
    