%% function to apply SH-ridge to the DWI data (select l2 tunning parameter by BIC)

function [f_est_all, dwi_est_all, df_all, RSS_all, BIC_all, index_selected] ...
    = SH_ridge(DWI, design_matrix, penalty_matrix, lambda_seq)
%%%%%%%%%% Output
%% f_est_all: estimated SH coefficient for all lambda
%% dwi_est_all: estimated dwi for all lambda
%% RSS_all: RSS for all lambda
%% df_all: degree of freedom of the model for all lambda
%% FOD_selected: selected estimated FOD in this lambda sequence using BIC 
%% User can utilize other selection criterion such as AIC using outputs
%% index_selected: the index of selected FOD in the lambda sequence

%%%%%%%%%% Input
%% DWI: DWI data
%% design_matrix: SH_vertex*R_matrix
%% penalty_matrix: ridge penalty matrix. (could use other penalty too)
%% lambda_seq: penalty parameter sequence (suggest provide lambda from large to small)
%% suggested lambda_seq: 10^-2 to 10^-6 with 100 log-equallyspaced lambda (or 10^-1 to 10^-5)
n_sample = size(DWI,1);

f_est_all = zeros(size(design_matrix,2),length(lambda_seq));
dwi_est_all = zeros(n_sample,length(lambda_seq));
df_all = zeros(length(lambda_seq),1);
RSS_all = zeros(length(lambda_seq),1);
BIC_all = zeros(length(lambda_seq),1);

for i =1 : size(lambda_seq,2)
    lambda_h=lambda_seq(i);
    temp=(design_matrix'*design_matrix + lambda_h*penalty_matrix)\design_matrix';
        f_est_all(:,i)=temp*DWI;
        %%% fitted dwi
        dwi_est_all(:,i)= design_matrix*f_est_all(:,i); %%fitted dwi values
        %%% BIC
        temp1=design_matrix*temp;
        df_all(i,1)=sum(diag(temp1));
        %%RSS
        RSS_all(i,1)=sum((dwi_est_all(:,i)-DWI).^2);
        %%%BIC
        BIC_all(i,1)=n_sample.*log(RSS_all(i,1)/n_sample)+df_all(i,1).*log(n_sample);
end
    
%%% choose lambda based on BIC 
index_selected = find(BIC_all==min(BIC_all));
index_selected = index_selected(1);    % in case of equal BIC, take the smallest model
    
