
clear all;
lmax = 8;
b = [1, 1]; 
ratio = [10, 10];
n_sample = 41;
sigma = 0.05;
separation = 90;
DTI_lmax = 8;

folder_path = strcat('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_review/','2fib','_sep',num2str(separation),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'_SH',num2str(DTI_lmax),'/');    
% folder_path = strcat('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_review/','1fib','_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'/');    
% folder_path = strcat('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_review/','1fib','_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'/');    

% folder_path = strcat('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_new_results/','3fib_sep',num2str(round(sep*180/pi,0)),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'/');    
% folder_path = strcat('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_new_results/','missRes_2fib_sep',num2str(round(sep*180/pi,0)),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_ratioR',num2str(10),'_n',num2str(n_sample),'_sig',num2str(sigma),'/');    
% folder_path = strcat('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_new_results/','3fib_plain_sep',num2str(round(sep*180/pi,0)),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'/');    

N_rep=100;  %% number of replicates     
  

for rep = 1:N_rep

	display(rep); 
%     temp_name = strcat('0fib','_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(1),'_n',num2str(n_sample),'_sig',num2str(sigma),'_rep',num2str(rep),'.mat');
%     temp_name = strcat('1fib','_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'_rep',num2str(rep),'.mat');
%     temp_name = strcat('3fib_sep',num2str(round(sep*180/pi,0)),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'_rep',num2str(rep),'.mat');
    temp_name = strcat('2fib_sep',num2str(separation),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'_rep',num2str(rep),'.mat');
    filename = strcat(folder_path, temp_name);
	load(filename);


%     if(weighted==1)
%         folder_path = strcat('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_new_results/','weighted_2fib_sep',num2str(round(sep*180/pi,0)),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'/');
%     else
%         folder_path = strcat('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_new_results/','2fib_sep',num2str(round(sep*180/pi,0)),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'/');    
%     end
%     folder_path = strcat('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_new_results/','missRes_2fib_sep',num2str(round(sep*180/pi,0)),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_ratioR',num2str(10),'_n',num2str(n_sample),'_sig',num2str(sigma),'/');    
%     folder_path = strcat('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_new_results/','3fib_sep',num2str(round(sep*180/pi,0)),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'/');    
%     folder_path = strcat('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_new_results/','3fib_plain_sep',num2str(round(sep*180/pi,0)),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'/');    
    clearvars C_trans C_trans_symm Constraint design_SH_lmax8 design_SH_lmax12 design_SH_lmax16 design_SN pos_plot SH_J5_all SH_vertex_lmax8 SH_vertex_lmax12 SH_vertex_lmax16 SN_vertex_symm X X_inv Penalty_matrix_lmax8 Penalty_matrix_lmax12 Penalty_matrix_lmax16;
%     newfilename = strcat(folder_path,'new_',temp_name);
%     folder_path = strcat('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_review/','1fib','_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'/');    
%     folder_path = strcat('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/simulation_review/','2fib_sep',num2str(round(sep*180/pi,0)),'_lmax',num2str(lmax),'_b',num2str(b(1)),'_ratio',num2str(ratio(1)),'_n',num2str(n_sample),'_sig',num2str(sigma),'/');    
%     filename = strcat(folder_path, temp_name);
    save(filename);
    
end
