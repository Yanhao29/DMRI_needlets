function [beta,thresh_now,iter] = superCSD(obs, A, L, f_ini, lambda,  thresh, tau, maxit,ridge,lam_r)

% obs:observations (DWI signal)
% A: design matrix of basis functions (SH or needlets)
% L: basis functions evaluated on a fine grid
% lambda: penalty parameter (usually 1)
% thresh: when to stop (very small change in FOD)
% tau: cutoff value of the mean(FOD) (smaller than tau*mean(FOD_ini) will set to be 0)
% maxit: maximum number of iterations
% ridge: logical true or false (in needlets case need ridge regression)
% lam_r: penalty parameter for ridge case (1e-4 seems working good)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Super csd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fmatrix_initial_csd = f_ini;
fmatrix_old_csd = fmatrix_initial_csd;
F_old_csd = L*fmatrix_old_csd;

%%% parameters control the non-negative constraints
tao = tau*mean(F_old_csd);%%% the non-negative threshold. 0.9 yields the best result
lambda_csd = lambda;
loop_csd = maxit;
thresh_now = 99;
thresh_csd = thresh;
count = 0;

while thresh_now>thresh_csd&count<loop_csd
    Lpenalty_csd = L(find(F_old_csd<tao),:);
    big_design_csd = [A;lambda_csd.*Lpenalty_csd];
    big_DWI = zeros(size(big_design_csd , 1) , 1);
    
    big_DWI(1:length(obs)) = obs;
    
    if(ridge==false) %(SH case)
        if rank(big_design_csd)<size(big_design_csd , 2)
            %fprintf('there is a problem with rank when loop is: ');
%             if rank(Lpenalty_csd) > 0
                fmatrix_old_csd = (big_design_csd'*big_design_csd+lam_r.*eye(size(big_design_csd,2)))\big_design_csd'*big_DWI;
%             else 
%                 beta = fmatrix_old_csd;
%                 thresh_now = thresh_now;
%                 iter = count;
%                 sprintf('no scsd performed\n')
%                 rank(Lpenalty_csd)
%                 break;
%             end
        else
            %fmatrix_estimation = (design_matrix'*design_matrix + lambda*Penalty_matrix)\design_matrix'*DWI_simulated;
            fmatrix_old_csd = (big_design_csd'*big_design_csd)\big_design_csd'*big_DWI;
        end
        F_now_csd = L*fmatrix_old_csd;
        thresh_now = max(abs(F_now_csd-F_old_csd));
        F_old_csd = F_now_csd;
        count = count+1;
    else  %(needlet case)
        fmatrix_old_csd = (big_design_csd'*big_design_csd+lam_r.*eye(size(big_design_csd,2)))\big_design_csd'*big_DWI;
        F_now_csd = L*fmatrix_old_csd;
        thresh_now = max(abs(F_now_csd-F_old_csd));
        F_old_csd = F_now_csd;
        count = count+1;
    end
end

beta = fmatrix_old_csd;
% thresh_now = thresh_now;
iter = count;

end
