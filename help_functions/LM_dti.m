function [beta, iter, fitted] = LM_dti(X,Y,beta_ini,thresh)
% levenberg marquardt algorithm

    n = length(Y);
    p = length(beta_ini);

    beta_c = beta_ini;
    beta_new = zeros(p,1); 
    iter = 0;

    while norm(beta_c-beta_new)>thresh

        beta_c = beta_new;
        iter = iter+1;

        temp = exp(-X*beta_c);
        y_c = Y-temp;
        design_matrix = -X.*temp(:,ones(1,p));
        diff = (design_matrix'*design_matrix)\(design_matrix')*y_c;
        beta_new = beta_c+diff;

    end

    beta = beta_new;
    fitted = design_matrix*beta;
end