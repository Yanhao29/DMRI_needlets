%% ridge penalty matrix construction
function penalty_matrix = penalty_matrix_ridge(lmax)
    penalty_matrix = zeros( (lmax+1)*(lmax+2)/2 , 1);
    for l = 0:2:lmax
        for m = (-l):l
            penalty_matrix((l+1)*(l+2)/2 - (l-m)) = l^2*(l+1)^2;  %% more penalty on higher order SH basis
        end
    end
    penalty_matrix = diag(penalty_matrix);
end
