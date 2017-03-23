function inn_prod = spneedlet_inn_prod( j1, k1, B )
%compute the inner product between psi_j1,k1 and all the other psi_j2,k2

%compute alm first
l_max = ceil(B^(j1+1));
bw = l_max+1;
alm = zeros(bw, 2*bw-1);
l_st = floor(B^(j1-1));
l_en = l_max;
Nside = 2^max((ceil(log2(fix(B^(j1+1)))-1)), 0);
tp = pix2ang(Nside, 'nest', false);
Npix = 12*Nside^2;
lambda = 4*pi/Npix;
for l = l_st:l_en
    for m = -l:l
        alm(l+1, m+bw) = sqrt(lambda)*fun_b(l/B^j1, B)*conj(spharmonic_eval( l, m, tp{k1}(1), tp{k1}(2) ));
    end
end

%needlet transform
inn_prod = spneedlet_tran_alm( alm, l_max, B );

end

