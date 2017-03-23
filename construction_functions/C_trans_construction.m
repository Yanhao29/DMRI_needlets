%% lmax give s the SH level

function C_trans = C_trans_construction(lmax)

%     path.cur='/Users/hao/Dropbox/stats_project/DTI/';
%     path.save = (strcat(path.cur,'matrix'));
%     addpath (strcat(path.cur,'construction_functions/'));
    
    B = 2; % bandwidth

    L = (lmax+1)*(lmax+2)/2;
    jmax = fix( log(lmax)/log(B)+1 );
    if abs(jmax-(log(lmax)/log(B)+1))<eps
        jmax = jmax-1;
    end
    N = 12*(1-4^(jmax+1))/(1-4);

    %do the spherical harmonic transform
    bw = lmax+1;

    C_trans = zeros(L,N+1);
    for l = 0:2:lmax
        for m = (-l):l
            SHindex = (l+1)*(l+2)/2 - (l-m);
            coef = zeros(bw, 2*bw-1);
            if m<0
                coef(l+1, m+bw) = (-1)^m/sqrt(2);
                coef(l+1, -m+bw) = 1/sqrt(2);
            elseif m>0
                coef(l+1, m+bw) = 1i*(-1)^(m+1)/sqrt(2);
                coef(l+1, -m+bw) = 1i/sqrt(2);
            else
                coef(l+1,bw) = 1;
            end
            sp_coef = spneedlet_tran_alm(coef,lmax,B);
            if(l==0&&m==0)
                sp_coef_vec = 1;
            else
                sp_coef_vec = 0;
            end

            for i =1:(jmax+1)
                sp_coef_vec = [sp_coef_vec,sp_coef{i}];
            end

            C_trans(SHindex,:) = sp_coef_vec;
        end
    end

%     filename = fullfile(path.save, strcat('C_lmax',num2str(lmax),'.mat' ));
%     save(filename,'C_trans');
end