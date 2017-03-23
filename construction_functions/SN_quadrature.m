function psi = SN_quadrature(J, jmax, BW)

Nside = 2^(J-1);
pix = cell2mat((pix2ang(Nside,'nest',false)));
phi = pix(2,:);   %%phi
theta = pix(1,:);                    %% theta
%BW=2
psi = cell(jmax+1,1);
    %psi_j = zeros(12*(j+1)^2,size(pos,2))
    
for j = 1:(jmax+1)
    psi{j} = zeros(12*(4)^(j-1),size(phi,2));
end

for k = 1:size(phi,2)
    current_ND = spneedlet_eval(theta(k),phi(k),BW,jmax);
    for j = 1:(jmax+1)
        psi{j}(:,k) = current_ND{j};
    end
end

