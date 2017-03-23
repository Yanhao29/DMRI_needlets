function Rmatrix = Response_Rmatrix_construction(b, ratio, J, lmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%specify response function to use when estimating FOD
%%% maybe misspecified 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%parameter seting for specified response function
w_response = 1;
b_response = b;   %%could be misspecified 
ratio_response = ratio; %%could be misspecified 

theta0_response = 0;  %%align with z-axis: this is always the case 
phi0_response = 0;
%% get wavelet grid corresponds to J
options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 1;

[vertex_sampling,~] = compute_semiregular_sphere(J,options);
pos = vertex_sampling{end};
%options.spherical = 1;
%plot_mesh(pos,face_sampling{end}); axis tight;  %plot the generated mesh
%lighting flat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%from perform_spherical_interpolation: (x,y,z) --> spehrical coordinates (r,\theta, \phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = atan2(pos(2,:),pos(1,:))/(2*pi);   %%phi
phi = phi+(phi<0);
theta = acos(pos(3,:))/(pi);                    %% theta

%%% evaluation of the specified R on the wavelet grid J
R = zeros(size(theta,2) , 1);

for at = 1:size(theta,2) %% the vertex

R(at) = myresponse_crossing(b_response,ratio_response, w_response, theta0_response , phi0_response, theta(at)*pi , phi(at)*2*pi  );

end

% path.cur = '/Users/hao/Dropbox/stats_project/FOD_codes_simulation/';
% addpath(strcat(path.cur,'matrix'));
% load(strcat('SH_vertex_J', num2str(J), '_lmax', num2str(lmax), '.mat'));
% Bmatrix = SH_vertex;
% clearvars SH_vertex;
% Bmatrix = sharmonic_vertex_all(J,lmax);
Bmatrix = SH_vertex(J, lmax, 0);
Rmatrix = ((Bmatrix'*Bmatrix)\(Bmatrix'*R));   

%%% get the block diagonal R matrix for the response function:
Rmatrix_temp = zeros(lmax/2+1 , 1);   

for l = 0:2:lmax

    Rmatrix_temp(l/2+1) = Rmatrix((l+1)*(l+2)/2-l);  %%get m=0 coefficients 

end

%%%  each block corresponds to an order l
%%% and all elements in the block are the same, corresponding to m=0 phase 
for l = 0:2:lmax
    for m = (-l):l

        Rmatrix((l+1)*(l+2)/2 - (l-m)) = Rmatrix_temp(l/2+1)*sqrt(4*pi / (2*l+1));

    end
end

Rmatrix = diag(Rmatrix);