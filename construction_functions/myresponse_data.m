
%%%%%11-05-2011 Jun
%%%return the DWI for a single fiber.

function y = myresponse_data( b , ratio ,ev1, theta0 , phi0 , theta , phi)


%%% b control the scale , and around 0.6 is kind of optimal. 
%%% ratio  = e1 / e2;

%%%theta0  , phi0 control the rotation corresponding to theta and phi in
%%%the spherical coordinate system.

%%theta , phi are the function parameter.

%%%%e1>>e2
%%%%% e1 , e2 , e2 D = diag{e2 , e2 , e1}

%%%% The response function along z axis
%%%% R(theta) = exp(-u'Du)

e2 = ev1*1/ratio;
e1 = ev1*1;

D = [e2,0,0;0,e2,0;0,0,e1];


u = [sin(theta)*cos(phi) ; sin(theta)*sin(phi) ; cos(theta)];


%myresponse = @(theta , phi ) exp(-b.*( e2.*sin(theta).^2 + e1.*cos(theta).^2 )  )

%%the rotation matrix around x
%R_theta = [1,0,0;0,cos(theta0) , sin(theta0);0,-sin(theta0) , cos(theta0)];

%%the rotation matrix around z
%R_phi = [cos(phi0),sin(phi0),0;-sin(phi0) , cos(phi0),0;0,0,1];

%%% the rotation matrix

R = [cos(phi0)*cos(theta0),-sin(phi0),cos(phi0)*sin(theta0);sin(phi0)*cos(theta0),cos(phi0),sin(phi0)*sin(theta0);-sin(theta0),0,cos(theta0)];
%R_theta*R_phi;


%u'*R'*D*R*u

y = exp(-b*u'*R*D*R'*u);
%%%%%%%%%%plot


%theta = pi/(4*B) : 2*pi/(4*B) : (4*B-1)*pi/(4*B); 
%phi = 0 : pi/B : (2*B - 1)*2*pi/(2*B); 


%[theta , phi] = meshgrid(theta , phi);

%rho = myresponse(theta , phi);
%r = rho.*sin(theta);
%x = r.*cos(phi);    % spherical coordinate equations
%y = r.*sin(phi);
%z = rho.*cos(theta);


%surf(x,y,z);

%%%% SH coefficients of the z-aligned response function
%%%% for Jun to upadte: write this as .m file 