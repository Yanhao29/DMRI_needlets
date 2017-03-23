%%%11_15_2011 Jun
%%% Obtain the SH coefficients for the Dirac function.

function delta = Dirac_SH_coe(lmax, theta0 , phi0)
%%% Bmatrix is the (4B^2 * ((l_max+1)(l_max+2)/2)) matrix , the row
%%% corresponding to the evaluation of each SH basis on a particular point
%%% and column corresponding to the evaluation of a particular SH basis on
%%% the equal angle grid.

%%% theta0 and phi0 control the direction of the delta function.

%%% function returns delta , the SH coefficients of the Dirac delta
%%% function with the direction specified by theta0 and phi0

%% 
%%%determinet the largest possible SH basis order: depends on B
%if mod(B,2) == 1
%	lmax = B - 1;
%	else
%	lmax = B - 2;
%end


% theta0 = theta0 ; %why?
%phi0 = pi/2 + phi0;

%%


%%%
R = (lmax+1)*(lmax+2)/2;  %% number of real SH basis 
delta = zeros(R,1);


%%% use the function "myspharm.m" written by Jun
%%j theta , k phi 

for l = 0:2:lmax     %% even order SH
for m = (-l):l       %% SH phase      

Rindex = (l+1)*(l+2)/2 - (l-m); %% column/SH basis  index

if (m<0 && m>= (-l))  %% m<0: see Descoteaux 2007 equation 3

delta( Rindex) = ( (-1)^(m)*spharmonic_eval(l,m,theta0, phi0) + spharmonic_eval(l,-m,theta0, phi0))/sqrt(2);

elseif(m == 0)  %% m=0

delta(Rindex) = spharmonic_eval(l,m,theta0, phi0);

else  %%m>0

delta(Rindex) = 1i*( (-1)^(m+1)*spharmonic_eval(l,m,theta0, phi0) + spharmonic_eval(l,-m,theta0, phi0))/sqrt(2);

end

end
end


delta = real(delta);


end