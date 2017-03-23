%% Evaluate sphere harmoics on mesh corresponding to J an lmax (on sphere)
function Bmatrix = sharmonic(pos,lmax)


phi = atan2(pos(2,:),pos(1,:))/(2*pi);   %%phi
phi = phi+(phi<0);
theta = acos(pos(3,:))/(pi);                    %% theta

%%%
L = (lmax+1)*(lmax+2)/2;  %% number of real SH basis 
Bmatrix = zeros(size(theta,2) , L);

%%% use the function "spharmonic_eval.m" written by Jun
for at = 1:size(theta,2) %% the vertex

for l = 0:2:lmax     %% even order SH
for m = (-l):l       %% SH phase      

Rindex = (l+1)*(l+2)/2 - (l-m); %% column/SH basis  index

if (m<0 && m>= (-l))  %% m<0: see Descoteaux 2007 equation 3

Bmatrix(at , Rindex) = ( (-1)^(m)*spharmonic_eval(l,m,theta(at)*pi , phi(at)*2*pi ) + spharmonic_eval(l,-m,theta(at)*pi, phi(at)*2*pi ))/sqrt(2);

elseif(m == 0)  %% m=0

Bmatrix(at , Rindex) = spharmonic_eval(l,m,theta(at)*pi , phi(at)*2*pi );

else  %%m>0

Bmatrix(at , Rindex) = 1i*( (-1)^(m+1)*spharmonic_eval(l,m,theta(at)*pi , phi(at)*2*pi ) + spharmonic_eval(l,-m,theta(at)*pi , phi(at)*2*pi ))/sqrt(2);

end

end
end

end


Bmatrix= real(Bmatrix);