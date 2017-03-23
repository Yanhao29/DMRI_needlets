

%%%%11-04-2011 Jun
%%%% add Rician noise on each points of the orignal grid (2B)*(2B)

function y = add_Rician_noise(dwi , sigma)

%%%sigma is the std of the error
%%% dwi is a vector or matrix of true dwis without noise , (4B^2)*1 or (2B)*(2B)
%%% We add Rician noise on each of it.

[a,b] = size(dwi);

errors = sigma*randn(2,a*b);

y = dwi;

for k = 1:a*b
	y(k) = sqrt( ( y(k) + errors(1,k) )^2 + errors(2,k)^2 );
end
