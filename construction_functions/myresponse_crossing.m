

%%%%11-05-2011 Jun
%%%%%return the DWI for the crossing fiber cases.


function y = myresponse_crossing( b , ratio ,w, theta0 , phi0 , theta , phi)

%%%%The number of the fiber.

num_fiber = length(b);
y = 0;
for i = 1:num_fiber
    y  =y+ w(i)*myresponse(b(i) , ratio(i) , theta0(i) , phi0(i) , theta , phi);
end

