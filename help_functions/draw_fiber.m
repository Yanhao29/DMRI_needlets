
%%%add fiber direction onto the plot
function draw_fiber(theta0 , phi0, scale, rho)

    n = numel(theta0);
    theta0 = reshape(theta0,numel(theta0),1);
    phi0 = reshape(phi0,numel(phi0),1);
    
    xx = 0:scale/100:scale;
    for i = 1:n
        x = rho*1.1/scale*xx*sin(theta0(i))*cos(phi0(i));
        y = rho*1.1/scale*xx*sin(theta0(i))*sin(phi0(i));
        z = rho*1.1/scale*xx*cos(theta0(i));
        plot3(x,y,z);
        plot3(-x,-y,-z);
    end

end