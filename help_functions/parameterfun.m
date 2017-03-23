function z = parameterfun(x,y)
    sigma2 = x(1);
    s0 = x(2);
%     val = 5*log(sigma2)+sum(y.^2./(2*sigma2))+5*s0^2/(2*sigma2)-sum(log(besseli(0,y.*s0/sigma2)));
%     z = 5*log(sigma2)+sum(y.^2./(2*sigma2))+5*s0^2/(2*sigma2)-sum(y.*s0/sigma2-log(sqrt(2*pi*y.*s0/sigma2)));
%     z = double(z);

    temp = 0;
    for i=1:length(y)
        xx = y(i);
        ttemp = -log(xx)+log(sigma2)+(xx^2+s0^2)/(2*sigma2)-xx*s0/sigma2+log(sqrt(2*pi*xx*s0/sigma2));
        temp = temp+ttemp;
    end
    
    z = double(temp);
    
    