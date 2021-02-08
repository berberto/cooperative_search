function r = rhoinit (x, y)

    P = Param_maze;
    
%    if ((x - P.x)^2+(y - P.y)^2 < P.sigma^2)
%        r = 2/(pi*P.sigma^2)*(1 - ((x - P.x)^2 + (y - P.y)^2)/(P.sigma^2));
%    else
%       r = 0;
%   end
    
%     r = exp(-((x - P.x)^2+(y - P.y)^2)/(2*P.sigma^2));
    r = exp(-((x - P.x)^2+(y - P.y)^2)^2/(P.sigma)^4)/(2.78416*P.sigma^2);
end
