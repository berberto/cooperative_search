function uinit = initcond(region)

    load('params');
    nr = numel(region.x);
    f = @rhoinit;
    
    uinit = zeros(P.N,nr);
    
    uinit(1,:) = P.v0;
    for i = 1:nr
        uinit(2,i) = f(region.x(i),region.y(i));
    end
        
end
