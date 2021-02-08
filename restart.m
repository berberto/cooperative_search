function uinit = restart(region)

    load('params');
    nr = numel(region.x);
    f = @rhoinit;
    
    uinit = zeros(P.N,nr);

    load('start');
    sol = u.NodalSolution;

    for i = 1:nr
        uinit(1,i) = sol(i,1,2);
        uinit(2,i) = f(region.x(i),region.y(i));
    end
    
end
