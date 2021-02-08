% Define d coefficients -- Only when u(1) is the value!
function dc = dcoeff_val(region,state)
    
    P = Param_maze;
    nr = numel(region.x);
    
    % Allocate d
    dc = ones(P.N,nr);
    
    % Set d entries (u(1) is the value)
    dc(1,:) = P.eps*exp(state.u(1,:)/(2*P.D*P.gammat));
    
end
