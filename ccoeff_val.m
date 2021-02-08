% Define c coefficients
function cc = ccoeff_val(region,state)

    P = Param_maze;
    n1 = 16;                % # entries needed in c tensor
    nr = numel(region.x);   % # points in region

    % Allocate c
    cc = zeros(n1,nr);
    
    % Set c entries (u(1) is the value)
    cc(1,:) = P.D*exp(state.u(1,:)/(2*P.D*P.gammat));
    cc(4,:) = cc(1,:);
    cc(5,:) = - state.u(2,:)/P.gamma;
    cc(8,:) = cc(5,:);
    cc(13,:) = P.D*ones(1,nr);
    cc(16,:) = cc(13,:);

end
