% Define f coefficients -- Only if u(1) is the value!
function fc = fcoeff_val(region,state)
    
    P = Param_maze;
    nr = numel(region.x);   % # points in region
    % Allocate f
    fc = zeros(2,nr);
    
    % Set f entries (u(1) is the value)
%     fc(1,:) = P.q + P.gt*state.u(2,:) -...          % interaction and time
%         .5*(state.ux(1,:).^2 + state.uy(1,:).^2);   % 1/2*grad(phi)^2
    fc(1,:) = - (P.q + P.gt*state.u(2,:)).*exp(state.u(1,:)/(2*P.D*P.gammat));
    
    % source of particles
    fc(2,:) = P.pumping*exp(-(region.x.^2 + region.y.^2)/(2*P.sigma ^2));
    
end
