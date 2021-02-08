function pars = Param_maze()    % gives a structure with fields (.)

    %Parameters
    pars.N = 2;       % number of equations
    pars.D = 1;       % diffusion constant
    pars.gt = 0;    % interaction constant
    pars.q = 10;      % time cost
    pars.eps = 1e-1;  % time constant (fast timescale)
    pars.phdep = 1.0; % prop.const. between outgone density and pheromone deposited
    pars.reg = 1e-3;  % parameter to regularize logarithm
    pars.pumping = 0; % maximum value density source
    
    pars.gamma = 1.0;   % control cost parameter
    pars.alpha = 0.0;     % risk sensitivity parameter
    
    pars.gammat = pars.gamma/(1.0 - 2.0 * pars.D * pars.alpha * pars.gamma);
    
    
    % position and width of density at time 0
    pars.x = 1.;
    pars.y = 1.;
    pars.sigma = .1;
    % normalization factor
    pars.norm = 1.;
    
    % value at time 0
    pars.v0 = 0;

end
