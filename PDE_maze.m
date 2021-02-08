%% 
clear
clc

dir = 'output_maze';

% Set parameters -- Parameters will not be global, do it in all functions!
P = Param_maze;

if (2.0 * P.D * P.alpha * P.gamma) >= 1
    fprintf('Invalid risk-sensitivity parmeter\n');
    exit(1);
end

% Create a system of 2 PDE
N = P.N;
model = createpde(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   GENERATE THE GEOMETRY, FROM FILE IMPORT, PDETOOL OR MANUALLY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Creation of the geometry - From PDETOOL (ONLY FOR NEW GEOMETRY)
    
% Open pdetool
% run `pdetool' in the command window

% Draw and then export data (graphics data, formula, names)
% Draw > Export...
% save('geometry_M3','gd','sf','ns');

% import geometry
%load('geometry_M3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   VILLA PISANI

%% Create geometry from polygons in files
 
% columns are x and y coordinates of the paths
P1 = importdata('centered_1.dat');
P2 = importdata('centered_2.dat');
P3 = importdata('centered_3.dat');
P4 = importdata('centered_4.dat');
P5 = importdata('centered_5.dat');
P6 = importdata('centered_6.dat');
P7 = importdata('centered_7.dat');
P8 = importdata('centered_8.dat');
P9 = importdata('centered_9.dat');
    side = max(P1(:,1)) - min(P1(:,1)); % size of the largest (outer) polygon
    sideY = max(P1(:,2)) - min(P1(:,2));

% create csg polygons
P1 = [2; length(P1); P1(:,1)/side; P1(:,2)/side];
P2 = [2; length(P2); P2(:,1)/side; P2(:,2)/side];
P3 = [2; length(P3); P3(:,1)/side; P3(:,2)/side];
P4 = [2; length(P4); P4(:,1)/side; P4(:,2)/side];
P5 = [2; length(P5); P5(:,1)/side; P5(:,2)/side];
P6 = [2; length(P6); P6(:,1)/side; P6(:,2)/side];
P7 = [2; length(P7); P7(:,1)/side; P7(:,2)/side];
P8 = [2; length(P8); P8(:,1)/side; P8(:,2)/side];
P9 = [2; length(P9); P9(:,1)/side; P9(:,2)/side];

% columns are polygons (unneeded entries are 0)
gd = zeros(numel(P1),9);
    gd(1:numel(P1),1) = P1;
    gd(1:numel(P2),2) = P2;
    gd(1:numel(P3),3) = P3;
    gd(1:numel(P4),4) = P4;
    gd(1:numel(P5),5) = P5;
    gd(1:numel(P6),6) = P6;
    gd(1:numel(P7),7) = P7;
    gd(1:numel(P8),8) = P8;
    gd(1:numel(P9),9) = P9;

% define formula for combining polygons (outer - inner(s))
sf = 'P1 - (P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9)';

% define correspondence between polygons and their names in the formula
ns = zeros(2,9);
    ns(1,:) = double('P');
    ns(2,:) = double('1'):double('9');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   CIRCLE

%%  Create geometry of square with round target
% xmin=-5;
% ymin=-5;
% xmax=5;
% ymax=5;
% rect = [3; 4; xmin; xmax; xmax; xmin; ymin; ymin; ymax; ymax];
%     side = max(rect(3:6))-min(rect(3:6));
%     sideY = max(rect(7:10))-min(rect(7:10));
% circ = zeros(length(rect),1);
%     rad = 0.2;
%     circ(1:4) = [1; 0; 0; rad];

% gd = [rect circ];
% sf = 'R - C';
% ns = [double('R') double('C')];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Assign the geometry to the model and plot it

% create geometry of the model
g = decsg(gd,sf,ns);                % decompose geometry into minimal regions
geom = geometryFromEdges(model,g);  % assign Geometry property to model

%Generate  mesh
generateMesh(model,'Hmax',.1);  % set mesh with specified maximum edge length
[p,e,t] = meshToPet(model.Mesh); % saves p (points) e (edges) t (triangles) of the mesh
p = jigglemesh(p,e,t);           % improve quality of the mesh

% save mesh variables for usage in external functions
save('mesh','p','t','e');

%% Plot the geometry (leave commented in simulations)
% pdeplot(model)
% grid on
% axis equal
%     
% hold on
% pdegplot(model,'EdgeLabels','on');
% axis equal


%% Set coefficients of the PDE

% set target edge(s)
target = 5:8;
others = setxor(1:geom.NumEdges, target);

%% set mesh on the target (straight element, exit of Villa pisani)
et = e(:,e(5,:)==target);
pt = p(:,et(1,:));
npt = 20;
dpt = pt(:,2) - pt(:,1);
delta = dpt*size(pt,2)/npt;
delta = delta.*ones(1,npt);
for i = 0:(npt-1)
    xt = pt(1,1)*ones(1,npt) + i*delta(1,i+1);
    yt = pt(2,1)*ones(1,npt) + i*delta(2,i+1); 
end

%% set mesh on circular target (for the circle problem)
% npt = 240;
% dtheta = 2*pi/npt;
% theta = 0:dtheta:2*pi;
% delta = zeros(2,npt);
% xt = zeros(1,npt);
% yt = zeros(1,npt);
% for i = 1:npt
%     xt(i) = 1.01*rad*cos((i-.5)*dtheta);
%     yt(i) = 1.01*rad*sin((i-.5)*dtheta);
%     delta(:,i) = 1.01*rad*[cos(i*dtheta)-cos((i-1)*dtheta);
%                       sin(i*dtheta)-sin((i-1)*dtheta)];
% end


% Set boundary conditions (separately for target and others)
Q = [0 0; 0 0];
G = [0;0];
BCtarget = applyBoundaryCondition(model,'mixed','edge', target,...
                    'u',0,'EquationIndex',2, 'q',Q, 'g',G);
BCothers = applyBoundaryCondition(model,'neumann','edge', others);
BCs = model.BoundaryConditions;

% Calculate normalization of density
Nin0 = 0;
f = @rhoinit;
for triangle = 1:length(t)
   % vertices numbers (1 -> 2 -> 3 goes counterclockwise):
   v1n=t(1,triangle);
   v2n=t(2,triangle);
   v3n=t(3,triangle);
   % vertices coordinates:
   v1=p(:,v1n);
   v2=p(:,v2n);
   v3=p(:,v3n);
   % calculate area of the triangle
   d12 = v2 - v1;
   d13 = v3 - v1;
   Atriangle = .5*(d12(1)*d13(2) - d12(2)*d13(1));
   % calculate approximate integral
   Nin0 = Nin0 + Atriangle*(f(v1(1),v1(2)) + f(v2(1),v2(2)) + f(v3(1),v3(2)))/3;
end
P.norm = Nin0;
save('params','P');
disp(P.norm);


%Set coefficients and initial conditions
u0 = @initcond;
m = 0;
d = @dcoeff_val;
a = 0;
c = @ccoeff_val;
f = @fcoeff_val;

specifyCoefficients(model,'m',m,...
                          'd',d,...
                          'c',c,...
                          'a',a,...
                          'f',f);
Cs = model.EquationCoefficients;


%% Solve and save the solution of the PDE

% Rectangular mesh (for interpolation)
n1 = 250;
n2 = floor(n1*sideY/side);
X = linspace(min(p(1,:)),max(p(1,:)),n1); %min(p(1,:)):0.04:max(p(1,:));
Y = linspace(min(p(2,:)),max(p(2,:)),n2); %min(p(2,:)):0.04:max(p(2,:));
[xx,yy] = meshgrid(X,Y);

% Time steps
dt = .002;
tf = 2.5;
times = 0:dt:tf;
subst = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   for circle only (for stochastic simulations in C)
%
% % Output directory
% dir = sprintf('%s/circle_long',dir);
% mkdir(dir);
%
% % Write parameters on C header file
% outname = sprintf('%s/parameters.h', dir);
% outpars = fopen(outname, 'w');
% fprintf(outpars, '#ifndef PARAMETERS_H\n#define PARAMETERS_H\n\n');
% fprintf(outpars, '#define dt\t%f\n', dt);      % time step for the field
% fprintf(outpars, '#define D\t%f\n', P.D);      % diffusion constant
% fprintf(outpars, '#define gamma\t%f\n', P.gamma);      % control cost rate
% fprintf(outpars, '#define q\t%f\n', P.q);      % time cost rate
% fprintf(outpars, '#define R\t%f\n', rad);      % radius of the target
% fprintf(outpars, '#define Lx\t%d\n', n1);      % lattice horizontal side
% fprintf(outpars, '#define Ly\t%d\n', n2);      % lattice vertical side
% fprintf(outpars, '#define x0\t%f\n', P.x);     % initial center of particle distribution -- x
% fprintf(outpars, '#define y0\t%f\n', P.y);     % initial center of particle distribution -- y
% fprintf(outpars, '#define s0\t%f\n', P.sigma); % initial size of particle distribution
% fprintf(outpars, '#define xmin\t%f\n', xmin);  
% fprintf(outpars, '#define xmax\t%f\n', xmax);
% fprintf(outpars, '#define ymin\t%f\n', ymin);
% fprintf(outpars, '#define ymax\t%f\n', ymax);
% fprintf(outpars, '\n#endif\n');
% fclose(outpars);


% Open output files
details = sprintf('nosource_rs%0.2f_D%0.2f_q%0.1f_gt%0.1f',...
                                     P.alpha, P.D, P.q, P.gt);

outname = sprintf('%s/%s_sol.dat', dir, details);
fprintf('Output:  %s\n\n', outname);
outsol = fopen(outname,'w');

outname = sprintf('%s/%s_out.dat', dir, details);
outout = fopen(outname,'w');

outname = sprintf('%s/%s_cos.dat', dir, details);
outcos = fopen(outname,'w');

outname = sprintf('%s/%s_phe.dat', dir, details);
outphe = fopen(outname,'w');

outname = sprintf('%s/lattice.dat', dir);
outlat = fopen(outname,'w');

outname = sprintf('%s/%s_con.dat', dir, details);
outcon = fopen(outname,'w');


% Print lattice coordinates
for i = 1:length(X)
    fprintf(outlat,'%0.6e\n', X(i));
end
for j = 1:length(Y)
    fprintf(outlat,'%0.6e\n', Y(j));
end
fclose(outlat);


%Solve the PDE
u = u0;
stepcount = 0;
phrmn=0;
time = 0;
source = P.pumping*sqrt(2*pi)*P.sigma;  % total # particles pumping rate
check = true;
while time < tf  %check
    
    % last solution is the new initial condition
    setInitialConditions(model,u);
    
    % solve in [t,t+dt]
    u = solvepde(model, time:dt/subst:time+dt);
    
    % store the solution and the gradient of the value
    solt = u.NodalSolution;
    conx = u.XGradients(:,1);
    cony = u.YGradients(:,1);
    
    % local rates of cost
    ccon = 1/P.gamma*(conx.^2 + cony.^2).*solt(:,2);  % control cost
    ctim = P.q*solt(:,2);                             % time cost
    ccol = P.gt/2*(solt(:,2).^2);                     % collision cost
    
    % Calculate fraction of particles inside
    % (Desirability is the number of particles which exited the domain)
    % and space integral of control cost, time cost and collision cost
    Nin = 0;
    Ccon = 0;
    Ctim = 0;
    Ccol = 0;
    for triangle = 1:length(t)
       % vertices numbers (1 -> 2 -> 3 goes counterclockwise):
       v1n=t(1,triangle);
       v2n=t(2,triangle);
       v3n=t(3,triangle);
       % vertices coordinates:
       v1=p(:,v1n);
       v2=p(:,v2n);
       v3=p(:,v3n);
       % calculate area of the triangle
       d12 = v2 - v1;
       d13 = v3 - v1;
       Atriangle = .5*(d12(1)*d13(2) - d12(2)*d13(1));
       % calculate approximate integrals
       Nin  = Nin + Atriangle*(solt(v1n,2) + solt(v2n,2) + solt(v3n,2))/3;
       Ccon = Ccon + Atriangle*(ccon(v1n) + ccon(v2n) + ccon(v3n))/3;
       Ctim = Ctim + Atriangle*(ctim(v1n) + ctim(v2n) + ctim(v3n))/3;
       Ccol = Ccol + Atriangle*(ccol(v1n) + ccol(v2n) + ccol(v3n))/3;
    end
    fprintf(outout, '%0.6e\t%0.6e\n', time, 1-Nin);
    fprintf(outcos, '%0.6e\t%0.6e\t%0.6e\t%0.6e\n', time, Ccon, Ctim, Ccol);
    
    % calculate outgoing particle flux at from the target
    [Jxt, Jyt] = evaluateCGradient(u,xt,yt,2,1);
    Jxt = reshape(-Jxt,size(xt));
    Jyt = reshape(-Jyt,size(yt));
    outfl = 0;
    for i = 1:length(xt)
        outfl = outfl + Jyt(i)*delta(1,i) - Jxt(i)*delta(2,i);
    end
    % pheromone deposited is the outgoing flux integrated in time
    phrmn = phrmn + P.phdep*outfl*dt;
    for i = 1:length(xt)
        fprintf(outphe, '%0.6e\t%0.6e\t%0.6e\n', xt(i), yt(i), phrmn);
    end
    fprintf(outphe, '\n\n');
   
    % update Neumann BC at the target
    G = [phrmn; 0];
    BCtarget.g = G;
    
    % iteration goes on until check is true:
    % - if no pumping, number of particles inside is below threshold;
    % - otherwise, difference between injected and outgoing particles per
    % time step (relative to pumping) below threshold (true at st. state)
    if P.pumping == 0
        fprintf('t = %0.3f\tphe = %0.5f\tNout/N = %0.5f\n',...
                        time, phrmn, 1-Nin);
        check = Nin > .03;
    else
        nst = abs(source-outfl)/source;
        fprintf('t = %0.3f\tphe = %0.5f\tNout/N = %0.5f\tn-st = %0.5f\n',...
                        time, phrmn, 1-Nin, nst);
        check = nst > .001;
    end
    
    % evaluate interpolating function on the mesh
    M1 = interpolateSolution(u,xx,yy,1,1);
    M2 = interpolateSolution(u,xx,yy,2,1);
    M1 = reshape(M1,size(xx))';
    M2 = reshape(M2,size(xx))';
    
    % evaluate the control (grad(value)/gamma)
    [Ux,Uy] = evaluateGradient(u,xx,yy,1,1);
    Ux = reshape(1/P.gamma*Ux,size(xx))';
    Uy = reshape(1/P.gamma*Uy,size(xx))';
    Ux(isnan(Ux)) = 1000;   Uy(isnan(Uy)) = 1000;   

    % evaluate the flux of particle density
    [Jx,Jy] = evaluateCGradient(u,xx,yy,2,1);
    Jx = reshape(-Jx,size(xx))';
    Jy = reshape(-Jy,size(xx))';

    % print value, density and current density
    % print drift field
    fprintf(outsol, '#\t%d\n', stepcount);
    for j = 1:length(Y)
        for i = 1:length(X)
            fprintf(outsol,'%0.6e\t%0.6e\t%0.6e\t%0.6e\t%0.6e\t%0.6e\n',...
                    X(i), Y(j), M1(i,j), M2(i,j),...
                    Jx(i,j), Jy(i,j) );
            fprintf(outcon, '%0.6e\t%0.6e\n', Ux(i,j), Uy(i,j));
        end
        fprintf(outsol,'\n');
    end
    fprintf(outsol,'\n\n');
    fprintf(outcon,'\n');

    stepcount = stepcount + 1;
    time = time + dt;

end
fclose(outsol);
fclose(outcon);
fclose(outcos);
fclose(outphe);
fclose(outout);


%%
% Plot the solution on separate figures
% for i = 1:30:steps
%     figure                                  % plot on new figure
%     
%     % plot solution 1 (value)
%     valplot = subplot(2,1,1);
%     colormap(valplot, winter)
%     title('Value');
%     pdeplot(model, 'xydata', sol(:,1,i));
% 
%     % plot solution 2 (1-particle density)
%     rhoplot = subplot(2,1,2);
%     colormap(rhoplot, jet)
%     title('Density');
%     pdeplot(model,'xydata', sol(:,2,i)); %, 'zdata', sol(:,2,i));
% end
