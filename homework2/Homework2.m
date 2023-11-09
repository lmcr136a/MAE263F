% PLEASE COMMENT A LOT
clear all;
close all;
clc

%% Global variables
global Fg mMat dt
global kappaBar EI voronoiLength % Bending
global GJ % Twisting
global EA refLen % Stretching

%% Inputs
nv = 50; % number of nodes
ne = nv - 1; % number of edges
ndof = 3*nv + ne; % number of DOF = 4*nv - 1
dt = 0.01; % Time step size

% Geometry
RodLength = 0.2;
natR = 0.02; % Natural radius
r0 = 0.001; % Cross-sectional radius
rho = 1000; % Density

% Material parameters
Y = 10e6; % Young's Modulus
nu = 0.5; % Poisson's ratio
G = Y / (2 * (1 + nu)); % Shear modulus

% Gravity
g = [0; 0; -9.81];
totalTime = 5; % Simulation time

%% Stiffness variables
EI = Y * pi * r0^4/4; % Bending stiffness
GJ = G * pi * r0^4/2; % Shearing stiffness
EA = Y * pi * r0^2; % Stretching stiffness

%% Tolerance
tol = EI / RodLength^2 * 1e-6;

%% Mass
totalM = pi * r0^2 * RodLength * rho;
dm = totalM / ne; % mass per edge
massVector = zeros(ndof, 1);

for c = 1:nv % Loop over nodes
    ind = [4*c-3; 4*c-2; 4*c-1];
    
    if c == 1
        massVector(ind) = dm/2;
    elseif c == nv
        massVector(ind) = dm/2;
    else
        massVector(ind) = dm;
    end
end

for c=1:ne % Loop over edges
    massVector(4*c) = 1/2 * dm * r0^2;
end

mMat = diag(massVector); % ndof x ndof sized mass matrix

%% Geometry
nodes = zeros(nv, 3);
dTheta = (RodLength / natR) * (1/ne);
for c = 1 :nv
    nodes(c,1) = natR * cos( (c-1) * dTheta);
    nodes(c,2) = natR * sin( (c-1) * dTheta);
    nodes(c,3) = 0;
end

%% Initial DOF vector
q0 = zeros(ndof, 1); % ndof = 4N-1
for c=1:nv
    ind = [4*c-3; 4*c-2; 4*c-1];
    q0(ind) = nodes(c,:);
end

%% Reference length (edge length)
refLen = zeros(ne, 1);
for c = 1:ne % loop over the edges
    dx = nodes(c+1,:) - nodes(c,:);
    refLen(c) = norm( dx );
end

%% Voronoi length (Length associated with each node)
voronoiLength = zeros(nv, 1);
for c=1:nv % loop over the nodes
    if c==1
        voronoiLength(c) = 0.5 * refLen(c);
    elseif c==nv
        voronoiLength(c) = 0.5 * refLen(c-1);
    else
        voronoiLength(c) = 0.5 * (refLen(c-1) + refLen(c));
    end
end

%% Gravity
Fg = zeros(ndof, 1);
for c=1:nv % loop over the nodes
    ind = [4*c-3; 4*c-2; 4*c-1];
    Fg(ind) = massVector(ind) .* g;
end

%% Reference frame (Space parallel transport at t=0)
a1 = zeros(ne, 3); % First reference director for all the edges
a2 = zeros(ne, 3); % Second reference director for all the edges
tangent = computeTangent( q0 ); % Tangent
t0 = tangent(1,:); % Tangent on first edge
t1 = [0;0;-1]; % "arbitrary" vector
a1Tmp = cross(t0, t1); % is perpendicular to both t0 and t1

if abs(a1Tmp) < 1e-6 % ==0?
    t1 = [0;1;0];
    a1Tmp = cross(t0, t1);
end
a1(1,:) = a1Tmp / norm( a1Tmp ); % Plug into my big a1 matrix
a2(1,:) = cross(tangent(1,:), a1(1,:));

% Done with the first edge
% Space parallel transport to construct the reference frame
for c=2:ne
    t0 = tangent(c-1,:); % tanget on (c-1)-th edge
    t1 = tangent(c,:); % tanget on c-th edge
    a1_0 = a1(c-1, :);
    a1_l = parallel_transport( a1_0, t0, t1);
    a1(c,:) = a1_l / norm( a1_l );
    a2(c,:) = cross( t1, a1(c,:) );
end

%% Material frame
theta = q0(4:4:end);
[m1, m2] = computeMaterialDirectors(a1, a2, theta);
%% Reference twist
refTwist = zeros(nv, 1); % Reference twist
% It is zero because I initialized a1, a2 using space parallel transport
%% Natural curvature
kappaBar = getkappa(q0, m1, m2);
%% Fixed and free DOFs
fixedIndex = 1:7;
freeIndex = 8:ndof;

%% Time stepping scheme
Nsteps = round(totalTime/dt);
ctime = 0; % current time (utility variable)
endZ = zeros(Nsteps, 1); % z-coordinate of last node
% Initialize position and velocity
q = q0; % Position
u = zeros(size(q)); % Velocity

plotname = "ElasticRod_t="+mat2str(0);
plotrod(q, a1, a2, m1, m2, ctime, plotname);

for timeStep = 1:Nsteps
    fprintf('Current time=%f\n', ctime);
    [q, u, a1, a2] = objfun(q0, u, a1, a2, freeIndex, ...
    tol, refTwist);
    ctime = ctime + dt;
    % Update q
    q0 = q;
    % Store
    endZ(timeStep) = q(end);
    if mod(timeStep, 1) == 0
        theta = q(4:4:end);
        [m1, m2] = computeMaterialDirectors(a1, a2, theta);

        plotname = "ElasticRod_t="+mat2str(timeStep);
        plotrod(q, a1, a2, m1, m2, ctime, plotname);

    end
end

%% Visualization
figure(2);
timearray = (1:1:Nsteps) * dt;
plot( timearray, endZ, 'ro-');
box on
xlabel('Time, t [sec]');
ylabel('z-coord of last node, \delta_z [m]');

