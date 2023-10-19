%% 
function [terminal_vel, position_R2_y, vel_R2_y, turning_ang] = ...
    simul_visc_spheres(N, prob_num, ques_num, Rs, totalTime, dt)
%% Given Parameters
R1 = Rs(1);
R2 = Rs(2);
R3 = Rs(3);
rho_metal = 7000;
rho_f = 1000;
rho = rho_metal - rho_f;

RodLength = 0.1;
deltaL = RodLength / (N-1);

r0 = 0.001;
Y = 1e9; % Using Y instead of E to avoid ambiguity
g = 9.8; % m/s^2
visc = 1000; % Pa-s

EI = Y * pi * r0^4 / 4;
EA = Y * pi * r0^2;

nodes = zeros(N, 2);
for c = 1:N
    nodes(c,1) = (c-1) * deltaL;
end

% Tolerance
tol = EI / RodLength^2 * 1e-3;

%% Matrices
% Mass matrix
M = zeros(2*N,2*N);
M(1,1) = 4/3*pi*R1^3*rho_metal;
M(2,2) = 4/3*pi*R1^3*rho_metal;
M(3,3) = 4/3*pi*R2^3*rho_metal;
M(4,4) = 4/3*pi*R2^3*rho_metal;
M(5,5) = 4/3*pi*R3^3*rho_metal;
M(6,6) = 4/3*pi*R3^3*rho_metal;
% Viscous damping matrix
C = zeros(2*N,2*N);
C1 = 6*pi*visc*R1;
C2 = 6*pi*visc*R2;
C3 = 6*pi*visc*R3;
C(1,1) = C1;
C(2,2) = C1;
C(3,3) = C2;
C(4,4) = C2;
C(5,5) = C3;
C(6,6) = C3;
% Gravity
W = zeros(2*N,1);
W(2) = -4/3*pi*R1^3*rho*g;
W(4) = -4/3*pi*R2^3*rho*g;
W(6) = -4/3*pi*R3^3*rho*g;

%% About DOF
% Initial DOF vector
q0 = zeros(2*N,1);
for c=1:N
    q0 ( 2*c - 1 ) = nodes(c,1); % x coordinate
    q0 ( 2*c ) = nodes(c,2); % y coordinate
end

% New position and velocity
q = q0; % DOF vector
u = (q - q0) / dt; % Velocity vector

%% Simul on time steps
Nsteps = round( totalTime / dt );

% Answer setting
showing_ts = [0, 0.01, 0.05, 0.1, 1.0, 10];
showing_steps = showing_ts / dt +1;

turning_ang = zeros(Nsteps, 1);
turning_ang(1) = 0;
position_R2_y = zeros( Nsteps, 1); % y-position of R2
position_R2_y(1) = q(4);
vel_R2_y = zeros( Nsteps, 1); % y-velocity of R2
vel_R2_y(1) = u(4);
%
if R1 == R2
    file_postfix = "_R="+mat2str(R2);
else
    file_postfix = "";
end

%
f1 = figure(1);
plot( q(1:2:end), q(2:2:end), 'ro-');
axis equal
drawnow
saveas(f1, "Figures/Problem"+prob_num+"_Q"+ques_num+file_postfix+"_t=0.0.png")

% Time marching scheme
for c=2:Nsteps
    q = q0; % Guess
    % Newton Raphson
    err = 10 * tol;
    while err > tol
        % Inertia
        f = M / dt * ( (q-q0) / dt - u );
        J = M / dt^2;
        %
        % Elastic forces
        %
        % Linear spring 1 between
        xk = q(1);
        yk = q(2);
        xkp1 = q(3);
        ykp1 = q(4);
        dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
        dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
        f(1:4) = f(1:4) + dF;
        J(1:4,1:4) = J(1:4,1:4) + dJ;

        % Linear spring 2 between nodes 2 and 3
        xk = q(3);
        yk = q(4);
        xkp1 = q(5);
        ykp1 = q(6);

        dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
        dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
        f(3:6) = f(3:6) + dF;
        J(3:6,3:6) = J(3:6,3:6) + dJ;

        % Bending spring between nodes 1, 2, and 3
        xkm1 = q(1);
        ykm1 = q(2);
        xk = q(3);
        yk = q(4);
        xkp1 = q(5);
        ykp1 = q(6);
        curvature0 = 0;
        dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
        curvature0, deltaL, EI);
        dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
        curvature0, deltaL, EI);
        f(1:6) = f(1:6) + dF;
        J(1:6,1:6) = J(1:6,1:6) + dJ;

        % Viscous force
        f = f + C * ( q - q0 ) / dt;
        J = J + C / dt;

        % Weight
        f = f - W;
        
        % Update
        q = q - J \ f;
        err = sum( abs(f) );
    end
    % Update
    u = (q - q0) / dt; % Velocity
    
    fprintf('Time: %f | Y-vel R1=%f, R2=%f, R3=%f\n',...
       (c-1) * dt, u(2), u(4), u(6) );

    q0 = q; % Old position
    % Q1
    % c is a index num, when c=1 -> q0, when c=2 -> t=0.01
    if ismember(c, showing_steps)
        f1 = figure(1);
        plot( q(1:2:end), q(2:2:end), 'ro-');
        axis equal
        drawnow
        saveas(f1, "Figures/Problem"+prob_num+"_Q"+ques_num+file_postfix+"_t="+(c-1)*dt+".png")
    end
    % Store
    position_R2_y(c) = q(4);
    vel_R2_y(c) = u(4);

    theta1 = atan2(q(6)-q(4), q(5)-q(3));
    theta2 = atan2(q(4)-q(2), q(3)-q(1));
    turning_ang(c) = (theta1+theta2)*180/pi;
end

f1 = figure(1);
plot( q(1:2:end), q(2:2:end), 'ro-');
axis equal
drawnow
saveas(f1, "Figures/Problem"+prob_num+"_Q"+ques_num+file_postfix+"_t="+totalTime+".png")
terminal_vel = u;

end