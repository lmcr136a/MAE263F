%% 
function [terminal_vel, position_R2_y, vel_R2_y, turning_ang] = ...
    simul_visc_spheres_explicit(N, prob_num, ques_num, Rs, totalTime, dt)
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
visc = 1; % Pa-s

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
C = zeros(6,6);
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


%{
f1 = figure(1);
plot( q(1:2:end), q(2:2:end), 'ro-');
axis equal
drawnow
saveas(f1, "Figures/Problem"+prob_num+"_Q"+ques_num+file_postfix+"_t=0.0_explicit.png")
%}
% Time marching scheme
for c=2:Nsteps
    q = q0; % Guess
    
    dE = zeros(2*N,1);

    % Es 1
    xk = q(1);
    yk = q(2);
    xkp1 = q(3);
    ykp1 = q(4);
    dE(1:4) = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
    % Es 2
    xk = q(3);
    yk = q(4);
    xkp1 = q(5);
    ykp1 = q(6);
    dE(3:6) = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
    % Eb
    xkm1 = q(1);
    ykm1 = q(2);
    xk = q(3);
    yk = q(4);
    xkp1 = q(5);
    ykp1 = q(6);
    curvature0 = 0;
    dEb = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, ...
    curvature0, deltaL, EI);
    dE = dE + dEb;
 
    q = q0 + dt *(u - (dt * inv(M))*(-W + C * u + dE));
    % Update
    u = (q - q0) / dt; % Velocity
    
    %fprintf('Time: %f | Y-vel R1=%f, R2=%f, R3=%f\n',...
    %    (c-1) * dt, u(2), u(4), u(6) );


    q0 = q;

    % Q1
    % c is a index num, when c=1 -> q0, when c=2 -> t=0.01
    if ismember(c, showing_steps)
        fprintf('Time: %f | Y-vel R1=%f, R2=%f, R3=%f\n',...
            (c-1) * dt, u(2), u(4), u(6) );
        f1 = figure(1);
        plot( q(1:2:end), q(2:2:end), 'ro-');
        axis equal
        drawnow
        saveas(f1, "Figures/Problem"+prob_num+"_Q"+ques_num+"_t="+(c-1)*dt+"_explicit.png")
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
filename = "Figures/Problem"+prob_num+"_Q"+ques_num+"_t="+totalTime+"_explicit.png";
saveas(f1, filename);

terminal_vel = u;

end