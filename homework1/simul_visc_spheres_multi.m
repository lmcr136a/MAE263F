%% 
function [terminal_vel, position_R2_y, vel_R2_y, turning_ang] = ...
    simul_visc_spheres_multi(N, prob_num, ques_num, totalTime, dt)
%% Given Parameters
ndof = N*2;
rho_metal = 7000;
rho_f = 1000;
rho = rho_metal - rho_f;

RodLength = 0.1; 
deltaL = RodLength / (N-1);

% Radii of spheres
R = zeros(N,1);
R(:) = deltaL/10;
midNode = (N+1)/2;
R(midNode) = 0.025;

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

%% 
% Mass, M
M = zeros(ndof, ndof);
for k=1:N
    M(2*k-1, 2*k-1) = 4/3*pi*R(k)^3*rho_metal; % Mass for x_k
    M(2*k, 2*k) = M(2*k-1, 2*k-1); % Mass for y_k
end

C = zeros(ndof,ndof);
for k=1:N
    C(2*k-1, 2*k-1) = 6 * pi * visc * R(k);
    C(2*k, 2*k) = C(2*k-1, 2*k-1);
end

% Weight vector, W
W = zeros(ndof, 1);
for k=1:N
    W(2*k-1) = 0; % weight along x is zero
    W(2*k) = -4/3*pi*R(k)^3*rho*g;
end

% Initial DOF
q0 = zeros(ndof, 1);
for c=1:N % loop over nodes
    q0( 2*c-1 ) = nodes(c,1); % x1, x2, x3
    q0( 2*c ) = nodes(c,2); % y1, y2, y3
end

u0 = zeros(ndof, 1); % old velocity (initial velocity)

% New position and velocity
q = q0; % DOF vector
u = (q - q0) / dt; % Velocity vector

% tolerance
tol = EI/RodLength^2 * 1e-3; % small enouch force that can be neglected

%% Simul on time steps
Nsteps = round( totalTime / dt );

% Answer setting
showing_ts = [0, 1, 10, 50];
showing_steps = showing_ts / dt +1;

turning_ang = zeros(Nsteps, 1);
turning_ang(1) = 0;
position_R2_y = zeros( Nsteps, 1); % y-position of R2
position_R2_y(1) = q(4);
vel_R2_y = zeros( Nsteps, 1); % y-velocity of R2
vel_R2_y(1) = u(4);

%
f1 = figure(1);
plot( q(1:2:end), q(2:2:end), 'ro-');
axis equal
drawnow
saveas(f1, "Figures/Problem"+prob_num+"_Q"+ques_num+"_t=0.0.png")

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
        for k=1:N-1
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            l_k = deltaL;
            dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
            ind = [2*k-1, 2*k, 2*k+1, 2*k+2];
            f(ind) = f(ind) + dF;
            J(ind,ind) = J(ind,ind) +dJ;
        end
        % Bending spring
        for k=2:N-1
            xkm1 = q(2*k-3);
            ykm1 = q(2*k-2);
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            curvature0 = 0;
            l_k = deltaL;
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
            ind = [2*k-3, 2*k-2, 2*k-1, 2*k, 2*k+1, 2*k+2];
            f(ind) = f(ind) + dF;
            J(ind, ind) = J(ind, ind) +dJ;
        end
        
        % Viscous force
        f = f + C * (q-q0) / dt;
        J = J + C / dt;

        % Weight
        f = f - W;

        % Update
        q = q - J \ f;
        err = sum ( abs(f) );
    end
    % Update
    u = (q - q0) / dt;
    q0 = q;
    u0 = u;
    %
    % Q1
    % c is a index num, when c=1 -> q0, when c=2 -> t=0.01
    fprintf('Time: %f | Y-vel R1=%f, R2=%f, R3=%f\n',...
        (c-1) * dt, u(2), u(4), u(6) );
    f1 = figure(1);
    plot( q(1:2:end), q(2:2:end), 'ro-');
    axis equal
    drawnow
    if ismember(c, showing_steps)
        saveas(f1, "Figures/Problem"+prob_num+"_Q"+ques_num+"_t="+(c-1)*dt+".png")
    end
    %
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
saveas(f1, "Figures/Problem"+prob_num+"_t="+totalTime+".png");

terminal_vel = u;

end