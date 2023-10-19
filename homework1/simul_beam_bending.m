%% 
function [ymax, vmax] = ...
    simul_beam_bending(N, prob_num, P,R,Y,r,RodLength,totalTime, dt)
%% Given Parameters
ndof = N*2;
rho = 2700;
deltaL = RodLength / (N-1);

g = 9.8; % m/s^2

EI = Y * pi * (R^4-r^4) / 4;
EA = Y * pi * (R^2-r^2);

nodes = zeros(N, 2);
for c = 1:N
    nodes(c,1) = (c-1) * deltaL;
end

%% 
% Mass, M
M = zeros(ndof, ndof);
for k=1:N
    M(2*k-1, 2*k-1) = pi*(R^2-r^2)*RodLength/(N-1)*rho; % Mass for x_k
    M(2*k, 2*k) = M(2*k-1, 2*k-1); % Mass for y_k
end

% Weight vector, W
W = zeros(ndof, 1);
for k=1:N
    W(2*k-1) = 0; % weight along x is zero
    W(2*k) = -4/3*pi*R^3*rho*g;
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
showing_ts = [0, 0.0001, 0.001, 0.03,];
showing_steps = showing_ts / dt +1;

ymax = zeros(Nsteps, 1);
ymax(1) = 0;
vmax = zeros(Nsteps, 1);
vmax(1) = 0;

file_postfix = "_P="+mat2str(min(P));


f1 = figure(1);
plot( q(1:2:end), q(2:2:end), 'ro-');
axis equal
drawnow
saveas(f1, "Figures/Problem"+prob_num+file_postfix+"_t=0.0.png")

% Time marching scheme
for c=2:Nsteps

    % Fixed and free DOFs
    fixedDOF = [1;2; ndof];
    freeDOF = 3:ndof-1;
    bcv = [0;0;0]; 

    q = q0; % Guess
    q(fixedDOF) = bcv;

    % Newton Raphson
    err = 10 * tol;
    while err > tol
        q_free = q(freeDOF);

        f = M / dt * ( (q-q0)/dt - u0 );
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

        % Weight
        f = f - P;
        
        % At this point, we have f and J
        f_free = f(freeDOF);
        J_free = J(freeDOF, freeDOF);

        % Update
        dq_free = J_free \ f_free;
        q_free = q_free - dq_free;

        err = sum ( abs(f_free) );
        q(freeDOF) = q_free;
    end
    % Update
    u = (q - q0) / dt;
    q0 = q;
    u0 = u;
    %
    
    fprintf('Time: %f | Y-vel R1=%f, R2=%f, R3=%f\n',...
        (c-1) * dt, u(2), u(4), u(6) );
    f1 = figure(1);
    plot( q(1:2:end), q(2:2:end), 'ro-');
    axis equal
    drawnow
    % Q1
    % c is a index num, when c=1 -> q0, when c=2 -> t=0.01
    if ismember(c, showing_steps)
        saveas(f1, "Figures/Problem"+prob_num+file_postfix+"_t="+(c-1)*dt+".png")
    end
    
    % Store
    
    ymax(c) = min(q(2:2:end));
    vmax(c) = min(u(2:2:end));
    
end

f1 = figure(1);
plot( q(1:2:end), q(2:2:end), 'ro-');
axis equal
drawnow
saveas(f1, "Figures/Problem"+prob_num+file_postfix+"_t="+totalTime+".png")

end