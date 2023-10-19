clear all;
close all;

N = 21;
totalTime = 50;
prob_num = 2;
ques_num = 1;
dt = 0.01;
[u_i, r2pos_i, r2vel_i, theta_i]= ...
    simul_visc_spheres_multi(N, prob_num, ques_num, totalTime, dt);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 1 and 2

% Question 1

Nsteps = round( totalTime / dt );
timeArray = (1:Nsteps) * dt;

f2=figure(2);
plot(timeArray, r2pos_i, 'k-');
xlabel('Time, t [sec]');
ylabel('Position of the middle node, [meter]');
saveas(f2, "Figures/Problem2_R2pos.png");
f3=figure(3);
plot(timeArray, r2vel_i, 'r-');
xlabel('Time, t [sec]');
ylabel('Velocity of the middle node, [meter/sec]');
saveas(f3, "Figures/Problem2_R2vel.png");


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 3

% Question 3
% terminal vel vs N
%
Ns = [3, 5, 9, 19, 29, 39, 49,59,69,79,89,99];
terminal_velocities = [];
for iter=1:length(Ns)
    N = Ns(iter);
    [u, q,q,q]= ...
        simul_visc_spheres_multi(N, prob_num, ques_num, totalTime, dt);
    terminal_velocities = [terminal_velocities, u(end)];
end

f4=figure(4);
plot(Ns, terminal_velocities, 'k-');
xlabel('The number of nodes');
ylabel('Terminal_velocities [meter/sec]');
saveas(f4, "Figures/Problem2_terminalV_N.png")

% terminal vel vs dt
clear all;
close all;

totalTime = 50;
N = 21;
dts = [5e-4, 1e-3,5e-3, 1e-2,5e-2,1e-1,5e-1, 1e0];
terminal_velocities = [];
for iter=1:length(dts)
    dt = dts(iter);
    fprintf(mat2str(dt)+"\n");
    [u, q,q,q]= ...
        simul_visc_spheres_multi(N, 2, 3, totalTime, dt);
    terminal_velocities = [terminal_velocities, u(end)];
end

f4=figure(4);
plot(dts,terminal_velocities, 'k-');
xlabel('The time step size');
ylabel('Terminal_velocities [meter/sec]');
saveas(f4, "Figures/Problem2_terminalV_dt.png")

