
clear all;
close all;
%% Implicit Simulation
N = 3;
totalTime = 10;
prob_num = 1;
ques_num = 1;
Rs = [0.005, 0.025, 0.005];
dt = 0.01;

[u_i, r2pos_i, r2vel_i, theta_i]= ...
    simul_visc_spheres(N, prob_num, ques_num, Rs, totalTime, dt);

%% Explicit Simulation

totalTime_e = 10;
dt_e = 0.00001;

[u_e, r2pos_e, r2vel_e, theta_e]= ...
    simul_visc_spheres_explicit(N, prob_num, ques_num, Rs, totalTime_e, dt_e);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 1 and 2

% Question 1

Nsteps = round( totalTime / dt );
timeArray = (1:Nsteps) * dt;

f2=figure(2);
plot(timeArray, r2pos_i, 'k-');
xlabel('Time, t [sec]');
ylabel('Position of R2, [meter]');
saveas(f2, "Figures/Problem1_Q1_R2_position_Implicit.png")
f3=figure(3);
plot(timeArray, r2vel_i, 'r-');
xlabel('Time, t [sec]');
ylabel('Velocity of R2, [meter/sec]');
saveas(f3, "Figures/Problem1_Q1_R2_velocity_Implicit.png")

Nsteps = round( totalTime_e / dt_e );
timeArray = (1:Nsteps) * dt_e;

f4=figure(4);
plot(timeArray, r2pos_e, 'k-');
xlabel('Time, t [sec]');
ylabel('Position of R2, [meter]');
saveas(f4, "Figures/Problem1_Q1_R2_position_Explicit.png")
f5=figure(5);
plot(timeArray, r2vel_e, 'r-');
xlabel('Time, t [sec]');
ylabel('Velocity of R2, [meter/sec]');
saveas(f5, "Figures/Problem1_Q1_R2_velocity_Explicit.png")

% Question 2
fprintf("Q2) Implicit Sim Terminal Velocity => | R1: x=%f, y=%f | R2: x=%f, y=%f |R3: x=%f, y=%f |", ...
    u_i(1),u_i(2),u_i(3),u_i(4),u_i(5),u_i(6))
%
writematrix(mat2str(u_i), "results/P1_Q1_TerminalVelocity_"+dt_e+".txt");
fprintf("Q2) Explicit Sim Terminal Velocity => | R1: x=%f, y=%f | R2: x=%f, y=%f |R3: x=%f, y=%f |", ...
    u_e(1),u_e(2),u_e(3),u_e(4),u_e(5),u_e(6));
writematrix(mat2str(u_e), "results/P1_Q1_TerminalVelocity_"+dt_e+"_explicit.txt");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 3
ques_num = 3;
dt = 0.01;

% R1=R2=R3=0.005 
Rs = [0.005, 0.005, 0.005];
[q, q, q, theta_005]= ...
    simul_visc_spheres(N, prob_num, ques_num, Rs, totalTime, dt);

% R1=R2=R3=0.025 
Rs = [0.025, 0.025, 0.025];
[q, q, q, theta_025]= ...
    simul_visc_spheres(N, prob_num, ques_num, Rs, totalTime, dt);
%
%
Nsteps = round( totalTime / dt );
timeArray = (1:Nsteps) * dt;
Nsteps_e = round( totalTime_e / dt_e );
timeArray_e = (1:Nsteps_e) * dt_e;

f6=figure(6);
plot(timeArray, theta_i, 'r-');
hold on
%plot(timeArray_e, theta_e, 'b-');
plot(timeArray, theta_005, 'g-');
plot(timeArray, theta_025, 'b-');
hold off
xlabel('Time, t [sec]');
ylabel('Turning Angles, [°]');
legend('R1&R3=0.005, R2=0.025',  ...
    'All 0.005', 'All 0.025' )
saveas(f6, "Figures/Problem1_Q3_turning_angle.png")

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 4
% Reset parameters
N = 3;
totalTime_e = 10;
prob_num = 1;
ques_num = 4;
Rs = [0.005, 0.025, 0.005];

dts = [0.0005, 0.01];

for i = 1:length(dts)
    dt_e = dts(i);
    [u_e, r2pos_e, r2vel_e, theta_e]= ...
        simul_visc_spheres_explicit(N, prob_num, ques_num, Rs, totalTime, dt_e);

    Nsteps = round( totalTime_e / dt_e );
    timeArray = (1:Nsteps) * dt_e;
    
    f4=figure(4);
    plot(timeArray, r2pos_e, 'k-');
    xlabel('Time, t [sec]');
    ylabel('Position of R2, [meter]');
    saveas(f4, "Figures/Problem1_Q4_R2_position_Explicit_"+mat2str(dt_e)+".png")
    f5=figure(5);
    plot(timeArray, r2vel_e, 'r-');
    xlabel('Time, t [sec]');
    ylabel('Velocity of R2, [meter/sec]');
    saveas(f5, "Figures/Problem1_Q4_R2_velocity_Explicit_"+mat2str(dt_e)+".png")
    
    fprintf("Q2) Explicit Sim Terminal Velocity => | R1: x=%f, y=%f | R2: x=%f, y=%f |R3: x=%f, y=%f |", ...
        u_e(1),u_e(2),u_e(3),u_e(4),u_e(5),u_e(6));
    writematrix(mat2str(u_e), "results/P1_Q4_TerminalVelocity_"+mat2str(dt_e)+"_explicit.txt");
end
