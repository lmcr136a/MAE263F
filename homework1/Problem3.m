clear all;
close all;

nv = 50;
totalTime = 1;
prob_num = 3;
ques_num = 1;
dt = 0.01;

d = 0.75;
R = 0.013;
r = 0.011;
I = pi*(R^4 - r^4)/4;
RodLength = 1; 
Y = 70e9; % Using Y instead of E to avoid ambiguity
node_P_idx = round(d*(nv-1));

Pval = 2000;
P = zeros(2*nv, 1);
P(2*node_P_idx) = -Pval;

[ymax, vmax]= ...
    simul_beam_bending(nv, prob_num, P,R,Y,r,RodLength,...
    totalTime, dt);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 1 and 2

Nsteps = round( totalTime / dt );
timeArray = (1:Nsteps) * dt;
ymax_th = ones(Nsteps, 1);
c_dld = min([d, RodLength-d]);
ymax_th = ymax_th*(Pval*c_dld*((RodLength)^2 - c_dld^2)^1.5)/(9*sqrt(3)*Y*I*RodLength);


f2=figure(2);
plot(timeArray, ymax, 'k-');
hold on
plot(timeArray, -ymax_th, 'r-');
hold off
xlabel('Time, t [sec]');
ylabel('ymax, [meter]');
saveas(f2, "Figures/Problem3_Q1_ymax_experiment.png")

f3=figure(3);
plot(timeArray, ymax_th, 'k-');
xlabel('Time, t [sec]');
ylabel('ymax in theroy, [meter]');
saveas(f3, "Figures/Problem3_Q1_ymax_theory.png")

f4=figure(4);
plot(timeArray, vmax, 'r-');
xlabel('Time, t [sec]');
ylabel('vmax, [meter/sec]');
saveas(f4, "Figures/Problem3_Q1_vmax.png")

%% When P = 20000


Pval = 20000;
P = zeros(2*nv, 1);
P(2*node_P_idx) = -Pval;

[ymax, vmax]= ...
    simul_beam_bending(nv, prob_num, P,R,Y,r,RodLength,...
    totalTime, dt);

Nsteps = round( totalTime / dt );
timeArray = (1:Nsteps) * dt;

ymax_th = ones(Nsteps, 1);
c_dld = min([d, RodLength-d]);
ymax_th = ymax_th*(Pval*c_dld*((RodLength)^2 - c_dld^2)^1.5)/(9*sqrt(3)*Y*I*RodLength);


f2=figure(2);
plot(timeArray, ymax, 'k-');
hold on
plot(timeArray, -ymax_th, 'r-');
hold off
xlabel('Time, t [sec]');
ylabel('ymax, [meter]');
saveas(f2, "Figures/Problem3_Q1_ymax_experiment_P=20000.png")


f3=figure(3);
plot(timeArray, ymax_th, 'k-');
xlabel('Time, t [sec]');
ylabel('ymax in theroy, [meter]');
saveas(f3, "Figures/Problem3_Q1_ymax_theory_P=20000.png")
f4=figure(4);
plot(timeArray, vmax, 'r-');
xlabel('Time, t [sec]');
ylabel('vmax, [meter/sec]');
saveas(f4, "Figures/Problem3_Q1_vmax_P=20000.png")




%% Q2  P vs ymax

Ps = [1, 10, 100, 500, 1000, 2000, 5000, 1e4, 2e4, 5e4, 1e5, 5e5];
ymax_exs = [];
ymax_ths = [];

for iter=1:length(Ps)
    Pval = Ps(iter);
    P = zeros(2*nv, 1);
    P(2*node_P_idx) = -Pval;

    c_dld = min([d, RodLength-d]);
    ymax_th = (Pval*c_dld*((RodLength)^2 - c_dld^2)^1.5)/(9*sqrt(3)*Y*I*RodLength);

    [ymax, vmax]= ...
        simul_beam_bending(nv, prob_num, P,R,Y,r,RodLength,...
        totalTime, dt);
    ymax_exs = [ymax_exs, min(ymax)];
    ymax_ths = [ymax_ths, -ymax_th];
    fprintf("\n");
    fprintf(mat2str(Pval));
    fprintf("\n");
    fprintf(mat2str(min(ymax)));
    fprintf("\n");
    fprintf(mat2str(-ymax_th));
end

f5=figure(5);
semilogx(Ps, ymax_exs, 'k-');
hold on
semilogx(Ps, ymax_ths, 'r-');
hold off
ylim([-0.5, 0]);
xlabel('P, [N]');
ylabel('ymax, [meter]');
legend("Experiment", "Theory");
saveas(f5, "Figures/Problem3_Q2_P_vs_y.png")


