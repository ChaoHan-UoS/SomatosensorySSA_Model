%% Initializing model and protocal setup
% Replace the default value of a certain parameter with different values par 
% to test the sensitivity of SSA to this parameter

plot_TC = 0; % Plots the shape of a typical tuning curve

% Plotting specification
AXES_FONTSIZE = 10;
LineWidth = 1;
MarkerSize = 10;
COLORBAR_FONTSIZE = 10;


%% Network Parameters
M = 5; % Number of rows
N = 4; % Number of columns

%% ---- L4 ---- 
% Time constants (in seconds)
tau_h = 0.001; % Membrane time constant
tau_rec = 0.500; % Recovery time constant of L4 synapses 
tau_rec_s = 0.800; % Recovery time constant of sensory input (thalamocortical) synapses

% Utilization rates
U = 0.5; % Portion of available fraction of resources that is utilized in response to an action potential
U_s = 0.8; % Same as U, only for the thalamo-cotical synapses that convey the sensory input

% Connection strengths
J = 2.2; % Recurrent connections
J_1 = 0.05; % Lateral Connections from vertical and horizontal neighboring barrels 
J_2 = 0.001; % Lateral connections from diagonal neighboring barrels

% Gain function
Gain_thres = 5; % Threshold of gain 
Gain_slope = 1; % Slope of gain 
Gain_max = 1000;

%% ---- L6 ---- 
% Time constants (in seconds)
tau_h_L6 = 0.001;
tau_rec_L6 = 1;

% Utilization rate
U_L6 = 0.5; 

% Connection strengths
J_L6 = 2.5; 
J_1_L6 = 0.03; 
J_2_L6 = 0.001;

% Gain function:
Gain_thres_L6 = 3; 
Gain_slope_L6 = 1; 
Gain_max_L6 = 1000;

%% ----- Thalamus ----- :
Ntc = 100; % Number of TC cells in each barreloid
Nre = 100; % Number of RE cells in each barreloid

% Izhikevich model parameters
mu_a = 0; sigma_a = 0.00;
a = [0.005*ones(4*Ntc, 1)+mu_a+sigma_a*randn(4*Ntc, 1); 0.02*ones(4*Nre, 1)]; 
b = [0.26*ones(4*Ntc, 1); 0.2*ones(4*Nre, 1)];
p_th = 0.6; % Proportion of bursting TC cells to all TC cells (bursting + non-bursting)
b([Ntc*p_th+(1:Ntc*(1-p_th)) Ntc+Ntc*p_th+(1:Ntc*(1-p_th)) 2*Ntc+Ntc*p_th+(1:Ntc*(1-p_th)) 3*Ntc+Ntc*p_th+(1:Ntc*(1-p_th))]) = 0.25;
mu_c = 0; sigma_c = 0;
c = [-52*ones(4*Ntc, 1)+mu_c+sigma_c*randn(4*Ntc, 1); -55*ones(4*Nre, 1)];
d = [2*ones(4*Ntc, 1); 4*ones(4*Nre, 1)];

tau_ampa = 5; tau_gabaa = 6; tau_gabab = 150; % (ms)
v_ampa = 0; v_gabaa = -75; v_gabab = -90; % (mV)

g = zeros(4*Ntc+4*Nre); % Maximum thalamic synaptic conductance
g= sparse(g); % Sparse matrix
G = g; % Real-time thalamic synaptic conductance
I_in = g; % Thalamic synaptic current

p_tc2re = 0.6; p_re2tc = 0.6; p_re2re_diag = 0.6; p_re2re = 0.2; 
g_tc2re_ampa = zeros(Nre, Ntc);
g_re2tc_gabaa = zeros(Ntc, Nre); 
g_re2re_gabaa = 0.5/(p_re2re*Nre)*(rand(Nre) < p_re2re); 

g_tc2re_ampa(1:Nre*p_th, 1:Ntc*p_th)= 2/(p_tc2re*Ntc*p_th)*(rand(Nre*p_th, Ntc*p_th) < p_tc2re); 
g_tc2re_ampa(Nre*p_th+1:Nre, Ntc*p_th+1:Ntc)= 2/(p_tc2re*Ntc*(1-p_th))*(rand(Nre*(1-p_th), Ntc*(1-p_th)) < p_tc2re); 
g_re2tc_gabaa(1:Ntc*p_th, 1:Nre*p_th) = 0.01/(p_re2tc*Nre*p_th)*(rand(Ntc*p_th, Nre*p_th) < p_re2tc); 
g_re2tc_gabaa(Ntc*p_th+1:Ntc, Nre*p_th+1:Nre) = 0.01/(p_re2tc*Nre*(1-p_th))*(rand(Ntc*(1-p_th), Nre*(1-p_th)) < p_re2tc); 
g_re2re_gabaa(1:Nre*p_th, 1:Nre*p_th) = 0.5/(p_re2re_diag*Nre*p_th)*(rand(Nre*p_th) < p_re2re_diag); 
g_re2re_gabaa(Nre*p_th+1:Nre, Nre*p_th+1:Nre) = 0.5/(p_re2re_diag*Nre*(1-p_th))*(rand(Nre*(1-p_th)) < p_re2re_diag); 

for ii = 0:3
    g(4*Ntc+(1:Nre)+ii*Nre, (1:Ntc)+ii*Ntc) = g_tc2re_ampa;
    g((1:Ntc)+ii*Ntc, 4*Ntc+(1:Nre)+ii*Nre) = g_re2tc_gabaa;
    g(4*Ntc+(1:Nre)+ii*Nre, 4*Ntc+(1:Nre)+ii*Nre) = g_re2re_gabaa;
end

v = [-62.5*ones(4*Ntc, 1); -70*ones(4*Nre, 1)]; % Initial values of v
u = b.*v; % Initial values of u
vv = zeros(4*Ntc+4*Nre, 0); % Temporal collections of v
firings = []; % spike timings
del_t = 2; % time bin (ms)
n_act = zeros(4,1); % Spike count of all TC cells in each time bin 
n_act1 = zeros(4,1); % Spike count of bursting TC cells
n_act2 = zeros(4,1); % Spike count of non-bursting TC cells
A_tc = zeros(4,1); % Population activity in each time bin
A_tc1 = zeros(4,1); 
A_tc2 = zeros(4,1); 

%% -----  Co-th-co Loop ----- 
% Thalamus to L4
J_thL41 = 1; % Connection strength from bursting TCs to L4
J_thL42 = 0.05; % Connection strength from non-bursting TCs to L4
thres_thL4 = 0; 

% L4 to L6
J_L46 = 0.24; % Connection strength from L4 to L6
tau_rec_L46 = 1.2;
U_L46 = 0.5;

% Bottom-up stimulus inputs to thalamus
k_sth1 = 5; % Proportion of bursting TCs (p_th*Ntc) that receive sensory inputs
k_sth2 = 20; % Proportion of non-bursting TCs ((1-p_th)*Ntc) that receive sensory inputs

rng(666);
Conn_sth = zeros(4*Ntc+4*Nre, 1); 
Conn_sth([randi(Ntc*p_th,1,k_sth1) Ntc+randi(Ntc*p_th,1,k_sth1) ...
    2*Ntc+randi(Ntc*p_th,1,k_sth1) 3*Ntc+randi(Ntc*p_th,1,k_sth1)]) = 1;
Conn_sth([randi([Ntc*p_th+1 Ntc],1,k_sth2) Ntc+randi([Ntc*p_th+1 Ntc],1,k_sth2) ...
    2*Ntc+randi([Ntc*p_th+1 Ntc],1,k_sth2) 3*Ntc+randi([Ntc*p_th+1 Ntc],1,k_sth2)]) = 1; 

% Top-down L6 inputs to thalamus
w_L6tc = 0.001; % Connection strength from L6 to TCs
w_L6re = 0.4; % Connection strength from L6 to REs
thres_L62th = 0; 
k_L6th = 50; % Proportion of TC/RE cells (Ntc/Nre) that receive L6 feedback

Conn_L6th = zeros(4*Ntc+4*Nre, 1); 
% Subpopulation of bursting TCs 
Conn_L6th([randi(Ntc*p_th,1,k_L6th*p_th) Ntc+randi(Ntc*p_th,1,k_L6th*p_th) ...
    2*Ntc+randi(Ntc*p_th,1,k_L6th*p_th) 3*Ntc+randi(Ntc*p_th,1,k_L6th*p_th)]) = w_L6tc; 
% Subpopulation of REs coupled with bursting TCs 
Conn_L6th(4*Ntc+[randi(Nre*p_th,1,k_L6th*p_th) Nre+randi(Nre*p_th,1,k_L6th*p_th) ...
    2*Nre+randi(Nre*p_th,1,k_L6th*p_th) 3*Nre+randi(Nre*p_th,1,k_L6th*p_th)]) = w_L6re;
% Subpopulation of non-bursting TCs 
Conn_L6th([randi([Ntc*p_th+1 Ntc],1,k_L6th*(1-p_th)) Ntc+randi([Ntc*p_th+1 Ntc],1,k_L6th*(1-p_th)) ...
    2*Ntc+randi([Ntc*p_th+1 Ntc],1,k_L6th*(1-p_th)) 3*Ntc+randi([Ntc*p_th+1 Ntc],1,k_L6th*(1-p_th))]) = w_L6tc; 
% Subpopulation of REs coupled with non-bursting TCs 
Conn_L6th(4*Ntc+[randi([Nre*p_th+1 Nre],1,k_L6th*(1-p_th)) Nre+randi([Nre*p_th+1 Nre],1,k_L6th*(1-p_th)) ...
    2*Nre+randi([Nre*p_th+1 Nre],1,k_L6th*(1-p_th)) 3*Nre+randi([Nre*p_th+1 Nre],1,k_L6th*(1-p_th))]) = w_L6re;

%% Stimulus Parameters:
dt = 0.0001; % Time-step (in seconds)
dtt = dt*10^3; % Time-step (in miliseconds) for thalamic spiking network
duration = 0.010; % Total duration of each stimulus (in seconds) 
ramp_dur = 0.002; % durations of the ramps at the beginning and end of each stimulus (in seconds) 

t_eq = 1; % Equilibrium time given before the first stimuli (in seconds)
num_steps_eq = floor(t_eq/dt); % Total number of steps to reach the equilibrium
post_stim = 0.500; % Simulation keeps runnng for 2*post_stim after the last trial of stimulus

%% Turning Curve of Sensory Input:
lambda = 1.6; % lambda < 2 
r = [0 1 sqrt(2) 2]; % Distance between two barrels; 0 1 sqrt(2) 2 sqrt(5) sqrt(8) 3 ...
curve_type = 'linear'; % Type of the tuning curve

switch curve_type
    case 'exponential' 
        s_f = exp(-r/lambda);
    case 'linear' % Triangular TC
        s_f = -(1/lambda)*(r-lambda);
        s_f(s_f < 0) = 0;
    case 'Gaussian' 
        p1 = -0.5*(r/lambda).^2;
        p2 = sqrt(2*pi)*lambda;
        s_f = exp(p1)/p2;
end

% Plotting the Tuning Curve:
if plot_TC 
    f = figure;
    subplot(2,2,1)
    plot([0 lambda],[1 0],'-b',[lambda 2],[0 0],'-b',r,s_f,'bo','LineWidth',LineWidth,'MarkerSize',MarkerSize);
    set(gca,'FontSize', AXES_FONTSIZE);
    title('Tuning Curve');
    xlabel('d'); 
    ylabel('T_{p,q}^{m,n}');
    xlim([0 2]);
    saveas(f,'Figure/tuning_curve.pdf');
end


%% Saving:
% E_act = []; % just to save memory
% x_act = [];
% z_act = [];
% E_act_L6 = [];
% x_act_L6 = [];
save(['Simulation Results/run_par' num2str(par) '_initialization.mat']);
disp('Initialization done');
