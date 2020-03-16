%%  Model Initialization

plot_TC = 1; % Plots the shape of a typical tuning curve
plot_init = 1; % Plots the equilibria of E, x and z

Cond{1} = 'Low';
Cond{2} = 'High';
Cond{3} = 'Equal';
Cond{4} = 'Many Standard';
Cond{5} = 'Deviant Alone 1';
Cond{6} = 'Deviant Alone 2';

Cond_Code{1} = 'L';
Cond_Code{2} = 'H';
Cond_Code{3} = 'E';
Cond_Code{4} = 'MS';
Cond_Code{5} = 'DA1';
Cond_Code{6} = 'DA2';

PW= [3,2];
AW1 = [4,2];
AW2 = [4,1];
AW3 = [4,3];

%% Parameters of the Model:
%% ---- L4 ---- :
M = 5; % Number of rows
N = 4; % Number of columns

% Time constants:
tau_h = 0.001;
tau_rec = 0.200; 
% tau_rec = par; % recovery time constant of intracortical synapses (in seconds) %0.700
tau_rec_s = 0.800; % recovery time constant of sensory input (thalamocortical) synapses (in seconds)

U = 0.5; % Portion of available fraction of resources that is utilized in response to an action potential 0.5
% U = par;
U_s = 0.7; % Same as U, only for the thalamo-cotical synapses that convey the sensory input

% Connection strengths:
J = 2.5; % self-connection %2.5
J_1 = 0.1045 % Connections from vertical and horizontal neighboring barrels %0.1 0.092
% J_1 = par;
J_2 = 0.0; % Connections from diagonal neighboring barrels %0 

%Gain function:
Gain_thres = 3; % Threshold of gain 
% Gain_thres = par;
Gain_slope = 1; % Slope of gain 
Gain_max = 1000;

%% ---- L6 ---- :
% Time constants:
tau_h_L6 = 0.001;
tau_rec_L6 = 0.700; 
U_L6 = 0.55; % Portion of available fraction of resources that is utilized in response to an action potential 0.5

% Connection strengths:
J_L6 = 2.3; % self-connection %2.1 2.4
J_1_L6 = 0.02 % Connections from vertical and horizontal neighboring barrels %0.01 0.075
J_2_L6 = 0.00; % Connections from diagonal neighboring barrels %0.003 

%Gain function:
Gain_thres_L6 = 2; % Threshold of gain 
Gain_slope_L6 = 1; % Slope of gain 
Gain_max_L6 = 1000;

%% ----- Thalamus ----- :
% Time constants:
tau_h_Th = 0.001; 

%Gain function:
Gain_thres_Th = 3; % Threshold of gain 
Gain_slope_Th = 1; % Slope of gain 
Gain_max_Th = 1000;

%% -----  Cortico-thalamo-cortical Loop ----- :
w_L426 = 0.06; % L4 to L6 connection weight 0.06
w_L6Th = 0.15; % L6 to thalamus weight
tau_ThL4 = 0.100; % Latency of thalamus to L4 transmission

%% Stimulus Parameters
t_eq = 20; % The time given to reach equilibrium (in seconds). It is important to allow enough time, otherwise the response to the first stimulus will be distorted. %5
dt = 0.0001; % Time-step (in seconds)
num_steps_eq = floor(t_eq/dt); % Total number of steps to reach the equilibrium
post_stim = 0.100; % This is defined here so the simulations keeps runnng for 2*post_stim after the last stimulus offset

duration = 0.010; % Total duration of each stimulus (offset to onset) (in seconds) 0.050 0.010
ramp_dur = 0.001; % durations of the ramps at the beginning and end of each stimulus (in seconds) 0.005 0.002

%% Turning Curve of Sensory Input
lambda = 1.2 % lambda < 2 %1.2 1.16
% lambda = par;
r = [0 1 sqrt(2) 2]; % Distance between two barrels; 0 1 sqrt(2) 2 sqrt(5) sqrt(8) 3 ...
curve_type = 'linear'; % Type of the tuning curve

switch curve_type
    case 'exponential' 
        s_f = exp(-r/lambda);
    case 'linear' % This type is triangular
        s_f = -(1/lambda)*(r-lambda);
        s_f(s_f < 0) = 0;
    case 'Gaussian' 
        p1 = -0.5*(r/lambda).^2;
        p2 = sqrt(2*pi)*lambda;
        s_f = exp(p1)/p2;
end

%% Plotting the Tuning Curve 
if plot_TC 
    f = figure('Visible','off');
    plot([0 lambda], [1 0], '-b');
    hold on;
    plot([lambda 2], [0 0], '-b');
    plot(r, s_f, 'bo');
    title('Tuning Curve');
    xlabel('d');
    ylabel('T_{p,q}^{m,n}');
    xlim([0 2]);
    saveas(f,'Figure/tuning_curve.pdf');
end

%% Initializing Variables:
%% ---- L4 ---- :
E = zeros(M,N); % Acitivity(firing rate) of all barrels (in Hz) 
h = E; % Sypnatic inputs of all barrels
x = zeros(M,N); % Mean fractions of resources available for synaptic transmission in all barrels
z = 1; % Only one synapse is needed here; zeros(P,NE,M_aug); % Same as x, only for the thalamo-cortical synapses

H = zeros(M,N); % Sum of synaptic inputs

E_act = zeros(M,N,floor(t_eq/dt)); % Activity of all barrelds during the time allowed for reaching equilibrium
x_act = zeros(M,N,floor(t_eq/dt));
z_act = zeros(1,floor(t_eq/dt));

EUx = zeros(M,N);
EUx_1 = EUx;
EUx_2 = EUx;

%% ---- L6 ---- :
E_L6= zeros(M,N); % Acitivity(firing rate) of all barrels (in Hz) 
h_L6 = E_L6; % Sypnatic inputs of all barrels
x_L6 = zeros(M,N); % Mean fractions of resources available for synaptic transmission in all barrels
H_L6 = zeros(M,N); % Sum of synaptic inputs

E_act_L6 = zeros(M,N,floor(t_eq/dt)); % Activity of all barrelds during the time allowed for reaching equilibrium
x_act_L6 = zeros(M,N,floor(t_eq/dt));

EUx_L6 = zeros(M,N);
EUx_1_L6 = EUx_L6;
EUx_2_L6 = EUx_L6;

%% ----- Thalamus ----- :
E_Th= zeros(M,N); % Acitivity(firing rate) of all barreloids (in Hz) 
h_Th = E_Th; % Sypnatic inputs of all barrels
E_act_Th = zeros(M,N,floor(t_eq/dt)); % Activity of all barrelds during the time allowed for reaching equilibrium

%% The Dynamic Loop: Initial run to obtain steady-state and find initial values
% The network is allowed to reach a steady-state with no sensory input. 
for i = 1:floor(t_eq/dt)
    %% ---- L4 ---- :
    % Total synaptic inputs each L4 barrel receives from itself, 4 vertical and horizontal neighbors and 4 diagonal neighbors
    EUx = E*U.*x;
    Temp1 = zeros(M, N);
    Temp1(:, 2:N) = EUx(:, 1:N-1); % Input from left barrel
    Temp2 = zeros(M, N);
    Temp2(:, 1:N-1) = EUx(:, 2:N); % Input from right barrel
    Temp3 = zeros(M, N);
    Temp3(2:M, :) = EUx(1:M-1, :); % Input from top barrel
    Temp4 = zeros(M, N);
    Temp4(1:M-1, :) = EUx(2:M, :); % Input from bottom barrel
    EUx_1 = Temp1 + Temp2 + Temp3 + Temp4; % Input from vertical and horizontal neighboring barrels
    Temp1 = zeros(M, N);
    Temp1(2:M, 2:N) = EUx(1:M-1, 1:N-1); % Input from left-top barrel
    Temp2 = zeros(M, N);
    Temp2(2:M, 1:N-1) = EUx(1:M-1, 2:N); % Input from right-top barrel
    Temp3 = zeros(M, N);
    Temp3(1:M-1, 2:N) = EUx(2:M, 1:N-1); % Input from left-bottom barrel
    Temp4 = zeros(M, N);
    Temp4(1:M-1, 1:N-1) = EUx(2:M, 2:N); % Input from right-bottom barrel
    EUx_2 = Temp1 + Temp2 + Temp3 + Temp4; % Input from diagonal neighboring barrels
    H = J*EUx + J_1*EUx_1 + J_2*EUx_2; % The total synaptic input each L4 barrel receives
    
    % Adding thalamus-to-L4 latency:
    if i-floor(tau_ThL4/dt) > 0
        H = H + E_act_L6(:, :, i-floor(tau_ThL4/dt))*U_s*z; % One-to-one
    end
    
    % The dynamics of synaptic input:
    h = h + (dt/tau_h)*(-h + H);
    
    % Implementing the non-linearity:
    E = h;  
    E = E - Gain_thres; 
    E = E * Gain_slope;
    E(E <  0) = 0;   
    E(E >  Gain_max) = Gain_max;   

   % The variables' dynamics:
    x = x + dt*((1 - x)./tau_rec - U.*E.*x);         
    z = z + dt*((1 - z)./tau_rec_s);

    % Tracking the activities of all barrels
    E_act(:,:,i) = E; 
    x_act(:,:,i) = x; 
    z_act(i) = z;
    
    %% ---- L6 ---- :
    % Total synaptic inputs each L6 infrabarrel receives from itself, 4 vertical and horizontal neighbors, 4 diagonal neighbors as well as intra-column L4 barrel
    EUx_L6 = E_L6*U_L6.*x_L6;
    Temp1 = zeros(M, N);
    Temp1(:, 2:N) = EUx_L6(:, 1:N-1); % Input from left barrel
    Temp2 = zeros(M, N);
    Temp2(:, 1:N-1) = EUx_L6(:, 2:N); % Input from right barrel
    Temp3 = zeros(M, N);
    Temp3(2:M, :) = EUx_L6(1:M-1, :); % Input from top barrel
    Temp4 = zeros(M, N);
    Temp4(1:M-1, :) = EUx_L6(2:M, :); % Input from bottom barrel
    EUx_1_L6 = Temp1 + Temp2 + Temp3 + Temp4; % Input from vertical and horizontal neighboring barrels
    Temp1 = zeros(M, N);
    Temp1(2:M, 2:N) = EUx_L6(1:M-1, 1:N-1); % Input from left-top barrel
    Temp2 = zeros(M, N);
    Temp2(2:M, 1:N-1) = EUx_L6(1:M-1, 2:N); % Input from right-top barrel
    Temp3 = zeros(M, N);
    Temp3(1:M-1, 2:N) = EUx_L6(2:M, 1:N-1); % Input from left-bottom barrel
    Temp4 = zeros(M, N);
    Temp4(1:M-1, 1:N-1) = EUx_L6(2:M, 2:N); % Input from right-bottom barrel
    EUx_2_L6 = Temp1 + Temp2 + Temp3 + Temp4; % Input from diagonal neighboring barrels
    H_L6 = J_L6*EUx_L6 + J_1_L6*EUx_1_L6 + J_2_L6*EUx_2_L6 + w_L426*E; % The total synaptic input each L6 infrabarrel receives
    
    % The dynamics of synaptic input:
    h_L6 = h_L6 + (dt/tau_h_L6)*(-h_L6 + H_L6);
    
    % Implementing the non-linearity:
    E_L6 = h_L6;  
    E_L6 = E_L6 - Gain_thres_L6; 
    E_L6 = E_L6 * Gain_slope_L6;
    E_L6(E_L6 <  0) = 0;   
    E_L6(E_L6 >  Gain_max_L6) = Gain_max_L6;   

    % The variables' dynamics:
    x_L6 = x_L6 + dt*((1 - x_L6)./tau_rec_L6 - U_L6.*E_L6.*x_L6);         

    % Tracking the activities of all barrels
    E_act_L6(:,:,i) = E_L6; 
    x_act_L6(:,:,i) = x_L6;
   
    %% ----- Thalamus ----- :
    % The dynamics of synaptic input:
    h_Th = h_Th + (dt/tau_h_Th)*(-h_Th + w_L6Th*E_L6);
    
    % Implementing the non-linearity:
    E_Th = h_Th;  
    E_Th = E_Th - Gain_thres_Th; 
    E_Th = E_Th * Gain_slope_Th;
    E_Th(E_Th <  0) = 0;   
    E_Th(E_Th >  Gain_max_Th) = Gain_max_Th;   

    % Tracking the activities of all barrels
    E_act_Th(:,:,i) = E_Th; 
    
end

%% Recording Equilibrium Conditions:
%L4 
E_eq = E; 
h_eq = h;
x_eq = x; 
z_eq = z; 

%L6
E_eq_L6 = E_L6; 
h_eq_L6 = h_L6;
x_eq_L6 = x_L6;  

%Thalamus
E_eq_Th = E_Th; 
h_eq_Th = h_Th;

%% Plot Equilibrium:
if plot_init 
    figure
%     f = figure('Visible','off');
    plot(dt:dt:t_eq, reshape(E_act(ceil(M/2),ceil(N/2),:), [1, floor(t_eq/dt)]), 'b'); % Firing rate of L4 and L6 of barrel (3,2) as well as its corresponding barreloid
    hold on
    plot(dt:dt:t_eq, reshape(E_act_L6(ceil(M/2),ceil(N/2),:), [1, floor(t_eq/dt)]), 'r');
    plot(dt:dt:t_eq, reshape(E_act_Th(ceil(M/2),ceil(N/2),:), [1, floor(t_eq/dt)]), 'g');
    title('Mean Firing Rate of barrel (3,2) during Initial Run');
    xlabel('time(s)');
    ylabel('Spikes/s');
    legend('L4','L6','Th');
%     saveas(f,'Figure/init_firing.pdf');

    figure
%     f=figure('Visible','off');
    plot(dt:dt:t_eq, reshape(x_act(ceil(M/2),ceil(N/2),:), [1,floor(t_eq/dt)]), 'b'); % Resources of L4 and L6 of barrel (3,2) 
    hold on
    plot(dt:dt:t_eq, reshape(x_act_L6(ceil(M/2),ceil(N/2),:), [1, floor(t_eq/dt)]), 'r');
    title('Mean Resources of barrel (3,2) during Initial Run');
    xlabel('time(s)');
    legend('L4','L6');
%     saveas(f,'Figure/init_resource.pdf');

    figure
%     f = figure('Visible','off');
    plot(dt:dt:t_eq, z_act(:)); % Thalamocortical resource, 1 at equilibrium
    title('ThC Resources during Initial Run');
    xlabel('time(s)');
%     saveas(f,'Figure/init_ThCresource.pdf');
    
end

%% Saving
% E_act = [];% just to save memory
% x_act = [];
% z_act = [];
% E_act_L6 = [];
% x_act_L6 = [];
save(['Simulation Results/run_par' num2str(par) '_initialization.mat']);
disp('Initialization done');
