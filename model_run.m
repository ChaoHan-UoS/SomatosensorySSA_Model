%% Model Run
% This code runs the network saved and brought to equilibrium by
% model_init, according to the specified condition defined by cond

% cond = 1;
% A = 1;
% ISI = 1;
% stim = 1;
% prob = 1;

tic
n_stim = 120; % Total no. of stimuli (Best take a product of 4) %120
t_prot = n_stim*(duration + ISI) + 2*post_stim; % Total time of the protocol (in seconds) 
tmax_tot = t_eq + t_prot; % Maximum time that the simulation will reach (in seconds) 
num_steps = floor(tmax_tot/dt); % Total number of steps in the simulation

Prob.L = prob;
Barr_Stim.L = [AW1; PW]; % Location of stimulated barrels
    
Prob.H = 1 - prob;
Barr_Stim.H = [PW; AW1]; 

Prob.E = 0.5;
Barr_Stim.E = [AW1; PW]; 
        
Prob.MS = prob;
Barr_Stim.MS = [AW1; PW; AW2; AW3];

Prob.DA1 = prob;
Barr_Stim.DA1 = AW1; 
        
Prob.DA2 = prob;
Barr_Stim.DA2 = PW; 

switch cond
    case('Low')
        Probs = Prob.L;
        Barr_Stims = Barr_Stim.L;
    case('High')
        Probs = Prob.H;
        Barr_Stims = Barr_Stim.H;
    case('Equal')
        Probs = Prob.E;
        Barr_Stims = Barr_Stim.E;
    case('Many Standard')
        Probs = Prob.MS;
        Barr_Stims = Barr_Stim.MS;
    case('Deviant Alone 1')
        Probs = Prob.DA1;
        Barr_Stims = Barr_Stim.DA1;
    case('Deviant Alone 2')
        Probs = Prob.DA2;
        Barr_Stims = Barr_Stim.DA2;    
end

%% Oddball Sequence Generation
switch cond
    case('Low') % AW1 as standard and PW as deviant 
        rng(666);
        n_stim_temp = 120;
        Oddball = [ones(ceil((1-Probs)*n_stim_temp),1)*Barr_Stims(1,:); ones(ceil(Probs*n_stim_temp),1)*Barr_Stims(2,:)]; % [AW1 Pw]
        Oddball = Oddball(randperm(n_stim_temp),:); % n_stim*2
        Oddball = Oddball(1:n_stim,:);
    case('High') % AW1 as deviant and PW as standard  
        Oddball = [ones(ceil(Probs*n_stim),1)*Barr_Stims(1,:); ones(ceil((1-Probs)*n_stim),1)*Barr_Stims(2,:)];
        Oddball = Oddball(randperm(n_stim),:); 
    case('Equal')
        Oddball = [ones(ceil(Probs*n_stim),1)*Barr_Stims(1,:); ones(ceil((1-Probs)*n_stim),1)*Barr_Stims(2,:)]; % Here Probs=0.5
        Oddball = Oddball(randperm(n_stim),:); 
    case('Many Standard')  % AW1 PW AW2 AW3, each stimulated with probability of 0.25
        rng(666);
        n_stim_temp = 120;
        whis_reps = n_stim_temp/size(Barr_Stims,1); % Repetitons of each stimulated whisker 
        Oddball = [ones(whis_reps,1)*Barr_Stims(1,:); ones(whis_reps,1)*Barr_Stims(3,:);...
                   ones(whis_reps,1)*Barr_Stims(4,:); ones(whis_reps,1)*Barr_Stims(2,:)];
        Oddball = Oddball(randperm(n_stim_temp),:);
        Oddball = Oddball(1:n_stim,:);
    case('Deviant Alone 1') % AW1 as deviant 
        Oddball = [Barr_Stims*ones(1,ceil(Probs*n_stim)) zeros(1,ceil((1-Probs)*n_stim))];
        Oddball = Oddball(randperm(n_stim));
    case('Deviant Alone 2') % PW as deviant 
        Oddball = [Barr_Stims*ones(1,ceil(Probs*n_stim)) zeros(1,ceil((1-Probs)*n_stim))];
        Oddball = Oddball(randperm(n_stim));
end

Spa_Temp = zeros(1,size(Barr_Stims,1),num_steps); % Spatial-temporal structure of the oddball stimulus sequence 
ramp_step = dt/ramp_dur;
Single_Stim = [ramp_step:ramp_step:1, ones(1,floor((duration - 2*ramp_dur)/dt)), wrev(ramp_step:ramp_step:1)]; % Waveform of each single stimulus
Stim_Onsets = zeros(1,n_stim); % onset time indices of each stimulus 
Time_Ind = zeros(n_stim,2); % onset & offset time indices of each stimulus

for ns = 1:n_stim    
    Time_Ind(ns,:) = floor(t_eq/dt) + floor((ISI+duration)/dt)*(ns - 1) + [0, (length(Single_Stim)-1)]; % Indices in which there is a stimulus; 
    if Oddball(ns,1) % For deviant alone, index of barrel is zero
        Spa_Temp(1, sum(Barr_Stims == Oddball(ns,:), 2) == 2, Time_Ind(ns,1):Time_Ind(ns,2)) = Single_Stim; % Amplitude of each stimulus along time indices
    end
    Stim_Onsets(ns) = Time_Ind(ns,1);
end

%% Initializing variables and sensory inputs:
% ---- L4 ---- :
E = zeros(M,N); % Acitivity(firing rate) of all barrels (in Hz) 
h = zeros(M,N); % Sypnatic inputs of all barrels
x = ones(M,N); % Mean fractions of resources available for synaptic transmission in all barrels
z1 = ones(M,N,size(Barr_Stims,1)); 
z2 = ones(M,N,size(Barr_Stims,1)); 
E_act_overall = zeros(M,N,num_steps); % Tracking activity of all barrels 
x_act_overall = zeros(M,N,num_steps); 
z_act_overall = zeros(M,N,num_steps); % AW1
z_act_overall2 = zeros(M,N,num_steps); % PW
s_E1 = zeros(M,N,size(Barr_Stims,1)); % Sensory inputs to L4 (first half TC cells)
s_E2 = zeros(M,N,size(Barr_Stims,1)); % (second half TC cells)
% s_E_plot = zeros(M,N,num_steps);

% ---- L6 ---- :
E_L6= zeros(M,N); 
h_L6 = zeros(M,N);
x_L6 = ones(M,N); 
E_act_overall_L6 = zeros(M,N,num_steps); 
x_act_overall_L6 = zeros(M,N,num_steps); 
x_L46 = ones(M,N); 
E_act_overall_L46 = zeros(M,N,num_steps); % L4 input to L6 
x_act_overall_L46 = zeros(M,N,num_steps); 

% ----- Thalamus ----- :
I_s = zeros(4*Ntc+4*Nre,1);% Sensory input currents to TC neurons
% I_s_overall = zeros(4*Ntc+4*Nre,num_steps);
% I_L6th_overall = zeros(4*Ntc+4*Nre,num_steps);
% I_in_overall = zeros(4*Ntc+4*Nre,num_steps);
E_act_overall_tc = zeros(4,num_steps); 

%% Dynamic Loop
for i = 1:num_steps
    %% ----- Thalamus ----- :    
    fired = find(v >= 30); % indices of spikes
    firings = [firings; i*dt+0*fired,fired];
    
    % Firing rate of TC cells in corresponding stimulated barreloid
    if mod(i*dtt, del_t) < 10^-8 % Floating-point number tolerance
        A_tc = n_act/(Ntc*del_t)*1000; % Convert unit of A_tc from KHz to Hz; 4*1 matrix
        A_tc1 = n_act1/(Ntc*del_t)*1000; 
        A_tc2 = n_act2/(Ntc*del_t)*1000; 
        
        E_act_overall_tc(1, i : i+del_t/dtt-1) = A_tc(1); % Hz
        E_act_overall_tc(2, i : i+del_t/dtt-1) = A_tc(2);
        E_act_overall_tc(3, i : i+del_t/dtt-1) = A_tc(3);
        E_act_overall_tc(4, i : i+del_t/dtt-1) = A_tc(4);
        
        n_act(:) = 0;
        n_act1(:) = 0;
        n_act2(:) = 0;
    else
        n_act(1) = n_act(1) + length(fired(fired <= Ntc)); % AW1
        n_act(2) = n_act(2) + length(fired(fired > Ntc & fired <= 2*Ntc)); % PW
        n_act(3) = n_act(3) + length(fired(fired > 2*Ntc & fired <= 3*Ntc)); % AW2
        n_act(4) = n_act(4) + length(fired(fired > 3*Ntc & fired <= 4*Ntc)); % AW3
        
        n_act1(1) = n_act1(1) + length(fired(fired <= Ntc*p_th)); % AW1
        n_act1(2) = n_act1(2) + length(fired(fired > Ntc & fired <= Ntc+Ntc*p_th)); % PW
        n_act1(3) = n_act1(3) + length(fired(fired > 2*Ntc & fired <= 2*Ntc+Ntc*p_th)); % AW2
        n_act1(4) = n_act1(4) + length(fired(fired > 3*Ntc & fired <= 3*Ntc+Ntc*p_th)); % AW3
        
        n_act2(1) = n_act2(1) + length(fired(fired > Ntc*p_th & fired <= Ntc)); % AW1
        n_act2(2) = n_act2(2) + length(fired(fired > Ntc+Ntc*p_th & fired <= 2*Ntc)); % PW
        n_act2(3) = n_act2(3) + length(fired(fired > 2*Ntc+Ntc*p_th & fired <= 3*Ntc)); % AW2
        n_act2(4) = n_act2(4) + length(fired(fired > 3*Ntc+Ntc*p_th & fired <= 4*Ntc)); % AW3
    end
        
    vv(fired, end + 1) = 30;
    unfired = find(v < 30); 
    vv(unfired, end) = v(unfired);
    v(fired) = c(fired);
    u(fired) = u(fired) + d(fired);  
      
    G(:,fired) = G(:,fired) + g(:,fired);
    G(:, unfired(unfired <= 4*Ntc)) = G(:, unfired(unfired <= 4*Ntc)) - dtt*G(:, unfired(unfired <= 4*Ntc))/tau_ampa;
    G(:, unfired(unfired > 4*Ntc)) = G(:, unfired(unfired > 4*Ntc)) - dtt*G(:, unfired(unfired > 4*Ntc))/tau_gabaa;
    
    I_in(:, 1:4*Ntc) = G(:, 1:4*Ntc).*repmat(v-v_ampa, [1 4*Ntc]);
    I_in(:, 4*Ntc+1:end) = G(:, 4*Ntc+1:end).*repmat(v-v_gabaa, [1 4*Nre]);
    
    % Sensory input 
    for j = 1:size(Barr_Stims,1)
        I_s((j-1)*Ntc+1:j*Ntc) = -A*Spa_Temp(:,j,i);
%         I_s((j-1)*Nre+1+4*Ntc:j*Nre+4*Ntc) = -A*Spa_Temp(:,j,i); % RE cells in each barreloid receive sensory input
    end
    if sum(i == Stim_Onsets(2:end)) % TC cells receiving sensory input vary from trial to trial 
        Conn_sth = zeros(4*Ntc+4*Nre, 1);
        Conn_sth([randi(Ntc*p_th,1,k_sth1) Ntc+randi(Ntc*p_th,1,k_sth1) ...
            2*Ntc+randi(Ntc*p_th,1,k_sth1) 3*Ntc+randi(Ntc*p_th,1,k_sth1)]) = 1;
        Conn_sth([randi([Ntc*p_th+1 Ntc],1,k_sth2) Ntc+randi([Ntc*p_th+1 Ntc],1,k_sth2) ...
            2*Ntc+randi([Ntc*p_th+1 Ntc],1,k_sth2) 3*Ntc+randi([Ntc*p_th+1 Ntc],1,k_sth2)]) = 1; 
    end
    I_s = I_s.*Conn_sth; % (pA)
    
    
    % Cortico-thalamic input
    E_L6th = [E_L6(AW1(1), AW1(2))*ones(Ntc,1); E_L6(PW(1), PW(2))*ones(Ntc,1); E_L6(AW2(1), AW2(2))*ones(Ntc,1); E_L6(AW3(1), AW3(2))*ones(Ntc,1);...
              E_L6(AW1(1), AW1(2))*ones(Nre,1); E_L6(PW(1), PW(2))*ones(Nre,1); E_L6(AW2(1), AW2(2))*ones(Nre,1); E_L6(AW3(1), AW3(2))*ones(Nre,1)];
    I_L6th = -(E_L6th-thres_L62th) .* Conn_L6th; % (pA)
    I_L6th(I_L6th > 0) = 0;
    
    % background noisy input
    I_noise = 0.05*(rand(4*Ntc+4*Nre, 1)-0); % Uniform distributed network noise %0.1
%     I_noise = 0.25*(rand(4*Ntc+4*Nre, 1)-0.5); % Uniform distributed network noise %0.1
%     I_noise = [0.35*(rand(Ntc, 1)-0.5); 0.35*(rand(Nre, 1)-0.5)]; % Uniform distributed network noise
%     I_noise = zeros(4*Ntc+4*Nre, 1);
    
    I = I_s + I_L6th + I_noise + sum(I_in, 2);
    
    v = v + dtt*(0.04*v.^2 + 5*v + 140 - u - I);
    u = u + dtt*a.*(b.*v - u);  
    
    % Tracking sensory input
%     I_s_overall(:,i) = I_s;
%     I_L6th_overall(:,i) = I_L6th; 
%     I_in_overall(:,i) = sum(I_in, 2);
    
     
    %% ---- L4 ---- :
    A_tc_s = A_tc-thres_thL4;
    A_tc_s(A_tc_s < 0) = 0;
    A_tc_s1 = A_tc1-thres_thL4;
    A_tc_s1(A_tc_s1 < 0) = 0;
    A_tc_s2 = A_tc2-thres_thL4;
    A_tc_s2(A_tc_s2 < 0) = 0;
    % Calculating the sensory (thalamocortical) input to all barrels at the current time-step i:
    for j = 1:size(Barr_Stims,1)
        for m = 1:M
            for n = 1:N
                switch (m-Barr_Stims(j,1))^2+(n-Barr_Stims(j,2))^2 % Distance between the stimulated barrel and barrel under consideration
                    case 0 % Distance = 0
                        s_E1(m,n,j) = U_s*z1(m,n,j)*A_tc_s1(j)*s_f(1);
                        s_E2(m,n,j) = U_s*z2(m,n,j)*A_tc_s2(j)*s_f(1);
                    case 1 % Distance = 1
                        s_E1(m,n,j) = U_s*z1(m,n,j)*A_tc_s1(j)*s_f(2);
                        s_E2(m,n,j) = U_s*z2(m,n,j)*A_tc_s2(j)*s_f(2);
                    case 2 % Distance = sqrt(2)
                        s_E1(m,n,j) = U_s*z1(m,n,j)*A_tc_s1(j)*s_f(2);
                        s_E2(m,n,j) = U_s*z2(m,n,j)*A_tc_s2(j)*s_f(2);
                    case 4 % Distance = 2
                        s_E1(m,n,j) = U_s*z1(m,n,j)*A_tc_s1(j)*s_f(4);
                        s_E2(m,n,j) = U_s*z2(m,n,j)*A_tc_s2(j)*s_f(4);
                    otherwise % Distance > 2
                        s_E1(m,n,j) = 0;
                        s_E2(m,n,j) = 0;
                end
            end
        end
    end
    

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
    H = J*EUx + J_1*EUx_1 + J_2*EUx_2; % The total synaptic input each barrel receives

    % Adding thalamocortical input:
    H = H + sum(J_thL41*s_E1 + J_thL42*s_E2, 3);
%     s_E_thres = s_E - thres_thL4;
%     s_E_thres(s_E_thres <  0) = 0;
%     H = H + sum(J_thL4*s_E,3);
    

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
    z1 = z1 + dt*((1 - z1)./tau_rec_s - s_E1);
    z2 = z2 + dt*((1 - z2)./tau_rec_s - s_E2);

    % Tracking the activities of all barrels:
    E_act_overall(:,:,i) = E; 
    x_act_overall(:,:,i) = x; 
    z_act_overall(:,:,i) = z1(:,:,1); 
    z_act_overall2(:,:,i) = z1(:,:,2);
%     s_E_plot(:,:,i) = sum(J_thL4*s_E,3);
    
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
    H_L6 = J_L6*EUx_L6 + J_1_L6*EUx_1_L6 + J_2_L6*EUx_2_L6 + J_L46*U_L46.*E.*x_L46;  % The total synaptic input each L6 infrabarrel receives

    % Tracking L4 input to L6
    E_act_overall_L46(:,:,i) = J_L46*U_L46.*E.*x_L46;
    x_act_overall_L46(:,:,i) = x_L46;
    
    % The dynamics of synaptic input:
    h_L6 = h_L6 + (dt/tau_h_L6)*(-h_L6 + H_L6);
    
    % Implementing the non-linearity:
    E_L6 = h_L6;  
    E_L6 = E_L6 - Gain_thres_L6; 
    E_L6 = E_L6 * Gain_slope_L6;
    E_L6(E_L6 <  0) = 0;   
    E_L6(E_L6 > Gain_max_L6) = Gain_max_L6;   

    % The variables' dynamics:
    x_L6 = x_L6 + dt*((1 - x_L6)./tau_rec_L6 - U_L6.*E_L6.*x_L6);
    x_L46 = x_L46 + dt*((1 - x_L46)./tau_rec_L46 - U_L46.*E.*x_L46);

    % Tracking the activities of all barrels
    E_act_overall_L6(:,:,i) = E_L6; 
    x_act_overall_L6(:,:,i) = x_L6; 
    
end

% vv_samp = vv([50 150 400+50 400+150], :); % Temporal profiles of membrane potential of 2 TC and 2 RE neurons
vv_samp = vv([5 10], :); 
vv = []; % to save memory

% Recording Equilibrium Conditions:
E_eq_tc = E_act_overall_tc(:,num_steps_eq); % TC cells
z_eq = z_act_overall(:,:,num_steps_eq); % AW1
z_eq2 = z_act_overall2(:,:,num_steps_eq); % PW

E_eq = E_act_overall(:,:,num_steps_eq); % L4
x_eq = x_act_overall(:,:,num_steps_eq); 

E_eq_L6 = E_act_overall_L6(:,:,num_steps_eq); % L6
x_eq_L6 = x_act_overall_L6(:,:,num_steps_eq);  




%{
% Activities of L4 barrels surrounding PW 
E_plot = zeros(9,num_steps);
E_plot(1,:) = reshape(E_act_overall(PW(1)-1,PW(2)-1,:), [1, num_steps]);
E_plot(2,:) = reshape(E_act_overall(PW(1)-1,PW(2)+1,:), [1, num_steps]);
E_plot(3,:) = reshape(E_act_overall(PW(1)-1,PW(2),:), [1, num_steps]);
E_plot(4,:) = reshape(E_act_overall(PW(1),PW(2)-1,:), [1, num_steps]);
E_plot(5,:) = reshape(E_act_overall(PW(1),PW(2),:), [1, num_steps]);
E_plot(6,:) = reshape(E_act_overall(PW(1),PW(2)+1,:), [1, num_steps]);
E_plot(7,:) = reshape(E_act_overall(PW(1)+1,PW(2),:), [1, num_steps]);
E_plot(8,:) = reshape(E_act_overall(PW(1)+1,PW(2)+1,:), [1, num_steps]);
E_plot(9,:) = reshape(E_act_overall(PW(1)+1,PW(2)-1,:), [1, num_steps]);

% Resource of L4 barrels surrounding PW 
x_plot = zeros(9,num_steps);
x_plot(1,:) = reshape(x_act_overall(PW(1)-1,PW(2)-1,:), [1, num_steps]);
x_plot(2,:) = reshape(x_act_overall(PW(1)-1,PW(2)+1,:), [1, num_steps]);
x_plot(3,:) = reshape(x_act_overall(PW(1)-1,PW(2),:), [1, num_steps]);
x_plot(4,:) = reshape(x_act_overall(PW(1),PW(2)-1,:), [1, num_steps]);
x_plot(5,:) = reshape(x_act_overall(PW(1),PW(2),:), [1, num_steps]);
x_plot(6,:) = reshape(x_act_overall(PW(1),PW(2)+1,:), [1, num_steps]);
x_plot(7,:) = reshape(x_act_overall(PW(1)+1,PW(2),:), [1, num_steps]);
x_plot(8,:) = reshape(x_act_overall(PW(1)+1,PW(2)+1,:), [1, num_steps]);
x_plot(9,:) = reshape(x_act_overall(PW(1)+1,PW(2)-1,:), [1, num_steps]);


E_plot1 = zeros(9,num_steps);
E_plot1(1,:) = reshape(E_act_overall(AW1(1)-1,AW1(2)-1,:), [1, num_steps]);
E_plot1(2,:) = reshape(E_act_overall(AW1(1)-1,AW1(2)+1,:), [1, num_steps]);
E_plot1(3,:) = reshape(E_act_overall(AW1(1)-1,AW1(2),:), [1, num_steps]);
E_plot1(4,:) = reshape(E_act_overall(AW1(1),AW1(2)-1,:), [1, num_steps]);
E_plot1(5,:) = reshape(E_act_overall(AW1(1),AW1(2),:), [1, num_steps]);
E_plot1(6,:) = reshape(E_act_overall(AW1(1),AW1(2)+1,:), [1, num_steps]);
E_plot1(7,:) = reshape(E_act_overall(AW1(1)+1,AW1(2),:), [1, num_steps]);
E_plot1(8,:) = reshape(E_act_overall(AW1(1)+1,AW1(2)+1,:), [1, num_steps]);
E_plot1(9,:) = reshape(E_act_overall(AW1(1)+1,AW1(2)-1,:), [1, num_steps]);

x_plot1 = zeros(9,num_steps);
x_plot1(1,:) = reshape(x_act_overall(AW1(1)-1,AW1(2)-1,:), [1, num_steps]);
x_plot1(2,:) = reshape(x_act_overall(AW1(1)-1,AW1(2)+1,:), [1, num_steps]);
x_plot1(3,:) = reshape(x_act_overall(AW1(1)-1,AW1(2),:), [1, num_steps]);
x_plot1(4,:) = reshape(x_act_overall(AW1(1),AW1(2)-1,:), [1, num_steps]);
x_plot1(5,:) = reshape(x_act_overall(AW1(1),AW1(2),:), [1, num_steps]);
x_plot1(6,:) = reshape(x_act_overall(AW1(1),AW1(2)+1,:), [1, num_steps]);
x_plot1(7,:) = reshape(x_act_overall(AW1(1)+1,AW1(2),:), [1, num_steps]);
x_plot1(8,:) = reshape(x_act_overall(AW1(1)+1,AW1(2)+1,:), [1, num_steps]);
x_plot1(9,:) = reshape(x_act_overall(AW1(1)+1,AW1(2)-1,:), [1, num_steps]);
%}


%% Saving
if save_results    
    filename2 = ['Simulation Results/run_par' num2str(par) '_' Cond_Code{co}];
    
    if length(Prob_Arr) > 1
      filename2 = [filename2 '_P' num2str(prob)];
    end  
    if length(ISI_Arr) > 1
       filename2 = [filename2 '_ISI' num2str(ISI)];
    end 
    if length(A_Arr) > 1
       filename2 = [filename2 '_A' num2str(A)];
    end
    if num_trials > 1
       filename2 = [filename2 '_Tr' num2str(tr)];
    end
    
    filename2 = [filename2 '.mat'];
    save(filename2)
end

curr_time2 = clock;
disp([num2str(curr_time2(4)) ':' num2str(curr_time2(5)) ' ' cond ' condition, ISI = ' num2str(ISI) 's, A = ' num2str(A) 'Hz, Trial ' num2str(tr) ' done']); 

toc
