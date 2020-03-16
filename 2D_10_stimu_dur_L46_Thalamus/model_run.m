%% Model Run
% This code runs the network saved and brought to equilibrium by
% model_init, according to the specified condition defined by cond

% cond = 1;
% A = 1;
% ISI = 1;
% stim = 1;
% prob = 1;

tic
n_stim = 80; % Total no. of stimuli (Best take a product of 4) %120
t_prot = n_stim*(duration + ISI) + 2*post_stim; % Total time of the protocol
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
         rng(66);
        Oddball = [ones(ceil((1-Probs)*n_stim),1)*Barr_Stims(1,:); ones(ceil(Probs*n_stim),1)*Barr_Stims(2,:)];
        Oddball = Oddball(randperm(n_stim),:); % n_stim*2
%         Oddball = [3     2;4     2;4     2;4     2;4     2;3     2; 4     2; 4     2;  4     2;3     2; 4     2;   4     2;4     2; 4     2;3     2; 4     2; 4     2;...
%             3     2; 4     2;  4     2;4     2; 3     2; 3     2;4     2;4     2;4     2; 4     2; 3     2; 4     2;4     2;3     2; 4     2; 4     2; 4     2;4     2;4     2;...
%      3     2; 4     2; 4     2;4     2];
    case('High') % AW1 as deviant and PW as standard  
        Oddball = [ones(ceil(Probs*n_stim),1)*Barr_Stims(1,:); ones(ceil((1-Probs)*n_stim),1)*Barr_Stims(2,:)];
        Oddball = Oddball(randperm(n_stim),:); 
    case('Equal')
        Oddball = [ones(ceil(Probs*n_stim),1)*Barr_Stims(1,:); ones(ceil((1-Probs)*n_stim),1)*Barr_Stims(2,:)]; % Here Probs=0.5
        Oddball = Oddball(randperm(n_stim),:); 
    case('Many Standard')  % AW1 PW AW2 AW3, each stimulated with probability of 0.25
        rng(66);
        whis_reps = n_stim/size(Barr_Stims,1); % Repetitons of each stimulated whisker 
%         Oddball = repmat(Barr_Stims,whis_reps,1);
        Oddball = [ones(whis_reps,1)*Barr_Stims(1,:); ones(whis_reps,1)*Barr_Stims(3,:);...
            ones(whis_reps,1)*Barr_Stims(4,:); ones(whis_reps,1)*Barr_Stims(2,:)];
        Oddball = Oddball(randperm(n_stim),:);
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
        Spa_Temp(1, sum(Barr_Stims == Oddball(ns,:),2) == 2, Time_Ind(ns,1):Time_Ind(ns,2)) = Single_Stim; 
        % Amplitude of each stimulus along time indices
    end
    Stim_Onsets(ns) = Time_Ind(ns,1);
end

%% Returning to equilibrium conditions and initializing sensory inputs:
%% ---- L4 ---- :
E = E_eq;
x = x_eq; 
% z = z_eq*ones(M,N);
z = z_eq*ones(M,N,size(Barr_Stims,1));

% The following lines are for tracking activity of all barrels 
E_act_overall = zeros(M,N,num_steps); 
E_act_overall(:,:,1:floor(t_eq/dt)) = E_act;
x_act_overall = zeros(M,N,num_steps); 
x_act_overall(:,:,1:floor(t_eq/dt)) = x_act;
z_act_overall = zeros(M,N,num_steps); 
z_act_overall(:,:,1:floor(t_eq/dt)) = z_eq*ones(M,N,floor(t_eq/dt));

% Sensory inputs to each barrel
s_E = zeros(M,N,size(Barr_Stims,1)); 
s_e = zeros(M,N); 


%% ---- L6 ---- :
E_L6 = E_eq_L6;
x_L6 = x_eq_L6; 

% The following lines are for tracking activity of all barrels 
E_act_overall_L6 = zeros(M,N,num_steps); 
E_act_overall_L6(:,:,1:floor(t_eq/dt)) = E_act_L6;
x_act_overall_L6 = zeros(M,N,num_steps); 
x_act_overall_L6(:,:,1:floor(t_eq/dt)) = x_act_L6;

%% ----- Thalamus ----- :
E_Th = E_eq_Th;

% The following lines are for tracking activity of all barrels 
E_act_overall_Th = zeros(M,N,num_steps); 
E_act_overall_Th(:,:,1:floor(t_eq/dt)) = E_act_Th;

%% Dynamic Loop
i = floor(t_eq/dt) + 1;
while i < num_steps
    %% ---- L4 ---- :
%     % Calculating the sensory (thalamocortical) input to all barrels at the current time-step i:
%     for j = 1:size(Barr_Stims,1)
%         for m = 1:M
%             for n = 1:N
%                 switch (m-Barr_Stims(j,1))^2+(n-Barr_Stims(j,2))^2 % Distance between the stimulated barrel and barrel under consideration
%                     case 0 % Distance = 0
%                         s_E(m,n,j) = Spa_Temp(:,j,i)*s_f(1);
%                     case 1 % Distance = 1
%                         s_E(m,n,j) = Spa_Temp(:,j,i)*s_f(2);
%                     case 2 % Distance = sqrt(2)
%                         s_E(m,n,j) = Spa_Temp(:,j,i)*s_f(3);
%                     case 4 % Distance = 2
%                         s_E(m,n,j) = Spa_Temp(:,j,i)*s_f(4);
%                     otherwise % Distance > 2
%                         s_E(m,n,j) = 0;
%                 end
%             end
%         end
%     end
%     s_e = U_s*A*z.*sum(s_E,3);

    for j = 1:size(Barr_Stims,1)
        for m = 1:M
            for n = 1:N
                switch (m-Barr_Stims(j,1))^2+(n-Barr_Stims(j,2))^2 % Distance between the stimulated barrel and barrel under consideration
                    case 0 % Distance = 0
                        s_E(m,n,j) = U_s*A*z(m,n,j)*Spa_Temp(:,j,i)*s_f(1);
                    case 1 % Distance = 1
                        s_E(m,n,j) = U_s*A*z(m,n,j)*Spa_Temp(:,j,i)*s_f(2);
                    case 2 % Distance = sqrt(2)
                        s_E(m,n,j) = U_s*A*z(m,n,j)*Spa_Temp(:,j,i)*s_f(3);
                    case 4 % Distance = 2
                        s_E(m,n,j) = U_s*A*z(m,n,j)*Spa_Temp(:,j,i)*s_f(4);
                    otherwise % Distance > 2
                        s_E(m,n,j) = 0;
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

%     % Adding the sensory input to synaptic input H:
%     H = H + s_e;

    % Adding the sensory input to synaptic input H:
    H = H + sum(s_E,3);
    
    % Adding thalamus-to-L4 latency:
    if i-floor(tau_ThL4/dt) > 0
        H = H + E_act_overall_Th(:, :, i-floor(tau_ThL4/dt))*U_s; % One-to-one
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
    z = z + dt*((1 - z)./tau_rec_s - s_E);

    % Tracking the activities of all barrels:
    E_act_overall(:,:,i) = E; 
    x_act_overall(:,:,i) = x; 
    z_act_overall(:,:,i) = z(:,:,1); 
    
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
    E_L6(E_L6 > Gain_max_L6) = Gain_max_L6;   

    % The variables' dynamics:
    x_L6 = x_L6 + dt*((1 - x_L6)./tau_rec_L6 - U_L6.*E_L6.*x_L6);         

    % Tracking the activities of all barrels
    E_act_overall_L6(:,:,i) = E_L6; 
    x_act_overall_L6(:,:,i) = x_L6;
    
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
    E_act_overall_Th(:,:,i) = E_Th; 
    
    i = i + 1;
end

% Activities of L4 barrels surrounding PW to plot
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

% Resource of L4 barrels surrounding PW to plot
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

%{
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
    
    if p == 1  && cc == 1 && bb == 1 && aa == 1 && tr == 1 
        save(filename2)  
    else
        save(filename2, 'E_act_overall', 'Oddball') 
    end 
end

curr_time2 = clock;
disp([num2str(curr_time2(4)) ':' num2str(curr_time2(5)) ' ' cond ' condition, ISI = ' num2str(ISI) 's, A = ' num2str(A) 'Hz, Trial ' num2str(tr) ' done']); 

toc
