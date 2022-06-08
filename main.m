%% Run Oddball vs Many-standard conditions 
% The code is adapted from Yarden TS, Nelken I (2017)

close all;
clear;

%% Network Parameters:

% tau_rec
% Par_Arr = [0.1:0.05:1.5]; 
Par_Arr = 0.03:0.01:0.05;
% Par_Arr = [0.18];
% % U
% Par_Arr = [0.4];
% Par_Arr = [0.52];

% lambda
% Par_Arr = [1.32];
% Par_Arr = [1.2];

% J_1
% Par_Arr = [0.11];
% Par_Arr = [0.09];
% Par_Arr = [0.1];

% % Gain_thres
% Par_Arr = [7:0.5:8];

%% Dynamic Loop
for p = 1:length(Par_Arr)
    par = Par_Arr(p);

    Init = 1; % Init = 1 means that a new network will be initialized; set this to 0 to be on previously saved networks  
    if Init
        run model_init.m;
    else
        load(['Simulation Results/run_par' num2str(par) '_initialization.mat']);
    end

    % Stimulation Parameters:
    Prob_Arr = [0.25]; % Sets the probabilitis of stimulation in the Low Condition, 0.1
    A_Arr = [5]; % 6.5
    ISI_Arr = [1-0.010]; % Interval between stimuli (offset to onset) %0.2

    num_trials = 1; 
    for co = [1 4]  % There are 6 basic conditions 
        for cc = 1:length(Prob_Arr)
            for bb = 1:length(ISI_Arr) 
                for aa = 1:length(A_Arr)
                    for tr = 1:num_trials
                        save('Simulation Results/meta_data.mat','Par_Arr','p','par','Prob_Arr','A_Arr','ISI_Arr','num_trials','co', ...
                          'tr','cc','bb','aa','dt','Cond_Code');
                        clear;
                        load('Simulation Results/meta_data.mat');
                        load(['Simulation Results/run_par' num2str(par) '_initialization.mat']);
                        load('Simulation Results/meta_data.mat'); 
                        % Overload some variable like co which is always 4 saved in initialization after running to the outer loop

                        cond = Cond{co};
                        prob = Prob_Arr(cc);
                        ISI = ISI_Arr(bb);
                        A = A_Arr(aa);

                        save_results = 1; % Whether save the workspace genenrated by running model_run.m
                        run model_run.m;  
                    end
                end
            end
        end
        disp(['Calculation for network with parameter = ' num2str(par) ', ' Cond{co} ' condition completes']);
    end
    disp('--------------------------------------------------------------------------------------');
end
 

