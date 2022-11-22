%% Run simulations with different network and/or stimulus parameters
% The code implementation is inspired and partially adapted from Yarden TS, Nelken I (2017).
% SSA protocols for rat's whisker identity are based on Musall et al. (2017).

close all;
clear;

% SSA Conditions
Cond{1} = 'Low'; % Oddball low
Cond{2} = 'High'; % Oddball high
Cond{3} = 'Equal';
Cond{4} = 'Many Standard';
Cond{5} = 'Deviant Alone 1';
Cond{6} = 'Deviant Alone 2';
Cond{7} = 'Sequenced';
Cond{8} = 'Randomized';

Cond_Code{1} = 'L';
Cond_Code{2} = 'H';
Cond_Code{3} = 'E';
Cond_Code{4} = 'MS';
Cond_Code{5} = 'DA1';
Cond_Code{6} = 'DA2';
Cond_Code{7} = 'S';
Cond_Code{8} = 'R';

% Network Parameters
Par_Arr = [0];

% Stimulus Parameters
Cond_Arr = [1 4]; % Protocol index
Prob_Arr = [0.25]; % Probability of deviant occurrence
A_Arr = [5]; % Stimulus intensity
ISI_Arr = [1-0.010]; % inter-stimulus interval (stimulus offset to onset)

% Number of repeated trials for each parameters setup 
num_trials = 1; 


%% Loops over parameters
for p = 1:length(Par_Arr)

    par = Par_Arr(p);
    run model_init.m;

    for co = 1:length(Cond_Arr) 
        for cc = 1:length(Prob_Arr)
            for bb = 1:length(ISI_Arr) 
                for aa = 1:length(A_Arr)
                    for tr = 1:num_trials
                        save('Simulation Results/meta_data.mat','Cond','Cond_Code','Par_Arr','Cond_Arr', ...
                            'Prob_Arr','ISI_Arr','A_Arr','num_trials','p','par','co','cc','bb','aa','tr','dt');
                        clear;
                        load('Simulation Results/meta_data.mat');

                        load(['Simulation Results/run_par' num2str(par) '_initialization.mat']);
                        load('Simulation Results/meta_data.mat'); 

                        cond = Cond{Cond_Arr(co)};
                        prob = Prob_Arr(cc);
                        ISI = ISI_Arr(bb);
                        A = A_Arr(aa);

                        save_results = 1; % Whether save the workspace genenrated by running model_run.m
                        run model_run.m;  
                    end
                end
            end
        end
        disp(['Calculation for network with parameter = ' num2str(par) ', ' cond ' condition completes']);
    end
    disp('--------------------------------------------------------------------------------------');
end
 

