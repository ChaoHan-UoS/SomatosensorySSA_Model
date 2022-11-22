%% Many-standard condition
clear
k = 1;
kk = 4;
load('Simulation Results/meta_data.mat'); 
load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat']);

n_stim_plot = 50;
stim_ind = 71;
t = Stim_Onsets(stim_ind);
t_marg = 0.5; % in seconds
time_win = n_stim_plot*(duration+ISI)+t_marg; 
t_onset = t - floor(t_marg/dt);
t_offset = t + floor(time_win/dt);
tim = -t_marg:dt:time_win; 
t_step = 1;

E_plot_L4_ms = reshape(E_act_overall(PW(1),PW(2),:), [1, num_steps]);
x_plot_L4_ms = reshape(x_act_overall(PW(1),PW(2),:), [1, num_steps]);
E_plot_L46_ms = reshape(E_act_overall_L46(PW(1),PW(2),:), [1, num_steps]);
E_plot_L6_ms = reshape(E_act_overall_L6(PW(1),PW(2),:), [1, num_steps]);
x_plot_L6_ms = reshape(x_act_overall_L6(PW(1),PW(2),:), [1, num_steps]);
x_plot_L6_ms_D2 = reshape(x_act_overall_L6(AW1(1),AW1(2),:), [1, num_steps]);
x_plot_L6_ms_D1 = reshape(x_act_overall_L6(AW2(1),AW2(2),:), [1, num_steps]);
x_plot_L6_ms_D3 = reshape(x_act_overall_L6(AW3(1),AW3(2),:), [1, num_steps]);
kk = 1;
S = load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat'],'E_act_overall','x_act_overall','E_act_overall_L46',...
    'E_act_overall_L6','x_act_overall_L6','Spa_Temp','E_act_overall_tc');
E_plot_L4_odd = reshape(S.E_act_overall(PW(1),PW(2),:), [1, num_steps]);
x_plot_L4_odd = reshape(S.x_act_overall(PW(1),PW(2),:), [1, num_steps]);
E_plot_L46_odd = reshape(S.E_act_overall_L46(PW(1),PW(2),:), [1, num_steps]);
E_plot_L6_odd = reshape(S.E_act_overall_L6(PW(1),PW(2),:), [1, num_steps]);
x_plot_L6_odd = reshape(S.x_act_overall_L6(PW(1),PW(2),:), [1, num_steps]);
x_plot_L6_odd_D2 = reshape(S.x_act_overall_L6(AW1(1),AW1(2),:), [1, num_steps]);
Spa_Temp_odd = S.Spa_Temp;
E_act_overall_tc_odd = S.E_act_overall_tc;

f = figure;
% ---- stimulus sequence ---- 
subplot(8,1,1);
plot(tim, reshape(Spa_Temp_odd(1,1,t_onset:t_offset), [1 length(t_onset:t_offset)])+1, 'b',...
     tim, reshape(Spa_Temp_odd(1,2,t_onset:t_offset), [1 length(t_onset:t_offset)])+2.5, 'r', LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca,'FontSize', AXES_FONTSIZE, 'XTickLabel', [], 'YTick', [1 2.5], 'YTickLabel', {'D2','C2'},'TickDir','out','box','off');
axis([-t_marg time_win 0 7]);
axis off

subplot(8,1,2);
plot(tim, reshape(Spa_Temp(1,1,t_onset:t_offset), [1 length(t_onset:t_offset)])+1, 'b',...
     tim, reshape(Spa_Temp(1,3,t_onset:t_offset), [1 length(t_onset:t_offset)])+1+1.5, 'b',...
     tim, reshape(Spa_Temp(1,4,t_onset:t_offset), [1 length(t_onset:t_offset)])+1+2*1.5, 'b',... 
     tim, reshape(Spa_Temp(1,2,t_onset:t_offset), [1 length(t_onset:t_offset)])+1+3*1.5, 'g',LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca,'FontSize', AXES_FONTSIZE, 'XTickLabel', [], 'YTick', [1 1+1.5 1+2*1.5 1+3*1.5], 'YTickLabel', {'D2','D1','D3','C2'},'TickDir','out','box','off');
axis([-t_marg time_win 0 7]);
axis off

% ---- Odd deviant in TCs ---- 
subplot(4,1,2) 
% plot(tim, E_act_overall_tc_odd(1,t_onset:t_offset),'b',tim, E_act_overall_tc_odd(2,t_onset:t_offset),'r',...
%     'LineWidth',LineWidth,'MarkerSize',MarkerSize); % Rate of TC in AW1/PW
plot(tim, E_act_overall_tc_odd(2,t_onset:t_offset),'r','LineWidth',LineWidth,'MarkerSize',MarkerSize); % Rate of TC in oddball PW
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'TickDir','out','box','off','FontSize', AXES_FONTSIZE)  
ylabel('A_{TC} (Hz)')
% legend('odd D2','odd C2','Box','off')

% ---- MS deviant in TCs ---- 
subplot(4,1,3) 
% plot(tim, E_act_overall_tc(1,t_onset:t_offset),'b',tim, E_act_overall_tc(2,t_onset:t_offset),'g',...
%     'LineWidth',LineWidth,'MarkerSize',MarkerSize); % Rate of TC in AW1/PW
plot(tim, E_act_overall_tc(2,t_onset:t_offset),'g','LineWidth',LineWidth,'MarkerSize',MarkerSize); % Rate of TC in many-standards PW
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[stim_ind-1:t_step:stim_ind-1+time_win],'XLim',[-t_marg time_win],'TickDir','out','box','off','FontSize', AXES_FONTSIZE)  
xlabel('time (s)')
ylabel('A_{TC} (Hz)')
% legend('ms D2','ms C2','Box','off')

set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition'))


%% Mean thalamic activity
clear;
par_index = 1;
load('Simulation Results/meta_data.mat'); 
t_marg = 0.020; % in seconds; marginal time left before stimulus onset to show baseline activity
time_win = 1-t_marg; % in seconds; time window used to calculate SI and CSI, starting from t_marg before each stimulus onset to time_win 
n_odddev= 0; % Number of oddball deviant 
n_oddstd = 0; % Number of oddball standard 
n_msdev = 0; % Number of many-standard deviant 
E_sum_odddev = zeros(3, floor(time_win/dt)+floor(t_marg/dt)+1); 
E_sum_oddstd = E_sum_odddev;
E_sum_msdev = E_sum_odddev;
Spcount_odddev = zeros(3,1);
Spcount_oddstd = zeros(3,1);
Spcount_msdev = zeros(3,1);
CSI = zeros(3,1); % Quantify true deviance detection in L4, L6 and TCs
SI = zeros(3,1); % Quantify SSA in L4, L6 and TCs
for kk = [1 4]
    load(['Simulation Results/run_par' num2str(Par_Arr(par_index)) '_' Cond_Code{kk} '.mat']);
    for ns = 1:n_stim
        if kk == 1
            if Oddball(ns,:) == PW
                E_sum_odddev(1,:) = E_sum_odddev(1,:) + reshape(E_act_overall(PW(1), PW(2), Stim_Onsets(ns)-floor(t_marg/dt):Stim_Onsets(ns)+floor(time_win/dt)),[1 floor(time_win/dt)+floor(t_marg/dt)+1]); % L4
                E_sum_odddev(2,:) = E_sum_odddev(2,:) + reshape(E_act_overall_L6(PW(1), PW(2), Stim_Onsets(ns)-floor(t_marg/dt):Stim_Onsets(ns)+floor(time_win/dt)),[1 floor(time_win/dt)+floor(t_marg/dt)+1]); % L6
                E_sum_odddev(3,:) = E_sum_odddev(3,:) + E_act_overall_tc(2, Stim_Onsets(ns)-floor(t_marg/dt):Stim_Onsets(ns)+floor(time_win/dt)); % TC
                n_odddev = n_odddev + 1;
            else
                E_sum_oddstd(1,:) = E_sum_oddstd(1,:) + reshape(E_act_overall(AW1(1), AW1(2), Stim_Onsets(ns)-floor(t_marg/dt):Stim_Onsets(ns)+floor(time_win/dt)),[1 floor(time_win/dt)+floor(t_marg/dt)+1]);
                E_sum_oddstd(2,:) = E_sum_oddstd(2,:) + reshape(E_act_overall_L6(AW1(1), AW1(2), Stim_Onsets(ns)-floor(t_marg/dt):Stim_Onsets(ns)+floor(time_win/dt)),[1 floor(time_win/dt)+floor(t_marg/dt)+1]);
                E_sum_oddstd(3,:) = E_sum_oddstd(3,:) + E_act_overall_tc(1, Stim_Onsets(ns)-floor(t_marg/dt):Stim_Onsets(ns)+floor(time_win/dt));
                n_oddstd = n_oddstd + 1;
            end
        else
            if Oddball(ns,:) == PW
                E_sum_msdev(1,:) = E_sum_msdev(1,:) + reshape(E_act_overall(PW(1), PW(2), Stim_Onsets(ns)-floor(t_marg/dt):Stim_Onsets(ns)+floor(time_win/dt)),[1 floor(time_win/dt)+floor(t_marg/dt)+1]);
                E_sum_msdev(2,:) = E_sum_msdev(2,:) + reshape(E_act_overall_L6(PW(1), PW(2), Stim_Onsets(ns)-floor(t_marg/dt):Stim_Onsets(ns)+floor(time_win/dt)),[1 floor(time_win/dt)+floor(t_marg/dt)+1]);
                E_sum_msdev(3,:) = E_sum_msdev(3,:) + E_act_overall_tc(2, Stim_Onsets(ns)-floor(t_marg/dt):Stim_Onsets(ns)+floor(time_win/dt));           
                n_msdev = n_msdev + 1;
            end
        end
    end
end
E_sum_odddev = E_sum_odddev/n_odddev;
E_sum_oddstd = E_sum_oddstd/n_oddstd;
E_sum_msdev = E_sum_msdev/n_msdev;
Spcount_odddev = sum(E_sum_odddev*dt,2); % Area under E_sum
Spcount_oddstd= sum(E_sum_oddstd*dt,2);
Spcount_msdev = sum(E_sum_msdev*dt,2);
SI = (Spcount_odddev - Spcount_oddstd)./(Spcount_odddev + Spcount_oddstd)
CSI = (Spcount_odddev - Spcount_msdev)./(Spcount_odddev + Spcount_msdev)

% Plotting the average LATE repsonses of cells
t_step = 0.100; % in seconds
tmax = 1; % in seconds; max value of time-axis (no more than time_win)
tim = (-t_marg:dt:time_win)*10^3; % in ms

f = figure;
ymax = 20;
subplot(3,1,1)
plot(tim, E_sum_odddev(3,:), '-r', tim, reshape(E_sum_msdev(3,:), [1 floor(time_win/dt)+floor(t_marg/dt)+1]), '-g','LineWidth',LineWidth,'MarkerSize',MarkerSize);
set(gca,'XTick',[0:t_step:time_win]*10^3,'XLim',[-t_marg tmax]*10^3,'YLim',[0 ymax],'XTickLabelRotation',0,'TickDir','out','box','off','FontSize', AXES_FONTSIZE);
xlabel('time (ms)');
ylabel('A_{TC} (Hz)');
legend('odd dev','ms dev');
legend('boxoff');
title('Thalamus');
set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition'))

