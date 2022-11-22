%% Parameter dependence of SSA 
clear;
load('Simulation Results/meta_data.mat'); 
SI = zeros(num_trials,length(Par_Arr)); % Quantify SSA in L4
CSI = zeros(num_trials, length(Par_Arr)); % Quantifying true deviance detection in L4
time_win = 0.400; % in seconds; time window used to calculate SI and CSI, starting from each stimulus onset

for k = 1:length(Par_Arr)
    for kkk = 1:num_trials
        n_odddev = 0; % Number of oddball deviant 
        n_oddstd = 0; % Number of oddball standard 
        n_msdev = 0; % Number of many-standard deviant
        E_sum_odddev = zeros(1, time_win/dt);
        E_sum_oddstd = E_sum_odddev;
        E_sum_msdev = E_sum_odddev;
        Spcount_odddev = 0;
        Spcount_oddstd = 0;
        Spcount_msdev = 0;
        for kk = [1 4]
            if num_trials > 1
                load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '_Tr' num2str(kkk) '.mat']);
            else
                load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat']);
            end
            for ns = 1:n_stim
                if kk == 1
                    if Oddball(ns, :) == PW   
                        E_sum_odddev = E_sum_odddev + reshape(E_act_overall(PW(1), PW(2), Stim_Onsets(ns):Stim_Onsets(ns)+time_win/dt-1),1,time_win/dt); % L4
                        n_odddev = n_odddev + 1;
                    else
                        E_sum_oddstd = E_sum_oddstd + reshape(E_act_overall(AW1(1), AW1(2), Stim_Onsets(ns):Stim_Onsets(ns)+time_win/dt-1),1,time_win/dt);
                        n_oddstd = n_oddstd + 1;
                    end
                else
                    if Oddball(ns, :) == PW
                        E_sum_msdev = E_sum_msdev + reshape(E_act_overall(PW(1), PW(2), Stim_Onsets(ns):Stim_Onsets(ns)+time_win/dt-1),1,time_win/dt);
                        n_msdev = n_msdev + 1;
                    end
                end
            end
        end
        E_sum_odddev = E_sum_odddev/n_odddev;
        E_sum_oddstd = E_sum_oddstd/n_oddstd;
        E_sum_msdev = E_sum_msdev/n_msdev;
        Spcount_odddev = sum(E_sum_odddev*dt); % Area under E_sum
        Spcount_oddstd = sum(E_sum_oddstd*dt);
        Spcount_msdev = sum(E_sum_msdev*dt);
        SI(kkk, k) = (Spcount_odddev - Spcount_oddstd)/(Spcount_odddev + Spcount_oddstd);
        CSI(kkk, k) = (Spcount_odddev - Spcount_msdev)/(Spcount_odddev + Spcount_msdev);
    end
end
SI_mean = mean(SI,1);
SI_std = std(SI,0,1);
CSI_mean = mean(CSI,1);
CSI_std = std(CSI,0,1);

% f = figure('Visible','off');
f = figure;
hold on;
errorbar(Par_Arr, SI_mean, SI_std, '-or');
errorbar(Par_Arr, CSI_mean, CSI_std, '-og');
plot(Par_Arr, zeros(size(Par_Arr)), ':k');
% xlabel('\tau_{rec}(s)');
% xlim([0.1 1.5]);
% xticks([0.1:0.1:1.5]);
% ylim([-0.1 0.5]);
legend('SI','CSI','Location','northwest');
saveas(f,'Figures/ParCurve.pdf');
 

%% SSA and TDD
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

% Plotting the average LATE repsonses of L4, L6 and TC cells
t_step = 0.100; % in seconds
tmax = 0.450; % in seconds; max value of time-axis (no more than time_win)
tim = (-t_marg:dt:time_win)*10^3; % in ms
f = figure;

% ---- L4 ---- 
ymax = 30;
subplot(3,2,1)
plot(tim, E_sum_odddev(1,:), '-r', tim, E_sum_oddstd(1,:), '-b',...
     tim, E_sum_msdev(1,:), '-g','LineWidth',LineWidth,'MarkerSize',MarkerSize);
set(gca,'XTick',[0:t_step:time_win]*10^3,'XLim',[-t_marg tmax]*10^3,'YLim',[0 ymax],'YTick',0:15:ymax,'XTickLabelRotation',0,'TickDir','out','box','off','FontSize', AXES_FONTSIZE);
xlabel('time (ms)');
ylabel('A_{L4} (Hz)');
% axis square
legend('odd dev','odd std','ms dev');
legend('boxoff');
% title(['L4, SI = ' num2str(SI(1)) ', CSI = ' num2str(CSI(1))]);
title('L4');

% ---- L6 ---- 
ymax = 10;
subplot(3,2,3)
plot(tim, E_sum_odddev(2,:), '-r', tim, E_sum_oddstd(2,:), '-b',...
     tim, reshape(E_sum_msdev(2,:), [1 floor(time_win/dt)+floor(t_marg/dt)+1]), '-g','LineWidth',LineWidth,'MarkerSize',MarkerSize);
set(gca,'XTick',[0:t_step:time_win]*10^3,'XLim',[-t_marg tmax]*10^3,'YLim',[0 ymax],'XTickLabelRotation',0,'TickDir','out','box','off','FontSize', AXES_FONTSIZE);
xlabel('time (ms)');
ylabel('A_{L6} (Hz)');
% axis square
% title(['L6, SI = ' num2str(SI(2)) ', CSI = ' num2str(CSI(2))]);
title('L6');

% ---- Thalamus ---- 
ymax = 20;
subplot(3,2,5)
plot(tim, E_sum_odddev(3,:), '-r', tim, E_sum_oddstd(3,:), '-b',...
     tim, reshape(E_sum_msdev(3,:), [1 floor(time_win/dt)+floor(t_marg/dt)+1]), '-g','LineWidth',LineWidth,'MarkerSize',MarkerSize);
set(gca,'XTick',[0:t_step:time_win]*10^3,'XLim',[-t_marg tmax]*10^3,'YLim',[0 ymax],'XTickLabelRotation',0,'TickDir','out','box','off','FontSize', AXES_FONTSIZE);
xlabel('time (ms)');
ylabel('A_{TC} (Hz)');
% axis square
% title(['Thalamus, SI = ' num2str(SI(3)) ', CSI = ' num2str(CSI(3))]);
title('Thalamus');
set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition'))
saveas(f,'Figures/SSA_late.pdf');


% Plotting the average EARLY repsonses of L4, L6 and TC cells:
t_step = 0.010; % in seconds
tmin = 0.005;
tmax = 0.040;
f = figure;

% ---- L4 ---- :
subplot(3,2,1)
plot(tim, reshape(E_sum_odddev(1,:), [1 floor(time_win/dt)+floor(t_marg/dt)+1]), '-r',...
     tim, reshape(E_sum_oddstd(1,:), [1 floor(time_win/dt)+floor(t_marg/dt)+1]), '-b',...
     tim, reshape(E_sum_msdev(1,:), [1 floor(time_win/dt)+floor(t_marg/dt)+1]), '-g','LineWidth',LineWidth,'MarkerSize',MarkerSize);
set(gca,'XTick',[0:t_step:time_win]*10^3,'XLim',[-tmin tmax]*10^3,'XTickLabelRotation',0,'TickDir','out','box','off','FontSize', AXES_FONTSIZE);
xlabel('time (ms)');
ylabel('A_{L4} (Hz)');
% axis square
title('L4');

% ---- L6 ---- :
subplot(3,2,3)
plot(tim, reshape(E_sum_odddev(2,:), [1 floor(time_win/dt)+floor(t_marg/dt)+1]), '-r',...
     tim, reshape(E_sum_oddstd(2,:), [1 floor(time_win/dt)+floor(t_marg/dt)+1]), '-b',...
     tim, reshape(E_sum_msdev(2,:), [1 floor(time_win/dt)+floor(t_marg/dt)+1]), '-g','LineWidth',LineWidth,'MarkerSize',MarkerSize);
set(gca,'XTick',[0:t_step:time_win]*10^3,'XLim',[-tmin tmax]*10^3,'XTickLabelRotation',0,'TickDir','out','box','off','FontSize', AXES_FONTSIZE);
xlabel('time (ms)');
ylabel('A_{L6} (Hz)');
% axis square
title('L6');

% ---- Thalamus ---- :
subplot(3,2,5)
plot(tim, reshape(E_sum_odddev(3,:), [1 floor(time_win/dt)+floor(t_marg/dt)+1]), '-r',...
     tim, reshape(E_sum_oddstd(3,:), [1 floor(time_win/dt)+floor(t_marg/dt)+1]), '-b',...
     tim, reshape(E_sum_msdev(3,:), [1 floor(time_win/dt)+floor(t_marg/dt)+1]), '-g','LineWidth',LineWidth,'MarkerSize',MarkerSize);
set(gca,'XTick',[0:t_step:time_win]*10^3,'XLim',[-tmin tmax]*10^3,'XTickLabelRotation',0,'TickDir','out','box','off','FontSize', AXES_FONTSIZE);
xlabel('time (ms)');
ylabel('A_{TC} (Hz)');
% axis square
title('Thalamus');
set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition'))
saveas(f,'Figures/SSA_early.pdf');


%% Spatiotemporal distribution of population spikes
clear
k = 1;
kk = 1;
load('Simulation Results/meta_data.mat'); 
load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat']);
t = Stim_Onsets(1);

% ---- Snapshot of spatial propagation of PSs in L4 ---- 
f = figure; % fig1
for i = 1:11
    subplot(5,5,i)
    image([1 4], [1 5], E_act_overall(:,:,t+floor(i*0.005/dt)));
    set(gca,'XTick',1:4,'YTick',1:5,'YTickLabel',{'A','B','C','D','E'},'TickDir','out','box','off','FontSize',AXES_FONTSIZE);
    title([num2str(i*5),' ms']);
    colormap(hot(150))
    axis equal
    axis tight
    c = colorbar;
    c.TickDirection = 'out';
%     c.TickLength  = 0.03;
    c.Box = 'off';
    c.FontSize = COLORBAR_FONTSIZE;
    c.Label.String = 'A_{L4} (Hz)';
    c.Limits = [0 150];
    c.Ticks = [0 50 100 150];
    if i < 11
        c.Visible = 'off';
    end
end
set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition')) 
saveas(f,'Figures/Prop_L41.pdf');

% ---- Snapshot of spatial propagation of PSs in L6 ---- 
f = figure; % fig2
for i = 1:11
    subplot(5,5,i)
    image([1 4], [1 5], E_act_overall_L6(:,:,t+floor(i*0.005/dt)));
    set(gca,'XTick',1:4,'YTick',1:5,'YTickLabel',{'A','B','C','D','E'},'TickDir','out','box','off','FontSize',AXES_FONTSIZE);
    title([num2str(i*5),' ms']);
    colormap(hot(151))
    axis equal
    axis tight
    c = colorbar;
    c.TickDirection = 'out';
%     c.TickLength  = 0.03;
    c.Box = 'off';
    c.FontSize = COLORBAR_FONTSIZE;
    c.Label.String = 'A_{L6} (Hz)';
    c.Limits = [0 150];
    c.Ticks = [0 50 100 150];
    if i < 11
        c.Visible = 'off';
    end
end
set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition')) 
saveas(f,'Figures/Spat_temp_PS.pdf');

%% Temporal profiles of propogation of PSs
clear
k = 1;
kk = 1;
load('Simulation Results/meta_data.mat'); 
load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat']);
t = Stim_Onsets(1);

% ---- Temporal profiles of PS generated in certain barrels ---- :
E_plot = zeros(9,num_steps); % Activities of L4 barrels surrounding AW1 to plot
E_plot(1,:) = reshape(E_act_overall(AW1(1)-1,AW1(2)-1,:), [1, num_steps]);
E_plot(2,:) = reshape(E_act_overall(AW1(1)-1,AW1(2)+1,:), [1, num_steps]);
E_plot(3,:) = reshape(E_act_overall(AW1(1)-1,AW1(2),:), [1, num_steps]);
E_plot(4,:) = reshape(E_act_overall(AW1(1),AW1(2)-1,:), [1, num_steps]);
E_plot(5,:) = reshape(E_act_overall(AW1(1),AW1(2),:), [1, num_steps]);
E_plot(6,:) = reshape(E_act_overall(AW1(1),AW1(2)+1,:), [1, num_steps]);
E_plot(7,:) = reshape(E_act_overall(AW1(1)+1,AW1(2),:), [1, num_steps]);
E_plot(8,:) = reshape(E_act_overall(AW1(1)+1,AW1(2)+1,:), [1, num_steps]);
E_plot(9,:) = reshape(E_act_overall(AW1(1)+1,AW1(2)-1,:), [1, num_steps]);

E_plot1 = zeros(9,num_steps); % Activities of L6 barrels surrounding AW1 to plot
E_plot1(1,:) = reshape(E_act_overall_L6(AW1(1)-1,AW1(2)-1,:), [1, num_steps]);
E_plot1(2,:) = reshape(E_act_overall_L6(AW1(1)-1,AW1(2)+1,:), [1, num_steps]);
E_plot1(3,:) = reshape(E_act_overall_L6(AW1(1)-1,AW1(2),:), [1, num_steps]);
E_plot1(4,:) = reshape(E_act_overall_L6(AW1(1),AW1(2)-1,:), [1, num_steps]);
E_plot1(5,:) = reshape(E_act_overall_L6(AW1(1),AW1(2),:), [1, num_steps]);
E_plot1(6,:) = reshape(E_act_overall_L6(AW1(1),AW1(2)+1,:), [1, num_steps]);
E_plot1(7,:) = reshape(E_act_overall_L6(AW1(1)+1,AW1(2),:), [1, num_steps]);
E_plot1(8,:) = reshape(E_act_overall_L6(AW1(1)+1,AW1(2)+1,:), [1, num_steps]);
E_plot1(9,:) = reshape(E_act_overall_L6(AW1(1)+1,AW1(2)-1,:), [1, num_steps]);

% Long timescales 
time_win = 0.450; % in seconds
t_marg = 0.020;
t_pre = t - floor(t_marg/dt);
t_post = t + floor(time_win/dt);
tim = (-t_marg:dt:time_win)*10^3; % in milliseconds
t_step = 0.100; % in seconds

f = figure; 
ymax = 60;
subplot(3,2,1);
plot(tim, E_plot(5,t_pre:t_post),...
     tim, E_plot(7,t_pre:t_post),...
     tim, E_plot(8,t_pre:t_post))
set(gca,'XTick',[0:t_step:time_win]*10^3,'XLim',[-t_marg time_win]*10^3,'YTick',0:ymax/2:ymax,'YLim',[0 ymax],'TickDir','out','box','off','FontSize',AXES_FONTSIZE)
% set(gca,'XTick',[0:t_step:50],'XLim',[-t_marg time_win]*10^3,'YTick',0:50:150,'YLim',[0 150],'FontSize', AXES_FONTSIZE)
c = colorbar;
c.Visible = 'off';
ylabel('A_{L4} (Hz)')
legend('D2','E2','E3','Box','off');

subplot(3,2,3);
ymax = 60;
plot(tim, E_plot1(5,t_pre:t_post),...
     tim, E_plot1(7,t_pre:t_post),...
     tim, E_plot1(8,t_pre:t_post))
set(gca,'XTick',[0:t_step:time_win]*10^3,'XLim',[-t_marg time_win]*10^3,'YTick',0:ymax/2:ymax,'YLim',[0 ymax],'TickDir','out','box','off','FontSize', AXES_FONTSIZE)
c = colorbar;
c.Visible = 'off';
xlabel('time (ms)')
ylabel('A_{L6} (Hz)')
legend('D2','E2','E3','Box','off');
set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition')) 
saveas(f,'Figures/Temp_Long_PS.pdf');

% Short timescales： 
time_win = 0.040; % in seconds
t_marg = 0.005;
t_pre = t - floor(t_marg/dt);
t_post = t + floor(time_win/dt);
tim = (-t_marg:dt:time_win)*10^3; % in milliseconds
t_step = 0.010; % in seconds

f = figure; 
subplot(3,2,1);
plot(tim, E_plot(5,t_pre:t_post),...
     tim, E_plot(7,t_pre:t_post),...
     tim, E_plot(8,t_pre:t_post))
set(gca,'XTick',[0:t_step:time_win]*10^3,'XLim',[-t_marg time_win]*10^3,'YTick',0:50:150,'YLim',[0 150],'TickDir','out','box','off','FontSize',AXES_FONTSIZE)
% set(gca,'XTick',[0:t_step:50],'XLim',[-t_marg time_win]*10^3,'YTick',0:50:150,'YLim',[0 150],'FontSize', AXES_FONTSIZE)
c = colorbar;
c.Visible = 'off';
ylabel('A_{L4} (Hz)')
legend('D2','E2','E3','Box','off');

subplot(3,2,3);
plot(tim, E_plot1(5,t_pre:t_post),...
     tim, E_plot1(7,t_pre:t_post),...
     tim, E_plot1(8,t_pre:t_post))
set(gca,'XTick',[0:t_step:time_win]*10^3,'XLim',[-t_marg time_win]*10^3,'YTick',0:50:150,'YLim',[0 150],'TickDir','out','box','off','FontSize', AXES_FONTSIZE)
c = colorbar;
c.Visible = 'off';
xlabel('time (ms)')
ylabel('A_{L6} (Hz)')
legend('D2','E2','E3','Box','off');
set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition')) 
saveas(f,'Figures/Temp_short_PS.pdf');


%% Temporal raster and psth of TCs  
clear
k = 1;
kk = 1;
load('Simulation Results/meta_data.mat'); 
load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat']);
t = Stim_Onsets(1);

% Long timescales：
n_stim_plot = 1; 
time_win = n_stim_plot*(duration+ISI-0.550)*10^3; % (ms)
t_marg = 20; % (ms)
t_pre = t - ceil(t_marg*10^-3/dt);
t_post = t + ceil(time_win*10^-3/dt);
tim = -t_marg:(dt*10^3):time_win; % (in milliseconds)
t_step = 100; % ms
t_eq = t_eq*10^3; % ms

f = figure; 
subplot(3,2,1)
hold on
plot(firings(firings(:,2)<=Ntc, 1)*10^3,...
     firings(firings(:,2)<=Ntc, 2),'b.','MarkerSize',5); % tc in AW1
plot(firings(firings(:,2)>4*Ntc & firings(:,2)<=4*Ntc+Nre, 1)*10^3,...
     firings(firings(:,2)>4*Ntc & firings(:,2)<=4*Ntc+Nre, 2)-3*Ntc,'Color','#CCCCFF',...
     'Marker','.','MarkerSize',5,'LineStyle','none'); % re in AW1
xlabel('time (ms)')
ylabel('Neuron Index');
set(gca,'XTick',t_eq+[0:t_step:time_win],'XTickLabel',[0:t_step:time_win],'XLim',t_eq+[-t_marg time_win],...
    'YTick',0:50:Ntc+Nre,'YLim',[0 Ntc+Nre],'TickDir','out','box','off','FontSize',AXES_FONTSIZE) 
legend('D2 TC','D2 RE');

ymax = 40;
subplot(3,2,3)
plot(tim,E_act_overall_tc(1,t_pre:t_post),'b',...
    'LineWidth',LineWidth,'MarkerSize',MarkerSize); % Rate of TC in AW1
set(gca,'XTick',[0:t_step:time_win],'XLim',[-t_marg time_win],'YTick',0:ymax/2:ymax,'YLim',[0 ymax],'TickDir','out','box','off','FontSize',AXES_FONTSIZE);
xlabel('time (ms)')
ylabel('A_{TC} (Hz)')

subplot(3,2,5);
area(tim,reshape(Spa_Temp(1,1,t_pre:t_post), [1 length(t_pre:t_post)]),'LineStyle','none','FaceColor','b');
set(gca,'FontSize', AXES_FONTSIZE, 'XTickLabel', [], 'YTick', [0], 'YTickLabel', {'D2'},'TickDir','out','box','off');
axis([-t_marg time_win 0 4]);
set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition')) 
saveas(f,'Figures/Temp_long_TC.pdf');

% Short timescales： 
n_stim_plot = 1; 
time_win = n_stim_plot*40; % (ms)
t_marg = 5; % (ms)
t_pre = t - ceil(t_marg*10^-3/dt);
t_post = t + ceil(time_win*10^-3/dt);
tim = -t_marg:(dt*10^3):time_win; % (in milliseconds)
t_step = 10; % ms

f = figure; 
subplot(3,2,1)
hold on
plot(firings(firings(:,2)<=Ntc, 1)*10^3,...
     firings(firings(:,2)<=Ntc, 2),'b.','MarkerSize',5); % tc in AW1
plot(firings(firings(:,2)>4*Ntc & firings(:,2)<=4*Ntc+Nre, 1)*10^3,...
     firings(firings(:,2)>4*Ntc & firings(:,2)<=4*Ntc+Nre, 2)-3*Ntc,'Color','#CCCCFF',...
     'Marker','.','MarkerSize',5,'LineStyle','none'); % re in AW1
xlabel('time (ms)')
ylabel('Neuron Index');
set(gca,'XTick',t_eq+[0:t_step:time_win],'XTickLabel',[0:t_step:time_win],'XLim',t_eq+[-t_marg time_win],...
    'YTick',0:50:Ntc+Nre,'YLim',[0 Ntc+Nre],'TickDir','out','box','off','FontSize',AXES_FONTSIZE) 
legend('D2 TC','D2 RE');

subplot(3,2,3)
plot(tim,E_act_overall_tc(1,t_pre:t_post),'b',...
    'LineWidth',LineWidth,'MarkerSize',MarkerSize); % Rate of TC in AW1
set(gca,'XTick',[0:t_step:time_win],'XLim',[-t_marg time_win],'YTick',0:50:100,'TickDir','out','box','off','FontSize',AXES_FONTSIZE);
xlabel('time (ms)')
ylabel('A_{TC} (Hz)')

subplot(3,2,5);
area(tim,reshape(Spa_Temp(1,1,t_pre:t_post), [1 length(t_pre:t_post)]),'LineStyle','none','FaceColor','b');
set(gca,'FontSize', AXES_FONTSIZE, 'XTickLabel', [], 'YTick', [0], 'YTickLabel', {'D2'},'TickDir','out','box','off');
axis([-t_marg time_win 0 4]);
set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition')) 
saveas(f,'Figures/Temp_short_TC.pdf');


%% Oddvall Low condition, cortical activity 
clear
k = 1;
kk = 1;
load('Simulation Results/meta_data.mat'); 
load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat']);

n_stim_plot = 15; 
t = Stim_Onsets(1);
t_marg = 0.500; 
time_win = n_stim_plot*(duration+ISI)+t_marg; 
t_onset = t - floor(t_marg/dt);
t_offset = t + floor(time_win/dt);
tim = -t_marg:dt:time_win; 
t_step = 5;

% --- Thalamocortical synapse dynamics ---
f = figure; 
% subplot(10,1,1) % AW1
% plot(tim, reshape(z_act_overall(4,2, t_onset:t_offset),[1,length(t_onset:t_offset)]), '-b',...
%      tim, reshape(z_act_overall(3,2, t_onset:t_offset),[1,length(t_onset:t_offset)]), '-.b', 'LineWidth',LineWidth,'MarkerSize',MarkerSize)
% set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'TickDir','out','FontSize',AXES_FONTSIZE)  
% ylabel('ThC Res.')
% legend('z_{D,2}^{D,2}','z_{D,2}^{C,2}','Box','off');
% 
% subplot(10,1,2) % PW
% plot(tim, reshape(z_act_overall2(3,2, t_onset:t_offset),[1,length(t_onset:t_offset)]), '-r',...
%      tim, reshape(z_act_overall2(4,2, t_onset:t_offset),[1,length(t_onset:t_offset)]), '-.r', 'LineWidth',LineWidth,'MarkerSize',MarkerSize)
% set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'TickDir','out','FontSize',AXES_FONTSIZE)
% ylabel('ThC Res.')
% legend('z_{C,2}^{C,2}','z_{C,2}^{D,2}','Box','off');

% ---- L4 ---- 
E_plot1 = zeros(9,num_steps); % Activities of L6 barrels surrounding AW1 in L4
E_plot1(1,:) = reshape(E_act_overall(AW1(1)-1,AW1(2)-1,:), [1, num_steps]);
E_plot1(2,:) = reshape(E_act_overall(AW1(1)-1,AW1(2)+1,:), [1, num_steps]);
E_plot1(3,:) = reshape(E_act_overall(AW1(1)-1,AW1(2),:), [1, num_steps]);
E_plot1(4,:) = reshape(E_act_overall(AW1(1),AW1(2)-1,:), [1, num_steps]);
E_plot1(5,:) = reshape(E_act_overall(AW1(1),AW1(2),:), [1, num_steps]);
E_plot1(6,:) = reshape(E_act_overall(AW1(1),AW1(2)+1,:), [1, num_steps]);
E_plot1(7,:) = reshape(E_act_overall(AW1(1)+1,AW1(2),:), [1, num_steps]);
E_plot1(8,:) = reshape(E_act_overall(AW1(1)+1,AW1(2)+1,:), [1, num_steps]);
E_plot1(9,:) = reshape(E_act_overall(AW1(1)+1,AW1(2)-1,:), [1, num_steps]);

x_plot1 = zeros(9,num_steps); % Resource of L6 barrels surrounding AW1 in L4
x_plot1(1,:) = reshape(x_act_overall(AW1(1)-1,AW1(2)-1,:), [1, num_steps]);
x_plot1(2,:) = reshape(x_act_overall(AW1(1)-1,AW1(2)+1,:), [1, num_steps]);
x_plot1(3,:) = reshape(x_act_overall(AW1(1)-1,AW1(2),:), [1, num_steps]);
x_plot1(4,:) = reshape(x_act_overall(AW1(1),AW1(2)-1,:), [1, num_steps]);
x_plot1(5,:) = reshape(x_act_overall(AW1(1),AW1(2),:), [1, num_steps]);
x_plot1(6,:) = reshape(x_act_overall(AW1(1),AW1(2)+1,:), [1, num_steps]);
x_plot1(7,:) = reshape(x_act_overall(AW1(1)+1,AW1(2),:), [1, num_steps]);
x_plot1(8,:) = reshape(x_act_overall(AW1(1)+1,AW1(2)+1,:), [1, num_steps]);
x_plot1(9,:) = reshape(x_act_overall(AW1(1)+1,AW1(2)-1,:), [1, num_steps]);

subplot(10,1,3) 
plot(tim,E_plot1(5,t_onset:t_offset),'b',tim,E_plot1(3,t_onset:t_offset),'r','LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'TickDir','out','box','off','FontSize', AXES_FONTSIZE)  
ylabel('A_{L4} (Hz)')
legend('D2','C2','Box','off')

subplot(10,1,4) 
plot(tim,x_plot1(5,t_onset:t_offset),'b',tim,x_plot1(3,t_onset:t_offset),'r','LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'TickDir','out','FontSize',AXES_FONTSIZE)  
ylabel('x_{L4}')

% ---- L6 ---- :
E_plot = zeros(9,num_steps); % Activities of L6 barrels surrounding AW1 in L6
E_plot(1,:) = reshape(E_act_overall_L6(AW1(1)-1,AW1(2)-1,:), [1, num_steps]);
E_plot(2,:) = reshape(E_act_overall_L6(AW1(1)-1,AW1(2)+1,:), [1, num_steps]);
E_plot(3,:) = reshape(E_act_overall_L6(AW1(1)-1,AW1(2),:), [1, num_steps]);
E_plot(4,:) = reshape(E_act_overall_L6(AW1(1),AW1(2)-1,:), [1, num_steps]);
E_plot(5,:) = reshape(E_act_overall_L6(AW1(1),AW1(2),:), [1, num_steps]);
E_plot(6,:) = reshape(E_act_overall_L6(AW1(1),AW1(2)+1,:), [1, num_steps]);
E_plot(7,:) = reshape(E_act_overall_L6(AW1(1)+1,AW1(2),:), [1, num_steps]);
E_plot(8,:) = reshape(E_act_overall_L6(AW1(1)+1,AW1(2)+1,:), [1, num_steps]);
E_plot(9,:) = reshape(E_act_overall_L6(AW1(1)+1,AW1(2)-1,:), [1, num_steps]);

x_plot = zeros(9,num_steps); % Resource of L6 barrels surrounding AW1 in L6
x_plot(1,:) = reshape(x_act_overall_L6(AW1(1)-1,AW1(2)-1,:), [1, num_steps]);
x_plot(2,:) = reshape(x_act_overall_L6(AW1(1)-1,AW1(2)+1,:), [1, num_steps]);
x_plot(3,:) = reshape(x_act_overall_L6(AW1(1)-1,AW1(2),:), [1, num_steps]);
x_plot(4,:) = reshape(x_act_overall_L6(AW1(1),AW1(2)-1,:), [1, num_steps]);
x_plot(5,:) = reshape(x_act_overall_L6(AW1(1),AW1(2),:), [1, num_steps]);
x_plot(6,:) = reshape(x_act_overall_L6(AW1(1),AW1(2)+1,:), [1, num_steps]);
x_plot(7,:) = reshape(x_act_overall_L6(AW1(1)+1,AW1(2),:), [1, num_steps]);
x_plot(8,:) = reshape(x_act_overall_L6(AW1(1)+1,AW1(2)+1,:), [1, num_steps]);
x_plot(9,:) = reshape(x_act_overall_L6(AW1(1)+1,AW1(2)-1,:), [1, num_steps]);

subplot(10,1,5) 
plot(tim,E_plot(5,t_onset:t_offset),'b',tim,E_plot(3,t_onset:t_offset),'r','LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'TickDir','out','box','off','FontSize',AXES_FONTSIZE)  
ylabel('A_{L6} (Hz)')
legend('D2','C2','Box','off')

subplot(10,1,6) 
plot(tim,x_plot(5,t_onset:t_offset),'b',tim,x_plot(3,t_onset:t_offset),'r','LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'TickDir','out','FontSize',AXES_FONTSIZE)  
ylabel('x_{L6}')

% ---- stimulus sequence ---- 
subplot(12,1,10); 
plot(tim, reshape(Spa_Temp(1,1,t_onset:t_offset), [1 length(t_onset:t_offset)])+1, 'b',...
     tim, reshape(Spa_Temp(1,2,t_onset:t_offset), [1 length(t_onset:t_offset)])+2.5, 'r', 'LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca,'FontSize', AXES_FONTSIZE, 'XTickLabel', [0:t_step:time_win], 'YTick', [1 2.5], 'YTickLabel', {'D2','C2'},'TickDir','out','box','off');
axis([-t_marg time_win 0 4]);
xlabel('time (s)')

set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition')) 
saveas(f,'Figures/Odd_Low_corti.pdf');


%% Oddvall Low condition, thalamic activity 
clear
k = 1;
kk = 1;
load('Simulation Results/meta_data.mat'); 
load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat']);
t = Stim_Onsets(1);

n_stim_plot = 15; 
t_marg = 0.500; % s
time_win = n_stim_plot*(duration+ISI)+t_marg; % s
t_pre = t - floor(t_marg/dt);
t_post = t + floor(time_win/dt);
tim = -t_marg:dt:time_win; 
t_step = 5; % s

f = figure; 
subplot(2,1,1)
hold on
plot(firings(firings(:,2)<=Ntc, 1),...
     firings(firings(:,2)<=Ntc, 2), 'b.','MarkerSize',5); % tc in AW1

plot(firings(firings(:,2)>4*Ntc & firings(:,2)<=4*Ntc+Nre, 1),...
     firings(firings(:,2)>4*Ntc & firings(:,2)<=4*Ntc+Nre, 2)-3*Ntc,'Color','#CCCCFF',...
     'Marker','.','MarkerSize',5,'LineStyle','none'); % re in AW1

plot(firings(firings(:,2)>Ntc & firings(:,2)<=2*Ntc, 1),...
     firings(firings(:,2)>Ntc & firings(:,2)<=2*Ntc, 2)+Nre, 'r.','MarkerSize',5); % tc in PW

plot(firings(firings(:,2)>4*Ntc+Nre & firings(:,2)<=4*Ntc+2*Nre, 1),...
     firings(firings(:,2)>4*Ntc+Nre & firings(:,2)<=4*Ntc+2*Nre, 2)-2*Ntc,'Color','#FFCCCC',...
     'Marker','.','MarkerSize',5,'LineStyle','none'); % re in PW
% axis square
ylabel('Neuron Index');
set(gca,'XTick',t_eq+[0:t_step:time_win],'XTickLabel',[0:t_step:time_win],'YTick',0:100:400,'XLim',t_eq+[-t_marg time_win],'YLim',[0 2*(Ntc+Nre)],...
    'TickDir','out','box','off','FontSize',AXES_FONTSIZE) 
legend('D2 TC','D2 RE','C2 TC','C2 RE');

subplot(4,1,3)
plot(tim, E_act_overall_tc(1,t_pre:t_post),'b',tim, E_act_overall_tc(2,t_pre:t_post),'r',...
    'LineWidth',LineWidth,'MarkerSize',MarkerSize); % Rate of TC in AW1/PW
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'TickDir','out','FontSize', AXES_FONTSIZE,'box','off');
ylabel('A_{TC} (Hz)')

subplot(12,1,10);
plot(tim, reshape(Spa_Temp(1,1,t_pre:t_post), [1 length(t_pre:t_post)])+1, 'b',...
     tim, reshape(Spa_Temp(1,2,t_pre:t_post), [1 length(t_pre:t_post)])+2.5, 'r', 'LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca,'FontSize', AXES_FONTSIZE, 'XTickLabel', [0:t_step:time_win], 'YTick', [1 2.5], 'YTickLabel', {'D2','C2'},'TickDir','out','box','off');
axis([-t_marg time_win 0 4]);
xlabel('time (s)')

set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition')) 
saveas(f,'Figures/Odd_Low_thalam.pdf');


%% Average number of TC neurons that burst in reponse to a deviant/standard stimulus in a trial
clear
k = 1;
kk = 1;
load('Simulation Results/meta_data.mat'); 
load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat']);

t_early = 0.040; % in seconds
fir_dev_count = zeros(n_stim*prob, 1);
fir_std_count = zeros(n_stim*(1-prob), 1);
n_dev= 0; 
n_std = 0; 

for ns = 1:n_stim
    if Oddball(ns,:) == PW
        % Seleting fired TC cells in PW (deviant whisker) during late phase of each trial
        firings_dev_trial = firings(:,2)>Ntc & firings(:,2)<=2*Ntc & firings(:,1)>=(ns+t_early) & firings(:,1)<=(ns+1);

        n_dev = n_dev + 1;
        fir_dev_count(n_dev) = length(unique(firings(firings_dev_trial, 2)));
    else
        % Seleting fired TC cells in AW1 (standard whisker) during late phase of each trial
        firings_std_trial = firings(:,2)<=Ntc & firings(:,1)>=(ns+t_early) & firings(:,1)<=(ns+1); 

        n_std = n_std + 1;
        fir_std_count(n_std) = length(unique(firings(firings_std_trial, 2)));     
    end
end

fir_dev_mean = mean(fir_dev_count)
fir_std_mean = mean(fir_std_count)


%% Percentage of deviant/standard trials in which a given TC neuron bursts 
clear
k = 1;
kk = 1;
load('Simulation Results/meta_data.mat'); 
load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat']);

t_early = 0.040; % in seconds
fir_dev_count = zeros(n_stim*prob, 1);
fir_std_count = zeros(n_stim*(1-prob), 1);
n_dev= 0; 
n_std = 0; 

neu_PW_ind = Ntc+10;  % TC cell index in PW
neu_AW1_ind = 10;  % TC cell index in AW1

for ns = 1:n_stim
    if Oddball(ns,:) == PW
        % Spikes emitted by a given TC cell in PW (deviant whisker) during late phase of each trial
        firings_dev_trial = firings(:,2)==neu_PW_ind & firings(:,1)>=(ns+t_early) & firings(:,1)<=(ns+1);

        n_dev = n_dev + 1;
        fir_dev_count(n_dev) = length(unique(firings(firings_dev_trial, 2)));
    else
        % Spikes emmitted by a given TC cells in AW1 (standard whisker) during late phase of each trial
        firings_std_trial = firings(:,2)==neu_AW1_ind & firings(:,1)>=(ns+t_early) & firings(:,1)<=(ns+1); 

        n_std = n_std + 1;
        fir_std_count(n_std) = length(unique(firings(firings_std_trial, 2)));     
    end
end

fir_dev_mean = mean(fir_dev_count)
fir_std_mean = mean(fir_std_count)


%% Many-standard condition
clear
k = 1;
kk = 4;
load('Simulation Results/meta_data.mat'); 
load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat']);

n_stim_plot = 10;
stim_ind = 98;
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
subplot(11,1,1);
plot(tim, reshape(Spa_Temp_odd(1,1,t_onset:t_offset), [1 length(t_onset:t_offset)])+1, 'b',...
     tim, reshape(Spa_Temp_odd(1,2,t_onset:t_offset), [1 length(t_onset:t_offset)])+2.5, 'r', LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca,'FontSize', AXES_FONTSIZE, 'XTickLabel', [], 'YTick', [1 2.5], 'YTickLabel', {'D2','C2'},'TickDir','out','box','off');
axis([-t_marg time_win 0 7]);
% axis off

subplot(11,1,2);
plot(tim, reshape(Spa_Temp(1,1,t_onset:t_offset), [1 length(t_onset:t_offset)])+1, 'b',...
     tim, reshape(Spa_Temp(1,3,t_onset:t_offset), [1 length(t_onset:t_offset)])+1+1.5, 'b',...
     tim, reshape(Spa_Temp(1,4,t_onset:t_offset), [1 length(t_onset:t_offset)])+1+2*1.5, 'b',... 
     tim, reshape(Spa_Temp(1,2,t_onset:t_offset), [1 length(t_onset:t_offset)])+1+3*1.5, 'g',LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca,'FontSize', AXES_FONTSIZE, 'XTickLabel', [], 'YTick', [1 1+1.5 1+2*1.5 1+3*1.5], 'YTickLabel', {'D2','D1','D3','C2'},'TickDir','out','box','off');
axis([-t_marg time_win 0 7]);
% axis off

% ---- Odd vs MS deviant in L4 ---- 
subplot(11,1,3) % Firing rate 
plot(tim,E_plot_L4_odd(t_onset:t_offset),'r',tim,E_plot_L4_ms(t_onset:t_offset),'g','LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'TickDir','out','box','off','FontSize', AXES_FONTSIZE)  
ylabel('A_{L4} (Hz)')
legend('odd C2','ms C2','Box','off')

subplot(11,1,4) % Resource
plot(tim,x_plot_L4_odd(t_onset:t_offset),'r',tim,x_plot_L4_ms(t_onset:t_offset),'g','LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'TickDir','out','FontSize', AXES_FONTSIZE)  
ylabel('x_{L4}')

% ---- L4-to-L6 input in Odd vs MS deviant ---- 
subplot(11,1,5) % Firing rate 
plot(tim,E_plot_L46_odd(t_onset:t_offset),'r',tim,E_plot_L46_ms(t_onset:t_offset),'g','LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'TickDir','out','box','off','FontSize', AXES_FONTSIZE)  
ylabel('A_{L46} (Hz/s)')

% ---- Odd vs MS deviant in L6 ---- 
subplot(11,1,6) % Firing rate 
plot(tim,E_plot_L6_odd(t_onset:t_offset),'r',tim,E_plot_L6_ms(t_onset:t_offset),'g','LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'TickDir','out','box','off','FontSize', AXES_FONTSIZE)  
ylabel('A_{L6} (Hz)')

subplot(11,1,7) % Resource
plot(tim,x_plot_L6_odd(t_onset:t_offset),'r',tim,x_plot_L6_ms(t_onset:t_offset),'g','LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'TickDir','out','FontSize', AXES_FONTSIZE)  
ylabel('x_{L6}')

% ---- L6 Resource in ODD vs MS standard ---- 
subplot(11,1,8) 
plot(tim,x_plot_L6_odd_D2(t_onset:t_offset),'b',tim,x_plot_L6_ms_D2(t_onset:t_offset),'b--',...
    'LineWidth',LineWidth,'MarkerSize',MarkerSize)
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'TickDir','out','FontSize', AXES_FONTSIZE) 
legend('odd D2','ms D2','Box','off')
ylabel('x_{L6}')

% ---- Odd deviant in TCs ---- 
subplot(11,1,10) 
% plot(tim, E_act_overall_tc_odd(1,t_onset:t_offset),'b',tim, E_act_overall_tc_odd(2,t_onset:t_offset),'r',...
%     'LineWidth',LineWidth,'MarkerSize',MarkerSize); % Rate of TC in AW1/PW
plot(tim, E_act_overall_tc_odd(2,t_onset:t_offset),'r','LineWidth',LineWidth,'MarkerSize',MarkerSize); % Rate of TC in oddball PW
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[],'XLim',[-t_marg time_win],'TickDir','out','box','off','FontSize', AXES_FONTSIZE)  
ylabel('A_{TC} (Hz)')
% legend('odd D2','odd C2','Box','off')

% ---- MS deviant in TCs ---- 
subplot(11,1,11) 
% plot(tim, E_act_overall_tc(1,t_onset:t_offset),'b',tim, E_act_overall_tc(2,t_onset:t_offset),'g',...
%     'LineWidth',LineWidth,'MarkerSize',MarkerSize); % Rate of TC in AW1/PW
plot(tim, E_act_overall_tc(2,t_onset:t_offset),'g','LineWidth',LineWidth,'MarkerSize',MarkerSize); % Rate of TC in many-standards PW
set(gca,'XTick',[0:t_step:time_win],'XTickLabel',[stim_ind-1:t_step:stim_ind-1+time_win],'XLim',[-t_marg time_win],'TickDir','out','box','off','FontSize', AXES_FONTSIZE)  
xlabel('time (s)')
ylabel('A_{TC} (Hz)')
% legend('ms D2','ms C2','Box','off')

set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition'))
saveas(f,'Figures/MS.pdf');


%% Firing patterns of TC and RE cells
AXES_FONTSIZE = 10;
LineWidth = 1;
MarkerSize = 10;
COLORBAR_FONTSIZE = 10;
f = figure;

% Initial bursting thalamic reticular (RE) Cells 
subplot(3,3,1) 
a=0.02; b=0.2; c=-55; d=4;
V=-70; u=b*V;
VV=[];  uu=[];
tau = 0.25; tspan = 0:tau:190;
T1=tspan(end)/10;
for t=tspan 
    if (t>T1)
        I=10;
    else
        I=0;
    end
    V = V + tau*(0.04*V^2+5*V+140-u+I);
    u = u + tau*a*(b*V-u);
    if V > 30
        VV(end+1)=30;
        V = c;
        u = u + d;
    else
        VV(end+1)=V;
    end
    uu(end+1)=u;
end
plot(tspan,VV,[0 T1 T1 max(tspan)],-90+[0 0 10 10],'LineWidth',LineWidth,'MarkerSize',MarkerSize);
set(gca,'FontSize', AXES_FONTSIZE) 
axis([0 max(tspan) -100 30])
hold on
plot([140 190],[-90 -90],'k',[tspan(end) tspan(end)],[-60 -10],'k','LineWidth',2.5)
axis off;
axis square
title('Bursting RE Neuron ');

% Rebound Burst Thalamocortical (TC) Relay Cells 
subplot(3,3,4) 
a=0.005; b=0.26; c=-52; d=2;
V=-62.5; u=b*V;
VV=[];  uu=[];
tau = 0.25; tspan = 0:tau:190;
T1=tspan(end)/10;
for t=tspan
    if t > T1 && t < T1 + 120 
        I=-10;
    else
        I=0;
    end
    V = V + tau*(0.04*V^2+5*V+140-u+I);
    u = u + tau*a*(b*V-u);
    if V > 30
        VV(end+1)=30;
        V = c;
        u = u + d;
    else
        VV(end+1)=V;
    end
    uu(end+1)=u;
end
plot(tspan,VV,[0 T1 T1 T1 + 120 T1 + 120 max(tspan)],-90+[0 0 -10 -10 0 0],'LineWidth',LineWidth,'MarkerSize',MarkerSize);
set(gca,'FontSize', AXES_FONTSIZE) 
axis([0 max(tspan) -100 30])
axis off;
axis square
title('Rebound Burst TC Neuron');
set(gcf,'PaperUnits','normalized','PaperPosition',[0.05 0.1 1 0.7]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition')) 
saveas(f,'Figures/Fir_Pattern.pdf');


%% PS and resource depletion in a single population with depressing synapses

% L4 network parameters:
tau_h = 0.001; % in seconds
tau_rec = 0.500; 
J = 2.2;
U = 0.5;
Gain_thres = 5; 
Gain_slope = 1; 

% Time steps (in seconds)
dt = 0.0001; 
t0 = 0.010;
t1 = 0.010;
t2 = 0.030;
num_steps0 = ceil(t0/dt); 
num_steps1 = ceil(t1/dt); 
num_steps2 = ceil(t2/dt); 
num_steps = num_steps0 + num_steps1 + num_steps2 + 1;

% Initial values:
x = 1:-0.2:0.4; 
x = x';
M = length(x);
E = zeros(M,1); % Acitivity of the population (in Hz) 
h = zeros(M,1); % Sypnatic current

% Tracking variables
E_seq = zeros(M,num_steps); 
x_seq = zeros(M,num_steps); 
I_ext_seq = zeros(M,num_steps); 


for i = 1:num_steps
    
    if i>num_steps0 && i<num_steps0+num_steps1
        I_ext = 50*ones(M,1);
    else
        I_ext = zeros(M,1);
    end

    % Implementing the non-linearity:
    E = h;
    E = E - Gain_thres;
    E = E * Gain_slope;
    E(E <  0) = 0;
    
    % The variables' dynamics:
    h = h + (dt/tau_h)*(-h + J*U*x.*E + I_ext);
    x = x + dt*((1-x)/tau_rec - U*x.*E);
    
    E_seq(:,i) = E;
    x_seq(:,i) = x;
    I_ext_seq(:,i) = I_ext;
end


% Plotting
tim = (-t0:dt:t1+t2)*1000; % (in ms)

f = figure;
subplot(3,3,1)
plot(tim,E_seq,'LineWidth',LineWidth,'MarkerSize',MarkerSize);
set(gca,'XTick',[0:0.010:t1+t2]*1000,'XLim',[-t0 t1+t2]*1000,'TickDir','out','box','off','FontSize',AXES_FONTSIZE);
axis square
box off
xlabel('time (ms)');
ylabel('Population activity (Hz)');
legend('Init Res = 1','0.8','0.6','0.4')
legend('boxoff')

subplot(3,3,4)
plot(tim,x_seq,'LineWidth',LineWidth,'MarkerSize',MarkerSize);
set(gca,'XTick',[0:0.010:t1+t2]*1000,'YTick',0.2:0.2:1,'XLim',[-t0 t1+t2]*1000,'YLim',[0.2 1],'TickDir','out','box','off','FontSize',AXES_FONTSIZE);
axis square
xlabel('time (ms)');
ylabel('Synaptic resource');
set(gcf,'PaperUnits','normalized','PaperPosition',[0.05 0.1 1 0.7]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition')) 
saveas(f,'Figures/PS_resource.pdf');


%% Barplot for oddball and many-standards protocols
AXES_FONTSIZE = 10;
f = figure;
subplot(11,6,1)
y = [0.25 0.75 0 0];
b = bar(y, 'FaceColor', 'flat', 'EdgeColor','flat');
b.CData(1,:) = [1 0 0];
b.CData(2,:) = [0 0 1];
xticklabels({'C2','D2'})
xtickangle(60)
box off
title('Oddball deviant');
set(gca,'FontSize',AXES_FONTSIZE,'YLim',[0 1]);

subplot(11,6,3)
y = [0.25 0.25 0.25 0.25];
b = bar(y, 'FaceColor', 'flat', 'EdgeColor','flat');
b.CData(1,:) = [0 1 0];
b.CData(2,:) = [0 0 1];
b.CData(3,:) = [0 0 1];
b.CData(4,:) = [0 0 1];
xticklabels({'C2','D2','D1','D3'})
xtickangle(60)
box off
title('Many-standards deviant');
set(gca,'FontSize',AXES_FONTSIZE,'YLim',[0 1]);


set(gcf,'PaperUnits','normalized','PaperPosition',[0.05 0.1 1 0.7]) 
set(gcf,'Units','normalized','position',get(gcf,'PaperPosition')) 
saveas(f,'Figures/Bar_OddvsMS.pdf');











