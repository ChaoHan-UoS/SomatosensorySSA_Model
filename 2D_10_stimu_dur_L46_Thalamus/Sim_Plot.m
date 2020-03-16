%% 1. Plot parameter dependence on SI and CSI 
clear;
load('Simulation Results/meta_data.mat'); 
SI = zeros(num_trials,length(Par_Arr)); % Index quantifying SSA
CSI = zeros(num_trials, length(Par_Arr)); % Index quantifying true deviance detection
time_win = 0.020; % (in seconds) Spiking response within a window from 0 to 20ms after stimulation onset is used
for k = 1:length(Par_Arr)
    for kkk = 1:num_trials
        n_odddev = 0; % Number of oddball deviant p=0.25
        n_oddstd = 0; % Number of oddball standard p=0.75
        n_msdev = 0; % Number of many-standard deviant p=0.25
        E_sum_odddev = zeros(1, 1, time_win/dt);
        E_sum_oddstd = E_sum_odddev;
        E_sum_msdev = E_sum_odddev;
        Spcount_odddev = 0;
        Spcount_oddstd = 0;
        Spcount_msdev = 0;
        for kk = [1 4]
            load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '_Tr' num2str(kkk) '.mat']);
            for ns = 1:n_stim
                if kk == 1
                    if Oddball(ns, :) == PW
                        E_sum_odddev = E_sum_odddev + E_act_overall(PW(1), PW(2), Stim_Onsets(ns):Stim_Onsets(ns)+time_win/dt-1);
                        n_odddev = n_odddev + 1;
                    else
                        E_sum_oddstd = E_sum_oddstd + E_act_overall(AW1(1), AW1(2), Stim_Onsets(ns):Stim_Onsets(ns)+time_win/dt-1);
                        n_oddstd = n_oddstd + 1;
                    end
                else
                    if Oddball(ns, :) == PW
                        E_sum_msdev = E_sum_msdev + E_act_overall(PW(1), PW(2), Stim_Onsets(ns):Stim_Onsets(ns)+time_win/dt-1);
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
SI_mean = mean(SI);
SI_std = std(SI);
CSI_mean = mean(CSI);
CSI_std = std(CSI);

% f = figure('Visible','off');
figure
errorbar(Par_Arr, SI_mean, SI_std, '-or');
hold on;
errorbar(Par_Arr, CSI_mean, CSI_std, '-og');
plot(Par_Arr, zeros(size(Par_Arr)), ':k');
% xlabel('\tau_{rec}(s)');
% xlim([0.1 1.5]);
% xticks([0.1:0.1:1.5]);
% ylim([-0.1 0.5]);
legend('SI','CSI','Location','northwest');
% saveas(f,'Figure/CSIvsTau_rec.pdf');

 
%% 2. SSA and TDD
% ---- L4 ---- :
clear;
load('Simulation Results/meta_data.mat'); 
CSI = zeros(size(Par_Arr)); % Index quantifying true deviance detection
SI = zeros(size(Par_Arr)); % Index quantifying SSA
time_win = 0.500; % (in seconds) Spiking response within a window from 0 to time_win s after the onset of stimulation 
for k = 1:length(Par_Arr)
    n_odddev = 0; % Number of oddball deviant p=0.25
    n_oddstd = 0; % Number of oddball standard p=0.75
    n_msdev = 0; % Number of many-standard deviant p=0.25
    E_sum_odddev = zeros(1, 1, time_win/dt);
    E_sum_oddstd = E_sum_odddev;
    E_sum_msdev = E_sum_odddev;
    Spcount_odddev = 0;
    Spcount_oddstd = 0;
    Spcount_msdev = 0;
    for kk = [1 4]
        load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat']);
        for ns = 1:n_stim  
            if kk == 1
                if Oddball(ns,:) == PW
                    E_sum_odddev = E_sum_odddev + E_act_overall(PW(1), PW(2), Stim_Onsets(ns):Stim_Onsets(ns)+time_win/dt-1);
                    n_odddev = n_odddev + 1;
                else
                    E_sum_oddstd = E_sum_oddstd + E_act_overall(AW1(1), AW1(2), Stim_Onsets(ns):Stim_Onsets(ns)+time_win/dt-1);
                    n_oddstd = n_oddstd + 1;
                end
            else
                if Oddball(ns,:) == PW
                    E_sum_msdev = E_sum_msdev + E_act_overall(PW(1), PW(2), Stim_Onsets(ns):Stim_Onsets(ns)+time_win/dt-1);
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
    SI(k) = (Spcount_odddev - Spcount_oddstd)/(Spcount_odddev + Spcount_oddstd);
    CSI(k) = (Spcount_odddev - Spcount_msdev)/(Spcount_odddev + Spcount_msdev);
end
disp(['SI = ' num2str(SI) ', CSI = ' num2str(CSI)]);

figure
% f = figure('Visible','off');
tim = dt:dt:time_win; % (in seconds)
tim = tim*1000; % (in milliseconds)
plot(tim, reshape(E_sum_odddev, [1, time_win/dt]), '-r');
hold on;
plot(tim, reshape(E_sum_oddstd, [1, time_win/dt]), '-b');
plot(tim, reshape(E_sum_msdev, [1, time_win/dt]), '-g');
xlabel('time(ms)');
ylabel('E(Spikes/s)');
legend('odd dev','odd std','ms dev');
legend('boxoff');
title(['L4, SI = ' num2str(SI) ', CSI = ' num2str(CSI) ', \tau_{rec} = 0.7 s']);
% saveas(f,'Figure/SSA&TDD.pdf');

% ---- L6 ---- :
clear;
load('Simulation Results/meta_data.mat'); 
CSI = zeros(size(Par_Arr)); % Index quantifying true deviance detection
SI = zeros(size(Par_Arr)); % Index quantifying SSA
time_win = 0.500;  % (in seconds) Spiking response within a window from 0 to time_win s after the onset of stimulation 
for k = 1:length(Par_Arr)
    n_odddev = 0; % Number of oddball deviant p=0.25
    n_oddstd = 0; % Number of oddball standard p=0.75
    n_msdev = 0; % Number of many-standard deviant p=0.25
    E_sum_odddev = zeros(1, 1, time_win/dt);
    E_sum_oddstd = E_sum_odddev;
    E_sum_msdev = E_sum_odddev;
    Spcount_odddev = 0;
    Spcount_oddstd = 0;
    Spcount_msdev = 0;
    for kk = [1 4]
        load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat']);
        for ns = 1:n_stim  
            if kk == 1
                if Oddball(ns,:) == PW
                    E_sum_odddev = E_sum_odddev + E_act_overall_L6(PW(1), PW(2), Stim_Onsets(ns):Stim_Onsets(ns)+time_win/dt-1);
                    n_odddev = n_odddev + 1;
                else
                    E_sum_oddstd = E_sum_oddstd + E_act_overall_L6(AW1(1), AW1(2), Stim_Onsets(ns):Stim_Onsets(ns)+time_win/dt-1);
                    n_oddstd = n_oddstd + 1;
                end
            else
                if Oddball(ns,:) == PW
                    E_sum_msdev = E_sum_msdev + E_act_overall_L6(PW(1), PW(2), Stim_Onsets(ns):Stim_Onsets(ns)+time_win/dt-1);
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
    SI(k) = (Spcount_odddev - Spcount_oddstd)/(Spcount_odddev + Spcount_oddstd);
    CSI(k) = (Spcount_odddev - Spcount_msdev)/(Spcount_odddev + Spcount_msdev);
end
disp(['SI = ' num2str(SI) ', CSI = ' num2str(CSI)]);

figure
% f = figure('Visible','off');
tim = dt:dt:time_win; % (in seconds)
tim = tim*1000; % (in milliseconds)
plot(tim, reshape(E_sum_odddev, [1, time_win/dt]), '-r');
hold on;
plot(tim, reshape(E_sum_oddstd, [1, time_win/dt]), '-b');
plot(tim, reshape(E_sum_msdev, [1, time_win/dt]), '-g');
xlabel('time(ms)');
ylabel('E(Spikes/s)');
legend('odd dev','odd std','ms dev');
legend('boxoff');
title(['L6, SI = ' num2str(SI) ', CSI = ' num2str(CSI) ', \tau_{rec} = 0.7 s']);
% saveas(f,'Figure/SSA&TDD.pdf');

% ---- thalamus ---- :
clear;
load('Simulation Results/meta_data.mat'); 
CSI = zeros(size(Par_Arr)); % Index quantifying true deviance detection
SI = zeros(size(Par_Arr)); % Index quantifying SSA
time_win = 0.500;  % (in seconds) Spiking response within a window from 0 to time_win s after the onset of stimulation 
for k = 1:length(Par_Arr)
    n_odddev = 0; % Number of oddball deviant p=0.25
    n_oddstd = 0; % Number of oddball standard p=0.75
    n_msdev = 0; % Number of many-standard deviant p=0.25
    E_sum_odddev = zeros(1, 1, time_win/dt);
    E_sum_oddstd = E_sum_odddev;
    E_sum_msdev = E_sum_odddev;
    Spcount_odddev = 0;
    Spcount_oddstd = 0;
    Spcount_msdev = 0;
    for kk = [1 4]
        load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat']);
        for ns = 1:n_stim  
            if kk == 1
                if Oddball(ns,:) == PW
                    E_sum_odddev = E_sum_odddev + E_act_overall_Th(PW(1), PW(2), Stim_Onsets(ns):Stim_Onsets(ns)+time_win/dt-1);
                    n_odddev = n_odddev + 1;
                else
                    E_sum_oddstd = E_sum_oddstd + E_act_overall_Th(AW1(1), AW1(2), Stim_Onsets(ns):Stim_Onsets(ns)+time_win/dt-1);
                    n_oddstd = n_oddstd + 1;
                end
            else
                if Oddball(ns,:) == PW
                    E_sum_msdev = E_sum_msdev + E_act_overall_Th(PW(1), PW(2), Stim_Onsets(ns):Stim_Onsets(ns)+time_win/dt-1);
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
    SI(k) = (Spcount_odddev - Spcount_oddstd)/(Spcount_odddev + Spcount_oddstd);
    CSI(k) = (Spcount_odddev - Spcount_msdev)/(Spcount_odddev + Spcount_msdev);
end
disp(['SI = ' num2str(SI) ', CSI = ' num2str(CSI)]);

figure
% f = figure('Visible','off');
tim = dt:dt:time_win; % (in seconds)
tim = tim*1000; % (in milliseconds)
plot(tim, reshape(E_sum_odddev, [1, time_win/dt]), '-r');
hold on;
plot(tim, reshape(E_sum_oddstd, [1, time_win/dt]), '-b');
plot(tim, reshape(E_sum_msdev, [1, time_win/dt]), '-g');
xlabel('time(ms)');
ylabel('E(Spikes/s)');
legend('odd dev','odd std','ms dev');
legend('boxoff');
title(['thalamus, SI = ' num2str(SI) ', CSI = ' num2str(CSI) ', \tau_{rec} = 0.7 s']);
% saveas(f,'Figure/SSA&TDD.pdf');

%% 3. Propogation of PS
clear
k = 1;
kk = 1;
load('Simulation Results/meta_data.mat'); 
load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat']);
t = Stim_Onsets(1);
time_win = 0.300; 

% ---- L4 ---- :
figure
% f = figure('Visible','off');
for i = 1:11
    subplot(4,3,i)
    image([1 4], [1 5], E_act_overall(:,:,t+(i-1)*0.005/dt));
    set(gca, 'XTick', [1:4], 'YTick', [1:5]);
    title(['t = ',num2str((i-1)*5),'ms']);
    colormap hot
    axis equal
    axis tight
    c = colorbar;
    c.Label.String = 'E (Spikes/s)';
    c.Limits = [0 90];
    c.Ticks = [0 30 60 90];
    c.TickLabels = [0 30 60 90];
    if i > 1
        c.Visible = 'off';
    end
end
% saveas(f,'Figure/Prop_PS1.pdf');

figure
% f = figure('Visible','off');
subplot(2,1,1);
imagesc([0 time_win*10^3], [1 9], E_plot(:,t:t+time_win/dt));
set(gca,'YDir','normal','YTick',1:9)
xlabel('time(ms)')
ylabel('Barrel')  
colormap jet
c = colorbar;
c.Label.String = 'E (Spikes/s)';
title('L4');
subplot(2,1,2);
tim = (0:dt:time_win)*10^3; % (in milliseconds)
hold on
for i = 1:9
    plot(tim, E_plot(i,t:t+time_win/dt))
end
c = colorbar;
c.Visible = 'off';
xlabel('time(ms)')
ylabel('E (Spikes/s)')
% saveas(f,'Figure/Prop_PS2.pdf');

% ---- L6 ---- :
% Activities of L4 barrels surrounding PW to plot
E_plot = zeros(9,num_steps);
E_plot(1,:) = reshape(E_act_overall_L6(PW(1)-1,PW(2)-1,:), [1, num_steps]);
E_plot(2,:) = reshape(E_act_overall_L6(PW(1)-1,PW(2)+1,:), [1, num_steps]);
E_plot(3,:) = reshape(E_act_overall_L6(PW(1)-1,PW(2),:), [1, num_steps]);
E_plot(4,:) = reshape(E_act_overall_L6(PW(1),PW(2)-1,:), [1, num_steps]);
E_plot(5,:) = reshape(E_act_overall_L6(PW(1),PW(2),:), [1, num_steps]);
E_plot(6,:) = reshape(E_act_overall_L6(PW(1),PW(2)+1,:), [1, num_steps]);
E_plot(7,:) = reshape(E_act_overall_L6(PW(1)+1,PW(2),:), [1, num_steps]);
E_plot(8,:) = reshape(E_act_overall_L6(PW(1)+1,PW(2)+1,:), [1, num_steps]);
E_plot(9,:) = reshape(E_act_overall_L6(PW(1)+1,PW(2)-1,:), [1, num_steps]);

figure
% f = figure('Visible','off');
for i = 1:11
    subplot(4,3,i)
    image([1 4], [1 5], E_act_overall_L6(:,:,t+(i-1)*0.005/dt));
    set(gca, 'XTick', [1:4], 'YTick', [1:5]);
    title(['t = ',num2str((i-1)*5),'ms']);
    colormap hot
    axis equal
    axis tight
    c = colorbar;
    c.Label.String = 'E (Spikes/s)';
    c.Limits = [0 90];
    c.Ticks = [0 30 60 90];
    c.TickLabels = [0 30 60 90];
    if i > 1
        c.Visible = 'off';
    end
end
% saveas(f,'Figure/Prop_PS1.pdf');

figure
% f = figure('Visible','off');
subplot(2,1,1);
imagesc([0 time_win*10^3], [1 9], E_plot(:,t:t+time_win/dt));
set(gca,'YDir','normal','YTick',1:9)
xlabel('time(ms)')
ylabel('Barrel')  
colormap jet
c = colorbar;
c.Label.String = 'E (Spikes/s)';
title('L6');
subplot(2,1,2);
tim = (0:dt:time_win)*10^3; % (in milliseconds)
hold on
for i = 1:9
    plot(tim, E_plot(i,t:t+time_win/dt))
end
c = colorbar;
c.Visible = 'off';
xlabel('time(ms)')
ylabel('E (Spikes/s)')
% saveas(f,'Figure/Prop_PS2.pdf');

% ---- thalamus ---- :
% Activities of L4 barrels surrounding PW to plot
E_plot = zeros(9,num_steps);
E_plot(1,:) = reshape(E_act_overall_Th(PW(1)-1,PW(2)-1,:), [1, num_steps]);
E_plot(2,:) = reshape(E_act_overall_Th(PW(1)-1,PW(2)+1,:), [1, num_steps]);
E_plot(3,:) = reshape(E_act_overall_Th(PW(1)-1,PW(2),:), [1, num_steps]);
E_plot(4,:) = reshape(E_act_overall_Th(PW(1),PW(2)-1,:), [1, num_steps]);
E_plot(5,:) = reshape(E_act_overall_Th(PW(1),PW(2),:), [1, num_steps]);
E_plot(6,:) = reshape(E_act_overall_Th(PW(1),PW(2)+1,:), [1, num_steps]);
E_plot(7,:) = reshape(E_act_overall_Th(PW(1)+1,PW(2),:), [1, num_steps]);
E_plot(8,:) = reshape(E_act_overall_Th(PW(1)+1,PW(2)+1,:), [1, num_steps]);
E_plot(9,:) = reshape(E_act_overall_Th(PW(1)+1,PW(2)-1,:), [1, num_steps]);

figure
% f = figure('Visible','off');
for i = 1:11
    subplot(4,3,i)
    image([1 4], [1 5], E_act_overall_Th(:,:,t+(i-1)*0.005/dt));
    set(gca, 'XTick', [1:4], 'YTick', [1:5]);
    title(['t = ',num2str((i-1)*5),'ms']);
    colormap hot
    axis equal
    axis tight
    c = colorbar;
    c.Label.String = 'E (Spikes/s)';
%     c.Limits = [0 90];
%     c.Ticks = [0 30 60 90];
%     c.TickLabels = [0 30 60 90];
    if i > 1
        c.Visible = 'off';
    end
end
% saveas(f,'Figure/Prop_PS1.pdf');

figure
% f = figure('Visible','off');
subplot(2,1,1);
imagesc([0 time_win*10^3], [1 9], E_plot(:,t:t+time_win/dt));
set(gca,'YDir','normal','YTick',1:9)
xlabel('time(ms)')
ylabel('Barrel')  
colormap jet
c = colorbar;
c.Label.String = 'E (Spikes/s)';
title('Thalamus');
subplot(2,1,2);
tim = (0:dt:time_win)*10^3; % (in milliseconds)
hold on
for i = 1:9
    plot(tim, E_plot(i,t:t+time_win/dt))
end
c = colorbar;
c.Visible = 'off';
xlabel('time(ms)')
ylabel('E (Spikes/s)')
% saveas(f,'Figure/Prop_PS2.pdf');

%% 4. Low condition 
clear
k = 1;
kk = 1;
load('Simulation Results/meta_data.mat'); 
load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat']);
n_stim_plot = 40;
t_onset = Stim_Onsets(1);
t_offset = t_onset + n_stim_plot*(duration+ISI)/dt - 1;
tim = dt:dt:n_stim_plot*(duration+ISI);

% --- thalamocortical synapse dynamics ---
figure
subplot(2,1,1) 
plot(tim, reshape(z_act_overall(4,2, t_onset:t_offset),[1,length(t_onset:t_offset)]), '-b')
xlim([0 n_stim_plot*(duration+ISI)])
ylabel('z')
subplot(2,1,2) 
plot(tim, reshape(Spa_Temp(1,1,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-b')
hold on;
plot(tim, reshape(Spa_Temp(1,2,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-r')
xlim([0 n_stim_plot*(duration+ISI)])
set(gca,'XTickLabel',[],'YTickLabel',[])

% ---- L4 ---- :
figure
% f = figure('Visible','off');
subplot(4,1,1)  
plot(tim, E_plot(5,t_onset:t_offset), '-b')
xlim([0 n_stim_plot*(duration+ISI)])
ylabel('E(Spikes/s)')
set(gca,'XTickLabel',[])  
c = colorbar;
c.Visible = 'off';
title('L4');
subplot(4,1,2) 
plot(tim, reshape(Spa_Temp(1,1,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-b')
hold on;
plot(tim, reshape(Spa_Temp(1,2,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-r')
xlim([0 n_stim_plot*(duration+ISI)])
set(gca,'XTickLabel',[],'YTickLabel',[])
c = colorbar;
c.Visible = 'off';
subplot(4,1,3)  
imagesc([0 n_stim_plot*(duration+ISI)], [1 9], E_plot(:,t_onset:t_offset));
set(gca,'YDir','normal','XTickLabel',[], 'YTick', [1:2:9]);
ylabel('Column')
colormap hot
c = colorbar;
c.Label.String = 'E (Spikes/s)';
subplot(4,1,4) 
imagesc([0 n_stim_plot*(duration+ISI)], [1 9], x_plot(:,t_onset:t_offset));
set(gca,'YDir','normal', 'YTick', [1:2:9]);
xlabel('time(s)')
ylabel('Column')
colormap hot
c = colorbar;
c.Label.String = 'x';
% saveas(f,'Figure/Low_condition1.pdf');

figure
% f = figure('Visible','off');
subplot(3,1,1) 
plot(tim, E_plot(7,t_onset:t_offset), '-b')
hold on
plot(tim, E_plot(5,t_onset:t_offset), '-r')
xlim([0 n_stim_plot*(duration+ISI)])
set(gca,'XTickLabel',[])  
ylabel('E(Spikes/s)')
title('L4');
subplot(3,1,2) 
plot(tim, x_plot(7,t_onset:t_offset), '-b')
hold on
plot(tim, x_plot(5,t_onset:t_offset), '-r')
xlim([0 n_stim_plot*(duration+ISI)])
set(gca,'XTickLabel',[])  
ylabel('x')
subplot(3,1,3) 
plot(tim, reshape(Spa_Temp(1,1,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-b')
hold on;
plot(tim, reshape(Spa_Temp(1,2,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-r')
set(gca,'YTickLabel',[])  
xlim([0 n_stim_plot*(duration+ISI)])
legend('AW1','PW')
xlabel('time(s)')
% saveas(f,'Figure/Low_condition2.pdf');

% ---- L6 ---- :
E_plot = zeros(9,num_steps);
E_plot(1,:) = reshape(E_act_overall_L6(PW(1)-1,PW(2)-1,:), [1, num_steps]);
E_plot(2,:) = reshape(E_act_overall_L6(PW(1)-1,PW(2)+1,:), [1, num_steps]);
E_plot(3,:) = reshape(E_act_overall_L6(PW(1)-1,PW(2),:), [1, num_steps]);
E_plot(4,:) = reshape(E_act_overall_L6(PW(1),PW(2)-1,:), [1, num_steps]);
E_plot(5,:) = reshape(E_act_overall_L6(PW(1),PW(2),:), [1, num_steps]);
E_plot(6,:) = reshape(E_act_overall_L6(PW(1),PW(2)+1,:), [1, num_steps]);
E_plot(7,:) = reshape(E_act_overall_L6(PW(1)+1,PW(2),:), [1, num_steps]);
E_plot(8,:) = reshape(E_act_overall_L6(PW(1)+1,PW(2)+1,:), [1, num_steps]);
E_plot(9,:) = reshape(E_act_overall_L6(PW(1)+1,PW(2)-1,:), [1, num_steps]);
% Resource of L4 barrels surrounding PW to plot
x_plot = zeros(9,num_steps);
x_plot(1,:) = reshape(x_act_overall_L6(PW(1)-1,PW(2)-1,:), [1, num_steps]);
x_plot(2,:) = reshape(x_act_overall_L6(PW(1)-1,PW(2)+1,:), [1, num_steps]);
x_plot(3,:) = reshape(x_act_overall_L6(PW(1)-1,PW(2),:), [1, num_steps]);
x_plot(4,:) = reshape(x_act_overall_L6(PW(1),PW(2)-1,:), [1, num_steps]);
x_plot(5,:) = reshape(x_act_overall_L6(PW(1),PW(2),:), [1, num_steps]);
x_plot(6,:) = reshape(x_act_overall_L6(PW(1),PW(2)+1,:), [1, num_steps]);
x_plot(7,:) = reshape(x_act_overall_L6(PW(1)+1,PW(2),:), [1, num_steps]);
x_plot(8,:) = reshape(x_act_overall_L6(PW(1)+1,PW(2)+1,:), [1, num_steps]);
x_plot(9,:) = reshape(x_act_overall_L6(PW(1)+1,PW(2)-1,:), [1, num_steps]);
figure
% f = figure('Visible','off');
subplot(4,1,1)  
plot(tim, E_plot(5,t_onset:t_offset), '-b')
xlim([0 n_stim_plot*(duration+ISI)])
ylabel('E(Spikes/s)')
set(gca,'XTickLabel',[])  
c = colorbar;
c.Visible = 'off';
title('L6');
subplot(4,1,2) 
plot(tim, reshape(Spa_Temp(1,1,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-b')
hold on;
plot(tim, reshape(Spa_Temp(1,2,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-r')
xlim([0 n_stim_plot*(duration+ISI)])
set(gca,'XTickLabel',[],'YTickLabel',[])
c = colorbar;
c.Visible = 'off';
subplot(4,1,3)  
imagesc([0 n_stim_plot*(duration+ISI)], [1 9], E_plot(:,t_onset:t_offset));
set(gca,'YDir','normal','XTickLabel',[], 'YTick', [1:2:9]);
ylabel('Column')
colormap hot
c = colorbar;
c.Label.String = 'E (Spikes/s)';
subplot(4,1,4) 
imagesc([0 n_stim_plot*(duration+ISI)], [1 9], x_plot(:,t_onset:t_offset));
set(gca,'YDir','normal', 'YTick', [1:2:9]);
xlabel('time(s)')
ylabel('Column')
colormap hot
c = colorbar;
c.Label.String = 'x';
% saveas(f,'Figure/Low_condition1.pdf');

figure
% f = figure('Visible','off');
subplot(3,1,1) 
plot(tim, E_plot(7,t_onset:t_offset), '-b')
hold on
plot(tim, E_plot(5,t_onset:t_offset), '-r')
xlim([0 n_stim_plot*(duration+ISI)])
set(gca,'XTickLabel',[])  
ylabel('E(Spikes/s)')
title('L6');
subplot(3,1,2) 
plot(tim, x_plot(7,t_onset:t_offset), '-b')
hold on
plot(tim, x_plot(5,t_onset:t_offset), '-r')
xlim([0 n_stim_plot*(duration+ISI)])
set(gca,'XTickLabel',[])  
ylabel('x')
subplot(3,1,3) 
plot(tim, reshape(Spa_Temp(1,1,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-b')
hold on;
plot(tim, reshape(Spa_Temp(1,2,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-r')
set(gca,'YTickLabel',[])  
xlim([0 n_stim_plot*(duration+ISI)])
legend('AW1','PW')
xlabel('time(s)')
% saveas(f,'Figure/Low_condition2.pdf');

%% 5. Many-standard condition
clear
k = 1;
kk = 4;
load('Simulation Results/meta_data.mat'); 
load(['Simulation Results/run_par' num2str(Par_Arr(k)) '_' Cond_Code{kk} '.mat']);
n_stim_plot = 40;
t_onset = Stim_Onsets(1);
t_offset = t_onset + n_stim_plot*(duration+ISI)/dt - 1;
tim = dt:dt:n_stim_plot*(duration+ISI);

% ---- L4---- :
figure
% f = figure('Visible','off');
subplot(4,1,1)  
plot(tim, E_plot(5,t_onset:t_offset), '-b')
xlim([0 n_stim_plot*(duration+ISI)])
ylabel('E(Spikes/s)')
set(gca,'XTickLabel',[])  
c = colorbar;
c.Visible = 'off';
title('L4');
subplot(4,1,2) 
plot(tim, reshape(Spa_Temp(1,1,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-b')
hold on
plot(tim, reshape(Spa_Temp(1,2,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-r')
plot(tim, reshape(Spa_Temp(1,3,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]),'-g')
plot(tim, reshape(Spa_Temp(1,4,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]),'-k')
xlim([0 n_stim_plot*(duration+ISI)])
set(gca,'XTickLabel',[],'YTickLabel',[])  
c = colorbar;
c.Visible = 'off';
subplot(4,1,3)  
imagesc([0 n_stim_plot*(duration+ISI)], [1 9], E_plot(:,t_onset:t_offset));
set(gca,'YDir','normal','XTickLabel',[], 'YTick', [1:2:9]);
ylabel('Column')
colormap hot
c = colorbar;
c.Label.String = 'E (Spikes/s)';
subplot(4,1,4) 
imagesc([0 n_stim_plot*(duration+ISI)], [1 9], x_plot(:,t_onset:t_offset));
set(gca,'YDir','normal', 'YTick', [1:2:9]);
ylabel('Column')
colormap hot
c = colorbar;
c.Label.String = 'x';
% saveas(f,'Figure/MS_condition1.pdf');

figure
% f = figure('Visible','off');
subplot(3,1,1) 
plot(tim, E_plot(7,t_onset:t_offset), '-b')
hold on
plot(tim, E_plot(5,t_onset:t_offset), '-r')
xlim([0 n_stim_plot*(duration+ISI)])
set(gca,'XTickLabel',[])  
ylabel('E(Spikes/s)')
title('L4');
subplot(3,1,2) 
plot(tim, x_plot(7,t_onset:t_offset), '-b')
hold on
plot(tim, x_plot(5,t_onset:t_offset), '-r')
xlim([0 n_stim_plot*(duration+ISI)])
set(gca,'XTickLabel',[])  
ylabel('x')
subplot(3,1,3) 
plot(tim, reshape(Spa_Temp(1,1,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-b')
hold on
plot(tim, reshape(Spa_Temp(1,2,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-r')
plot(tim, reshape(Spa_Temp(1,3,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]),'-g')
plot(tim, reshape(Spa_Temp(1,4,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]),'-k')
set(gca,'YTickLabel',[])  
xlim([0 n_stim_plot*(duration+ISI)])
xlabel('time(s)')
legend('AW1','PW','AW2','AW3')
% saveas(f,'Figure/MS_condition2.pdf');

% ---- L6 ---- :
E_plot = zeros(9,num_steps);
E_plot(1,:) = reshape(E_act_overall_L6(PW(1)-1,PW(2)-1,:), [1, num_steps]);
E_plot(2,:) = reshape(E_act_overall_L6(PW(1)-1,PW(2)+1,:), [1, num_steps]);
E_plot(3,:) = reshape(E_act_overall_L6(PW(1)-1,PW(2),:), [1, num_steps]);
E_plot(4,:) = reshape(E_act_overall_L6(PW(1),PW(2)-1,:), [1, num_steps]);
E_plot(5,:) = reshape(E_act_overall_L6(PW(1),PW(2),:), [1, num_steps]);
E_plot(6,:) = reshape(E_act_overall_L6(PW(1),PW(2)+1,:), [1, num_steps]);
E_plot(7,:) = reshape(E_act_overall_L6(PW(1)+1,PW(2),:), [1, num_steps]);
E_plot(8,:) = reshape(E_act_overall_L6(PW(1)+1,PW(2)+1,:), [1, num_steps]);
E_plot(9,:) = reshape(E_act_overall_L6(PW(1)+1,PW(2)-1,:), [1, num_steps]);
% Resource of L4 barrels surrounding PW to plot
x_plot = zeros(9,num_steps);
x_plot(1,:) = reshape(x_act_overall_L6(PW(1)-1,PW(2)-1,:), [1, num_steps]);
x_plot(2,:) = reshape(x_act_overall_L6(PW(1)-1,PW(2)+1,:), [1, num_steps]);
x_plot(3,:) = reshape(x_act_overall_L6(PW(1)-1,PW(2),:), [1, num_steps]);
x_plot(4,:) = reshape(x_act_overall_L6(PW(1),PW(2)-1,:), [1, num_steps]);
x_plot(5,:) = reshape(x_act_overall_L6(PW(1),PW(2),:), [1, num_steps]);
x_plot(6,:) = reshape(x_act_overall_L6(PW(1),PW(2)+1,:), [1, num_steps]);
x_plot(7,:) = reshape(x_act_overall_L6(PW(1)+1,PW(2),:), [1, num_steps]);
x_plot(8,:) = reshape(x_act_overall_L6(PW(1)+1,PW(2)+1,:), [1, num_steps]);
x_plot(9,:) = reshape(x_act_overall_L6(PW(1)+1,PW(2)-1,:), [1, num_steps]);
figure
% f = figure('Visible','off');
subplot(4,1,1)  
plot(tim, E_plot(5,t_onset:t_offset), '-b')
xlim([0 n_stim_plot*(duration+ISI)])
ylabel('E(Spikes/s)')
set(gca,'XTickLabel',[])  
c = colorbar;
c.Visible = 'off';
title('L6');
subplot(4,1,2) 
plot(tim, reshape(Spa_Temp(1,1,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-b')
hold on
plot(tim, reshape(Spa_Temp(1,2,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-r')
plot(tim, reshape(Spa_Temp(1,3,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]),'-g')
plot(tim, reshape(Spa_Temp(1,4,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]),'-k')
xlim([0 n_stim_plot*(duration+ISI)])
set(gca,'XTickLabel',[],'YTickLabel',[])  
c = colorbar;
c.Visible = 'off';
subplot(4,1,3)  
imagesc([0 n_stim_plot*(duration+ISI)], [1 9], E_plot(:,t_onset:t_offset));
set(gca,'YDir','normal','XTickLabel',[], 'YTick', [1:2:9]);
ylabel('Column')
colormap hot
c = colorbar;
c.Label.String = 'E (Spikes/s)';
subplot(4,1,4) 
imagesc([0 n_stim_plot*(duration+ISI)], [1 9], x_plot(:,t_onset:t_offset));
set(gca,'YDir','normal', 'YTick', [1:2:9]);
ylabel('Column')
colormap hot
c = colorbar;
c.Label.String = 'x';
% saveas(f,'Figure/MS_condition1.pdf');

figure
% f = figure('Visible','off');
subplot(3,1,1) 
plot(tim, E_plot(7,t_onset:t_offset), '-b')
hold on
plot(tim, E_plot(5,t_onset:t_offset), '-r')
xlim([0 n_stim_plot*(duration+ISI)])
set(gca,'XTickLabel',[])  
ylabel('E(Spikes/s)')
title('L6');
subplot(3,1,2) 
plot(tim, x_plot(7,t_onset:t_offset), '-b')
hold on
plot(tim, x_plot(5,t_onset:t_offset), '-r')
xlim([0 n_stim_plot*(duration+ISI)])
set(gca,'XTickLabel',[])  
ylabel('x')
subplot(3,1,3) 
plot(tim, reshape(Spa_Temp(1,1,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-b')
hold on
plot(tim, reshape(Spa_Temp(1,2,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]), '-r')
plot(tim, reshape(Spa_Temp(1,3,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]),'-g')
plot(tim, reshape(Spa_Temp(1,4,t_onset:t_offset), [1,int64(n_stim_plot*(duration+ISI)/dt)]),'-k')
set(gca,'YTickLabel',[])  
xlim([0 n_stim_plot*(duration+ISI)])
xlabel('time(s)')
legend('AW1','PW','AW2','AW3')
% saveas(f,'Figure/MS_condition2.pdf');


%}