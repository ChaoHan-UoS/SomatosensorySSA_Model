%%  Plot simulation results 

close;
clear;

plot_1_B = 0;
plot_1_D = 1; % Propogation of PS and mean dev and std response
plot_1_DD = 1; 
plot_1_E = 0;
plot_3_A_C = 1; % SSA in oddball protocol
plot_4_C = 1; % True deviance detection
plot_4_F_H = 1; % SSA in many-standard protocol
plot_3_F = 0;
plot_2_A = 0;
plot_2_B_C = 0;
plot_5_B_left = 0;
plot_5_B_right = 0;
plot_5_C_left = 0;
plot_5_C_right = 0;
plot_6_A_D = 0;
plot_7_A = 0;
plot_late_resp = 0;

AXES_FONTSIZE = 9;
LineWidth = 0.8;
COLORBAR_FONTSIZE = 9;

%% Fig1.B
if plot_1_B 
    load('Simulation Results/net1/run_net1_L.mat')
    load('Simulation Results/net1/meta_data.mat')
    t = find(Spec_Temp(1,:,2), floor(duration/dt));
    tim = dt:dt:duration;
    tim = tim*10^3;
    
    figure(1)
    subplot(2,1,1)
    plot(tim, E_mean(Rec_Column,t), '-b','LineWidth',LineWidth)
    hold on
    plot(tim, I_mean(Rec_Column,t), '-r','LineWidth',LineWidth)
    pos = get(gcf,'position');
    set(gcf,'position',[pos(1),pos(2),(0.6)*pos(3),pos(4)]);
    ylabel('Spikes/s')
    legend('E','I')
    legend('boxoff')
    set(gca,'FontSize',AXES_FONTSIZE,'XTick',[0 50],'XTickLabel',[0 50])
    xlim([0 50])
    axis square

    subplot(2,1,2)
    plot(tim, x_mean(Rec_Column,t), '-b','LineWidth',LineWidth)
    hold on;
    plot(tim, y_mean(Rec_Column,t), '-r','LineWidth',LineWidth)
    pos = get(gcf,'position');
    set(gcf,'position',[pos(1),pos(2),(0.6)*pos(3),pos(4)]);
    xlabel('time(ms)')
    ylabel('x')
    legend('x','y')
    legend('boxoff')
    set(gca,'FontSize',AXES_FONTSIZE,'XTick',[0 50],'XTickLabel',[0 50])
    ylim([0 1])
    xlim([0 50])
    axis square
    
end


%% Fig1.E
if plot_1_E 
    load('Simulation Results/net1/run_net1_L.mat')
    load('Simulation Results/net1/meta_data.mat')
    t_onset = find(Spec_Temp(1,:,2), 1);
    t_offset = t_onset + floor(5*(duration+ISI)/dt) - 1;
    tim = dt:dt:5*(duration+ISI);
    
    figure(4)
    subplot(2,1,1)
    plot(tim, E_mean(Rec_Column,t_onset:t_offset), '-b','LineWidth',LineWidth)
    hold on
    plot(tim, I_mean(Rec_Column,t_onset:t_offset), '-r','LineWidth',LineWidth)
    pos = get(gcf,'position');
    set(gcf,'position',[pos(1),pos(2),0.6*pos(3),pos(4)]);
    ylabel('E(Spikes/s)')
    legend('E','I')
    legend('boxoff')
    set(gca,'FontSize',AXES_FONTSIZE)
    xlim([0 5*(duration+ISI)])
    axis square

    subplot(2,1,2)
    plot(tim, x_mean(Rec_Column,t_onset:t_offset), '-b','LineWidth',LineWidth)
    hold on;
    plot(tim, y_mean(Rec_Column,t_onset:t_offset), '-r','LineWidth',LineWidth)
    pos = get(gcf,'position');
    set(gcf,'position',[pos(1),pos(2),0.6*pos(3),pos(4)]);
    xlabel('time(s)')
    legend('x','y')
    legend('boxoff')
    set(gca,'FontSize',AXES_FONTSIZE)
    xlim([0 5*(duration+ISI)])
    ylim([0 1])
    axis square
end


%% Low condition 
if plot_3_A_C
    load('Simulation Results/net1/run_net1_L.mat')
    load('Simulation Results/net1/meta_data.mat')
    n_stim_plot = 50;
    t_onset = Stim_Onsets(1);
    t_offset = t_onset + floor(n_stim_plot*(duration+ISI)/dt)-1;
%     t_offset = t_onset + floor(n_stim_plot*(duration+ISI)/dt);
    tim = dt:dt:n_stim_plot*(duration+ISI);
       
    figure
    subplot(4,1,1)  
    plot(tim, E_plot(5,t_onset:t_offset), '-b')
    xlim([0 n_stim_plot*(duration+ISI)])
    ylabel('E(Spikes/s)')
    set(gca,'XTickLabel',[])  
    
%     subplot(4,1,2) 
%     imagesc([0 n_stim_plot*(duration+ISI)], [Rec_Column - Freq_Diff Rec_Column], [Spec_Temp(1,t_onset:t_offset,1);Spec_Temp(1,t_onset:t_offset,2)]);
%     set(gca,'YDir','normal')   
%     map = [1 1 1; 0 0 1];
%     ax1 = subplot(4,1,2); 
%     colormap(ax1, map)
%     set(gca, 'YTick', [Rec_Column - Freq_Diff Rec_Column]);
%     ylabel('Column') 
    
    subplot(4,1,2) 
    plot(tim, reshape(Spec_Temp(1,1,t_onset:t_offset), [1,ceil(n_stim_plot*(duration+ISI)/dt)]), '-b')
    hold on;
    plot(tim, reshape(Spec_Temp(1,2,t_onset:t_offset), [1,ceil(n_stim_plot*(duration+ISI)/dt)]), '-r')
    set(gca,'YTickLabel',[])  
    xlim([0 n_stim_plot*(duration+ISI)])
    set(gca,'XTickLabel',[])  

    subplot(4,1,3)  
    imagesc([0 n_stim_plot*(duration+ISI)], [1 9], E_plot(:,t_onset:t_offset));
    set(gca,'YDir','normal')  
    ax2 = subplot(4,1,3); 
    c = hot;
%     c = flipud(c);
    colormap(ax2, c)
    ylabel('Column') 
    set(gca,'XTickLabel',[])  
    colorbar
    ylabel('Column') 
    set(gca, 'YTick', [1:9]);
    c = colorbar;
    c.Label.String = 'E (Spikes/s)';
%     c.Limits = [0 80];
%     c.TickLabels = [0 20 40 60 80]
    c.Position = [0.922 0.329 0.021 0.156]
   
    subplot(4,1,4) 
    imagesc([0 n_stim_plot*(duration+ISI)], [1 9], x_plot(:,t_onset:t_offset));
    set(gca,'YDir','normal')  
    ax3 = subplot(4,1,4); 
    colormap(ax3, hot)
    ylabel('Column') 
    colorbar
    xlabel('time(s)')
    ylabel('Column') 
    set(gca, 'YTick', [1:9]);
    c = colorbar;
    c.Label.String = 'x';
%     c.Limits = [0.35 0.7];
%     c.TickLabels = [0.4 0.5 0.6 0.7]
    c.Position = [0.922 0.111 0.021 0.156]
 
   
    figure
    subplot(3,1,1) 
    plot(tim, E_plot(7,t_onset:t_offset), '-b')
    hold on
    plot(tim, E_plot(5,t_onset:t_offset), '-r')
    xlim([0 n_stim_plot*(duration+ISI)])
    ylabel('E(Spikes/s)')

    subplot(3,1,2) 
    plot(tim, x_plot(7,t_onset:t_offset), '-b')
    hold on
    plot(tim, x_plot(5,t_onset:t_offset), '-r')
    xlim([0 n_stim_plot*(duration+ISI)])
    ylabel('x')
    
    subplot(3,1,3) 
    plot(tim, reshape(Spec_Temp(1,1,t_onset:t_offset), [1,ceil(n_stim_plot*(duration+ISI)/dt)]), '-b')
    hold on;
    plot(tim, reshape(Spec_Temp(1,2,t_onset:t_offset), [1,ceil(n_stim_plot*(duration+ISI)/dt)]), '-r')
    set(gca,'YTickLabel',[])  
    xlim([0 n_stim_plot*(duration+ISI)])
    set(gca,'XTickLabel',[])  
    legend('AW1','PW')
end


%% Many-standard
if plot_4_F_H
    load('Simulation Results/net1/run_net1_MS.mat')
    load('Simulation Results/net1/meta_data.mat')
    n_stim_plot = 50;
    t_onset = Stim_Onsets(1);
    t_offset = t_onset + floor(n_stim_plot*(duration+ISI)/dt) - 1;
%     t_offset = t_onset + floor(n_stim_plot*(duration+ISI)/dt);
    tim = dt:dt:n_stim_plot*(duration+ISI);
    
    figure
    subplot(4,1,1)  
    plot(tim, E_plot(5,t_onset:t_offset), '-b')
    xlim([0 n_stim_plot*(duration+ISI)])
    ylabel('E(Spikes/s)')   
    set(gca,'XTickLabel',[])  
    
    subplot(4,1,2) 
    plot(tim, reshape(Spec_Temp(1,1,t_onset:t_offset), [1,ceil(n_stim_plot*(duration+ISI)/dt)]),'-b');
    hold on
    plot(tim, reshape(Spec_Temp(1,2,t_onset:t_offset), [1,ceil(n_stim_plot*(duration+ISI)/dt)]),'-r');
    plot(tim, reshape(Spec_Temp(1,3,t_onset:t_offset), [1,ceil(n_stim_plot*(duration+ISI)/dt)]),'-g');
    plot(tim, reshape(Spec_Temp(1,4,t_onset:t_offset), [1,ceil(n_stim_plot*(duration+ISI)/dt)]),'-k');
    set(gca,'YTickLabel',[])  
    xlim([0 n_stim_plot*(duration+ISI)])
%     legend('AW1','PW','AW2','AW3')
    set(gca,'XTickLabel',[])    
    
%     subplot(4,1,2) 
%     imagesc([0 n_stim_plot*(duration+ISI)], [1 4], reshape(Spec_Temp(1,:,t_onset:t_offset), [4, length(t_onset:t_offset)]));
%     set(gca,'YDir','normal')   
%     map = [1 1 1; 0 0 1];
%     ax1 = subplot(4,1,2); 
%     colormap(ax1, map)  
%     set(gca, 'YTickLabels', {7 5 9 8});
%     ylabel('Barrel') 
 
    subplot(4,1,3)  
    imagesc([0 n_stim_plot*(duration+ISI)], [1 9], E_plot(:,t_onset:t_offset));
    set(gca,'YDir','normal')  
    ax2 = subplot(4,1,3); 
    colormap(ax2, hot)
    ylabel('Column') 
    set(gca,'XTickLabel',[])  
    colorbar
    ylabel('Column') 
    set(gca, 'YTick', [1:9]);
    c = colorbar;
    c.Label.String = 'E (Spikes/s)';
%     c.Limits = [0 80];
%     c.TickLabels = [0 20 40 60 80]
    c.Position = [0.922 0.329 0.021 0.156]
   
    subplot(4,1,4) 
    imagesc([0 n_stim_plot*(duration+ISI)], [1 9], x_plot(:,t_onset:t_offset));
    set(gca,'YDir','normal')  
    ax3 = subplot(4,1,4); 
    colormap(ax3, hot)
    ylabel('Column') 
    colorbar
    xlabel('time(s)')
    ylabel('Column')
    set(gca, 'YTick', [1:9]);
    c = colorbar;
    c.Label.String = 'x';
%     c.Limits = [0.35 0.7];
%     c.TickLabels = [0.4 0.5 0.6 0.7]
    c.Position = [0.922 0.111 0.021 0.156]
    
       
    figure
    subplot(3,1,1) 
%     plot(tim, E_plot(7,t_onset:t_offset))
    hold on
    plot(tim, E_plot(5,t_onset:t_offset),'-r')
%     plot(tim, E_plot(9,t_onset:t_offset))
%     plot(tim, E_plot(8,t_onset:t_offset))
    xlim([0 n_stim_plot*(duration+ISI)])
    ylabel('E(Spikes/s)')

    subplot(3,1,2) 
%     plot(tim, x_plot(7,t_onset:t_offset))
    hold on
    plot(tim, x_plot(5,t_onset:t_offset),'-r')
%     plot(tim, x_plot(9,t_onset:t_offset))
%     plot(tim, x_plot(8,t_onset:t_offset))
    xlim([0 n_stim_plot*(duration+ISI)])
    ylabel('x')
    
    subplot(3,1,3) 
    plot(tim, reshape(Spec_Temp(1,1,t_onset:t_offset), [1,ceil(n_stim_plot*(duration+ISI)/dt)]),'-b');
    hold on
    plot(tim, reshape(Spec_Temp(1,2,t_onset:t_offset), [1,ceil(n_stim_plot*(duration+ISI)/dt)]),'-r');
    plot(tim, reshape(Spec_Temp(1,3,t_onset:t_offset), [1,ceil(n_stim_plot*(duration+ISI)/dt)]),'-g');
    plot(tim, reshape(Spec_Temp(1,4,t_onset:t_offset), [1,ceil(n_stim_plot*(duration+ISI)/dt)]),'-k');
    set(gca,'YTickLabel',[])  
    xlim([0 n_stim_plot*(duration+ISI)])
    xlabel('time(s)')
    legend('AW1','PW','AW2','AW3')
    

end

%% 
if plot_3_F
    nev_cond_code{1} = 'L';
    nev_cond_code{2} = 'H';
    nev_cond_code{3} = 'E';
    for ii = 1:3 
        save('counter.mat','ii','nev_cond_code')
        clear all
        load('counter.mat')
        load(['Simulation Results/net1/run_net1_' nev_cond_code{ii} '.mat'])
        j = 0;
        k = 0;
        E_sum_f2 = zeros(1, length(Single_Stim));
        E_sum_f1 = zeros(1, length(Single_Stim));
        for ns = 1:n_stim  
            if Spec_Temp(1,Stim_Onsets(ns),2) > 0
                E_sum_f2 = E_sum_f2 + E_mean(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1);
                j = j + 1;
            else
                E_sum_f1 = E_sum_f1 + E_mean(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1);
                k = k + 1;
            end
        end
        E_mean_f2 = E_sum_f2/j;
        E_mean_f1 = E_sum_f1/k;
        figure(1)
        tim = dt:dt:length(Single_Stim)*dt;
        tim = tim*1000;
        if ii == 1
            subplot(1,3,1)
            plot(tim, E_mean_f2, '-r')
            hold on
            plot(tim, E_mean_f1, '-b')
            axis([0 50 0 60])
            xlabel('time(ms)')
            ylabel('E(Spikes/s)')
            legend('f2(deviant)','f1(standard)')
        elseif ii == 2
            subplot(1,3,3)
            plot(tim, E_mean_f2, '-r')
            hold on
            plot(tim, E_mean_f1, '-b')
            axis([0 50 0 60])
            legend('f1(deviant)','f2(standard)')
        else
            subplot(1,3,2)
            plot(tim, E_mean_f1, '-r')
            hold on
            plot(tim, E_mean_f2, '-b')
            axis([0 50 0 60])
            legend('f1(equal)','f2(equal)')
        end
    end
end
            


%% 
if plot_2_A
    load('Simulation Results/net1/run_net1_L.mat')
    load('Simulation Results/net1/meta_data.mat')
    n_stim_plot = 50;
    t_onset = find(Spec_Temp(1,:,1), 1) - 2*(duration+ISI)/dt;
    t_offset = find(Spec_Temp(1,:,1), 1) + floor(n_stim_plot*(duration+ISI)/dt) - 1;
    tim = dt:dt:(n_stim_plot+2)*(duration+ISI);
    E_act_sel = reshape(E_act_overall(1,5:10:95,:), 10, num_steps);
       
    figure(1)
    subplot(5,1,1) 
    imagesc([0 n_stim_plot*(duration+ISI)], [10 12], [Spec_Temp(1,t_onset:t_offset,1);Spec_Temp(1,t_onset:t_offset,2)]);
    set(gca,'YDir','normal')   
    map = [1 1 1; 0 0 1];
    ax1 = subplot(5,1,1); 
    colormap(ax1, map)
    axis off
    text(-1,10,'f1')
    text(-1,12,'f2')
     
    subplot(2,1,2) 
    imagesc([0 n_stim_plot*(duration+ISI)], [5 95], E_act_sel(:,t_onset:t_offset))
    set(gca,'YDir','normal')   
    ax2 = subplot(2,1,2); 
    colormap(ax2, hot)
    xlabel('time(s)')
    ylabel('Neuron') 
    set(gca, 'YTick', 5:10:95);         
    
end

%%

if plot_2_B_C
    sel1 = 23;
    sel2 = 25;
    sel_neu = 8; % no.75 neuron
    load('Simulation Results/net1/run_net1_L.mat')
    load('Simulation Results/net1/meta_data.mat')
    tim = dt:dt:duration;
    tim = tim*1000;
    E_act_sel = reshape(E_act_overall(1,5:10:95,:), 10, num_steps);
    E_sum_ST = zeros(1,duration/dt);
    figure(1)
    subplot(1,4,1) 
    for i = 1:10
        plot(tim, E_act_sel(i,Time_Ind(sel1,1):Time_Ind(sel1,2)), '-b', 'LineWidth',0.1)
        hold on
        E_sum_ST = E_sum_ST + E_act_sel(i,Time_Ind(sel1,1):Time_Ind(sel1,2));
    end
    plot(tim, E_sum_ST/10, '-k','LineWidth',1)
    axis([0 50 0 85])
    xlabel('time(ms)')
    ylabel('E(Spikes/s)')
    title('Single Trial')
    text(2,82,'No PS')
    box off
    
    E_sum_ST = zeros(1,duration/dt);
    subplot(1,4,2) 
    for i = 1:10
        plot(tim, E_act_sel(i,Time_Ind(sel2,1):Time_Ind(sel2,2)), '-r', 'LineWidth',0.1)
        hold on
        E_sum_ST = E_sum_ST + E_act_sel(i,Time_Ind(sel2,1):Time_Ind(sel2,2));
    end
    plot(tim, E_sum_ST/10, '-k','LineWidth',1)
    axis([0 50 0 85])
    text(2,82,'PS')
    box off
   
    j = 0;
    k = 0;
    E_sum_f2 = zeros(1, length(Single_Stim));
    E_sum_f1 = zeros(1, length(Single_Stim));
    for ns = 1:n_stim  
        if Spec_Temp(1,Stim_Onsets(ns),2) > 0
            E_sum_f2 = E_sum_f2 + E_act_sel(sel_neu, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1);
            subplot(1,4,4)
            plot(tim, E_act_sel(sel_neu, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1), '-r', 'LineWidth',0.1)
            hold on
            j = j + 1;
        else
            E_sum_f1 = E_sum_f1 + E_act_sel(sel_neu, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1);
            subplot(1,4,3)
            plot(tim, E_act_sel(sel_neu, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1), '-b', 'LineWidth',0.1)
            hold on
            k = k + 1;
        end
    end
    E_mean_f1 = E_sum_f1/k;
    E_mean_f2 = E_sum_f2/j;

    subplot(1,4,3)
    plot(tim, E_mean_f1, '-k','LineWidth',1)
    xlabel('time(ms)')
    ylabel('E(Spikes/s)')
    title('Single Neuron')
    axis([0 50 0 85])
    text(2,82,'f1')
    box off
    
    subplot(1,4,4)
    plot(tim, E_mean_f2, '-k','LineWidth',1)
    axis([0 50 0 85])
    text(2,82,'f2')
    box off

end
            

%% 
if plot_5_B_left
    load('Simulation Results/net1/run_net1_L_P0.1_A5.mat')
    color{1} = 'k';
    color{2} = 'b';
    color{3} = 'r';
    f1 = 10;
    f2 = 12;
    CSI = zeros(length(Probs_Arr),length(A_Arr)); %3*20
    for j = 1:length(Probs_Arr)
        for k = 1:length(A_Arr)
            sw = zeros(1,4);
            for i = 1:2
                load(['Simulation Results/net1/run_net1_' nev_cond_code{i} '_P' num2str(Probs_Arr(j)) '_A' num2str(A_Arr(k)) '.mat'], 'E_mean', 'Oddball')
                m = 0;
                n = 0;
                E_sum_a_f2 = zeros(1, length(Single_Stim)+0.045/dt);
                E_sum_b_f2 = zeros(1, 0.005/dt);
                E_sum_a_f1 = zeros(1, length(Single_Stim)+0.045/dt);
                E_sum_b_f1 = zeros(1, 0.005/dt);
                for ns = 1:n_stim
                    if Oddball(ns) == f2
                        E_sum_a_f2 = E_sum_a_f2 + E_mean(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1+0.045/dt); % integrating E_mean from stimulus onset to 45 ms after offeset
                        E_sum_b_f2 = E_sum_b_f2 + E_mean(Rec_Column, Stim_Onsets(ns)-0.005/dt:Stim_Onsets(ns)-1); % baseline 5 ms before stimulus onset
                        m = m + 1;
                    else
                        E_sum_a_f1 = E_sum_a_f1 + E_mean(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1+0.045/dt); % integrating E_mean from stimulus onset to 45 ms after offeset
                        E_sum_b_f1 = E_sum_b_f1 + E_mean(Rec_Column, Stim_Onsets(ns)-0.005/dt:Stim_Onsets(ns)-1); % baseline 5 ms before stimulus onset
                        n = n + 1;
                    end
                end
                E_mean_a_f2 = E_sum_a_f2/m;
                E_mean_b_f2 = E_sum_b_f2/m;
                E_mean_a_f1 = E_sum_a_f1/n;
                E_mean_b_f1 = E_sum_b_f1/n;
                Spcount_f2 = sum(E_mean_a_f2(:)*dt) - sum(E_mean_b_f2(:)*dt);
                Spcount_f1 = sum(E_mean_a_f1(:)*dt) - sum(E_mean_b_f1(:)*dt);
                if i == 1
                    sw(1) = Spcount_f2; % sw = [d(f2) s(f1) s(f2) d(f1)]
                    sw(2) = Spcount_f1;
                else
                    sw(3) = Spcount_f2; % sw = [d(f2) s(f1) s(f2) d(f1)]
                    sw(4) = Spcount_f1;
                end
            end
            CSI(j,k) = (sw(4)+sw(1)-sw(2)-sw(3))/(sw(4)+sw(1)+sw(2)+sw(3));
        end
    end
    figure(1)
    subplot(1,2,1)
    for p = 1:length(Probs_Arr)
        plot(A_Arr, CSI(p,:), color{p})
        hold on
    end
    xlabel('A(Spikes/s)')
    ylabel('CSI')
    ylim([-0.05 1])
    legend('P=10%','P=30%','P=50%')
end


if plot_5_B_right
    clear 
    load('Simulation Results/net1/run_net1_L_P0.1_ISI300.mat')
    color{1} = 'k';
    color{2} = 'b';
    color{3} = 'r';
    f1 = 10;
    f2 = 12;
    CSI = zeros(length(Probs_Arr),length(ISI_Arr));
    for j = 1:length(Probs_Arr)
        for k = 1:length(ISI_Arr)
            sw = zeros(1,4);
            for i = 1:2
                load(['Simulation Results/net1/run_net1_' nev_cond_code{i} '_P' num2str(Probs_Arr(j)) '_ISI' num2str(ISI_Arr(k)*1000) '.mat'], 'E_mean', 'Oddball', 'Stim_Onsets')
                m = 0;
                n = 0;
                E_sum_a_f2 = zeros(1, length(Single_Stim)+0.045/dt);
                E_sum_b_f2 = zeros(1, 0.005/dt);
                E_sum_a_f1 = zeros(1, length(Single_Stim)+0.045/dt);
                E_sum_b_f1 = zeros(1, 0.005/dt);
                for ns = 1:n_stim
                    if Oddball(ns) == f2
                        E_sum_a_f2 = E_sum_a_f2 + E_mean(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1+0.045/dt); % integrating E_mean from stimulus onset to 45 ms after offeset
                        E_sum_b_f2 = E_sum_b_f2 + E_mean(Rec_Column, Stim_Onsets(ns)-0.005/dt:Stim_Onsets(ns)-1); % baseline 5 ms before stimulus onset
                        m = m + 1;
                    else
                        E_sum_a_f1 = E_sum_a_f1 + E_mean(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1+0.045/dt); % integrating E_mean from stimulus onset to 45 ms after offeset
                        E_sum_b_f1 = E_sum_b_f1 + E_mean(Rec_Column, Stim_Onsets(ns)-0.005/dt:Stim_Onsets(ns)-1); % baseline 5 ms before stimulus onset
                        n = n + 1;
                    end
                end
                E_mean_a_f2 = E_sum_a_f2/m;
                E_mean_b_f2 = E_sum_b_f2/m;
                E_mean_a_f1 = E_sum_a_f1/n;
                E_mean_b_f1 = E_sum_b_f1/n;
                Spcount_f2 = sum(E_mean_a_f2(:)*dt) - sum(E_mean_b_f2(:)*dt);
                Spcount_f1 = sum(E_mean_a_f1(:)*dt) - sum(E_mean_b_f1(:)*dt);
                if i == 1
                    sw(1) = Spcount_f2; % sw = [d(f2) s(f1) s(f2) d(f1)]
                    sw(2) = Spcount_f1;
                else
                    sw(3) = Spcount_f2; % sw = [d(f2) s(f1) s(f2) d(f1)]
                    sw(4) = Spcount_f1;
                end
            end
            CSI(j,k) = (sw(4)+sw(1)-sw(2)-sw(3))/(sw(4)+sw(1)+sw(2)+sw(3));
        end
    end
    subplot(1,2,2)
    for p = 1:length(Probs_Arr)
        plot(ISI_Arr, CSI(p,:), color{p})
        hold on
    end
    xlabel('ISI(s)')
    ylabel('CSI')
    ylim([-0.05 1])
end


%% 
if plot_5_C_left
    load('Simulation Results/net1/run_net1_L_F2_A5.mat')
    color{1} = 'k';
    color{2} = 'b';
    color{3} = 'r';
    CSI = zeros(length(Freq_Diff_Arr),length(A_Arr)); %3*20
    for j = 1:length(Freq_Diff_Arr)
        for k = 1:length(A_Arr)
            sw = zeros(1,4);
            for i = 1:2
                load(['Simulation Results/net1/run_net1_' nev_cond_code{i} '_F' num2str(Freq_Diff_Arr(j)) '_A' num2str(A_Arr(k)) '.mat'], 'E_mean', 'Oddball')
                m = 0;
                n = 0;
                E_sum_a_f2 = zeros(1, length(Single_Stim)+0.045/dt);
                E_sum_b_f2 = zeros(1, 0.005/dt);
                E_sum_a_f1 = zeros(1, length(Single_Stim)+0.045/dt);
                E_sum_b_f1 = zeros(1, 0.005/dt);
                for ns = 1:n_stim
                    if Oddball(ns) == Rec_Column + Freq_Diff_Arr(j)/2 % f2 = 12 13 14
                        E_sum_a_f2 = E_sum_a_f2 + E_mean(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1+0.045/dt); % integrating E_mean from stimulus onset to 45 ms after offeset
                        E_sum_b_f2 = E_sum_b_f2 + E_mean(Rec_Column, Stim_Onsets(ns)-0.005/dt:Stim_Onsets(ns)-1); % baseline 5 ms before stimulus onset
                        m = m + 1;
                    else
                        E_sum_a_f1 = E_sum_a_f1 + E_mean(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1+0.045/dt); % integrating E_mean from stimulus onset to 45 ms after offeset
                        E_sum_b_f1 = E_sum_b_f1 + E_mean(Rec_Column, Stim_Onsets(ns)-0.005/dt:Stim_Onsets(ns)-1); % baseline 5 ms before stimulus onset
                        n = n + 1;
                    end
                end
                E_mean_a_f2 = E_sum_a_f2/m;
                E_mean_b_f2 = E_sum_b_f2/m;
                E_mean_a_f1 = E_sum_a_f1/n;
                E_mean_b_f1 = E_sum_b_f1/n;
                Spcount_f2 = sum(E_mean_a_f2(:)*dt) - sum(E_mean_b_f2(:)*dt);
                Spcount_f1 = sum(E_mean_a_f1(:)*dt) - sum(E_mean_b_f1(:)*dt);
                if i == 1
                    sw(1) = Spcount_f2; % sw = [d(f2) s(f1) s(f2) d(f1)]
                    sw(2) = Spcount_f1;
                else
                    sw(3) = Spcount_f2; % sw = [d(f2) s(f1) s(f2) d(f1)]
                    sw(4) = Spcount_f1;
                end
            end
            CSI(j,k) = (sw(4)+sw(1)-sw(2)-sw(3))/(sw(4)+sw(1)+sw(2)+sw(3));
        end
    end
    figure(1)
    subplot(1,2,1)
    for p = 1:length(Freq_Diff_Arr)
        plot(A_Arr, CSI(p,:), color{p})
        hold on
    end
    xlabel('A(Spikes/s)')
    ylabel('CSI')
    ylim([-0.05 1])
    legend('\Deltaf=2','\Deltaf=4','\Deltaf=6')
end


if plot_5_C_right
    clear 
    load('Simulation Results/net1/run_net1_L_F2_A5.mat')
    load('Simulation Results/net1/meta_data.mat')
    color{1} = 'k';
    color{2} = 'b';
    color{3} = 'r';
    CSI = zeros(length(Freq_Diff_Arr),length(ISI_Arr));
    for j = 1:length(Freq_Diff_Arr)
        for k = 1:length(ISI_Arr)
            sw = zeros(1,4);
            for i = 1:2
                load(['Simulation Results/net1/run_net1_' nev_cond_code{i} '_F' num2str(Freq_Diff_Arr(j)) '_ISI' num2str(ISI_Arr(k)*1000) '.mat'], 'E_mean', 'Oddball', 'Stim_Onsets')
                m = 0;
                n = 0;
                E_sum_a_f2 = zeros(1, length(Single_Stim)+0.045/dt);
                E_sum_b_f2 = zeros(1, 0.005/dt);
                E_sum_a_f1 = zeros(1, length(Single_Stim)+0.045/dt);
                E_sum_b_f1 = zeros(1, 0.005/dt);
                for ns = 1:n_stim
                    if Oddball(ns) == Rec_Column + Freq_Diff_Arr(j)/2 % f2 = 12 13 14
                        E_sum_a_f2 = E_sum_a_f2 + E_mean(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1+0.045/dt); % integrating E_mean from stimulus onset to 45 ms after offeset
                        E_sum_b_f2 = E_sum_b_f2 + E_mean(Rec_Column, Stim_Onsets(ns)-0.005/dt:Stim_Onsets(ns)-1); % baseline 5 ms before stimulus onset
                        m = m + 1;
                    else
                        E_sum_a_f1 = E_sum_a_f1 + E_mean(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1+0.045/dt); % integrating E_mean from stimulus onset to 45 ms after offeset
                        E_sum_b_f1 = E_sum_b_f1 + E_mean(Rec_Column, Stim_Onsets(ns)-0.005/dt:Stim_Onsets(ns)-1); % baseline 5 ms before stimulus onset
                        n = n + 1;
                    end
                end
                E_mean_a_f2 = E_sum_a_f2/m;
                E_mean_b_f2 = E_sum_b_f2/m;
                E_mean_a_f1 = E_sum_a_f1/n;
                E_mean_b_f1 = E_sum_b_f1/n;
                Spcount_f2 = sum(E_mean_a_f2(:)*dt) - sum(E_mean_b_f2(:)*dt);
                Spcount_f1 = sum(E_mean_a_f1(:)*dt) - sum(E_mean_b_f1(:)*dt);
                if i == 1
                    sw(1) = Spcount_f2; % sw = [d(f2) s(f1) s(f2) d(f1)]
                    sw(2) = Spcount_f1;
                else
                    sw(3) = Spcount_f2; % sw = [d(f2) s(f1) s(f2) d(f1)]
                    sw(4) = Spcount_f1;
                end
            end
            CSI(j,k) = (sw(4)+sw(1)-sw(2)-sw(3))/(sw(4)+sw(1)+sw(2)+sw(3));
        end
    end
    subplot(1,2,2)
    for p = 1:length(Freq_Diff_Arr)
        plot(ISI_Arr, CSI(p,:), color{p})
        hold on
    end
    xlabel('ISI(s)')
    ylabel('CSI')
    ylim([-0.05 1])
end

%% 
if plot_6_A_D
    load('Simulation Results/net1/run_net1_L_A1_tr1.mat')
    f1 = 10;
    f2 = 12;
    CSI = zeros(1,length(A_Arr)); %1*4
    for k = 1:length(A_Arr)
        sw = zeros(1,4);
        for i = 1:2
            load(['Simulation Results/net1/run_net1_' nev_cond_code{i} '_A' num2str(A_Arr(k)) '_tr1.mat'])
            m = 0;
            n = 0;
            E_sum_a_f2 = zeros(1, length(Single_Stim)+0.045/dt);
            E_sum_b_f2 = zeros(1, 0.005/dt);
            E_sum_a_f1 = zeros(1, length(Single_Stim)+0.045/dt);
            E_sum_b_f1 = zeros(1, 0.005/dt);
            for ns = 1:n_stim
                if Oddball(ns) == f2
                    E_sum_a_f2 = E_sum_a_f2 + E_mean(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1+0.045/dt); % integrating E_mean from stimulus onset to 45 ms after offeset
                    E_sum_b_f2 = E_sum_b_f2 + E_mean(Rec_Column, Stim_Onsets(ns)-0.005/dt:Stim_Onsets(ns)-1); % baseline 5 ms before stimulus onset
                    m = m + 1;
                else
                    E_sum_a_f1 = E_sum_a_f1 + E_mean(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1+0.045/dt); % integrating E_mean from stimulus onset to 45 ms after offeset
                    E_sum_b_f1 = E_sum_b_f1 + E_mean(Rec_Column, Stim_Onsets(ns)-0.005/dt:Stim_Onsets(ns)-1); % baseline 5 ms before stimulus onset
                    n = n + 1;
                end
            end
            E_mean_a_f2 = E_sum_a_f2/m;
            E_mean_b_f2 = E_sum_b_f2/m;
            E_mean_a_f1 = E_sum_a_f1/n;
            E_mean_b_f1 = E_sum_b_f1/n;
            Spcount_f2 = sum(E_mean_a_f2(:)*dt) - sum(E_mean_b_f2(:)*dt);
            Spcount_f1 = sum(E_mean_a_f1(:)*dt) - sum(E_mean_b_f1(:)*dt);
            if i == 1
                sw(1) = Spcount_f2; % sw = [d(f2) s(f1) s(f2) d(f1)]
                sw(2) = Spcount_f1;
            else
                sw(3) = Spcount_f2; % sw = [d(f2) s(f1) s(f2) d(f1)]
                sw(4) = Spcount_f1;
            end
        end
        CSI(1,k) = (sw(4)+sw(1)-sw(2)-sw(3))/(sw(4)+sw(1)+sw(2)+sw(3));
    end
    
    n_stim_plot = 40;
    t_onset = Stim_Onsets(1);
    t_offset = t_onset + floor(n_stim_plot*(duration+ISI)/dt) - 1;
    tim = dt:dt:n_stim_plot*(duration+ISI);
    map = [1 1 1; 0 0 1];
    for j = 1:4
        load(['Simulation Results/net1/run_net1_L_A' num2str(A_Arr(j)) '_tr1.mat'])
        if j == 1
            figure(1)
            subplot(4,1,1)
            plot(tim, E_mean(Rec_Column,t_onset:t_offset), '-k')
            ylabel('E(Spikes/s)')
            ylim([0 8])
            title('No PS')
            c ={ [ 'A = ' num2str(A_Arr(j)) ' Spikes/s'],['CSI = ' num2str(round(CSI(j),3)) ]};
            text(12.5,6.5,c)
            box off
            set(gca, 'YTick', 0:2:8);    

            subplot(4,1,2)
            imagesc([0 n_stim_plot*(duration+ISI)], [10 12], [Spec_Temp(1,t_onset:t_offset,1);Spec_Temp(1,t_onset:t_offset,2)]);
            set(gca,'YDir','normal')   
            ax1 = subplot(4,1,2); 
            colormap(ax1, map)
            axis off
            text(-0.2,10,'f1')
            text(-0.2,12,'f2')
        elseif j == 2
            subplot(4,1,3)
            plot(tim, E_mean(Rec_Column,t_onset:t_offset), '-k')
            ylabel('E(Spikes/s)')
            title('Selective')
            ylim([0 80])
            c ={ [ 'A = ' num2str(A_Arr(j)) ' Spikes/s'],['CSI = ' num2str(round(CSI(j),3)) ]};
            text(12.5,65,c)
            box off
            set(gca, 'YTick', 0:20:80);    
            
            subplot(4,1,4)
            imagesc([0 n_stim_plot*(duration+ISI)], [10 12], [Spec_Temp(1,t_onset:t_offset,1);Spec_Temp(1,t_onset:t_offset,2)]);
            set(gca,'YDir','normal')  
            ax1 = subplot(4,1,4); 
            colormap(ax1, map)
            xlabel('time(s)')
            axis off
            text(-0.2,10,'f1')
            text(-0.2,12,'f2')
        elseif j == 3
            figure(2)
            subplot(4,1,1)
            plot(tim, E_mean(Rec_Column,t_onset:t_offset), '-k')
            ylabel('E(Spikes/s)')
            title('Periodic')
            ylim([0 80])
            c ={ [ 'A = ' num2str(A_Arr(j)) ' Spikes/s'],['CSI = ' num2str(round(CSI(j),3)) ]};
            text(12.5,65,c)
            box off
            set(gca, 'YTick', 0:20:80);    
            
            subplot(4,1,2)
            imagesc([0 n_stim_plot*(duration+ISI)], [10 12], [Spec_Temp(1,t_onset:t_offset,1);Spec_Temp(1,t_onset:t_offset,2)]);
            set(gca,'YDir','normal')  
            ax1 = subplot(4,1,2); 
            colormap(ax1, map)
            axis off
            text(-0.2,10,'f1')
            text(-0.2,12,'f2')
        else
            subplot(4,1,3)
            plot(tim, E_mean(Rec_Column,t_onset:t_offset), '-k')
            ylabel('E(Spikes/s)')
            title('Reliable')
            ylim([0 80])
            c ={ [ 'A = ' num2str(A_Arr(j)) ' Spikes/s'],['CSI = ' num2str(round(CSI(j),3)) ]};
            text(12.5,65,c)
            box off
            set(gca, 'YTick', 0:20:80);    
            
            subplot(4,1,4)
            imagesc([0 n_stim_plot*(duration+ISI)], [10 12], [Spec_Temp(1,t_onset:t_offset,1);Spec_Temp(1,t_onset:t_offset,2)]);
            set(gca,'YDir','normal')   
            ax1 = subplot(4,1,4); 
            colormap(ax1, map)
            xlabel('time(s)') 
            axis off
            text(-0.2,10,'f1')
            text(-0.2,12,'f2')
        end
    end  
end

            
%%           
if plot_7_A
    load('Simulation Results/net1/run_net1_L_P0.1_A5.mat')%%%%%%%%%%%%%%%%
    f1 = 10;
    f2 = 12;
    CSI = zeros(length(ISI_Arr),length(A_Arr)); %46*35
    for j = 1:length(ISI_Arr)
        for k = 1:length(A_Arr)
            sw = zeros(1,4);
            for i = 1:2
                load(['Simulation Results/net1/run_net1_' nev_cond_code{i} '_ISI' num2str(ISI_Arr(j)) '_A' num2str(A_Arr(k)) '.mat'], 'E_mean', 'Oddball')
                m = 0;
                n = 0;
                E_sum_a_f2 = zeros(1, length(Single_Stim)+0.045/dt);
                E_sum_b_f2 = zeros(1, 0.005/dt);
                E_sum_a_f1 = zeros(1, length(Single_Stim)+0.045/dt);
                E_sum_b_f1 = zeros(1, 0.005/dt);
                for ns = 1:n_stim
                    if Oddball(ns) == f2
                        E_sum_a_f2 = E_sum_a_f2 + E_mean(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1+0.045/dt); % integrating E_mean from stimulus onset to 45 ms after offeset
                        E_sum_b_f2 = E_sum_b_f2 + E_mean(Rec_Column, Stim_Onsets(ns)-0.005/dt:Stim_Onsets(ns)-1); % baseline 5 ms before stimulus onset
                        m = m + 1;
                    else
                        E_sum_a_f1 = E_sum_a_f1 + E_mean(Rec_Column, Stim_Onsets(ns):Stim_Onsets(ns)+length(Single_Stim)-1+0.045/dt); % integrating E_mean from stimulus onset to 45 ms after offeset
                        E_sum_b_f1 = E_sum_b_f1 + E_mean(Rec_Column, Stim_Onsets(ns)-0.005/dt:Stim_Onsets(ns)-1); % baseline 5 ms before stimulus onset
                        n = n + 1;
                    end
                end
                E_mean_a_f2 = E_sum_a_f2/m;
                E_mean_b_f2 = E_sum_b_f2/m;
                E_mean_a_f1 = E_sum_a_f1/n;
                E_mean_b_f1 = E_sum_b_f1/n;
                Spcount_f2 = sum(E_mean_a_f2(:)*dt) - sum(E_mean_b_f2(:)*dt);
                Spcount_f1 = sum(E_mean_a_f1(:)*dt) - sum(E_mean_b_f1(:)*dt);
                if i == 1
                    sw(1) = Spcount_f2; % sw = [d(f2) s(f1) s(f2) d(f1)]
                    sw(2) = Spcount_f1;
                else
                    sw(3) = Spcount_f2; % sw = [d(f2) s(f1) s(f2) d(f1)]
                    sw(4) = Spcount_f1;
                end
            end
            CSI(j,k) = (sw(4)+sw(1)-sw(2)-sw(3))/(sw(4)+sw(1)+sw(2)+sw(3));
        end
    end
    
    figure(1)
    imagesc([0 35], [0.1 1], CSI(:,:));
    set(gca,'YDir','normal')
    xlabel('A(Spikes/s)')
    ylabel('ISI(s)')  
    colormap jet
    c = colorbar;
    c.Label.String = 'CSI';
    
end


%% Fig1.D
if plot_1_D 
    load('Simulation Results/net1/run_net1_L.mat')
    load('Simulation Results/net1/meta_data.mat')
    
    nev_cond_code{1} = 'L';
    nev_cond_code{2} = 'H';
    for ii = 1 
        save('counter.mat','ii','nev_cond_code')
        clear all
        load('counter.mat')
        load(['Simulation Results/net1/run_net1_' nev_cond_code{ii} '.mat'])
        j = 0;
        k = 0;
        E_sum2 = zeros(1,1,floor((duration+ISI)/dt)+1);
        E_sum1 = E_sum2;
        for ns = 1:n_stim-1  
            if Spec_Temp(1,2,Stim_Onsets(ns)) > 0
                E_sum2 = E_sum2 + E_act_overall(PW(1), PW(2), Stim_Onsets(ns):Stim_Onsets(ns+1));
                j = j + 1;
            else
                E_sum1 = E_sum1 + E_act_overall(AW1(1), AW1(2), Stim_Onsets(ns):Stim_Onsets(ns+1));
                k = k + 1;
            end
        end
        E_mean2 = E_sum2/j;
        E_mean1 = E_sum1/k;

        if ii == 1
            AXES_FONTSIZE = 9;
            LineWidth = 0.8;
            COLORBAR_FONTSIZE = 9;
            t = Stim_Onsets(1);

            figure
            tim = (0:dt:5*duration)*10^3;
            hold on
            for i = 1:9
                plot(tim, E_plot(i,t:t+5*duration/dt))
            end
            xlabel('time(ms)')
            ylabel('Column')  

                       
            figure
            subplot(1,2,1)
            imagesc([0 5*(duration)*10^3], [1 9], E_plot(:,t:t+5*(duration)/dt));
            set(gca,'YDir','normal')
            xlabel('time(ms)')
            ylabel('Column')  
            colormap jet
            colorbar
            set(gca,'FontSize',AXES_FONTSIZE)
            c = colorbar;
            set(c,'FontSize',COLORBAR_FONTSIZE);
            c.Label.String = 'E (Spikes/s)';
            axis square
            
            subplot(1,2,2)     
%             tim = dt:dt:(duration+ISI);
            tim = 0:dt:(duration+ISI);
            tim = tim*1000;
            plot(tim, reshape(E_mean2, [1, floor((duration+ISI)/dt)+1]), '-r', 'LineWidth', LineWidth)
            hold on
            plot(tim, reshape(E_mean1, [1, floor((duration+ISI)/dt)+1]), '-b', 'LineWidth', LineWidth)
            set(gca,'FontSize',AXES_FONTSIZE)
            xlabel('time(ms)')
            ylabel('E(Spikes/s)')
            xlim([0 2*duration*10^3])
            legend('dev','std')
            legend('boxoff')
%             legend('f2(deviant)','f1(standard)')
            pos = get(gcf,'position');
            set(gcf,'position',[pos(1),pos(2),pos(3),(0.6)*pos(4)]);
            axis square
            
              
            figure
            for i = 1:11
                subplot(4,3,i)
                image([1 4], [1 5], E_act_overall(:,:,t+(i-1)*0.005/dt));
                set(gca, 'YTick', [1:5], 'XTick', [1:4]);
                title(['t = ',num2str((i-1)*5),'ms']);
                colormap hot
%                 set(gca,'FontSize',AXES_FONTSIZE)
                axis equal
                axis tight
                c = colorbar;
%                 set(c,'FontSize',COLORBAR_FONTSIZE);
                c.Label.String = 'E (Spikes/s)';
                c.Limits = [0 90];
                c.Ticks = [0 30 60 90];
                c.TickLabels = [0 30 60 90];
%                 if i < 12
%                     colorbar('off') 
%                 end
            end
            
      

            
            
        end
    end
end


%% True deviance detection
if 1
    jj = 1;
    barplot = zeros(1,6);
    nev_cond_code{1} = 'L';
    nev_cond_code{2} = 'H';
    nev_cond_code{3} = 'E';
    nev_cond_code{4} = 'MS';
    nev_cond_code{5} = 'DA1';
    nev_cond_code{6} = 'DA2';
    color{1} = 'r';
    color{2} = 'b';
    color{3} = 'k';
    color{4} = 'g';
    color{5} = 'c';
    color{6} = 'm';
    colorb{1} = [1 0 0];
    colorb{2} = [0 0 1];
    colorb{3} = [0 0 0];
    colorb{4} = [0 1 0];
    colorb{5} = [0 1 1];
    colorb{6} = [1 0 1];
    N_stim = 2;
    for ii = [1 4] % Choose specific protocols
        save('counter.mat','ii','jj','nev_cond_code','color','colorb','N_stim','barplot')
        clear all
        AXES_FONTSIZE = 10;
        LineWidth = 1;
        COLORBAR_FONTSIZE = 10;
        load('counter.mat')
        load(['Simulation Results/net1/meta_data.mat'])
        load(['Simulation Results/net1/run_net1_' nev_cond_code{ii} '.mat'])
        E_sum = zeros(1,1,N_stim*length(Single_Stim));
        j = 0;
        for trnum = 1:num_trials
            load(['Simulation Results/net1/run_net1_' nev_cond_code{ii} '.mat'], 'E_act_overall', 'Oddball')
            for ns = 1:n_stim  
                if Oddball(ns,:) == PW
                    E_sum = E_sum + E_act_overall(PW(1), PW(2), Stim_Onsets(ns):Stim_Onsets(ns)+N_stim*length(Single_Stim)-1);
                    j = j + 1;
                end
            end
        end
        E_mean =  E_sum/j;
        figure(100) % Plot oddball and many-standard deviant
        tim = dt:dt:N_stim*length(Single_Stim)*dt;
        tim = tim*1000;
        plot(tim, reshape(E_mean, [1, N_stim*length(Single_Stim)]), color{ii}, 'LineWidth', LineWidth)
        hold on
        E_mean_bar = sum(E_mean(:)*dt);
        barplot(jj) = E_mean_bar;
        jj = jj + 1;
    end
    pos = get(gcf,'position');
    set(gcf,'position',[pos(1),pos(2),(0.6)*pos(3),(0.6)*pos(4)]);
    axis square
    xlabel('time(ms)')
    ylabel('E(Spikes/s)')
    legend('Many-standards','Oddball')
    legend('boxoff')
    title('f2')
end

