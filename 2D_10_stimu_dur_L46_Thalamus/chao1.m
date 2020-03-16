close all
clear all

L = 1;
DB = 0;

AXES_FONTSIZE = 11;
LineWidth = 0.8;
COLORBAR_FONTSIZE = 10;

%% Low Condition
if L
    load('Simulation Results/net1/run_net1_L.mat')
    load('Simulation Results/net1/meta_data.mat')
    
    figure(1)
    n_stim_plot = 10;
    onset = Stim_Onsets(1);
    offset = Stim_Onsets(1) + (duration + ISI) * n_stim_plot / dt -1;
    xaxis = dt:dt:(duration + ISI) * n_stim_plot;

    % Subplot 1
    s1 = subplot(4,1,1); 
    plot(xaxis,E_mean(11,onset:offset),'b','LineWidth',LineWidth); %plot on subplot
    pos = get(s1,'position');
    set(s1,'position',[0.0602,pos(2),0.8448,pos(4)]);
    xlim([0 n_stim_plot*(duration+ISI)])
    ylim([0 100])
    set(s1, 'XTickLabel',[],'FontSize',AXES_FONTSIZE);
    ylabel('E(Spikes/s)')

    % Subplot 2
    s2 = subplot(4,1,2); 
    imagesc([0 n_stim_plot*(duration+ISI)], [Rec_Column - Freq_Diff/2 Rec_Column + Freq_Diff/2],...
        [Spec_Temp(1,onset:offset,1);Spec_Temp(1,onset:offset,2)]);
    pos = get(s2,'position');
    set(s2,'position',[0.0602,0.6349,0.8448,(1/2)*pos(4)]);
    set(s2,'YDir','normal')   
    map = [1 1 1; 0 0 1];
    colormap(s2, map)
    set(s2, 'XTickLabel',[],'FontSize',AXES_FONTSIZE);
    set(s2, 'YTickLabel',{'f1','f2'});

    % Subplot 3
    s3 = subplot(4,1,3);
    imagesc([0 n_stim_plot*(duration+ISI)], [7 15], E_mean(7:15,onset:offset));
    pos = get(s3,'position');
    set(s3,'position',[0.0602,0.3715,0.8448,pos(4)+(1/3)*pos(4)]);
    set(s3,'YDir','normal');
    set(s3, 'YTick', [7 10 11 12 15],'FontSize',AXES_FONTSIZE);
    ylabel('Column');
    set(s3, 'XTickLabel',[]);
    colormap(s3, hot);
    c = colorbar;
    set(c,'FontSize',COLORBAR_FONTSIZE);
    c.Label.String = 'E (Spikes/s)';
    c.Limits = [0 80];
    c.TickLabels = [0 20 40 60 80];
    c.Position = [0.922 0.3715 0.021 pos(4)+(1/3)*pos(4)];

    % Subplot 4
    s4 = subplot(4,1,4) ; 
    imagesc([0 n_stim_plot*(duration+ISI)], [7 15], x_mean(7:15,onset:offset));
    pos = get(s4,'position');
    set(s4,'position',[0.0602,pos(2),0.8448,pos(4)+(1/3)*pos(4)]);
    set(s4,'YDir','normal') 
    set(s4, 'YTick', [7 10 11 12 15],'FontSize',AXES_FONTSIZE);
    ylabel('Column');
    xlabel('time(s)');
    colormap(s4, hot);
    c = colorbar;
    set(c,'FontSize',COLORBAR_FONTSIZE);
    c.Label.String = 'x';
    c.Limits = [0.35 0.7];
    c.TickLabels = [0.4 0.5 0.6 0.7];
    c.Position = [0.922 pos(2) 0.021 pos(4)+(1/3)*pos(4)]; 
     
    
    figure(2)  
    s1 = subplot(4,1,1) 
    plot(xaxis, E_mean(Rec_Column-Freq_Diff/2, onset:offset), '-b','LineWidth',LineWidth)
    hold on
    plot(xaxis, E_mean(Rec_Column, onset:offset), '-k','LineWidth',LineWidth)
    hold on
    plot(xaxis, E_mean(Rec_Column+Freq_Diff/2, onset:offset), '-r','LineWidth',LineWidth)
    ylabel('E(Spikes/s)')
    ylim([0 100])
    xlim([0 n_stim_plot*(duration+ISI)])
    pos = get(s1,'position');
    set(s1,'position',[0.0602    pos(2)    0.8448    pos(4)],'FontSize',AXES_FONTSIZE,'XTickLabel',[]);
    
    s2 = subplot(4,1,2) 
    plot(xaxis, x_mean(Rec_Column-Freq_Diff/2, onset:offset), '-b','LineWidth',LineWidth)
    hold on
    plot(xaxis, x_mean(Rec_Column, onset:offset), '-k','LineWidth',LineWidth)
    hold on
    plot(xaxis, x_mean(Rec_Column+Freq_Diff/2, onset:offset), '-r','LineWidth',LineWidth)
    ylabel('x')
    xlim([0 n_stim_plot*(duration+ISI)])
    ylim([0.2 0.8])
    legend(['Col' num2str(Rec_Column-Freq_Diff/2)],['Col' num2str(Rec_Column)],['Col' num2str(Rec_Column+Freq_Diff/2)])
    pos = get(s2,'position');
    set(s2,'position',[0.0602    pos(2)    0.8448    pos(4)],'YTick', [0.2:0.2:0.8],'FontSize',AXES_FONTSIZE,'XTickLabel',[])
    
    s3 = subplot(4,1,3)
    plot(xaxis, reshape(E_act_overall(Rec_Column,1,onset:offset),[1,length(xaxis)]), 'Color',[0.66 0 0.94],'LineWidth',LineWidth) 
    hold on;
    plot(xaxis, reshape(E_act_overall(Rec_Column,2,onset:offset),[1,length(xaxis)]), 'Color', [1 0.56 0],'LineWidth',LineWidth) 
    hold on;
    xlim([0 n_stim_plot*(duration+ISI)])
    ylim([0 100])
    ylabel('E(Spikes/s)')
    pos = get(s3,'position');
    set(s3,'position',[0.0602    pos(2)    0.8448    pos(4)],'FontSize',AXES_FONTSIZE,'XTickLabel',[])

    s4 = subplot(4,1,4)
    plot(xaxis, reshape(x_act_overall(Rec_Column,1,onset:offset),[1,length(xaxis)]), 'Color', [0.66 0 0.94],'LineWidth',LineWidth) 
    hold on;
    plot(xaxis, reshape(x_act_overall(Rec_Column,2,onset:offset),[1,length(xaxis)]),'Color', [1 0.56 0],'LineWidth',LineWidth) 
    hold on;
    xlim([0 n_stim_plot*(duration+ISI)])
    xlabel('time(s)')
    ylabel('x')
    legend('Neus with Inh Input in Col 11','Neus with exc input in Col 11')
    pos = get(s4,'position');
    set(s4,'position',[0.0602    pos(2)    0.8448    pos(4)],'FontSize',AXES_FONTSIZE)
    
    
%     figure(3)  
%     s1 = subplot(3,1,1); 
%     imagesc([0 n_stim_plot*(duration+ISI)], [Rec_Column - Freq_Diff/2 Rec_Column + Freq_Diff/2],...
%     [Spec_Temp(1,onset:offset,1);Spec_Temp(1,onset:offset,2)]);
%     pos = get(s1,'position');
%     set(s1,'position',[0.0602    0.8329    0.8448    0.1079]);
%     set(s1,'YDir','normal')   
%     map = [1 1 1; 0 0 1];
%     colormap(s1, map)
%     set(s1, 'XTickLabel',[],'FontSize',AXES_FONTSIZE);
%     set(s1, 'YTickLabel',{'f1','f2'});
%     
%     s2 = subplot(3,1,2) 
%     plot(xaxis, E_mean(Rec_Column-Freq_Diff/2, onset:offset), '-b','LineWidth',LineWidth)
%     hold on
%     plot(xaxis, E_mean(Rec_Column, onset:offset), '-k','LineWidth',LineWidth)
%     hold on
%     plot(xaxis, E_mean(Rec_Column+Freq_Diff/2, onset:offset), '-r','LineWidth',LineWidth)
%     ylabel('E(Spikes/s)')
%     ylim([0 100])
%     xlim([0 n_stim_plot*(duration+ISI)])
%     pos = get(s2,'position');
%     set(s2,'position',[0.0602    0.4686    0.8448    0.2877],'FontSize',AXES_FONTSIZE,'XTickLabel',[]);
%     
%     s3 = subplot(3,1,3) 
%     plot(xaxis, x_mean(Rec_Column-Freq_Diff/2, onset:offset), '-b','LineWidth',LineWidth)
%     hold on
%     plot(xaxis, x_mean(Rec_Column, onset:offset), '-k','LineWidth',LineWidth)
%     hold on
%     plot(xaxis, x_mean(Rec_Column+Freq_Diff/2, onset:offset), '-r','LineWidth',LineWidth)
%     ylabel('x')
%     xlim([0 n_stim_plot*(duration+ISI)])
%     legend(['Col' num2str(Rec_Column-Freq_Diff/2)],['Col' num2str(Rec_Column)],['Col' num2str(Rec_Column+Freq_Diff/2)])
%     pos = get(s3,'position');
%     set(s3,'position',[0.0602    0.1100    0.8448    0.2877],'FontSize',AXES_FONTSIZE)
%     xlabel('time(s)');
        
%     tim = xaxis;
%     t_onset = onset;
%     t_offset = offset;
%     figure(3)
%     subplot(6,1,1)
%     plot(tim, E_mean(Rec_Column-1,t_onset:t_offset), '-b')
%     hold on
%     plot(tim, reshape(E_act_overall(Rec_Column-1,1,t_onset:t_offset),[1,length(tim)]), 'm') 
%     hold on;
%     plot(tim, reshape(E_act_overall(Rec_Column-1,2,t_onset:t_offset),[1,length(tim)]), 'c') 
%     hold on;
%     xlim([0 n_stim_plot*(duration+ISI)])
%     ylim([0 70])
%     ylabel('E(Spikes/s)')
%     subplot(6,1,2)
%     plot(tim, x_mean(Rec_Column-1,t_onset:t_offset), '-b')
%     hold on
%     plot(tim, reshape(x_act_overall(Rec_Column-1,1,t_onset:t_offset),[1,length(tim)]), 'm') 
%     hold on;
%     plot(tim, reshape(x_act_overall(Rec_Column-1,2,t_onset:t_offset),[1,length(tim)]), 'c') 
%     hold on;
%     xlim([0 n_stim_plot*(duration+ISI)])
%     xlabel('time(s)')
%     ylabel('x')
%     legend('Mean value','Neus with inh input','Neus with exc input')
%     
%     subplot(6,1,3)
%     plot(tim, E_mean(Rec_Column,t_onset:t_offset), '-b')
%     hold on
%     plot(tim, reshape(E_act_overall(Rec_Column,1,t_onset:t_offset),[1,length(tim)]), 'm') 
%     hold on;
%     plot(tim, reshape(E_act_overall(Rec_Column,2,t_onset:t_offset),[1,length(tim)]), 'c') 
%     hold on;
%     xlim([0 n_stim_plot*(duration+ISI)])
%     ylim([0 70])
%     ylabel('E(Spikes/s)')
%     subplot(6,1,4)
%     plot(tim, x_mean(Rec_Column,t_onset:t_offset), '-b')
%     hold on
%     plot(tim, reshape(x_act_overall(Rec_Column,1,t_onset:t_offset),[1,length(tim)]), 'm') 
%     hold on;
%     plot(tim, reshape(x_act_overall(Rec_Column,2,t_onset:t_offset),[1,length(tim)]), 'c') 
%     hold on;
%     xlim([0 n_stim_plot*(duration+ISI)])
%     xlabel('time(s)')
%     ylabel('x')
%     legend('Mean value','Neus with inh input','Neus with exc input')
%         
%     subplot(6,1,5)
%     plot(tim, E_mean(Rec_Column+1,t_onset:t_offset), '-b')
%     hold on
%     plot(tim, reshape(E_act_overall(Rec_Column+1,1,t_onset:t_offset),[1,length(tim)]), 'm') 
%     hold on;
%     plot(tim, reshape(E_act_overall(Rec_Column+1,2,t_onset:t_offset),[1,length(tim)]), 'c') 
%     hold on;
%     xlim([0 n_stim_plot*(duration+ISI)])
%     ylim([0 70])
%     ylabel('E(Spikes/s)')
%     subplot(6,1,6)
%     plot(tim, x_mean(Rec_Column+1,t_onset:t_offset), '-b')
%     hold on
%     plot(tim, reshape(x_act_overall(Rec_Column+1,1,t_onset:t_offset),[1,length(tim)]), 'm') 
%     hold on;
%     plot(tim, reshape(x_act_overall(Rec_Column+1,2,t_onset:t_offset),[1,length(tim)]), 'c') 
%     hold on;
%     xlim([0 n_stim_plot*(duration+ISI)])
%     xlabel('time(s)')
%     ylabel('x')
%     legend('Mean value','Neus with inh input','Neus with exc input')

end

%% Many-standards Condition
if DB
    load('Simulation Results/net1/run_net1_DB.mat')
    load('Simulation Results/net1/meta_data.mat')
    figure(3)
    n_stim_plot = 20;
    onset = Stim_Onsets(1);
    offset = Stim_Onsets(1) + (duration + ISI) * n_stim_plot / dt -1;
    xaxis = dt:dt:(duration + ISI) * n_stim_plot;

    % Subplot 1
    s1 = subplot(4,1,1); 
    plot(xaxis,E_mean(11,onset:offset),'b','LineWidth',LineWidth); %plot on subplot
    pos = get(s1,'position');
    set(s1,'position',[0.0602, 0.8610, 0.8448, 0.1046]);
    xlim([0 n_stim_plot*(duration+ISI)])
    ylim([0 100])
    set(s1, 'XTickLabel',[],'FontSize',AXES_FONTSIZE);
    ylabel('E(Spikes/s)')

    % Subplot 2
    s2 = subplot(4,1,2); 
    imagesc([0 n_stim_plot*(duration+ISI)], [2 20], reshape(Spec_Temp(1,onset:offset,:),length(onset:offset),10)');
    pos = get(s2,'position');
    set(s2,'position',[0.0602, 0.6267, 0.8448, 0.2103]);
    set(s2,'YDir','normal')   
    map = [1 1 1; 0 0 1];
    colormap(s2, map)
    set(s2, 'YTick', [2:2:20], 'XTickLabel', [], 'FontSize', AXES_FONTSIZE);
%     set(s2, 'YTickLabel',{'f1','f2'});

    % Subplot 3
    s3 = subplot(4,1,3);
    imagesc([0 n_stim_plot*(duration+ISI)], [1 21], E_mean(:,onset:offset));
    pos = get(s3,'position');
    set(s3,'position',[0.0602, 0.3542, 0.8448, 0.2460]);
    set(s3,'YDir','normal');
    set(s3, 'YTick', [3:4:19], 'FontSize', AXES_FONTSIZE);
    ylabel('Column');
    set(s3, 'XTickLabel',[]);
    colormap(s3, hot);
    c = colorbar;
    set(c,'FontSize',COLORBAR_FONTSIZE);
    c.Label.String = 'E (Spikes/s)';
    c.Limits = [0 80];
    c.TickLabels = [0 20 40 60 80];
    c.Position = [0.922 0.3542 0.021 0.2460];

    % Subplot 4
    s4 = subplot(4,1,4) ; 
    imagesc([0 n_stim_plot*(duration+ISI)], [1 21], x_mean(:,onset:offset));
    pos = get(s4,'position');
    set(s4,'position',[0.0602, 0.0787, 0.8448, 0.2460]);
    set(s4,'YDir','normal') 
    set(s4, 'YTick', [3:4:19],'FontSize',AXES_FONTSIZE);
    ylabel('Column');
    xlabel('time(s)');
    colormap(s4, hot);
    c = colorbar;
    set(c,'FontSize',COLORBAR_FONTSIZE);
    c.Label.String = 'x';
    c.Limits = [0.35 0.7];
    c.TickLabels = [0.4 0.5 0.6 0.7];
    c.Position = [0.922 0.0787 0.021 0.2460];
end


