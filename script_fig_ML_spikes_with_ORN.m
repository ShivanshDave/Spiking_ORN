%% R99.F2A
SpikeEN = flip([0;1]); 
PULSE.ton = [0;0];
PULSE.toff = [1;1];
PULSE.conc = [20;20];
PULSE.tspan = [-0.001 0.6];
DATA = simulate_ORN(PULSE,SpikeEN);

%%
plt.Lwd = 1.2;
plt.FTsz = 14;
plt.Xoff = 0.02;
plt.FGpos = [10 10 900 400];
plt.scale = [3 7 7 7 1];
plt.ytick = [-60,0,40];
plt.xtick = [0:0.1:0.7];
plt.fname = '.\Report\figs\v1\fig_ML_spikes_with_ORN.png';

plot_pulse_currents_overlap(plt,DATA)
%%
function plot_pulse_currents_overlap(plt,D)

    figure('Renderer', 'painters', 'Position', plt.FGpos);
    plt.t = tiledlayout(sum(plt.scale),1,'TileSpacing','compact','Padding','compact');
    plt.X = [D.PULSE.tspan(1)-plt.Xoff, D.PULSE.tspan(2)];
    plt.rot = 0;
    
    %-- Stim --
    nexttile([plt.scale(1) 1])
    TT = linspace(D.T(1),D.T(end),100);
    OD = simulate_pulse_train(TT,D.PULSE.ton,D.PULSE.toff,D.PULSE.conc);
    OD = OD(end,:);
    plt.ax1 = plot(TT,OD,'k-','LineWidth',plt.Lwd);
    ylabel({'O_{stim}','(uM)'})
    set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
        'YTick', [0 D.PULSE.conc(1)],'YTickLabel',{'  ','20'},'YTickLabelRotation',plt.rot,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')
    
    % Spikes
    nexttile([plt.scale(2) 1])
    hold on
    for k = 1:size(D.PULSE.ton,1)               
        plt.axl(k) = plot(D.T,real(D.PRED.spkV(:,k)),'-','LineWidth',plt.Lwd);
    end
    xlabel('Time (sec)')
    ylabel({'V_{ML}','(mV)'})
%     lgd = legend({'ORN+Spikes','ORN'},'Location','best');
%     title(lgd,'Conc. (uM)')    
    set(gca,'XLim',plt.X,'XColor','none','XTick',[],'XTickLabel',[],...
        'YLim',[plt.ytick(1) plt.ytick(end)],'YTick',plt.ytick,...
        'tickdir', 'out','FontSize',plt.FTsz,'YTickLabelRotation',plt.rot,...
        'color','none','box','off')
    lgd = legend(plt.axl(:),{'ORN+Spikes','ORN','Spk-ID','Spk-ID'},'Location','best');
    set(lgd,'Box','off')
    
    % ornV
    nexttile([plt.scale(3) 1])
    hold on
    for k = 1:size(D.PULSE.ton,1)               
        plt.axl(k) = plot(D.T,real(D.PRED.ornV(:,k)),'-','LineWidth',plt.Lwd);
    end
    xlabel('Time (sec)')
    ylabel({'V_{ORN}','(mV)'})
%     lgd = legend({num2str(D.PULSE.conc)},'Location','best');
%     title(lgd,'Conc. (uM)')    
    set(gca,'XLim',plt.X,'XColor','none','XTick',[],'XTickLabel',[],...
        'YLim',[-50,-30],'YTick',[-50,-40,-30],'YTickLabelRotation',plt.rot,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box','off')

    
    % Im
    nexttile([plt.scale(4) 1])
    hold on
    for k = 1:size(D.PULSE.ton,1)               
        plt.axl(k) = plot(D.T,real(D.PRED.Im(:,k)),'-','LineWidth',plt.Lwd);
    end
    xlabel('Time (sec)')
    ylabel({'I_{ORN}','(pA)'})
    spikes = script_spikes_ID(real(D.PRED.spkV),D.T,0);
    dind=8; set(gca,'ColorOrderIndex',1);
    i=1; plt.sc = scatter(D.T(spikes{i,1}-dind), 1*real(D.PRED.Im(spikes{i,1}-dind,i)), 45, 'om');
%     i=2; scatter(D.T(spikes{i,1}-dind), 1*real(D.PRED.Im(spikes{i,1}-dind,i)), 45, 'om')
    lgd = legend(plt.sc, {'Identifying spikes'},'Location','southwest');
    set(lgd,'Box','off')
    
%     lgd = legend({num2str(D.PULSE.conc)},'Location','best');
%     title(lgd,'Conc. (uM)')    
    set(gca,'XLim',plt.X,'XColor','none','XTick',[],'XTickLabel',[],...
        'YLim',[-60,20],'YTick',[-60,-20,20],'YTickLabelRotation',plt.rot,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box','off')

    % TIME axis
    nexttile([plt.scale(5) 1])
    axis();
    xlabel('Time (sec)')
    set(gca,'XLim',plt.X,'XTick',plt.xtick,...
        'YColor','none','YTick', [], 'YTickLabel', [],...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')
    
    exportgraphics(gcf,plt.fname,'Resolution',300)
end