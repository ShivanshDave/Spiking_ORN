%% R99.F2A
SpikeEN = 1; 
PULSE.ton = [0];
PULSE.toff = [1];
PULSE.conc = [20];
PULSE.tspan = [-0.001 0.6];
DATA = simulate_ORN(PULSE,SpikeEN);

%%
plt.Lwd = 1.2;
plt.FTsz = 16;
plt.Xoff = 0.02;
plt.FGpos = [10 10 900 500];
plt.scale = [3 12 4 4 1];
plt.ytick = [-55,-25,0,35];
plt.xtick = [0:0.1:0.7];
plt.fname = '.\Report\figs\fig_ML_spikes.png';

plot_pulse_currents_overlap(plt,DATA)
%%
function plot_pulse_currents_overlap(plt,D)

    figure('Renderer', 'painters', 'Position', plt.FGpos);
    plt.t = tiledlayout(sum(plt.scale),1,'TileSpacing','loose','Padding','compact');
    plt.X = [D.PULSE.tspan(1)-plt.Xoff, D.PULSE.tspan(2)];
    
    %-- Stim --
    nexttile([plt.scale(1) 1])
    TT = linspace(D.T(1),D.T(end),100);
    OD = simulate_pulse_train(TT,D.PULSE.ton,D.PULSE.toff,D.PULSE.conc);
    OD = OD(end,:);
    plt.ax1 = plot(TT,OD,'k-','LineWidth',plt.Lwd);
    ylabel('Stim')
    set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
        'YTick', [0 D.PULSE.conc(1)],'YTickLabel',{'','20'},'YTickLabelRotation',90,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')
    
    % Spikes
    nexttile([plt.scale(2) 1])
    hold on
    for k = 1:size(D.PULSE.ton,1)               
        plt.axl(k) = plot(D.T,real(D.PRED.spkV(:,k)),'-','LineWidth',plt.Lwd);
    end
    xlabel('Time (sec)')
    ylabel('ML Vspk(mV)')
%     lgd = legend({num2str(D.PULSE.conc)},'Location','best');
%     title(lgd,'Conc. (uM)')    
    set(gca,'XLim',plt.X,'XColor','none','XTick',[],'XTickLabel',[],...
        'YLim',[plt.ytick(1) plt.ytick(end)],'YTick',plt.ytick,...
        'tickdir', 'out','FontSize',plt.FTsz,'YTickLabelRotation',90,...
        'color','none','box','off','ColorOrder',[0.4940 0.1840 0.5560])
    
    % nK
    nexttile([plt.scale(3) 1])
    hold on
    for k = 1:size(D.PULSE.ton,1)               
        plt.axl(k) = plot(D.T,real(D.PRED.nK(:,k)),'-','LineWidth',plt.Lwd);
    end
    xlabel('Time (sec)')
    ylabel('nK')
%     lgd = legend({num2str(D.PULSE.conc)},'Location','best');
%     title(lgd,'Conc. (uM)')    
    set(gca,'XLim',plt.X,'XColor','none','XTick',[],'XTickLabel',[],...
        'YLim',[0,0.6],'YTick',[0,0.6],'YTickLabelRotation',90,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box','off','ColorOrder',hsv(8))
    
    % Ca / CaFR
    nexttile([plt.scale(4) 1])
    hold on
    yyaxis left
    for k = 1:size(D.PULSE.ton,1)               
        plt.axl(k) = plot(D.T,real(D.PRED.CaFR(:,k)),'-','LineWidth',plt.Lwd);
    end
    xlabel('Time (sec)')
    ylabel('CaFR')
%     lgd = legend({num2str(D.PULSE.conc)},'Location','best');
%     title(lgd,'Conc. (uM)')    
    set(gca,'XLim',plt.X,'XColor','none','XTick',[],'XTickLabel',[],...
        'YLim',[1 5],'YTick',[1 5],'YTickLabelRotation',90,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box','off')
    
    yyaxis right
    plot(D.T,real(D.PRED.Ca(:,k)),'--','LineWidth',plt.Lwd);
    ylabel('Ca')
    set(gca,'XLim',plt.X,'XColor','none','XTick',[],'XTickLabel',[],...
        'YLim',[0 7],'YTick',[0 7],'YTickLabelRotation',90,...
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