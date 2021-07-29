%% R99.X
SpikeEN = 0; plt.N = 8;
PULSE.ton = 0.000*ones(plt.N,1);
PULSE.toff = [10,5,2,1, 0.5,0.25,0.1,0.05]';
PULSE.conc = 20*ones(plt.N,1);
PULSE.tspan = [-1 12];
DATA = simulate_ORN(PULSE,SpikeEN);

%%
plt.Lwd = 1.2;
plt.FTsz = 16;
plt.Xoff = 0.25;
plt.FGpos = [10 10 900 500];
plt.scale = [3 12 1];
plt.ytick = [-50,-25,0];
plt.xtick = [-1,0,1,2,5,10];
plt.fname = '.\Report\figs\fig_txn_compare_dur.png';

plot_pulse_currents_overlap(plt,DATA)
%%
function plot_pulse_currents_overlap(plt,D)

    figure('Renderer', 'painters', 'Position', plt.FGpos);
    plt.t = tiledlayout(sum(plt.scale),1,'TileSpacing','none','Padding','compact');
    plt.X = [D.PULSE.tspan(1)-plt.Xoff, D.PULSE.tspan(2)];
    
    nexttile([plt.scale(1) 1])
    TT = linspace(D.T(1),D.T(end),100);
    OD = simulate_pulse_train(TT,D.PULSE.ton,D.PULSE.toff,D.PULSE.conc);
%     OD = OD(end,:);
    plt.ax1 = plot(TT,OD,'-','LineWidth',plt.Lwd);
    ylabel({'Conc.','(uM)'})
    set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
        'YTick', [0 D.PULSE.conc(1)],...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off','ColorOrder',turbo(8))

    nexttile([plt.scale(2) 1])
    hold on
    for k = 1:size(D.PULSE.ton,1)               
        plt.axl(k) = plot(D.T,real(D.PRED.Im(:,k)),'-','LineWidth',plt.Lwd);
    end
    xlabel('Time (sec)')
    ylabel({'Cell Current','(pA)'})
    lgd = legend({num2str(D.PULSE.toff)},'NumColumns',2,...
        'Location','best');
    title(lgd,'Duration. (Sec)')    
    set(gca,'XLim',plt.X,'XColor','none','XTick',[],'XTickLabel',[],...
        'YLim',[plt.ytick(1) plt.ytick(end)],'YTick',plt.ytick,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box','off','ColorOrder',turbo(8))

    nexttile([plt.scale(3) 1])
    axis();
    xlabel('Time (sec)')
    set(gca,'XLim',plt.X,'XTick',plt.xtick,...
        'YColor','none','YTick', [], 'YTickLabel', [],...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')
    
    exportgraphics(gcf,plt.fname,'Resolution',300)
end