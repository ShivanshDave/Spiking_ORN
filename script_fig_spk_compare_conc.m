%% Figure-1 (Reisert 1999) Pulse-1, Conc-8
SpikeEN = 1; plt.N = 8;
PULSE.ton = 0.000*ones(plt.N,1);
PULSE.toff = 1.00*ones(plt.N,1);
PULSE.conc = [300,100,50,20,10,5,2,1]';
PULSE.tspan = [-0.5 3];
DATA = simulate_ORN(PULSE,SpikeEN);

%%
plt.Lwd = 1.2;
plt.FTsz = 16;
plt.Xoff = 0.1;
plt.FGpos = [10 10 700 900];
plt.scale = [2 4 1];
plt.ytick = [-65,0,20];
plt.xtick = 0:3;
plt.fname = '.\Report\figs\fig_spk_compare_conc.png';

plot_pulse_currents(plt,DATA)

% plot pulse and each current
function plot_pulse_currents(plt, D)

    figure('Renderer', 'painters', 'Position', plt.FGpos);
    nc=size(D.PULSE.ton,1);
    plt.t = tiledlayout(sum([plt.scale,(nc-1)*plt.scale(2)]),1,...
        'TileSpacing','none','Padding','compact');
    plt.X = [D.PULSE.tspan(1)-plt.Xoff, D.PULSE.tspan(2)];
    
    nexttile([plt.scale(1) 1])
    TT = linspace(D.T(1),D.T(end),100);
    OD = simulate_pulse_train(TT,D.PULSE.ton,D.PULSE.toff,D.PULSE.conc);
    OD = OD(end,:);
    plt.ax1 = plot(TT,OD,'k-','LineWidth',plt.Lwd);
    ylabel({'Stim'})
    set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
        'YTick', [0 1], 'YTickLabel', {'0','Var'},'YTickLabelRotation',0,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')
    
    clr = turbo(nc).*0.75;
    for k = 1:nc
        nexttile([plt.scale(2) 1])       
        plt.axl(k) = plot(D.T,real(D.PRED.Im(:,k)),'-',...
            'LineWidth',plt.Lwd,'Color',clr(k,:));
        text(D.PULSE.tspan(2),20,[num2str(D.PULSE.conc(k)),' \muM'],...
            'HorizontalAlignment','right','FontSize',plt.FTsz,'Color',clr(k,:))
        set(gca,'XLim',plt.X,'XColor','none','XTick',[],'XTickLabel',[],...
            'YLim',[plt.ytick(1) plt.ytick(end)],...
            'YColor','none','YTick',[],'YTickLabel',[],...
            'tickdir','out','FontSize',plt.FTsz,...
            'color','none','box','off')
        
    end
    ylabel({'I_M (pA)'})
    set(gca,'XLim',plt.X,'XColor','none','XTick',[],'XTickLabel',[],...
        'YLim',[plt.ytick(1) plt.ytick(end)],...
        'YColor',[0.1500 0.1500 0.1500],'YTick',plt.ytick,...
        'YTickLabel',{plt.ytick(1),'',plt.ytick(end)},...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box','off')

    nexttile([plt.scale(3) 1])
    axis();
    xlabel('Time (sec)')
    set(gca,'XLim',plt.X,'XTick',plt.xtick,...
        'YColor','none','YTick', [], 'YTickLabel', [],...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')
    
    exportgraphics(gcf,plt.fname,'Resolution',300)
end
