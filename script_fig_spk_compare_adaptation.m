%% X
SpikeEN = 1; plt.N = 4;
PULSE.ton = [0 4].*ones(plt.N,1);
PULSE.toff = [4 5].*ones(plt.N,1);
PULSE.conc(:,1) = flip([10,5,1.75,0]');
PULSE.conc(:,2) = 20*ones(plt.N,1);
PULSE.tspan = [-0.25 6];
DATA = simulate_ORN(PULSE,SpikeEN);

%%
plt.Lwd = 1.2;
plt.FTsz = 18;
plt.Xoff = 0.1;
plt.FGpos = [10 10 1200 500];
plt.scale = [3 5 1];
plt.ytick = [-50,0,20];
plt.xtick = [0 2 4 5 6];
plt.tspan = PULSE.tspan;
plt.fname = '.\Report\figs\v1\fig_spk_compare_adaptation.png';

plot_pulse_currents(plt,DATA)

% plot pulse and each current
function plot_pulse_currents(plt, D)

    figure('Renderer', 'painters', 'Position', plt.FGpos);
    nc=size(D.PULSE.ton,1);
    plt.t = tiledlayout(sum([plt.scale,(nc-1)*plt.scale(2)]),1,...
        'TileSpacing','none','Padding','tight');
    plt.X = [plt.tspan(1)-plt.Xoff, plt.tspan(2)];
    
    nexttile([plt.scale(1) 1])
    TT = linspace(D.T(1),D.T(end),100);
    OD = simulate_pulse_train(TT,D.PULSE.ton,D.PULSE.toff,D.PULSE.conc);
    OD = OD(end,:);
    plt.ax1 = plot(TT,OD,'k-','LineWidth',plt.Lwd);
    ylabel({'Stim'})
    set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
        'YTick', [0 10 20], 'YTickLabel', {'0','Var','20'},'YTickLabelRotation',0,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')
    
    clr = flip(turbo(8).*0.75);
    for k = 1:nc
        nexttile([plt.scale(2) 1])       
        plt.axl(k) = plot(D.T,real(D.PRED.Im(:,k)),'-',...
            'LineWidth',plt.Lwd,'Color',clr(k,:));
        text(plt.tspan(1),-10,[num2str(D.PULSE.conc(k,1)),' \muM'],...
            'HorizontalAlignment','left','FontSize',plt.FTsz-1,'Color',clr(k,:))
        set(gca,'XLim',plt.X,'XColor','none','XTick',[],'XTickLabel',[],...
            'YLim',[plt.ytick(1) plt.ytick(end)],...
            'YColor','none','YTick',[],'YTickLabel',[],...
            'tickdir','out','FontSize',plt.FTsz,...
            'color','none','box','off')   
    end
    ylabel({'I_{ORN} (pA)'})
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
