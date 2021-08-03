%% X
SpikeEN = 1; plt.N = 5;
PULSE.ton = 0.000*ones(plt.N,1);
PULSE.toff = [10,0.5,0.25,0.1,0.05]';
PULSE.conc = 20*ones(plt.N,1);
PULSE.tspan = [-.01 12];
DATA = simulate_ORN(PULSE,SpikeEN);

%%
plt.Lwd = 1.2;
plt.FTsz = 16;
plt.Xoff = 0.05;
plt.FGpos = [10 10 500 500];
plt.scale = [3 4 1 1];
plt.ytick = [-55,0,20];
plt.xtick = 0:0.25:1;
plt.tspan = [0 .75];
plt.fname = '.\Report\figs\fig_spk_compare_dur.png';

plot_pulse_currents(plt,DATA)

% plot pulse and each current
function plot_pulse_currents(plt, D)

    figure('Renderer', 'painters', 'Position', plt.FGpos);
    nc=size(D.PULSE.ton,1);
    plt.t = tiledlayout(sum([plt.scale,(nc-1)*plt.scale(2)]),1,...
        'TileSpacing','none','Padding','compact');
    plt.X = [plt.tspan(1)-plt.Xoff, plt.tspan(2)];
    
    nexttile([plt.scale(1) 1])
    TT = linspace(D.T(1),D.T(end),100);
    OD = simulate_pulse_train(TT,D.PULSE.ton,D.PULSE.toff,D.PULSE.conc);
    OD = OD(3,:);
    plt.ax1 = plot(TT,OD,'k--','LineWidth',plt.Lwd);
    ylabel({'Stim'})
    set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
        'YTick', [0 20], 'YTickLabel', {'0','20'},'YTickLabelRotation',0,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')
    
    nexttile([1 1])
    axis off
%     set(gca,'XColor','none','XTick',[],'XTickLabel', [],...
%         'YColor','none','YTick', [], 'YTickLabel', [],...
%         'tickdir', 'out','FontSize',plt.FTsz,...
%         'color','none','box', 'off')
    
    clr = summer(nc);
    for k = 1:nc
        nexttile([plt.scale(2) 1])       
        plt.axl(k) = plot(D.T,real(D.PRED.Im(:,k)),'-',...
            'LineWidth',plt.Lwd,'Color',clr(2,:));
        text(plt.tspan(1),15,[num2str(D.PULSE.toff(k)),' s'],...
            'HorizontalAlignment','left','FontSize',plt.FTsz-2,'Color','black')
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
