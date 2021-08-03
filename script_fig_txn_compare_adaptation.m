%% R99.F4
SpikeEN = 0; plt.N = 8;
PULSE.ton = [-4 0].*ones(plt.N,1);
PULSE.toff = [0 1].*ones(plt.N,1);
PULSE.conc(:,1) = 5*ones(plt.N,1);
PULSE.conc(:,2) = [300,100,50,20,10,5,2,1]';
PULSE.tspan = [-5 4];
DATA = simulate_ORN(PULSE,SpikeEN);

%%
plt.Lwd = 1.2;
plt.FTsz = 16;
plt.Xoff = 0.1;
plt.FGpos = [10 10 900 500];
plt.scale = [3 12 1];
plt.ytick = [-50,-25,0];
plt.xtick = [-4,0,1,4];
plt.fname = '.\Report\figs\fig_txn_compare_adaptation.png';

plot_pulse_currents_overlap(plt,DATA)
%%
function plot_pulse_currents_overlap(plt,D)

    figure('Renderer', 'painters', 'Position', plt.FGpos);
    plt.t = tiledlayout(sum(plt.scale),1,'TileSpacing','tight','Padding','compact');
    plt.X = [D.PULSE.tspan(1)-plt.Xoff, D.PULSE.tspan(2)];
    
    nexttile([plt.scale(1) 1])
    TT = linspace(D.T(1),D.T(end),100);
    OD = simulate_pulse_train(TT,D.PULSE.ton,D.PULSE.toff,D.PULSE.conc);
    OD = OD(5,:);
    plt.ax1 = plot(TT,OD,'k-','LineWidth',plt.Lwd);
    ylabel({'Conc.','(uM)'})
    set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
        'YTick', [0 5 10], 'YTickLabel', {'0','5','Var'},...
        'YTickLabelRotation',0,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')

    nexttile([plt.scale(2) 1])
    hold on
    for k = 1:size(D.PULSE.ton,1)               
        plt.axl(k) = plot(D.T,real(D.PRED.Im(:,k)),'-','LineWidth',plt.Lwd);
    end
    xlabel('Time (sec)')
    ylabel({'Cell Current','(pA)'})
    lgd = legend({num2str(D.PULSE.conc(:,2))},'Location','best');
    title(lgd,'Conc. (uM)')    
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