%% R99.F2A
SpikeEN = 0;
PULSE.ton = 0.000*ones(8,1);
PULSE.toff = 1.00*ones(8,1);
PULSE.conc = [300,100,50,20,10,5,2,1]';
PULSE.tspan = [-1 4];
DATA = simulate_ORN(PULSE,SpikeEN);

%%
plt.Lwd = 1.2;
plt.FTsz = 14;
plt.Xoff = 0.1;
plt.FGpos = [10 10 900 500];
plt.scale = [1 4];
plt.ytick = [-50,-25,0];
plt.xtick = [0:4];
plt.fname = '.\Report\figs\fig_txn_compare.png';

plot_pulse_currents_overlap(plt,DATA)
%%
function plot_pulse_currents_overlap(plt,D)

    figure('Renderer', 'painters', 'Position', plt.FGpos);
    t = tiledlayout(sum(plt.scale),1,'TileSpacing','none','Padding','compact');
    
    nexttile([plt.scale(1) 1])
    TT = linspace(D.T(1),D.T(end),100);
    OD = simulate_pulse_train(TT,D.PULSE.ton,D.PULSE.toff,D.PULSE.conc);
    OD = OD(end,:);
    ax1 = plot(TT,OD,'k-','LineWidth',plt.Lwd);
    set(gca,'XColor','none','YTick', [0 1], 'YTickLabel', [],...
        'FontSize',plt.FTsz,...
        'color','none','box', 'off', 'tickdir', 'out')
    xlim([D.PULSE.tspan(1)-plt.Xoff, D.PULSE.tspan(2)])
    ylabel({'Conc.','(uM)'})

    nexttile([plt.scale(2) 1])
    hold on
    for k = 1:size(D.PULSE.ton,1)               
        axl(k) = plot(D.T,real(D.PRED.Im(:,k)),'-','LineWidth',plt.Lwd);
%         set(gca,'ColorOrderIndex',k)
    end
%     axis on
    ylim([-50 0])
    xlim([D.PULSE.tspan(1)-plt.Xoff, D.PULSE.tspan(2)])
    set(gca,'YTick', plt.ytick,'XTick',plt.xtick,'FontSize',plt.FTsz,...
        'color','none','box', 'off', 'tickdir', 'out')
    xlabel('Time (sec)')
    ylabel('Cell Current (pA)')
    legend(num2str(D.PULSE.conc))   
    
    exportgraphics(gcf,plt.fname,'Resolution',300)
end