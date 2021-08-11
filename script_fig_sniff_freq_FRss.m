% script_fig_sniff_freq_FRss
%% X
b_frq = 30; % Breaths/min
conc = [10,50,300]';
np = 5; 

T=60/b_frq;
SpikeEN = 1; 
plt.N = length(conc); 
PULSE.ton = T*(0:np-1).*ones(plt.N,1);
PULSE.toff = round(T/3,1) + T*(0:np-1).*ones(plt.N,1);
PULSE.conc = conc.*ones(1,np).*ones(plt.N,1);
PULSE.tspan = [-1 np*T+1];

%%
DATA = simulate_ORN(PULSE);

%% Plot
plt.Lwd = 1;
plt.FTsz = 14;
plt.FGpos = [10 10 1000 600];
plt.Xoff = 0.2;
plt.xtick = -1:4;
plt.scale = [1,2,4];
plt.fname = ['.\Report\figs\sniff\fig_spk_sniffing_' num2str(b_frq) 'bpm.png'];

plot_one_bfreq(plt,DATA);


%%


function plot_one_bfreq(plt,DATA)
    figure('Renderer', 'painters', 'Position', plt.FGpos);
    plt.t = tiledlayout(3*plt.scale(1),sum(plt.scale(2:3)),'TileSpacing','tight','Padding','tight');
    plt.X = [DATA.PULSE.tspan(1)-plt.Xoff, DATA.PULSE.tspan(2)];
    plt.ax = [];

    % PULSE
    plt.ax(end+1) = nexttile([plt.scale(1) plt.scale(2)]);
    TT = linspace(DATA.T(1),DATA.T(end),100);
    OD = simulate_pulse_train(TT,PULSE.ton,PULSE.toff,PULSE.conc);
    % OD = OD(end,:);
    plot(TT,OD,'LineWidth',plt.Lwd);
    % legend(num2str(PULSE.conc))
    % xlabel('Time (sec)')
    ylabel('Conc. (uM)')
    set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
            'YTickLabelRotation',0,...
            'tickdir', 'out','FontSize',plt.FTsz,...
            'color','none','box', 'off')

    % I ORN
    plt.ax(end+1) = nexttile([plt.scale(1) plt.scale(3)]);
    plot(DATA.T,real(DATA.PRED.Im),'LineWidth',plt.Lwd)
    hold on
    spikes = script_spikes_ID(real(DATA.PRED.spkV),DATA.T,0);
    dind=8;
    set(gca,'ColorOrderIndex',1)
    for i=1:plt.N
        scatter(DATA.T(spikes{i,1}-dind),...
            1.1*real(DATA.PRED.Im(spikes{i,1}-dind,i)), 25, '^','filled')
    end
    % xlabel('Time (sec)')
    ylabel('I_{ORN} (pA)')
    set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
            'YLim',[-55 20],'YTick', [-55 0 20], 'YTickLabelRotation',0,...
            'tickdir', 'out','FontSize',plt.FTsz,...
            'color','none','box', 'off')

    % Ca
    plt.ax(end+1) = nexttile([plt.scale(1) plt.scale(2)]);
    plot(DATA.T,real(DATA.PRED.Ca),'LineWidth',plt.Lwd)
    % xlabel('Time (sec)')
    ylabel('Calcium')
    % title('Ca')
    set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
            'YLim',[0 10],'YTick', [0 5 10], 'YTickLabelRotation',0,...
            'tickdir', 'out','FontSize',plt.FTsz,...
            'color','none','box', 'off')

    % V SPK
    plt.ax(end+1) = nexttile([plt.scale(1) plt.scale(3)]);
    plot(DATA.T,real(DATA.PRED.spkV),'LineWidth',plt.Lwd)
    % xlabel('Time (sec)')
    ylabel('V_{SPK} (mV)')
    set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
            'YLim',[-60 40],'YTick', [-60 0 40], 'YTickLabelRotation',0,...
            'tickdir', 'out','FontSize',plt.FTsz,...
            'color','none','box', 'off')

    % CaFR    
    plt.ax(end+1) = nexttile([plt.scale(1) plt.scale(2)]);
    plot(DATA.T,real(DATA.PRED.CaFR),'LineWidth',plt.Lwd)
    xlabel('Time (sec)')
    ylabel('CaFR')
    set(gca,'XLim',plt.X,...
            'YLim',[0,10],'YTick', [1 5 10], 'YTickLabelRotation',0,...
            'tickdir', 'out','FontSize',plt.FTsz,...
            'color','none','box', 'off')

    % nK
    plt.ax(end+1) = nexttile([plt.scale(1) plt.scale(3)]);
    plot(DATA.T,real(DATA.PRED.nK),'LineWidth',plt.Lwd)
    xlabel('Time (sec)')
    ylabel('nK')
    set(gca,'XLim',plt.X,...
            'YLim',[-.05,.6],'YTick', [0 0.3 0.6], 'YTickLabelRotation',0,...
            'tickdir', 'out','FontSize',plt.FTsz,...
            'color','none','box', 'off')

    linkaxes(plt.ax(:),'x')    
    exportgraphics(gcf,plt.fname,'Resolution',300)

end