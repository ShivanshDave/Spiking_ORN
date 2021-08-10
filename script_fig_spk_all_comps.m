% Odor-pulse setup
PULSE.ton = [ 0.2000  2 
              0.2000 2
              0.2000 2];
PULSE.toff = [1.2000  3 
              1.2000  3
              1.2000  3];
PULSE.conc = [30 5
              5  30
              0  0];
PULSE.tspan = [0 4];


%% RUN
SpikeEN=[1;1;1];
DATA = simulate_ORN(PULSE,SpikeEN);

%% Plot
plt.Lwd = 1;
plt.FTsz = 14;
plt.FGpos = [10 10 1000 600];
plt.Xoff = 0.2;
plt.xtick = -1:4;
plt.scale = [1,2,4];
plt.fname = '.\Report\figs\fig_spk_all_components.png';

figure('Renderer', 'painters', 'Position', plt.FGpos);
plt.t = tiledlayout(3*plt.scale(1),sum(plt.scale(2:3)),'TileSpacing','tight','Padding','tight');
plt.X = [DATA.PULSE.tspan(1)-plt.Xoff, DATA.PULSE.tspan(2)];

% PULSE
nexttile([plt.scale(1) plt.scale(2)])
TT = linspace(DATA.T(1),DATA.T(end),100);
OD = simulate_pulse_train(TT,PULSE.ton,PULSE.toff,PULSE.conc);
% OD = OD(end,:);
plot(TT,OD,'LineWidth',plt.Lwd);
% legend(num2str(PULSE.conc))
% xlabel('Time (sec)')
ylabel('Conc. (uM)')
set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
        'YTick', [0 15 30], 'YTickLabelRotation',0,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')

% I ORN
nexttile([plt.scale(1) plt.scale(3)])
plot(DATA.T,real(DATA.PRED.Im),'LineWidth',plt.Lwd)
% xlabel('Time (sec)')
ylabel('I_{ORN} (pA)')
set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
        'YLim',[-55 20],'YTick', [-55 0 20], 'YTickLabelRotation',0,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')

% Ca
nexttile([plt.scale(1) plt.scale(2)])
plot(DATA.T,real(DATA.PRED.Ca),'LineWidth',plt.Lwd)
% xlabel('Time (sec)')
ylabel('Calcium')
% title('Ca')
set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
        'YLim',[0 10],'YTick', [0 5 10], 'YTickLabelRotation',0,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')

% V SPK
nexttile([plt.scale(1) plt.scale(3)])
plot(DATA.T,real(DATA.PRED.spkV),'LineWidth',plt.Lwd)
% xlabel('Time (sec)')
ylabel('V_{SPK} (mV)')
set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
        'YLim',[-60 40],'YTick', [-60 0 40], 'YTickLabelRotation',0,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')

% CaFR    
nexttile([plt.scale(1) plt.scale(2)])
plot(DATA.T,real(DATA.PRED.CaFR),'LineWidth',plt.Lwd)
xlabel('Time (sec)')
ylabel('CaFR')
set(gca,'XLim',plt.X,'XTick', plt.xtick,...
        'YLim',[0,10],'YTick', [1 5 10], 'YTickLabelRotation',0,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')

% nK
nexttile([plt.scale(1) plt.scale(3)])
plot(DATA.T,real(DATA.PRED.nK),'LineWidth',plt.Lwd)
xlabel('Time (sec)')
ylabel('nK')
set(gca,'XLim',plt.X,'XTick', plt.xtick,...
        'YLim',[-.05,.6],'YTick', [0 0.3 0.6], 'YTickLabelRotation',0,...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')
    
exportgraphics(gcf,plt.fname,'Resolution',300)
    
% figure;
% nexttile; plot(DATA.T, DATA.IPREDn.IL); title('IL')
% nexttile; plot(DATA.T, DATA.IPREDn.ICACL); title('Icacl')
% nexttile; plot(DATA.T, DATA.IPREDn.ICNG); title('Icng')
% nexttile; plot(DATA.T, DATA.IPREDn.PRED_CURRENT); title('PC')

disp('---Done---')