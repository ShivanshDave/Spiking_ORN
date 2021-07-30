%% R99.X
plt.Lwd = 1.2;
plt.FTsz = 16;
plt.Xoff = 0.02;
plt.FGpos = [10 10 900 500];
plt.scale = [8 8 1];
plt.ytick = [-80,-40,0,40];
plt.xtick = [0:0.1:0.7];
plt.rot = 0;
plt.fname = '.\Report\figs\fig_ML_spikes_ID.png';

plt.f = figure('Renderer', 'painters', 'Position', plt.FGpos);
plt.t = tiledlayout(sum(plt.scale),1,'TileSpacing','compact','Padding','compact');

%% ----- ONR and ORN+ML --------
SpikeEN = flip([0;1]); 
PULSE.ton = [0;0];
PULSE.toff = [1;1];
PULSE.conc = [50;50];
PULSE.tspan = [-0.001 0.7];
DATA = simulate_ORN(PULSE,SpikeEN);
plt.X = [DATA.PULSE.tspan(1)-plt.Xoff, DATA.PULSE.tspan(2)];
%%
plt.ind = 1;
% plt.lgd = {'Stim:5uM ORN+Spikes','Stim:5uM ORN'};
plt.lgd = {'ORN+Spikes','ORN only','Spikes-1','Spikes-2'};
plt.ltitle = 'C = 50uM';
plot_pulse_spikes(plt,DATA)

%% ----- ORN+ML : C1 and C2 --------
SpikeEN = flip([1;1]); 
PULSE.ton = [0;0];
PULSE.toff = [1;1];
PULSE.conc = [300;20];
PULSE.tspan = [-0.001 0.7];
DATA = simulate_ORN(PULSE,SpikeEN);
plt.X = [DATA.PULSE.tspan(1)-plt.Xoff, DATA.PULSE.tspan(2)];
%%
plt.ind = 2;
plt.lgd = {'C=300uM','C=20uM','Spikes-1','Spikes-2'};
plt.ltitle = 'ORN+Spikes';
plot_pulse_spikes(plt,DATA)
%%

function plot_pulse_spikes(plt,D)

    % Spikes
    nexttile([plt.scale(plt.ind) 1])
    hold on
    for k = 1:size(D.PULSE.ton,1)               
        plt.axl(k) = plot(D.T,real(D.PRED.Im(:,k)),'-','LineWidth',plt.Lwd);
    end
    spikes = script_spikes_ID(real(D.PRED.spkV),D.T,0);
    dind=8; set(gca,'ColorOrderIndex',1);
    i=1; scatter(D.T(spikes{i,1}-dind), 1.1*real(D.PRED.Im(spikes{i,1}-dind,i)), 25, '^','filled')
    i=2; scatter(D.T(spikes{i,1}-dind), 1.1*real(D.PRED.Im(spikes{i,1}-dind,i)), 25, '^','filled')
    lgd = legend(plt.lgd,'Location','northeast','NumColumns',2);
    title(lgd,plt.ltitle)
    xlabel('Time (sec)')
    ylabel('I_{M} (pA)')

    set(gca,'XLim',plt.X,'XColor','none','XTick',[],'XTickLabel',[],...
        'YLim',[plt.ytick(1) plt.ytick(end)],'YTick',plt.ytick,...
        'tickdir', 'out','FontSize',plt.FTsz,'YTickLabelRotation',plt.rot,...
        'color','none','box','off')

    if plt.ind == 1; return; end
    
    % TIME axis
    nexttile([plt.scale(end) 1])
    axis();
    xlabel('Time (sec)')
    set(gca,'XLim',plt.X,'XTick',plt.xtick,...
        'YColor','none','YTick', [], 'YTickLabel', [],...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')
    
    exportgraphics(plt.f,plt.fname,'Resolution',300)
end