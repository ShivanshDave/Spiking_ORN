% script_fig_sniff_freq_FRss

%% SETUP
D.test_frq = [5,10,15,20,25,30,35,40,45,50];
D.conc = [10,20,50,300]';
D.np = 5; 
plt.raw_plots = 1;
plt.FR_plots = 1;
plt.pulseOnly = 1; % 1 for Pulse-Figure (methods section)
plt.plotsOnly = 1; % 1 for re-generate figures

if plt.pulseOnly; D.test_frq = 30; end 

%% COLLECT DATA
for i=1:length(D.test_frq)
    b_frq = D.test_frq(i); % Breaths/min    
    T=60/b_frq;
    SpikeEN = 1; 
    plt.N = length(D.conc); 
    PULSE.ton = T*(0:D.np-1).*ones(plt.N,1);
    PULSE.toff = round(T/3,1) + T*(0:D.np-1).*ones(plt.N,1);
    PULSE.conc = D.conc.*ones(1,D.np).*ones(plt.N,1);
    PULSE.tspan = [-1 D.np*T+1];
    if plt.pulseOnly; continue; end

    % %
    disp(['F=',num2str(b_frq)]);
    if ~plt.plotsOnly
        DATA = simulate_ORN(PULSE);
        DATA.spk = script_spikes_ID(real(DATA.PRED.spkV),DATA.T,0);
        D.(['f',num2str(b_frq)]) = DATA;
    else
        DATA = D.(['f',num2str(b_frq)]);
    end
    % %
    if plt.raw_plots
        % % Plot
        plt.Lwd = 1;
        plt.FTsz = 14;
        plt.FGpos = [10 10 1000 600];
        plt.Xoff = 0.2;
        plt.xtick = -1:4;
        plt.scale = [1,2,4];
        plt.fname = ['.\Report\figs\sniff\fig_spk_sniffing_' num2str(b_frq) 'bpm'];
        plot_one_bfreq(plt,DATA);
    end
end

%% FIND FR
D.FRss = nan(length(D.conc),length(D.test_frq));
D.FRmx = D.FRss;
D.num_spk_ss = D.FRss;
D.num_spk_mx = D.FRss;

for f=1:length(D.test_frq)
    Dspk = D.(['f',num2str(D.test_frq(f))]);    
    % LAST
    dur_stim = [Dspk.PULSE.ton(end),Dspk.PULSE.toff(end)];
    D.num_spk_ss(:,f) = cellfun(@(v) sum(v>=dur_stim(1) & v<=dur_stim(2)), Dspk.spk(:,3)); 
    D.FRss(:,f) = D.num_spk_ss(:,f)/diff(dur_stim);
    %FIRST
    dur_stim = [Dspk.PULSE.ton(1),Dspk.PULSE.toff(1)];
    D.num_spk_mx(:,f) = cellfun(@(v) sum(v>=dur_stim(1) & v<=dur_stim(2)), Dspk.spk(:,3)); 
    D.FRmx(:,f) = D.num_spk_mx(:,f)/diff(dur_stim);
end

%%
if plt.FR_plots
    % % Plot
    plt.Lwd = 2;
    plt.FTsz = 24;
    plt.FGpos = [10 10 900 700];
    plt.Xoff = 1;
    plt.xtick = -1:4;
    plt.scale = [7,1];
    plt.fname = ['.\Report\figs\v1\fig_sniff_freq_FR_tuning'];

    plot_freq_FRss(plt,D);
end
%%
if plt.pulseOnly
    % % Plot
    plt.Lwd = 1;
    plt.FTsz = 18;
    plt.FGpos = [10 10 1200 300];
    plt.Xoff = 0.25;
    plt.xtick = -1:12;
    plt.scale = [1,1];
    plt.fname = ['.\Report\figs\v1\fig_sniff_pulse'];

    plot_freq_pulse(plt,PULSE);
end

disp('---DONE---')
%%
function plot_freq_pulse(plt,PULSE)
    figure('Renderer', 'painters', 'Position', plt.FGpos);
    plt.t = tiledlayout(plt.scale(1),plt.scale(2),'TileSpacing','tight','Padding','tight');
    plt.X = [PULSE.tspan(1)-plt.Xoff, PULSE.tspan(2)];
    plt.ax = [];

    % PULSE
    plt.ax(end+1) = nexttile([plt.scale(1) plt.scale(2)]);
    TT = linspace(PULSE.tspan(1),PULSE.tspan(2),100);
    OD = simulate_pulse_train(TT,PULSE.ton,PULSE.toff,PULSE.conc);
    % OD = OD(end,:);
    plot(TT,OD,'LineWidth',plt.Lwd);
    lgd=legend(num2str(PULSE.conc(:,1)),'Box','off');
    title(lgd,'Conc. (uM)') 
    xlabel('Time (sec)')
    ylabel({'O_{stim} (uM)'})
    set(gca,'XLim',plt.X,'XTick', plt.xtick,...
            'YLim',[-30 300],'YTickLabelRotation',0,...
            'tickdir', 'out','FontSize',plt.FTsz,...
            'color','none','box', 'off')
        
    exportgraphics(gcf,[plt.fname '.png'],'Resolution',300)
end

function plot_freq_FRss(plt,D)
    %% INFORMATION GAIN
    figure('Renderer', 'painters', 'Position', plt.FGpos);
    plt.t = tiledlayout(2*plt.scale(1),1,'TileSpacing','tight','Padding','loose');
    plt.X = [D.test_frq(1)-plt.Xoff, D.test_frq(end)];
    plt.ax = [];

    % FR Max - First stim
    plt.ax(end+1) = nexttile([plt.scale(1) plt.scale(2)]);
    plot(D.test_frq, D.FRmx,'LineWidth',plt.Lwd)
    hold on
    set(gca,'ColorOrderIndex',1)
    scatter(D.test_frq, D.FRmx,75,'^')
    % xlabel('Time (sec)')
    ylabel({'(peak)'})
    title('First stimulus-pulse')
    set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
            'YTickLabelRotation',0,...
            'tickdir', 'out','FontSize',plt.FTsz,...
            'color','none','box', 'off')

    % FR ss - Last stim
    plt.ax(end+1) = nexttile([plt.scale(1) plt.scale(2)]);
    plot(D.test_frq, D.FRss,'LineWidth',plt.Lwd)
    hold on
    set(gca,'ColorOrderIndex',1)
    scatter(D.test_frq, D.FRss,75,'^')
    xlabel('Breathing frequency (bpm)')
    ylabel({'(steady-state)'})
    title('Last stimulus-pulse')
    lgd = legend({num2str(D.conc)},'Location','best','Box','off');
    title(lgd,'Conc. (uM)') 
    set(gca,'XLim',plt.X,'YLim',[-1 10],...
            'YTickLabelRotation',0,...
            'tickdir', 'out','FontSize',plt.FTsz,...
            'color','none','box', 'off')
        
    ylabel(plt.t, 'Sniffing sensitivity','FontSize',plt.FTsz+2)
        

    linkaxes(plt.ax(:),'x')    
    exportgraphics(gcf,[plt.fname '_gain.png'],'Resolution',300)
    
    
    %% COUNT
    figure('Renderer', 'painters', 'Position', plt.FGpos);
    plt.t = tiledlayout(2*plt.scale(1),1,'TileSpacing','tight','Padding','loose');
    plt.X = [D.test_frq(1)-plt.Xoff, D.test_frq(end)];
    plt.ax = [];

    % FR Max - First stim
    plt.ax(end+1) = nexttile([plt.scale(1) plt.scale(2)]);
    plot(D.test_frq, D.num_spk_mx,'LineWidth',plt.Lwd)
    hold on
    set(gca,'ColorOrderIndex',1)
    scatter(D.test_frq, D.num_spk_mx,75,'^')
    % xlabel('Time (sec)')
    ylabel({'(peak)'})
    title('First stimulus-pulse')
    set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
            'YTickLabelRotation',0,...
            'tickdir', 'out','FontSize',plt.FTsz,...
            'color','none','box', 'off')

    % FR ss - Last stim
    plt.ax(end+1) = nexttile([plt.scale(1) plt.scale(2)]);
    plot(D.test_frq, D.num_spk_ss,'LineWidth',plt.Lwd)
    hold on
    set(gca,'ColorOrderIndex',1)
    scatter(D.test_frq, D.num_spk_ss,75,'^')
    xlabel('Breathing frequency (bpm)')
    ylabel({'(steady-state)'})
    title('Last stimulus-pulse')
%     lgd = legend({num2str(D.conc)},'Location','east','Box','off');
%     title(lgd,'Conc. (uM)') 
    set(gca,'XLim',plt.X,'YLim',[-1 10],...
            'YTickLabelRotation',0,...
            'tickdir', 'out','FontSize',plt.FTsz,...
            'color','none','box', 'off')
        
    ylabel(plt.t, 'ORN Spike count','FontSize',plt.FTsz+2)
       
    linkaxes(plt.ax(:),'x')    
    exportgraphics(gcf,[plt.fname '_raw.png'],'Resolution',300)
    

end

function plot_one_bfreq(plt,DATA)
    figure('Renderer', 'painters', 'Position', plt.FGpos);
    plt.t = tiledlayout(3*plt.scale(1),sum(plt.scale(2:3)),'TileSpacing','tight','Padding','tight');
    plt.X = [DATA.PULSE.tspan(1)-plt.Xoff, DATA.PULSE.tspan(2)];
    plt.ax = [];

    % PULSE
    plt.ax(end+1) = nexttile([plt.scale(1) plt.scale(2)]);
    TT = linspace(DATA.T(1),DATA.T(end),100);
    OD = simulate_pulse_train(TT,DATA.PULSE.ton,DATA.PULSE.toff,DATA.PULSE.conc);
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
    spikes = DATA.spk; dind=8;
    set(gca,'ColorOrderIndex',1)
    for i=1:plt.N
        scatter(DATA.T(spikes{i,1}-dind),...
            1.1*real(DATA.PRED.Im(spikes{i,1}-dind,i)), 25, '^','filled')
    end
    % xlabel('Time (sec)')
    ylabel('I_{ORN} (pA)')
    set(gca,'XLim',plt.X,'XColor','none','XTick', [], 'XTickLabel', [],...
            'YTickLabelRotation',0,...
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
    ylabel('V_{ML} (mV)')
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
    exportgraphics(gcf,[plt.fname '.png'],'Resolution',300)
    
    xlim([DATA.PULSE.ton(end),DATA.PULSE.toff(end)])
    exportgraphics(gcf,[plt.fname '_last.png'],'Resolution',300)

end