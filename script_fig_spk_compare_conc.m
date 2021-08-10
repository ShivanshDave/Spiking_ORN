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

%% << F-3 >>
plt.Lwd = 1.2;
plt.FTsz = 14;
plt.Xoff = 0.1;
plt.FGpos = [10 10 600 900];
plt.scale = [0 5 0];
plt.ytick = [-50,-25,0];
plt.xtick = -1:4;
plt.fname = '.\Report\figs\fig_spk_compare_conc_quant.png';

plot_spk_quantify(plt,DATA)
%%
function plot_spk_quantify(plt,DATA)

    figure('Renderer', 'painters', 'Position', plt.FGpos);
    plt.t = tiledlayout(3*sum(plt.scale),1,'TileSpacing','tight','Padding','compact');
    
    
    %%
    spikes = script_spikes_ID(real(DATA.PRED.spkV),DATA.T,0);
    %%
    spk_t = spikes(:,3);
    idx = ~cellfun('isempty',spk_t);
    spk_freq = zeros(size(spk_t));
    spk_count = zeros(size(spk_t));
    spk_latency = nan(size(spk_t));
    
    spk_latency(idx) = cellfun(@(v)v(1),spk_t(idx));
    spk_count(idx) = cellfun(@(v) length(v),spk_t(idx));
    spk_freq(idx) = cellfun(@(v) v(end)-v(1),spk_t(idx));
    spk_freq(idx) = spk_count(idx)./spk_freq(idx);
    
    
    
    sz = 100;
    %% A - Freq
    nexttile([plt.scale(2) 1])
    hold on
%     plot(DATA.PULSE.conc, real(max(DATA.PRED.Im)),'k--')
        yyaxis right
    plot(DATA.PULSE.conc, real(min(DATA.PRED.Im)),'--')
    scatter(DATA.PULSE.conc, real(min(DATA.PRED.Im)),sz,0.75*turbo(8),'filled','s')
    ylabel({'Cell Current (pA)'})
    set(gca,'XColor','none','XTick', [], 'XTickLabel', [],...
        'xscale','log','tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box','off','YDir','reverse')
    
        yyaxis left
    plot(DATA.PULSE.conc, spk_freq,'-')
    scatter(DATA.PULSE.conc, spk_freq,sz,0.75*turbo(8),'s','LineWidth',1.5)
    ylabel({'Frequency (Hz)'})
%     xlabel('Concentration (uM)')
    
    set(gca,'XColor','none','XTick', [], 'XTickLabel', [],...
        'xscale','log','tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box','off')
    
    %% B - Spk count
    nexttile([plt.scale(2) 1])
    hold on
    plot(DATA.PULSE.conc, spk_count,'-','Color',[0 0.4470 0.7410])
    scatter(DATA.PULSE.conc, spk_count, sz,0.75*turbo(8),'o','LineWidth',1.5)
    ylabel({'Number of Spikes'})
%     xlabel('Concentration (uM)')
    
    set(gca,'XColor','none','XTick', [], 'XTickLabel', [],...
        'YColor',[0 0.4470 0.7410],'xscale','log','tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box','off')
    
    %% C - Spk latency
    nexttile([plt.scale(2) 1])
    hold on
        yyaxis right
    [~,id] = min(DATA.PRED.Im);
    plot(DATA.PULSE.conc, DATA.T(id),'--')
    scatter(DATA.PULSE.conc, DATA.T(id),sz,0.75*turbo(8),'filled','^')
    ylabel({'Time to paek (sec)'})
    ylim([0 1.2])
        yyaxis left
    plot(DATA.PULSE.conc, spk_latency,'k-')
    scatter(DATA.PULSE.conc, spk_latency, sz,0.75*turbo(8),'s','LineWidth',1.5)
    ylabel({'Latency (sec)'})
    ylim([0 1.2])
        
    
    xlabel('Concentration (uM)')
    set(gca,'xscale','log','tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box','off')

    %%
    exportgraphics(gcf,plt.fname,'Resolution',300)
end

% plot pulse and each current
function plot_pulse_currents(plt, D)

    figure('Renderer', 'painters', 'Position', plt.FGpos);
    nc=size(D.PULSE.ton,1);
    plt.t = tiledlayout(sum([plt.scale,(nc-1)*plt.scale(2)]),1,...
        'TileSpacing','none','Padding','compact');
    plt.X = [D.PULSE.tspan(1)-plt.Xoff, D.PULSE.tspan(2)];
    
    % STIM
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
    
    % SPK
    clr = turbo(nc).*0.75;
    for k = 1:nc
        nexttile([plt.scale(2) 1])       
        plt.axl(k) = plot(D.T,real(D.PRED.Im(:,k)),'-',...
            'LineWidth',plt.Lwd,'Color',clr(k,:));
        text(D.PULSE.tspan(2),-15,[num2str(D.PULSE.conc(k)),' \muM'],...
            'HorizontalAlignment','right','FontSize',plt.FTsz,'Color',clr(k,:))
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
    
    % X ax
    nexttile([plt.scale(3) 1])
    axis();
    xlabel('Time (sec)')
    set(gca,'XLim',plt.X,'XTick',plt.xtick,...
        'YColor','none','YTick', [], 'YTickLabel', [],...
        'tickdir', 'out','FontSize',plt.FTsz,...
        'color','none','box', 'off')
    
    exportgraphics(gcf,plt.fname,'Resolution',300)
end
