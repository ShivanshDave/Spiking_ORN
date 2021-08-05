% script_fig_sniff_freq_FRss
%% X
b_frq = 20; % Breaths/min
conc = [10,50,300]';

np = 10; 
T=60/b_frq;
SpikeEN = 1; 
plt.N = length(conc); 
PULSE.ton = T*(0:np-1).*ones(plt.N,1);
PULSE.toff = round(T/3,1) + T*(0:np-1).*ones(plt.N,1);
PULSE.conc = conc.*ones(1,np).*ones(plt.N,1);
PULSE.tspan = [-1:0.1:np*T+1];


DATA = simulate_ORN(PULSE);

%%
plt.Lwd = 1.2;
plt.FTsz = 16;
plt.Xoff = 0.1;
plt.FGpos = [10 10 900 600];
plt.scale = [3 4 1];
plt.ytick = [-50,0,20];
plt.xtick = [0 2 4 6];
plt.tspan = PULSE.tspan;
plt.fname = '.\Report\figs\fig_spk_compare_adaptation.png';

plot_pulse_currents(plt,DATA)

% plot pulse and each current
function plot_pulse_currents(plt, D)


end