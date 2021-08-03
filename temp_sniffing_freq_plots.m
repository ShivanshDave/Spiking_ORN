%%
np = 10;
freq = 0.5;
SpikeEN = 1; plt.N = 3;
T=1/freq;
% Ton = T/3;
%%
PULSE.ton = T*(0:np-1).*ones(plt.N,1);
PULSE.toff = round(T/3,1) + T*(0:np-1).*ones(plt.N,1);
PULSE.conc = [5,20,100]'.*ones(1,np).*ones(plt.N,1);
PULSE.tspan = [-T 11*T];
%%

DATA = simulate_ORN(PULSE,SpikeEN);

%%
% ;
figure(1);
clf
t = tiledlayout(3,2,'TileSpacing','none','Padding','compact');

nexttile
 TT = linspace(DATA.T(1),DATA.T(end),100);
 OD = simulate_pulse_train(TT,PULSE.ton,PULSE.toff,PULSE.conc);
% OD = OD(end,:);
plot(TT,OD);
% legend(num2str(PULSE.conc))
xlabel('Time (sec)')
ylabel('Conc. (uM)')
axis(axis + [0 0 -1 1])
set(gca,'FontSize',14);

nexttile; 
plot(DATA.T,real(DATA.PRED.Im))
xlabel('Time (sec)')
ylabel('Receptor (pA)')
title('ORN Current')
set(gca,'FontSize',14);
hold on
spikes = script_spikes_ID(real(DATA.PRED.spkV),DATA.T,0);
dind=8;
set(gca,'ColorOrderIndex',1)
for i=1:plt.N
    scatter(DATA.T(spikes{i,1}-dind),...
        1.1*real(DATA.PRED.Im(spikes{i,1}-dind,i)), 25, '^','filled')
end


nexttile; 
plot(DATA.T,real(DATA.PRED.Ca))
xlabel('Time (sec)')
ylabel('Calcium')
title('Ca')
set(gca,'FontSize',14);

nexttile
plot(DATA.T,real(DATA.PRED.spkV))
xlabel('Time (sec)')
ylabel('Spike (mV)')
title('ML Spikes')
set(gca,'FontSize',14);

nexttile
plot(DATA.T,real(DATA.PRED.CaFR))
xlabel('Time (sec)')
title('CaFR')
set(gca,'FontSize',14);

nexttile
plot(DATA.T,real(DATA.PRED.nK))
xlabel('Time (sec)')
title('nK')
set(gca,'FontSize',14);


