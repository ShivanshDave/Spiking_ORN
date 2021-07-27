% Odor-pulse setup
PULSE.ton = [ 0.2000  2 
              0.2000 2
              0.2000 2];
PULSE.toff = [1.2000  3 
              1.2000  3
              1.2000  3];
PULSE.conc = [20 5
              5  20
              0  0];
PULSE.tspan = [0 4];


%% RUN
DATA = simulate_ORN(PULSE);

%% Plot
figure(1);
clf
t = tiledlayout(3,2,'TileSpacing','none','Padding','compact');

nexttile
OD = simulate_pulse_train(TT,PULSE.ton,PULSE.toff,PULSE.conc);
% OD = OD(end,:);
plot(TT,OD);
legend(num2str(PULSE.conc))
xlabel('Time (sec)')
ylabel('Conc. (uM)')
axis(axis + [0 0 -1 1])

nexttile; 
plot(DATA.T,real(DATA.PRED.Im))
xlabel('Time (sec)')
ylabel('Receptor (pA)')
title('ORN Current')
set(gca,'FontSize',14);

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

% figure;
% nexttile; plot(DATA.T, DATA.IPREDn.IL); title('IL')
% nexttile; plot(DATA.T, DATA.IPREDn.ICACL); title('Icacl')
% nexttile; plot(DATA.T, DATA.IPREDn.ICNG); title('Icng')
% nexttile; plot(DATA.T, DATA.IPREDn.PRED_CURRENT); title('PC')

disp('---Done---')