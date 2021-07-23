%% Odor-pulse setup
PULSE.ton = [ 0.2000  2 
              0.2000 2
              0.2000 2];
PULSE.toff = [1.2000  3 
              1.2000  3
              1.2000  3];
PULSE.conc = [0  0
              5  50
              50 5];
PULSE.tspan = [0 4];

%% ORN System co-eff
P = struct('Sigma',0.0569, 'cap',0.0039, 'cc1lin',0.7750,...
        'cc2',26.3950,'ck1lin',8.5342,'ck2',0.3069,'clmax',0.9397,...
        'cnmax',0.9663,'cx1lin',1.2307,'cx2',10.9297,'ef',2.7583,...
        'gl',4.9195,'hmc1',1.4829,'hmc2',2.7678,'inf',1.7619,'inhmax',3.5697,...
        'k1',0.1143,'k2lin',12.9344,'kI',10.0453,'kinh',1.0018,'kinhcng',0.5181,...
        'n1',3.1844,'n2',3.1128,'nI',1.9848,'ninh',1.3081,'ninhcng',1.4511,...
        'pd',7.5749,'r1',3.1663,'r2',6.5597,'smax',45.5118,'vcl',-7.7902,...
        'vcng',0.0106,'vl',-44.0413);

% ML Spike sytem co-eff
S = struct;

    % Spike properties
    S.spkThr = -42; % (ORN_rest=-44) mV
    S.maxFR = 75; % Max firing rate Hz
    S.revCp = 0.3; % Reverse coupling from spkV to mem.Voltage

    % Membrane voltage parameters, adapted from (Anderson et. al., 2015)
    S.vCa = 120;                % Rev.Pot for Calcium channels
    S.gCa = 4.4;                % Calcium conductance
    S.vK  = -84;                % Rev.Pot for Potassium channels
    S.gK  =   8;                % Potassium conductance
    S.vL  = -44; % -60          % Rev.Pot for leak channels
    S.gL  =   2;                % Leak channels conductance
    S.Cm  =  20;                % Membrane Conductance
    % Ca2+ ion channel parameters
    S.va = -1.2; S.vb = 18; % phi_m = 2;
    % K+ ion channel parameters
    S.vc = 2 ;%+ S.burst*10; 
    S.vd = 30 ;%- S.burst*12.6; 
    S.phi_n = 0.04 ;%+ S.burst*0.19; % S.phi_n  = 0.02; 
    % Slow current feedback for bursting
    S.epsi_int = .01; S.v0_int = -40;


%% RUN
DATA = simulate_ORN(PULSE,P,S);

%% Plot
figure(1);
clf

% nexttile
% TT = linspace(DATA.T(1),DATA.T(end),100);
% OD = simulate_pulse_train(TT,PULSE.ton,PULSE.toff,PULSE.conc);
% OD = OD(end,:);
% plot(TT,OD,'r-');
% xlabel('Time (sec)')
% ylabel('Conc. (uM)')

nexttile
plot(DATA.T,real(DATA.PRED.ornV))
xlabel('Time (sec)')
ylabel('Olfaction (mV)')
title('ORN Voltage')

nexttile
plot(DATA.T,real(DATA.PRED.spkV))
xlabel('Time (sec)')
ylabel('Spike (mV)')
title('ML Spikes')

nexttile
plot(DATA.T,real(DATA.PRED.sFR))
xlabel('Time (sec)')
ylabel('sFR')

%%

function OUT = simulate_pulse_train(tnow,ton,toff,val,varargin)
%PULSE_TRAIN Generate a train of pulses.
%   V = PULSE_TRAIN(T,TON,TOFF,VAL) returns a pulse train
%   having pulses beginning at times TON and ending at times
%   TOFF with peak values equal to VAL.  Pulses may overlap.
% 
%   PULSE_TRAIN(T,TON,TOFF,VAL,SHARP) uses sharpness values
%   (default = 0.001) for the pulses.  Larger values make 
%   smoother pulse trains.  
%   
%   TON, TOFF, VAL and SHARP must be matrices of equal size.  
%   Each row of these matrices corresponds to a separate train.
%
%   Note:  To simulate pulse trains having different number
%   of pulses simply set the corresponding values in VAL to 
%   zero. 
%   
%   Example: Generate 2 pulse trains each with different 
%      characteristics.
%
%      t = linspace(0,10,100);
%      ton = [0.2 3 1;
%             0.2 5 6]
%      toff = [1.2 3.5 1;
%              2.2 5.5 6.5];
%      val = [1 2 0;
%             4 5 6];  
%
%      sharp = [0.001 0.01 0;
%               0.1 0.1 0.01];
%
%      T = pulse_train(t,ton,toff,val);
%      plot(t,T);
%
%

tnow = tnow(:);

if nargin == 5
	SHARPNESS = varargin{1};
else
	SHARPNESS = 0.001; 
end

NTP = length(tnow);
SIZ = size(ton);

if (NTP > 1)
	ton = repmat(ton,[1 1 length(tnow)]);
	toff = repmat(toff,[1 1 length(tnow)]);
	val = repmat(val,[1 1 length(tnow)]);
	tnow = repmat(reshape(tnow,[1 1 NTP]),[SIZ,1]);
	OUT = permute(sum(val.*(hv(tnow-ton,SHARPNESS)- hv(tnow-toff,SHARPNESS)),2),[1 3 2]);
 
else
	OUT = sum(val.*(hv(tnow-ton,SHARPNESS)- hv(tnow-toff,SHARPNESS)),2);
end

    function OUT = hv(x,SHARPNESS)

        OUT = 1./(1+exp(-x./SHARPNESS));

    end


end