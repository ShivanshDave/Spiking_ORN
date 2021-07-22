%% Odor-pulse setup
PULSE.ton = [ 0.2000  2 
              0.2000 2];
PULSE.toff = [1.2000  3 
              1.2000  3];
PULSE.conc = [0  0
              2  50];
PULSE.tspan = 0:0.001:4;

%% RUN
DATA = simulate_ORN(PULSE);

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
plot(DATA.T,real(DATA.PRED.V))
xlabel('Time (sec)')
ylabel('Mem. volt. (mV)')

nexttile
plot(DATA.T,real(DATA.PRED.spkV))
xlabel('Time (sec)')
ylabel('Mem. volt. (mV)')

nexttile
plot(DATA.T,real(DATA.PRED.nK))
xlabel('Time (sec)')
ylabel('nK [0:1]')

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