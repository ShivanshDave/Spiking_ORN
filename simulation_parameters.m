%% ORN System co-eff
P = struct('Sigma',0.0569, 'cap',0.0039, 'cc1lin',0.7750,...
        'cc2',26.3950,'ck1lin',8.5342,'ck2',0.3069,'clmax',0.9397,...
        'cnmax',0.9663,'cx1lin',1.2307,'cx2',10.9297,'ef',2.7583,...
        'gl',4.9195,'hmc1',1.4829,'hmc2',2.7678,'inf',1.7619,'inhmax',3.5697,...
        'k1',0.1143,'k2lin',12.9344,'kI',10.0453,'kinh',1.0018,'kinhcng',0.5181,...
        'n1',3.1844,'n2',3.1128,'nI',1.9848,'ninh',1.3081,'ninhcng',1.4511,...
        'pd',7.5749,'r1',3.1663,'r2',6.5597,'smax',45.5118,'vcl',-7.7902,...
        'vcng',0.0106,'vl',-44.0413);

%% ML Spike sytem co-eff
S = struct;

% Spike properties
S.spkThr = -43.25; % (ORN_rest=-44) mV
S.maxFR = 35; % Max firing rate Hz
S.revCp = 7; % Reverse coupling from spkV to mem.Voltage

% Ca2+ dependent firing rate modulation 
S.mCaFR = 1; S.pCaFR = 2; S.nCaFR = 100; %(min,pos,neg)
S.gIca = 10; % Ca@+ current gain
S.gIion = 25; % other ion channels activation

% ML Mem. volt. parameters, adapted from (Anderson et. al., 2015)
S.vCa = 120;                % Rev.Pot for Calcium channels
S.gCa = 4.4;                % Calcium conductance
S.vK  = -84;                % Rev.Pot for Potassium channels
S.gK  =   8;                % Potassium conductance
S.vL  = -44; % -60          % Rev.Pot for leak channels
S.gL  =   2;                % Leak channels conductance
S.Cm  =  20;                % Membrane Conductance
% Ca2+ ion channel parameters
S.va = -1.2; S.vb = 18;
% K+ ion channel parameters
S.vc = 2; S.vd = 30; S.phi_n = 0.04;
