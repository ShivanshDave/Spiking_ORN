function DATA = simulate_ORN(PULSE) 
% Initialize
simulation_parameters;

init_bLR    = 1.e-8; %1
init_aG     = 1.e-8; %2
init_cAMP   = 1.e-8; %3
init_Ca     = 1.e-8; %4
init_CAMK   = 1.e-8; %5
init_CaCAM  = 1.e-8; %6
init_IX     = 1.e-8; %7
init_ornV   = -44;   %8
init_spkV   = -47.14; %9 ML dV/dt
init_nk     = 0.04;  %10 ML dN/dt 
init_CaFR   = 1;     %11 ML d(slope)/dt


yinit = {init_bLR, init_aG, init_cAMP, init_Ca,...
        init_CaCAM,init_CAMK,init_IX,init_ornV,...
        init_spkV, init_nk, init_CaFR};
FN = {'bLR','aG','cAMP','Ca',...
    'CaCaM','aCaMK','IX','ornV',...
    'spkV','nK','CaFR'};
yinit = cell2struct(yinit,FN,2);
N = size(PULSE.ton,1);
var_names = fieldnames(yinit);
init_vals = struct2cell(yinit);
init_vals = [init_vals{:}];
init_vals = repmat(init_vals(:)',N,1);
init_vals = init_vals(:)';		

NVAR = length(var_names);
NCURVE = length(PULSE.ton(:,1));
NEQ = NVAR*NCURVE;
JP = spdiags(ones(NEQ,2*NVAR-1),[-(NEQ-NCURVE):NCURVE:(NEQ-NCURVE)],NEQ,NEQ);

% SImulate
%Solve for 6 seconds so that we may come to a steady-state from our initial conditions.
tspan = [PULSE.tspan'];
ODEOPTS = odeset('JPattern','on','MaxStep',0.9);

tic % start timer
[T,Y] = ode15s(@(t,y) SYSTEM(t,y,ODEOPTS,PULSE,P,S,N,JP), tspan, init_vals);

PRED = [];
for j = 1:length(var_names)
    PRED.(var_names{j}) = Y(:,((j-1)*N+1):(j*N));
end

% Compute currents
Incx = [];
inhcng = 1+(P.inhmax-1).*PRED.CaCaM.^P.ninhcng./(PRED.CaCaM.^P.ninhcng + P.kinhcng.^P.ninhcng);
%inhcng = 1;
Icng = (P.cnmax.*PRED.cAMP.^P.n1./(PRED.cAMP.^P.n1 + (inhcng.*P.hmc1).^P.n1)).*(P.vcng-PRED.ornV);
Icacl = (P.clmax.*PRED.Ca.^P.n2./(PRED.Ca.^P.n2 + P.hmc2.^P.n2)).*(P.vcl-PRED.ornV);
Il = P.gl.*(P.vl-PRED.ornV);        

% predicted currents
IPRED.PRED_CURRENT = (-Icng - Icacl);
IPRED.ICNG = -Icng;
IPRED.ICACL = -Icacl;
IPRED.INCX = -Incx;
IPRED.IL = -Il;

% normalized currents
NORMALIZING_FACTOR = 70;
IPREDn.ICNG = IPRED.ICNG/NORMALIZING_FACTOR;
IPREDn.ICACL = IPRED.ICACL/NORMALIZING_FACTOR;
IPREDn.INCX = IPRED.INCX/NORMALIZING_FACTOR;
IPREDn.IL = IPRED.IL/NORMALIZING_FACTOR;	
IPREDn.PRED_CURRENT = IPRED.PRED_CURRENT/NORMALIZING_FACTOR;
IPREDn.NORMALIZING_FACTOR = NORMALIZING_FACTOR;

DATA.PULSE = PULSE;
DATA.PRED = PRED;
% DATA.IPRED = IPRED;
DATA.IPREDn = IPREDn;
DATA.IPREDn.NF = NORMALIZING_FACTOR;
DATA.T = T;
DATA.P = P;
DATA.S = S;

disp(toc)

end

function dy = SYSTEM(t,y,ODEOPTS,PULSE,P,S,N,JP)
    
    if toc > 15
        disp(t)
        if ~strcmp(input('Run for 15 more sec? y/n : ','s'),'y')
            error('Timeout');
        else 
            tic;
        end
    end
    
    if strcmp(ODEOPTS,'jpattern')
		dy = JP;
    else
 
        %####--- ORN Transduction ----####%
        read = @(N,Y,i) y((i-1)*N+1:(i)*N,1);
        bLR     = read(N,y,1);
        aG      = read(N,y,2);
        cAMP    = read(N,y,3);
        Ca      = read(N,y,4);
        CaCAM   = read(N,y,5);
        CAMK    = read(N,y,6);
        IX      = read(N,y,7);
        ornV    = read(N,y,8);  
        
        % Define the tot variables to be 1.
        P.Rtot = 1;
        P.Gtot = 1;

        % Ligand-Receptor interaction
        k2 = P.k2lin.*bLR;
        % Transduction
        cc1 = P.cc1lin.*Ca;
        cx1 = P.cx1lin.*Ca;
        ck1 = P.ck1lin.*CaCAM;
        fca = (P.smax)./(1+(CAMK./(P.kinh)).^P.ninh);
        synth = aG.*fca;
        inhcng = 1+(P.inhmax-1).*CaCAM.^P.ninhcng./(CaCAM.^P.ninhcng + P.kinhcng.^P.ninhcng);	
        Icng = (P.cnmax.*cAMP.^P.n1./(cAMP.^P.n1 + (inhcng.*P.hmc1).^P.n1)).*(P.vcng-ornV);		
        Icacl = ((P.clmax.*Ca.^P.n2)./(Ca.^P.n2 + P.hmc2.^P.n2)).*(P.vcl-ornV);
        Il = P.gl.*(P.vl-ornV);

        % Odor stimulation & Ligand-Receptor interaction 
        hv = @(x) 1./(1+exp(-x./0.001));
        Ostim = sum(PULSE.conc.*(hv(t-PULSE.ton)- hv(t-PULSE.toff)),2);        
        D_bLR = P.k1*Ostim.*(P.Rtot-bLR) - P.r1.*bLR;
        D_aG = k2.*(P.Gtot-aG) - P.r2.*aG;
        % Transduction
        D_cAMP = synth - P.pd.*cAMP;
        D_Ca = P.inf.*Icng - (P.ef./(1 + (IX./P.kI).^P.nI)).*Ca + (-cc1 + P.cc2.*CaCAM);
        D_CaCAM = cc1 - P.cc2.*CaCAM;
        D_CAMK = ck1 - P.ck2.*CAMK;
        D_IX = cx1 - P.cx2.*IX;  %This has got to go back down in order for oscillations...
        D_txnV = (1./P.cap).*(Icng + Icacl + Il);
        
        %####--- ML Spike ----####%
        spkV    = read(N,y,9);
        nK      = read(N,y,10);
        CaFR    = read(N,y,11);

        % Match ML_SPK time with ORN_SYSTEM time
        ct = 1e3; % convert ms -> s 
        ct = ct*S.maxFR/10; % Default T=100ms,FR=10

        % Ca2+ dependent firing rate modulation
        ct = ct./(CaFR);
        D_CaFR = (1+Ca).*( S.pCaFR.*(D_Ca>0) ...
        - S.nCaFR.*(D_Ca<0 & CaFR>S.mCaFR) );
        % Ca2+ current
        Ica = S.gIca*D_Ca./(1+Ca); %nA

        % ML Ca2+ ion channel
        xi_m    = @(v) (v-S.va)/S.vb;                   % scaled argument for m-gate input
        minf    = @(v) 0.5*(1+tanh(xi_m(v)));       % m-gate activation function    
        % ML K+ ion channel
        xi_n    = @(v) (v-S.vc)./S.vd;                   % scaled argument for n-gate input
        ninf    = @(v) 0.5*(1+tanh(xi_n(v)));       % n-gate activation function
        tau_n   = @(v) 1./(S.phi_n.*cosh(xi_n(v)/2));  % n-gate activation t-const
        D_nK = ct.*(ninf(spkV)-nK)./tau_n(spkV);

        % Other ion channels activation
        Iion = S.gIion*(ornV>S.spkThr); % nA Thr@vL=(88,-60),(57,-44)    

        D_spkV = (ct/S.Cm).*( Iion + Ica ... 
            - S.gL*(spkV-S.vL) ... 
            - S.gK*nK.*(spkV-S.vK) ...
            - S.gCa*minf(spkV).*(spkV-S.vCa) );

        % spike reverse coupling
        D_ornV = D_txnV + (S.revCp).*D_spkV;
                
        %--- Passon derivatives ----%
        dy = [D_bLR;D_aG;D_cAMP;D_Ca;...
            D_CaCAM;D_CAMK;D_IX;D_ornV;...
            D_spkV;D_nK;D_CaFR]; 

    end
end
