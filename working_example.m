function working_example()

global X
X=0;

% Add labels and init-values for new variable
FN = {'y_fxn','y_track'}; % 
yinit = {0,1};

% Reduced complexity without the changing structure of the model
% Removed parallel computation for multiple ORNs
yinit = cell2struct(yinit,FN,2);
var_names = fieldnames(yinit);
init_vals = struct2cell(yinit);
init_vals = [init_vals{:}];
NVAR = length(var_names);
NCURVE = 1;
NEQ = NVAR*NCURVE;
JP = spdiags(ones(NEQ,2*NVAR-1),[-(NEQ-NCURVE):NCURVE:(NEQ-NCURVE)],NEQ,NEQ);

%% run simulation
tspan = [-pi pi];
ODEOPTS = odeset('JPattern','on','MaxStep',0.9);
tic % start timer
[T,Y] = ode15s(@(t,y) SYSTEM(t,y,ODEOPTS,JP), tspan, init_vals);
disp(['Simulation Time :',num2str(toc)])

%% plot
figure(1); clf;
nexttile;
plot(T,Y(:,1))
title('System')

nexttile;
plot(T,Y(:,2))
title('Expected : Triangular or square wave')
end


function dy = SYSTEM(t,y,ODEOPTS,JP) 
    global X


    if strcmp(ODEOPTS,'jpattern')
		dy = JP;
        return;
    end
    % Demo system
    D_y_fxn = -cos(t); % To get a sine-wave

    % Tracking logic
    base = 1; rise = 2;
    D_y_track = rise.*(D_y_fxn > 0.1) ...
        - 4*rise.*(D_y_fxn < 0.1 && y(2)>base);
    X = [X X(end)+1];
    


    dy = [D_y_fxn; D_y_track];
end
