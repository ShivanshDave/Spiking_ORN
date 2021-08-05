function spikes = script_spikes_ID(V,T,spikeThr)

sz_X = size(V,2);
spikes=cell(sz_X,3);

for i=1:sz_X
    % Find Peaks
    xval=V(:,i);
%     Fxval = bandpass_filter(xval,1/mean(diff(T)),10,50);
    [val,index]=LocalMaxima(xval);   
    
    % Filter Peaks
    refTh = 10;
    index(val<=spikeThr) = [];          % Apply Vm threshold
    val(val<=spikeThr) = [];
    ind=2; 
    while ind <= length(index)          % Apply refractory period      
       if abs(index(ind)-index(ind-1)) < refTh      
           index(ind)=[];
           val(ind)=[];
       else
           ind=ind+1;
       end
    end
    
    % Save Spikes
    spikes{i,1} = index;                % Save Spike index
    spikes{i,2} = xval(index);                  % Save Spike amplitude (mV)
    spikes{i,3} = T(index);             % Save Spike time (ms)
end

end


function [y] = bandpass_filter(x,Fs,fstart,fstop)
%
%
% Dinesh Natesan

% 
highpass = designfilt('highpassiir','FilterOrder',4, ...
    'PassbandFrequency',fstart,'PassbandRipple',0.1, ...
    'SampleRate',Fs);
lowpass = designfilt('lowpassiir','FilterOrder',4, ...
    'PassbandFrequency',fstop,'PassbandRipple',0.1, ...
    'SampleRate',Fs);

% x = to_column_matrix(x);    % convert to column matrix
y = nan(size(x));

for i=1:size(x,2)
    
    tx = filtfilt(highpass, x(:,i));
    y(:,i) = filtfilt(lowpass, tx);

end
end

function [lmval,indd]=LocalMinima(xx,filt)
%LMIN 	function [lmval,indd]=lmin(x,filt)
%	Find local minima in vector X, where LMVAL is the output
%	vector with minima values, INDD is the corresponding indeces
%	FILT is the number of passes of the small running average filter
%	in order to get rid of small peaks.  Default value FILT =0 (no
%	filtering). FILT in the range from 1 to 3 is usially sufficient to
%	remove most of a small peaks
%	Examples:
%	xx=0:0.01:35; y=sin(xx) + cos(xx ./3);
%	plot(xx,y); grid; hold on;
%	[a b]=lmin(y,2)
%	 plot(xx(a),y(a),'r+')
%	see also LMAX, MAX, MIN

%
%**************************************************|
% 	Serge Koptenko, Guigne International Ltd., |
%	phone (709)895-3819, fax (709)895-3822     |
%--------------06/03/97----------------------------|

x=xx;
len_x = length(x);
fltr=[1 1 1]/3;
if nargin <2, filt=0;
else
    x1=x(1); x2=x(len_x);
    
    for jj=1:filt,
        c=conv(fltr,x);
        x=c(2:len_x+1);
        x(1)=x1;
        x(len_x)=x2;
    end
end

lmval=[];
indd=[];
i=2;		% start at second data point in time series

while i < len_x-1,
    if x(i) < x(i-1)
        if x(i) < x(i+1)	% definite min
            lmval =[lmval x(i)];
            indd = [ indd i];
            
        elseif x(i)==x(i+1)&&x(i)==x(i+2)	% 'long' flat spot
            %lmval =[lmval x(i)];	%1   comment these two lines for strict case
            %indd = [ indd i];	%2 when only  definite min included
            i = i + 2;  		% skip 2 points
            
        elseif x(i)==x(i+1)	% 'short' flat spot
            %lmval =[lmval x(i)];	%1   comment these two lines for strict case
            %indd = [ indd i];	%2 when only  definite min included
            i = i + 1;		% skip one point
        end
    end
    i = i + 1;
end

if filt>0 && ~isempty(indd),
    if (indd(1)<= 3)||(indd(length(indd))+2>length(xx)),
        rng=1;	%check if index too close to the edge
    else rng=2;
    end
    
    for ii=1:length(indd),
        [val(ii) iind(ii)] = min(xx(indd(ii) -rng:indd(ii) +rng));
        iind(ii)=indd(ii) + iind(ii)  -rng-1;
    end
    indd=iind; lmval=val;
else
end

end

function [lmval,indd]=LocalMaxima(xx,filt)
%LMAX 	[lmval, indd]=lmax(xx,filt). Find local maxima in vector XX,where
%	LMVAL is the output vector with maxima values, INDD  is the
%	corresponding indexes, FILT is the number of passes of the small
%	running average filter in order to get rid of small peaks.  Default
%	value FILT =0 (no filtering). FILT in the range from 1 to 3 is
%	usially sufficient to remove most of a small peaks
%	For example:
%	xx=0:0.01:35; y=sin(xx) + cos(xx ./3);
%	plot(xx,y); grid; hold on;
%	[b,a]=lmax(y,2)
%	 plot(xx(a),y(a),'r+')
%	see also LMIN, MAX, MIN

%**************************************************|
% 	Serge Koptenko, Guigne International Ltd., |
%	phone (709)895-3819, fax (709)895-3822     |
%--------------06/03/97----------------------------|

x=xx;
len_x = length(x);
fltr=[1 1 1]/3;
if nargin <2, filt=0;
else
    x1=x(1); x2=x(len_x);
    for jj=1:filt,
        c=conv(fltr,x);
        x=c(2:len_x+1);
        x(1)=x1;
        x(len_x)=x2;
    end
end
lmval=[]; indd=[];
i=2;		% start at second data point in time series
while i < len_x-1,
    if x(i) > x(i-1)
        if x(i) > x(i+1)	% definite max
            lmval =[lmval x(i)];
            indd = [ indd i];
        elseif x(i)==x(i+1)&&x(i)==x(i+2)	% 'long' flat spot
            %lmval =[lmval x(i)];  	%1   comment these two lines for strict case
            %indd = [ indd i];	%2 when only  definite max included
            i = i + 2;  		% skip 2 points
        elseif x(i)==x(i+1)	% 'short' flat spot
            %lmval =[lmval x(i)];	%1   comment these two lines for strict case
            %indd = [ indd i];	%2 when only  definite max included
            i = i + 1;		% skip one point
        end
    end
    i = i + 1;
end
if filt>0 && ~isempty(indd),
    if (indd(1)<= 3)||(indd(length(indd))+2>length(xx)),
        rng=1;	%check if index too close to the edge
    else rng=2;
    end
    for ii=1:length(indd), 	% Find the real maximum value
        [val(ii) iind(ii)] = max(xx(indd(ii) -rng:indd(ii) +rng));
        iind(ii)=indd(ii) + iind(ii)  -rng-1;
    end
    indd=iind; lmval=val;
else
end

end