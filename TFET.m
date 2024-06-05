function [tfr_f, tfr_t, Ts1, Ts2, Ts, tfr1] = TFET(x,hlength,Chirp_rate,Type)
%   time-frequency extracting transform.
%   INPUT:
%       x       :	 Signal.
%       hlength :	 Window length.
%       Chirp_rate : Chirp rate estimation method 
%       Type :       TFET or TFET-2     

%   OUTPUT:
%       tfr1  : STFT result.
%       Ts    : TFET result.
%       Ts2   : the part of time-direction extracting result
%       Ts1   : the part of frequency-direction extracting result
%       tfr_t : the part of time-direction segmented result
%       tfr_f : the part of frequency-direction segmented result

% Written by J.L. Wu, Mar. 12, 2024

[xrow,xcol] = size(x);
N=xrow;

if (xcol~=1)
	error('X must be column vector');
end

if (nargin < 2)
	hlength=round(xrow/8);
end

t=1:N;
ft = 1:round(N/2);

hlength=hlength+1-rem(hlength,2);
ht = linspace(-0.5,0.5,hlength);
ht=ht';

BETA = 0.008;
% g
h = exp(-0.5*ht.^2/BETA);
% g'
dh = -ht .* h/BETA; 
% tg
th=ht.*h;

[hrow,~]=size(h); 
Lh=(hrow-1)/2;

% STFT with window of g, g', tg 
tfr1= zeros (N,N);
tfr2= zeros (N,N);
tfr3= zeros (N,N);

va=N/hlength;
    
for icol=1:N
    ti= t(icol); 
    tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
    indices= rem(N+tau,N)+1;
    rSig = x(ti+tau,1);
    tfr1(indices,icol)=rSig.*conj(h(Lh+1+tau));
    tfr2(indices,icol)=rSig.*conj(dh(Lh+1+tau));
    tfr3(indices,icol)=rSig.*conj(th(Lh+1+tau));
end

tfr1=fft(tfr1);
tfr2=fft(tfr2);
tfr3=fft(tfr3);

tfr1=tfr1(1:round(N/2),:);
tfr2=tfr2(1:round(N/2),:);
tfr3=tfr3(1:round(N/2),:);
 
% IF estimator of SET1
IF = (ft-1)'+real(va*1i*tfr2/2/pi./tfr1);

% GD estimator of TET1
GD = t+(hlength-1)*real(tfr3./tfr1);

dt_IF = zeros (round(N/2),N-1);
dt_GD = zeros (round(N/2),N-1);

dw_IF = zeros (round(N/2)-1,N);
dw_GD = zeros (round(N/2)-1,N);

for i=1:N
    dw_IF(:,i)=diff(IF(:,i));
    dw_GD(:,i)=diff(GD(:,i));
end
dw_IF(end+1,:)=dw_IF(end,:);
dw_GD(end+1,:)=dw_GD(end,:);

for i=1:round(N/2)
    dt_IF(i,:)=diff(IF(i,:));
    dt_GD(i,:)=diff(GD(i,:));
end
dt_IF(:,end+1)=dt_IF(:,end);
dt_GD(:,end+1)=dt_GD(:,end);

M=max(max(abs(tfr1)));

% Chirp rate estimate
c1 = dt_IF./dw_IF;   % 1/(β^2 * φ'')
c2 = dw_IF./dt_GD;   % β^2 * φ''^2
c3 = dw_GD./dw_IF;   % 1/φ''
c4 = dt_IF./dt_GD;   % φ''
c5 = dw_GD./dt_GD;   % β^2 * φ''
switch Chirp_rate
    case 1
        c = c1;
        Numerator = dw_IF;
    case 2
        c = sqrt(c2);
        Numerator = dt_GD;
    case 3
        c = c3;
        Numerator = dw_IF;
    case 4
        c = c4;
        Numerator = dt_GD;
    case 5
        c = c5;
        Numerator = dt_GD;
end

tfr_f=zeros(round(N/2),N);
tfr = tfr1;
tfr(find(abs(tfr1)<0.3*M))=0;

tfr_f(find(abs(c)<1 & abs(Numerator)>0.2))=tfr(find(abs(c)<1 & abs(Numerator)>0.2));

tfr_t=tfr-tfr_f;

if Chirp_rate == 1 || Chirp_rate == 3  
% in this case, 1/φ'' is estimated, so tfr_f and tfr_t need to be swapped
    tfr_t = tfr_f;
    tfr_f = tfr - tfr_t;
end

% IF and GD estimators of TFET2
if (strcmp(Type, 'TFET-2'))
    IF = IF-(GD-t).*c4;
    GD = c3.*((ft-1)'-IF) + GD;
end

[Ts1] = SET(x,tfr_f,round(IF));
[Ts2] = TET(x,tfr_t,round(GD));

Ts=Ts1+Ts2;

%the amplitude of tfr result has been pre-rectified for reconstruction
Ts1=Ts1/(sum(h)/2);
Ts2=Ts2/(N/2);
end


%% frequency-direction extracting
function [TFR] = SET(x,tfr,omega)
%   INPUT:
%       x       : Signal.
%       tfr     : STFT result for frequency-direction extracting transform.
%       omega   : IF estimator.

%   OUTPUT:
%       TFR   : SET result.

[nf,nt] = size(tfr);

TFR = zeros(nf,nt);
E=mean(abs(x));

for i=1:nf%frequency
    for j=1:nt%time
         if abs(tfr(i,j))>4*E
            if abs(omega(i,j)-i)<0.5 %default frequency resolution is 1Hz.
                TFR(i,j)=tfr(i,j);
            end
         end
    end
end
end


%% time-direction extracting
function [TFR] = TET(x,tfr,omega)

%   INPUT:
%       x       : Signal.
%       tfr     : STFT result for time-direction extracting transform.
%       omega   : GD estimator.

%   OUTPUT:
%       TFR   : TET result.

[nf,nt] = size(tfr);

E=mean(abs(x));

TFR=zeros(nf,nt);

for i=1:nf%frequency
    for j=1:nt%time
         if abs(tfr(i,j))>7*E
            if abs(omega(i,j)-j)<0.5 %default frequency resolution is 1Hz.
                TFR(i,j)=tfr(i,j);
            end
         end
    end
end
end