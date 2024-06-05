function [tfr_1, tfr_2, Ts1, Ts2, Ts,tfr1] = TFMST_Y(x,hlength)
%   time-frequency multisynchrosqueezing transform.
%	x       : Signal.
%	hlength : Window length.

%   Ts    : time-frequency multisynchrosqueezing transform.
%	tfr   : Time-Frequency Representation.

%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
%  Written by YuGang.

[xrow,xcol] = size(x);
N=xrow;

if (xcol~=1),
 error('X must be column vector');
end;

if (nargin < 2),
 hlength=round(xrow/8);
end;

t=1:N;
ft = 1:round(N/2);

hlength=hlength+1-rem(hlength,2);
ht = linspace(-0.5,0.5,hlength);ht=ht';

% Gaussian window
h = exp(-0.5*ht.^2/0.008);
% derivative of window
dh = -ht .* h/0.008;  % g'
%
th=h.*ht;

[hrow,~]=size(h); Lh=(hrow-1)/2;

tfr1= zeros (N,N);
tfr2= zeros (N,N);
tfr3= zeros (N,N);

va=N/hlength;
    
for icol=1:N
ti= t(icol); tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
indices= rem(N+tau,N)+1;
rSig = x(ti+tau,1);
tfr1(indices,icol)=rSig.*conj(h(Lh+1+tau));
tfr2(indices,icol)=rSig.*conj(dh(Lh+1+tau));
tfr3(indices,icol)=rSig.*conj(th(Lh+1+tau));
end;

tfr1=fft(tfr1);
tfr2=fft(tfr2);
tfr3=fft(tfr3);

tfr1=tfr1(1:round(N/2),:);
tfr2=tfr2(1:round(N/2),:);
tfr3=tfr3(1:round(N/2),:);

omega1 = zeros(round(N/2),N);
omega2= zeros(round(N/2),N);

for b=1:N
    omega1(:,b) = (ft-1)'+real(va*1i*tfr2(ft,b)/2/pi./tfr1(ft,b));
end

for a=1:round(N/2)
    omega2(a,:) = t+(hlength-1)*real(tfr3(a,t)./tfr1(a,t));
end

omega_t_dt = zeros (round(N/2),N-1);
omega_w_dt = zeros (round(N/2),N-1);

for i=1:round(N/2)
    omega_w_dt(i,:)=diff(omega1(i,:));
    omega_t_dt(i,:)=diff(omega2(i,:));
end
omega_t_dt(:,end+1)=omega_t_dt(:,end);
omega_w_dt(:,end+1)=omega_w_dt(:,end);

mm=max(max(abs(tfr1)));

c_rate_norm1=omega_w_dt./omega_t_dt;
 
tfr_1=zeros(round(N/2),N);
tfr1(find(abs(tfr1)<0.2*mm))=0;

tfr_1(find(abs(c_rate_norm1)<1 & abs(omega_t_dt)>0.2))=tfr1(find(abs(c_rate_norm1)<1 & abs(omega_t_dt)>0.2));% omega_t_dt should  be larger enough than zero.

tfr_2=tfr1-tfr_1;

[Ts1] = TFMST1_Y(tfr_1,round(omega1));
[Ts2] = TFMST2_Y(tfr_2,round(omega2));

Ts=Ts1+Ts2;
Ts1 = Ts1/(N/2);
Ts2 = Ts2/(sum(h)/2);
end


%% ---------------- % Computes the MSST.----------------
function [Ts] = TFMST1_Y(tfr,omega)

%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
%  Written by YuGang.

[~,N]=size(tfr);
Ts = zeros (round(N/2),N);
omega2=zeros(round(N/2),N);
[neta,nb]=size(tfr);
num=10;

if num>1
    for kk=1:num-1
        for b=1:nb
            for eta=1:neta
                k = omega(eta,b);
                if k>=1 && k<=neta
                    omega2(eta,b)=omega(k,b);
                end
            end
        end
        omega=omega2;
    end
else
    omega2=omega;
end

for b=1:nb%time
    % Reassignment step
    for eta=1:neta%frequency
         if abs(tfr(eta,b))>0.0001
            k = omega2(eta,b);
            if k>=1 && k<=neta
            Ts(k,b) = Ts(k,b) + tfr(eta,b);
            end
         end
    end
end

end


%% ---------------- % Computes the TMSST.----------------
function [Ts] = TFMST2_Y(tfr,omega)

%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

[~,N]=size(tfr);

omega2=omega;
omega22=zeros(round(N/2),N);
[neta,nb]=size(tfr);
num=10;

if num>1
    for kk=1:num-1
        for b=1:nb
            for eta=1:neta
                k2 = round(omega2(eta,b));
                if k2>=1 && k2<=nb
                    omega22(eta,b)=omega2(eta,k2);
                end
            end
        end
        omega2=omega22;
    end
else
    omega22=omega2;
end
omega22=round(round(omega22*2)/2);

t=1:N;
tfr2 = zeros(round(N/2),N);
    
for eta=1:round(N/2)%frequency
	tfr2(eta,:)=tfr(eta,:).*exp(-1j * 2 * pi*(eta-1)*(t/N));
end

Ts = zeros(round(N/2),N);
% Reassignment step
for b=1:N%time
    for eta=1:round(N/2)%frequency
         if abs(tfr(eta,b))>0.0001
            k2 = omega22(eta,b);
            if k2>=1 && k2<=N
                Ts(eta,k2) = Ts(eta,k2) + (tfr2(eta,b));
            end
         end
    end
end
end