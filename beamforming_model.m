close all
clear all
clc

load('rx_pattern.mat');
load('tx_pattern.mat');

fs = 375e6;
B = 375e6;
T = 20e-6;
alpha = B/T;
c = 3e8;
f0 = 9.6e9;
P = 1;
N = floor(fs*T);
Nd = N*(P+1);
Nant = 4;
t = (0:N-1)/fs;
tau = (0:P-1)*1e-6;
H = 550e3;
look_angle = deg2rad(20);
thb = deg2rad(2.4);
rp_step = 0.02;

rmin = H/cos(look_angle - thb/2);
r = rmin + c*tau/2;
th_rad = acos(H./r);
th_deg = rad2deg(th_rad);
th_deg = floor(th_deg*50)/50;

indx = round(th_deg/rp_step);

rx_db = [amp amp1 amp2 amp3];
a_db = rx_db + tx_amp;
a = db2mag(a_db);
w = zeros(length(a),Nant);

for i = 1:length(a)
    w(i,:) = a(i,:)/(a(i,:)*a(i,:)');
end

W = w(indx,:);
G = a(indx,:);

s_rx = exp(2*pi*1i*(-B/2*(t-tau') + 0.5*alpha*(t-tau').^2)+1i*f0*tau');

tndx = floor(fs*tau)+1;
s1 = zeros(P,Nd);
s2 = zeros(P,Nd);
s3 = zeros(P,Nd);
s4 = zeros(P,Nd);

for i = 1:P
    ind = tndx(i):N+tndx(i)-1;
    s1(i,ind) =  s_rx(i,:).*G(i,1);
    s2(i,ind) =  s_rx(i,:).*G(i,2);
    s3(i,ind) =  s_rx(i,:).*G(i,3);
    s4(i,ind) =  s_rx(i,:).*G(i,4);
end

s_ant1 = sum(s1,1);
s_ant2 = sum(s2,1);
s_ant3 = sum(s3,1);
s_ant4 = sum(s4,1);

% s_comb = s_ant1+s_ant2+s_ant3+s_ant4;
% plot(real(s_comb));

%Filter design
M = 9;
z_delay = T/M;
Nf = floor(z_delay*fs);
Ncoef = 31;
h = zeros(Ncoef,M);

Nfft = length(s_ant1);
ff = (-Nfft/2:Nfft/2-1)*fs/Nfft;

figure
for i = 0:M-1 
%hb(n+1,m+1) = (B/(M*Fs))*sinc((n-(Ncoef-1)/2)*B/(M*Fs))*exp(1i*(n-(Ncoef-1)/2)*pi*((2*m+1)/M-1)*B/Fs);
    h(:,i+1) = B/(M*fs)*sinc((-(Ncoef+1)/2+1:(Ncoef-1)/2)*B/M/fs).*exp(1i*pi*(-(Ncoef+1)/2+1:(Ncoef-1)/2)*((2*i+1)/M-1)*B/fs);
    s_filtered1(i+1,:) = conv(s_ant1,h(:,i+1), 'same');
    s_filtered2(i+1,:) = conv(s_ant2,h(:,i+1), 'same');
    s_filtered3(i+1,:) = conv(s_ant3,h(:,i+1), 'same');
    s_filtered4(i+1,:) = conv(s_ant4,h(:,i+1), 'same');
    S(i+1,:) = fftshift(fft(s_filtered4(i+1,:)));
    subplot(M,1,i+1);
    plot(ff,abs(S(i+1,:)));
%     plot(ff,real(s_filtered4(i+1,:)));
end


figure
plot(ff,abs(fftshift(fft(s_ant4))));
% plot(real(s_ant4));

figure
plot(ff,(abs(sum(S,1))));
% plot(real(sum(s_filtered4,1)));

Hf = sum(S,1);
aHf = 20*log10(abs(Hf));
figure
plot(ff,aHf, '-r');
hold on
plot(ff,20*log10(abs(S)));
grid on
hold off


