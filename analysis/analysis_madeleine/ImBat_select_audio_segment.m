% Investigating the mic data

load('audioConCat_1.mat');
figure(); plot(audioConCat);
Ge200305_small_seg = audioConCat(200000000:220000000);
Ge200305_small_seg = audioConCat(350000000:400000000);

figure(); plot(Ge200305_small_seg);
save('Ge200305_mic1_seg1','Ge200305_small_seg');

thresh=0.0013; %earthworks mics
thresh2=0.005; %voltage threshold for echolocation clicks for knowles mics
thresh3 = 0.0025;
calldist=0.015; %minimal distance between two clicks (in seconds)
errorval=0.005; 
fs=192000;
[b,a]=butter(8,2*[10e3 40e3]./fs,'bandpass');
convy=filter(b,a,Ge200305_small_seg);

[pks,indxc]=findpeaks(convy,'MinPeakHeight',thresh3);
pk_vect = NaN(size(Ge200305_small_seg,1),1);
pk_vect(indxc)=0.01;
figure(); hold on; plot(Ge200305_small_seg);
plot(pk_vect,'*r');

% Plot the FFT of the signal
y = fft(Ge200305_small_seg);
n = length(Ge200305_small_seg);          % number of samples
f = (0:n-1)*(fs/n);     % frequency range
power = abs(y).^2/n;    % power of the DFT

plot(f,power)
xlabel('Frequency')
ylabel('Power')

% Plot the circshifted FFT off the signal
y0 = fftshift(y);         % shift y values
f0 = (-n/2:n/2-1)*(fs/n); % 0-centered frequency range
power0 = abs(y0).^2/n;    % 0-centered power

plot(f0,power0)
xlabel('Frequency')
ylabel('Power')

% #1
m = length(Ge200305_small_seg);       % original sample length
n = pow2(nextpow2(m));  % transform length
y = fft(Ge200305_small_seg,n); 

f = (0:n-1)*(fs/n)/10;
power = abs(y).^2/n;      

plot(f(1:floor(n/2)),power(1:floor(n/2)))
xlabel('Frequency')
ylabel('Power')
