% G28 CPCP periodicity analysis
%
% Comparison between different z location. 
%

%% 

CPCPz3 = load('../CPCP_G28_100s_z=3.mat','CPCPt');
CPCPz3 = CPCPz3.CPCPt;

CPCPz2 = load('../CPCP_G28_100s_z=2.mat','CPCPt');
CPCPz2 = CPCPz2.CPCPt;

CPCPb = load('../CPCP_G28_100s_b=3.mat','CPCPt');
CPCPb = CPCPb.CPCPt;

time = 1:100;
figure;
plot(time,CPCPb,time,CPCPz2,time,CPCPz3);
legend({'b=3','z=2','z=3'})

CPCPbfft = fft(CPCPb);
CPCPbfft(1) = [];
n = length(CPCPbfft);
powerb = abs(CPCPbfft(1:floor(n/2))).^2; % power of first half of transform data
maxfreq = 1/2;                   % maximum frequency
freq = (1:n/2)/(n/2)*maxfreq;    % equally spaced frequency grid
periodb = 1./freq;

CPCPz2fft = fft(CPCPz2);
CPCPz2fft(1) = [];
n = length(CPCPz2fft);
powerz2 = abs(CPCPz2fft(1:floor(n/2))).^2; % power of first half of transform data
maxfreq = 1/2;                   % maximum frequency
freq = (1:n/2)/(n/2)*maxfreq;    % equally spaced frequency grid
periodz2 = 1./freq;

CPCPz3fft = fft(CPCPz3);
CPCPz3fft(1) = [];
n = length(CPCPz3fft);
powerz3 = abs(CPCPz3fft(1:floor(n/2))).^2; % power of first half of transform data
maxfreq = 1/2;                   % maximum frequency
freq = (1:n/2)/(n/2)*maxfreq;    % equally spaced frequency grid
periodz3 = 1./freq;


figure;
plot(periodb,powerb,periodz2,powerz2,periodz3,powerz3);
xlim([0 50]); %zoom in on max power
xlabel('Seconds/Cycle')
ylabel('Power')
legend({'b=3','z=2','z=3'})
