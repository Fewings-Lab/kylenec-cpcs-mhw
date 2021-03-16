% bandpass.m
% Kylene Cooley
% 16 Mar 2021
% bandpass filter a function with pl66

% Bandpass SST'
% 10-day low-pass filter
sstA = pl66tn(sstA,1,240); 
% 6-month (half of a year) low-pass filter
hrs = hours(years(0.5)); % define the cutoff frequency in hours
foo6 = pl66tn(sstA,1,hrs); % Evaluate the low-pass filtered signal

% Take high-pass part of signal
sstA = sstA-foo6;
% Replace one window-length with NaNs on each end
sstA(1:2*round(hrs))=NaN;
sstA(end-2*round(hrs):end)=NaN;
