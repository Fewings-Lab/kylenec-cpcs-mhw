% bandpassV3.m
% Kylene Cooley
% 31 Mar 2021
% bandpass filter a function with pl66

function [low1, high1, low2] = bandpassV3(dat,dt,TcL,TcH)

    % Bandpass SST'
    % low-pass filter with lower cutoff period in hours
    low1 = nan(size(dat)); % Initialize empty matrix to hold low-pass part
    for m=1:length(dat(:,1,1)) % loop along one spatial dimension since pl66 can only filter 2D arrays,
        % we can use either as long as time is the longer dimension
        low1(m,:,:) = pl66tn(squeeze(dat(m,:,:)),dt,TcL)'; % Evaluate the low-pass filtered signal  and assign to 2D longitude slice
    end

    % high pass part for variability maps
    high1 = dat-low1;
    % Replace one window-length with NaNs on each end
    high1(1:2*round(TcL))=NaN;
    high1(end-2*round(TcL):end)=NaN;

    % low-pass filter with higher cutoff period in
    % years
    hrs = hours(years(TcH)); % define the cutoff frequency in hours

    low2 = nan(size(dat)); % Initialize empty matrix to hold low-pass part
    for m=1:length(dat(:,1,1)) % loop along one spatial dimension since pl66 can only filter 2D arrays,
        % we can use either as long as time is the longer dimension
        low2(m,:,:) = pl66tn(squeeze(low1(m,:,:)),dt,hrs)'; % Evaluate the low-pass filtered signal  and assign to 2D longitude slice
    end 

%     % Take high-pass part of signal
%     datbp = low1-low2;
%     % Replace one window-length with NaNs on each end
%     datbp(1:2*round(hrs))=NaN;
%     datbp(end-2*round(hrs):end)=NaN;

end