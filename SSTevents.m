% SSTevents.m
% Kylene Cooley
% 29 Apr 2021
% ID the peak SST' events during the upwelling season for a time series and
% the maxima of the time derivative within 24 days before the event.

function [summerSSTa,summerDates,dsstM,tdTm] = SSTevents(t,sstA,dsstdt)
    % find all local maxima in SST'
    SSTlm = islocalmax(sstA);
    
    % limit max SST' to SST' above two standard deviations
    muA = mean(sstA,'omitnan');
    sigA = std(sstA,'omitnan');
    warmest = sstA>(muA+2*sigA);
    sstM = sstA(SSTlm & warmest);
    
    % find corresponding times of max SST' for plotting
    tTm = t(SSTlm & warmest);
    
    % limit SST' events to the upwelling season
    dv = datevec(t);
    monthnum = dv(SSTlm & warmest, 2);  % month number of all events
    summerSSTa = sstM(monthnum<3|monthnum==12);
    summerDates = tTm(monthnum<3|monthnum==12);
    
    % use summer dates to find the maxima dSST'/dt before each upwelling
    % season event
    [dsstM,tdTm] = peakwarm(t,dsstdt,summerDates);
end