% peakwarm.m
% Kylene Cooley
% 29 Apr 2021
% finds the local maxima in dSST'/dt preceding each event found in SST'

function [dsstM, tdTm] = peakwarm(t,dsstdt,summerDates)
    % find all local maxima in time derivative
    dSSTlm = islocalmax(dsstdt);
    
    % extract all dSST'/dt peak values and times
    dsstM = dsstdt(dSSTlm);
    tdTm = t(dSSTlm);
    
    % Reduce dSST'/dt points to within 24 days before the maximum SST'
    % dates
    maxdate = summerDates;              % latest allowed dates
    mindate = summerDates-days(24);     % earliest allowed dates
    tdTind = [];                        % holds indices found
    
    % loop for summerDate events
    for i=1:length(maxdate)
        k = find((tdTm>mindate(i) & tdTm<maxdate(i)),1,'last');
        tdTind = cat(1,tdTind,k);
    end
    
    % save only these times
    dsstM = dsstM(tdTind);
    tdTm = tdTm(tdTind);
end
