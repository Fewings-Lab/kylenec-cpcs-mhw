% acov.m
% Kylene Cooley
% Created 4 May 2020
% autocovariance

function R = acov(xx,yy)
    S = (xx-nanmean(xx)).*(yy-nanmean(yy)); 
    R = sum(S,'omitnan')/sum(~isnan(S)); % Normalized by N-1 excluding missing values
end