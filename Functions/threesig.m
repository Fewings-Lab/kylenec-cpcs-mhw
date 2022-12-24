% threesig.m

% 3-sigma filter
function signal = threesig(y)
    for h=1:3
        M = nanmean(y);
        S = nanstd(y,1); % Unbiased normalization
        y(abs(y-M)>3*S) = NaN; 
    end
    signal = y;
end
