function rhoyy = autocorrtl(Y,maxlag)
    rhoyy = zeros((2*maxlag+1),1);
    Ntot = length(Y);
    % compute and plot sample time-lagged autocorrelations:
    % Loop for time lags, tau = -maxlag to +maxlag
    for tau = 0:maxlag
        % indices of values to extract from Y(t)
        nn1 = 1:1:Ntot-tau;
        nn2 = 1+tau:1:Ntot;
        % use tau to extract xx and yy
        xx = Y(nn1);
        yy = Y(nn2);
        % find Ryy
        Ryy = acov(xx,yy);
        % get rhoyy 
        rhoyy(maxlag+1+tau) = acor(Ryy,xx,yy);
        if tau > 0
            % negative lag
            % Switch order of xx and yy in autocovariance and autocorrelation
            % find Ryy
            Ryy = acov(yy,xx);
            % get rhoyy
            rhoyy(maxlag+1-tau) = acor(Ryy,yy,xx);
        end
    end
end

% 3-sigma filter
function signal = threesig(y)
    for h=1:3
        M = nanmean(y);
        S = nanstd(y,1); % Unbiased normalization
        y(abs(y-M)>3*S) = NaN; 
    end
    signal = y;
end

% autocovariance
function R = acov(xx,yy)
    S = (xx-nanmean(xx)).*(yy-nanmean(yy)); 
    R = sum(S,'omitnan')/sum(~isnan(S)); % Normalized by N-1 excluding missing values
end

% autocorrelation
function rho = acor(Ryy,xx,yy)
    rho = Ryy/(nanstd(xx,1)*nanstd(yy,1));
end