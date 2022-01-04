% autocorrTL.m
% Kylene Cooley
% May 2020
% Revised 16 Mar 2021
% Computes the symmetric time-lagged autocorrelation given an array and the maximum
% time lag desired
function rhoyy = autocorrTL(Y,maxlag)
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
