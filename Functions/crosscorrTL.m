% crosscorrTL.m
% Kylene Cooley
% 16 Mar 2021
% Computes the symmetric time-lagged cross-correlation given two arrays and the maximum
% time lag desired. X and Y must be the same length??

function rhoyy = crosscorrTL(X,Y,maxlag)
    rhoyy = zeros((2*maxlag+1),1);
    Ntot = length(Y);
    % compute and plot sample time-lagged autocorrelations:
    % Loop for time lags, tau = -maxlag to +maxlag
    for tau = 0:maxlag
        % indices of values to extract from Y(t)
        nn1 = 1:1:Ntot-tau;
        nn2 = 1+tau:1:Ntot;
        % use tau to extract xx and yy
        xx = X(nn1);
        yy = Y(nn2);
        % find Ryy
        Ryy = acov(xx,yy);
        % get rhoyy 
        rhoyy(maxlag+1+tau) = acor(Ryy,xx,yy);
        if tau > 0
            % negative lag
            % Switch indices of xx and yy in autocovariance and autocorrelation
            xx = X(nn2);
            yy = Y(nn1);
            % find Ryy - keep the order the same between xx and yy
            Ryy = acov(xx,yy);
            % get rhoyy
            rhoyy(maxlag+1-tau) = acor(Ryy,xx,yy);
        end
    end
end