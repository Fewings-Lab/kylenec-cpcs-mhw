% pl66tn3D.m
% a function to do the pl66tn low-pass filter when your array has three
% dimensions to smooth in the time dim with cutoff period Tc in hours.
% Kylene Cooley
% 23 Aug 2021


function datLP = pl66tn3D(dat,dt,Tc)
    % Vector to store low-pass filtered signal
    datLP = NaN(size(dat));
    
    % Low-pass filter climatological annual cycle
    for l=1:length(dat(:,1,1)) % loop along first spatial dimension since pl66 can only filter 2D arrays,
        % we can use either as long as time is the longer dimension
        datLP(l,:,:) = (pl66tn(squeeze(dat(l,:,:)),dt,Tc))'; % apply filter and assign to 2D first-dimension slice
    end
end