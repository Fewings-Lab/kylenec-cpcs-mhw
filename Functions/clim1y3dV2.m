% clim1y3dV2.m 
% Kylene Cooley
% 31 Mar 2021
% Calculate the climatological annual cycle for a 3d data cube and apply a
% low-pass filter to smooth with cutoff period Tc in hours. Outputs are the
% low-pass filtered climatology and unfiltered 2D variability array of the
% unfiltered climatology

function [dat0, sigma, sigUP] = clim1y3dV2(dat, dn, dt, Tc)
    % Make a climatological annual cycle for each lat, lon pair
    dv = datevec(dn); % Convert to datevec
    yd = dn - datenum(dv(:,1),1,1) + 1; % Convert to yearday
    % Vector to use for matching times to 6-hourly data
    yhr = 1:(dt/24):(367-(dt/24));
    % Vector to store climatology
    dat0 = NaN(size(dat));

    % loop by yearday w/ time
    for i=1:366*(24/dt)
       times = ismembertol(yd,yhr(i)); % Logical of times in the 42-yr record that match the yearday and time
       mu = mean(dat(:,:,times),3,'omitnan'); % take the mean along time dimension
       k = find(times); % Indices of the nonzero points in times
       for j = 1:length(k) % Loop through these times
           dat0(:,:,k(j)) = mu; % Assign 2D mean array to the time slice
       end
    end 
 
    sigma = std(dat0,0,3,'omitnan'); % std of unfiltered climatology to plot outside function
    
    % standard deviation for map of only upwelling times
    load('tupwell.mat','tupwell')
    sigUP = std(dat0(:,:,tupwell),0,3,'omitnan'); % std of unfiltered climatology to plot outside function
    
    % Low-pass filter climatological annual cycle
    for l=1:length(dat(:,1,1)) % loop along first spatial dimension since pl66 can only filter 2D arrays,
        % we can use either as long as time is the longer dimension
        dat0(l,:,:) = (pl66tn(squeeze(dat0(l,:,:)),dt,Tc))'; % apply filter and assign to 2D first-dimension slice
    end
    
end