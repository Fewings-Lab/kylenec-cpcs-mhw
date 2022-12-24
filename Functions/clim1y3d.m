% clim1y3d.m 
% Kylene Cooley
% 17 Mar 2021
% Calculate the climatological annual cycle for a 3d data cube. Outputs are the
% low-pass filtered climatology and unfiltered 2D variability array of the
% unfiltered climatology
% Update 23 Aug 2021: no longer includes a low-pass filter within this part

function [dat0, sigma] = clim1y3d(dat, dn, dt)
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
end