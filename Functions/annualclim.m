% annualclim.m
% Kylene Cooley 
% 9 March 2021
% Returns annual climatology that is the same length of input time series
% and is filtered with Tc low-pass filter

function dat0 = annualclim(dat, dn, dt, Tc)
    % Convert datenums to date vector [yy mm dd hh mm ss]
    dv = datevec(dn); 

    % Calculate yearday by subtracting the datenum of Jan 1 of that year 
    % (and add 1 so the first day of the year is yearday 1 not 0)
    yd = dn - datenum(dv(:,1),1,1) + 1; 

    % Vectors to use for matching times and storing climatology
    yhr = 1:(dt/24):(367-(dt/24));
    foo = NaN(length(yd),1);

    % loop by values
    for i=1:366*(24/dt)
        values = ismembertol(yd,yhr(i));
        mu = mean(dat(values),'omitnan');
        foo(values) = mu;
    end 
    % Stitch 5 instances of foo together to filter
    foo2 = cat(1,foo,foo,foo,foo,foo);

    % Filter the working variable
    foo3 = pl66tn(foo2,1,Tc); % apply smoothing filter with cutoff period of Tc hours

    % take only the middle portion of filtered climatology
    dat0 = foo3(2*length(yd)+1:3*length(yd)); 

end
