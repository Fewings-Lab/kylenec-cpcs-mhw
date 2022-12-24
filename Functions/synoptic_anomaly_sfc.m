% synoptic_anomaly_sfc.m

function [out, info] = synoptic_anomaly_sfc(vars,latlim,lonlim,tstart,tend,anmly)

for k = 1:length(vars)
    var = vars{k};
    switch nargin
        case 6
            dat = read_era5_new('sfc','daily',tstart,tend,{var},lonlim,latlim,'anomaly',1);
            fields = fieldnames(dat);
            data = eval(join(['dat.',fields{8}],'.')); % outputs data vars retrieved from read_era5_new into a generically-named matrix
        case 5
            dat = read_era5_new('sfc','daily',tstart,tend,{var},lonlim,latlim);
            fields = fieldnames(dat);
            data = eval(join(['dat.',fields{7}],'.')); % outputs data vars retrieved from read_era5_new into a generically-named matrix
        otherwise
            print('Error: not enough or too many inputs')
            return
    end
    
    % filter to time scales of 10 days to 6 months with pl66
    
    % low-pass filter with a 10-day pl66 filter
    dt = dat.time(2)-dat.time(1); % sampling interval in matdate (number of days)  
    Tl = 10; % half-amplitude cutoff period in days to match matdate
    dat_low = NaN(size(data)); % this will be the 10-day low-pass filtered data cube
    for i=1:length(dat.lon) % looping through longitude slices since pl66 only handles 2D matrices
        szlo = size(dat_low(:,i,:));
        filtered = pl66tn(squeeze(data(:,i,:)),dt*24,Tl*24);
        if size(szlo) == size(filtered)
            dat_low(:,i,:) = filtered; % pl66 assumes time in hours
        else
            dat_low(:,i,:) = filtered';
        end
    end


    % high-pass filter with 6-month pl66 filter
    Th = hours(years(0.5));
    dat_lolo = NaN(size(dat_low));
    for j = 1:length(dat.lon)
        filtered = pl66tn(squeeze(dat_low(:,j,:)),dt*24,Th);
        szll = size(dat_lolo(:,j,:));
       if size(szll) == size(filtered)
            dat_lolo(:,j,:) = filtered; % pl66 assumes time in hours
        else
            dat_lolo(:,j,:) = filtered';
       end
    end

    dat_band = dat_low-dat_lolo;
    dat_band(:,:,1:2*round(Th/(dt*24))) = NaN;
    dat_band(:,:,end-2*round(Th/(dt*24)):end) = NaN;
    
    
    str1 = [join(['out.',vars{k}],'.'), '=dat_band'];
    eval(str1);
    switch nargin
        case 6
            dat_new = rmfield(dat,fields{8});
        case 5
            dat_new = rmfield(dat,fields{7});
        otherwise
            print('Error: something wrong with removing old data field from info!')
            return
    end
   
    str2 = [join(['info.',vars{k}],'.'), '= dat_new'];
    eval(str2);
end
return 
