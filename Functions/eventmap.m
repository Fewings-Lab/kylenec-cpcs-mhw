% eventmap.m
% A function for looping through some number of events and plotting a map
% view of the domain in a subplot grid of 4x4 figures.
% Kylene Cooley
% 26 Aug 2021

function eventmap(dat_cube,datefile,lat,lon,grouptitle,clabel)
latlim = [min(lat) max(lat)];
lonlim = [min(lon) max(lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = lat*ones(length(lon),1)';
Lon = ones(length(lat),1)*lon';

clim = [min(dat_cube,[],'all') -min(dat_cube,[],'all')]; % limits for the colorbar
lvls = min(clim,[],'omitnan'):.001:max(clim,[],'omitnan');

load(datefile,'t_dt')
t_dt = rmmissing(t_dt);

for i = 1:ceil(length(t_dt)/16)
    figure()
    for j = 1:16
        k = j+(i-1)*16;
        if k<=length(t_dt)
            subplot(4,4,j)
            h = worldmap(latlim,lonlim); % Map over Chile-Peru System
            setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
            fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds coastlines
            [~,~] = contourm(Lat,Lon,dat_cube(:,:,k),lvls,'Fill','on'); % Contour of SST'
            caxis(clim)
            cmocean('balance',length(lvls),'pivot',0)
            title(datestr(t_dt(k),24),'FontSize',14)
            daspect([1,1.75,1])
            tightmap
        end
        
    end
    if k>length(t_dt)
        sgtitle(grouptitle,'Interpreter','tex','FontSize',20)
        c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
        c.Label.Interpreter = 'latex';
        c.Label.String = clabel;
        c.Label.FontSize = 18;
    else
        sgtitle(grouptitle,'Interpreter','tex','FontSize',20)
        hp4 = get(subplot(4,4,16),'Position');
        c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
        c.Label.Interpreter = 'latex';
        c.Label.String = clabel;
        c.Label.FontSize = 18;
    end
end
end