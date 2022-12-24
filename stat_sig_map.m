% stat_sig_map.m

% show if the mean of the events is statistically significant from zero to
% 95% confidence level

load dSSTdt_cube.mat

dSST_avg = mean(dat_cube,3,'omitnan');

alpha = 0.05;

N = size(dat_cube,3);

sigma = std(dat_cube,1,3,'omitnan');

p = 1-(alpha/2);
q_t = tinv(p,N-1);

delta_mu = (q_t/sqrt(N)).*sigma;

mask_dSST = xor(dSST_avg>delta_mu,dSST_avg<-delta_mu);

dSST_avg_stat_sig = dSST_avg;
dSST_avg_stat_sig(~mask_dSST) = NaN;

save('stat_sig_dSST.mat','dSST_avg_stat_sig','mask_dSST')

% map masked dSST_avg

% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(dSST_info.lat) max(dSST_info.lat)];
lonlim = [min(dSST_info.lon) max(dSST_info.lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = dSST_info.lat*ones(length(dSST_info.lon),1)';
Lon = ones(length(dSST_info.lat),1)*dSST_info.lon';
clim = [min(dSST_avg_stat_sig,[],'all') max(dSST_avg_stat_sig,[],'all')]; % limits for the colorbar
lvls = [min(clim,[],'omitnan'):.001:max(clim,[],'omitnan')];

figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(Lat,Lon,dSST_avg_stat_sig,lvls,'Fill','on'); % Contour of SST'
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^\circ$C/day]";
c.Label.FontSize = 18;
sgtitle("Average Band-pass Filtered $$ \frac{ \partial SST'}{ \partial t} $$",'Interpreter','latex','FontSize',20)

% saveas(gcf,'avg_dSSTdt_stat_sig_v6.png')

%% extracting contour levels to find positve dSST'/dt polygon
% contourcs-stolen code
K = 0;
n0 = 1;
while n0<=size(C_dsst,2)
    K = K + 1;
    n0 = n0 + C_dsst(2,n0) + 1;
end
el = cell(K,1);
Cout = struct('Level',el,'Length',el,'X',el,'Y',el);
n0 = 1;
for k = 1:K
    Cout(k).Level = C_dsst(1,n0);
    idx = (n0+1):(n0+C_dsst(2,n0));
    Cout(k).Length = C_dsst(2,n0);
    Cout(k).X = C_dsst(1,idx);
    Cout(k).Y = C_dsst(2,idx);
    n0 = idx(end) + 1; % next starting index
end

level = -0.5;
n1 = 1;
len1 = Cout(1).Length;
while level<=0.05
    level = Cout(n1).Level;
    len1 = Cout(n1).Length;
    if len1 > 20
        n2 =1;
        x1 = Cout(n1).X;
        
        y1 = Cout(n1).Y;
        
        while n2 < n1
            len2 = Cout(n2).Length;
            if len2 > 20
                x2 = Cout(n2).X;
                y2 = Cout(n2).Y;
            end
            n2 = n2 + 1;
        end
    end
    n1 = n1 +1;
end
    
figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
plotm(y1,x1,'k')
plotm(y2,x2,'k')

%% plot the lines over contour map
load dSSTdt_cube.mat
load stat_sig_dSST.mat
load 0_05_level_lines.mat

% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(dSST_info.lat) max(dSST_info.lat)];
lonlim = [min(dSST_info.lon) max(dSST_info.lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = dSST_info.lat*ones(length(dSST_info.lon),1)';
Lon = ones(length(dSST_info.lat),1)*dSST_info.lon';
clim = [min(dSST_avg_stat_sig,[],'all') max(dSST_avg_stat_sig,[],'all')]; % limits for the colorbar
lvls = [min(clim,[],'omitnan'):.001:max(clim,[],'omitnan')];

figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(Lat,Lon,dSST_avg_stat_sig,lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'k','LineWidth',0.7)
plotm(y2,x2,'k','LineWidth',0.7)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^\circ$C/day]";
c.Label.FontSize = 8;
c.FontSize = 8;
set(gcf,'Units','centimeters','Position',[0 0 9.5 11],'PaperUnits','centimeters','PaperPosition',[0 0 9.5 11])
set(gca,'FontSize',8)
setm(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
t1 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);

% exportgraphics(gcf,'avg_dSSTdt_0_05_line_v6.png','Resolution',1600)

% saveas(gcf,'avg_dSSTdt_0_05_line_v4.png')
%% save the x, y coord of these lines
save('0_05_level_lines.mat','x1','y1','x2','y2')