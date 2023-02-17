% s06_4_R_maps.m

Qnet = load('qnet_ev_data.mat');
dSST = load('dSSTdt_cube.mat');
load event_dates.mat

% nix = find(~ismember(t_dt,dSST.dSST_info.time));

R_data = dSST.dat_cube - Qnet.dat_cube_dT;

%% Plot filled contours of dSST'/dt on world map for each event
% Coastline data set and coordinate limits around Chile-Peru System:

latlim = [min(Qnet.info.sshf.lat) max(Qnet.info.sshf.lat)];
lonlim = [min(Qnet.info.sshf.lon) max(Qnet.info.sshf.lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = Qnet.info.sshf.lat*ones(length(Qnet.info.sshf.lon),1)';
Lon = ones(length(Qnet.info.sshf.lat),1)*Qnet.info.sshf.lon';
clim = [min(R_data,[],'all') max(R_data,[],'all')]; % limits for the colorbar


load event_dates.mat

t_dt = rmmissing(t_dt);
% figure()
% subplot(4,4,16)
% hp4 = get(subplot(4,4,16),'Position');
% close

% parpool(20)
for i = 1:ceil(length(t_dt)/16)
   figure()
   for j = 1:16
       k = j+(i-1)*16;
      if k<=length(t_dt)
            subplot(4,4,j)
            h = worldmap(latlim,lonlim); % Map over Chile-Peru System
            setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
            plotm(coastlat,coastlon) % Adds coastlines
            
            [C,~] = contourm(Lat,Lon,R_data(:,:,k),100,'Fill','on'); % Contour of SST'
            caxis(clim)
            cmocean('balance','pivot',0)
            title(datestr(t_dt(k),24),'FontSize',14)
            daspect([1,1.75,1])
            tightmap
            % patchm([-36.5 -35.5 -35.5 -36.5],[-74.25 -74.25 -73.25 -73.25],'k','FaceAlpha',0.5,'LineStyle','--')
      end
      
   end
   if k>length(t_dt)
    sgtitle("Daily Band-pass Filtered Temperature Change from Q'_{net}",'Interpreter','tex','FontSize',20)
    c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
    c.Label.Interpreter = 'latex';
    c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ from $Q'_{net}$ [$^{\circ}$C/day]";
    c.Label.FontSize = 18;
   else
      sgtitle("Daily Band-pass Filtered Q'_{net}",'Interpreter','tex','FontSize',20)
      hp4 = get(subplot(4,4,16),'Position');
      c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
      c.Label.Interpreter = 'latex';
      c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^{\circ}$C/day]";
      c.Label.FontSize = 18;
   end
end

%% map of conditional average of events
load stat_sig_dSST.mat
load 0_05_level_lines.mat

R_avg = mean(R_data,3,'omitnan');

alpha = 0.05;

N = size(R_data,3);

sigma = std(R_data,1,3,'omitnan');

p = 1-(alpha/2);
q_t = tinv(p,N-1);

delta_mu = (q_t/sqrt(N)).*sigma;

% mask of statistically significant means
sigmask = xor(R_avg>delta_mu,R_avg<-delta_mu);

R_avg_stat_sig = R_avg;
R_avg_stat_sig(~sigmask) = NaN;

% mapping the average time rate of change of the anomaly
% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(Qnet.info.sshf.lat) max(Qnet.info.sshf.lat)];
lonlim = [min(Qnet.info.sshf.lon) max(Qnet.info.sshf.lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = Qnet.info.sshf.lat*ones(length(Qnet.info.sshf.lon),1)';
Lon = ones(length(Qnet.info.sshf.lat),1)*Qnet.info.sshf.lon';
clim = [min(dSST_avg_stat_sig,[],'all') max(dSST_avg_stat_sig,[],'all')]; % limits for the colorbar
lvls = [min(clim,[],'omitnan'):.001:max(clim,[],'omitnan')];

figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,R_avg_stat_sig,lvls,'Fill','on'); % Contour of SST'
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "$^\circ\textsf{C day}^{-1}$";
% sgtitle("Average band-pass filtered temperature change from residual",'Interpreter','tex','FontSize',20)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
plotm(y1,x1,'k','LineWidth',0.7)
plotm(y2,x2,'k','LineWidth',0.7)
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

exportgraphics(gcf,'avg_R_stat_sig_v5.png','Resolution',1600)



%% find and save outline of residual

% contourcs-stolen code
K = 0;
n0 = 1;
while n0<=size(C,2)
    K = K + 1;
    n0 = n0 + C(2,n0) + 1;
end
el = cell(K,1);
Cout = struct('Level',el,'Length',el,'X',el,'Y',el);
n0 = 1;
for k = 1:K
    Cout(k).Level = C(1,n0);
    idx = (n0+1):(n0+C(2,n0));
    Cout(k).Length = C(2,n0);
    Cout(k).X = C(1,idx);
    Cout(k).Y = C(2,idx);
    n0 = idx(end) + 1; % next starting index
end

level = -0.5;
n1 = 1;
len1 = Cout(1).Length;
while level<=0.035
    level = Cout(n1).Level;
    len1 = Cout(n1).Length;
    if len1 > 20
        n2 =1;
        p1 = Cout(n1).X;
        
        q1 = Cout(n1).Y;
        
        while n2 < n1
            len2 = Cout(n2).Length;
            if len2 > 20
                p2 = Cout(n2).X;
                q2 = Cout(n2).Y;
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
plotm(q1,p1,'k')
plotm(q2,p2,'k')

%% save residual lines
save('0_035_residual_lines.mat','p1','q1','p2','q2')
