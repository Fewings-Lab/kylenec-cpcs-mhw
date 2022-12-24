% Vector averaging practice

x = randn(4);
y = randn(4);
lat = [[15 15 15 15];[10 10 10 10];[5 5 5 5];[0 0 0 0]];
lon = [0 5 10 15; 0 5 10 15; 0 5 10 15; 0 5 10 15];

avgx = mean(x,'all');
avgy = mean(y,'all');
avglat = mean(lat,'all');
avglon = mean(lon,'all');

z = x +1i.*y;
mag = abs(z);
dir = angle(z);
avgz = mean(z,'all');

figure
load coastlines
worldmap([min(lat(:,1))-1 max(lat(:,1))+1],[min(lon(1,:))-1 max(lon(1,:))+1])
h = quivermc(lat,lon,x,y);
quivermc([0 0;avglat avglat],[0 avglon; 0 avglon],[0 0; 0 avgx],[0 0; 0 avgy],'color','b')
quivermc([0 0; avglat avglat],[0 avglon; 0 avglon],[0 0; 0 real(avgz)],[0 0; 0 imag(avgz)],'color','r')