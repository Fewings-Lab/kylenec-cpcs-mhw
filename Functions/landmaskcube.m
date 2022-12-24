% landmaskcube.m
% function outputs data cube the same size as input data cube but makes
% points over land equal to NaN
% Kylene Cooley
% 26 Aug 2021

function dat_cube = landmaskcube(dat_cube,lat,lon)
em = read_era5_new('other','daily','01/01/1979','12/31/2020',{'land_sea_mask'},[min(lon) max(lon)],[min(lat) max(lat)]);
boo = find(round(em.land_sea_mask));
land_mask = zeros(size(dat_cube,[1 2]));

land_mask(boo) = NaN;
for z = 1:length(squeeze(dat_cube(1,1,:)))
   dat_cube(:,:,z) = dat_cube(:,:,z)+land_mask;
end
end
