% anomalycube.m
% A function to output the anomaly data cube for inputs 3D time series, 
% climatology and date file.
% Kylene Cooley
% 26 Aug 2021

function datA_cube = anomalycube(ts,clim,time,datefile)
   load(datefile,'t_dt')
   
   datA = ts - clim;

   ind = ismember(time,t_dt);
   datA_cube = datA(:,:,ind);
   
end