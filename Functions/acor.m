% acor.m
% Kylene Cooley
% Created 4 May 2020
% autocorrelation

function rho = acor(Ryy,xx,yy)
    rho = Ryy/(nanstd(xx,1)*nanstd(yy,1));
end