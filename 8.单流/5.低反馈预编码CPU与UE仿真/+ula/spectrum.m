function [m,g]=spectrum(net,v)
% g(m+1) = f_m^H v, with deterministic first-index tie breaking.
g=fft(v,net.M)/sqrt(net.Nap);[~,ix]=max(abs(g).^2);m=ix-1;
end
