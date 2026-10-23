function [m,g]=spectrum(net,v)
% g(my+1,mz+1)=f(my,mz)'*v; my is the fast index.
g=fft2(reshape(v,net.Nside,net.Nside),net.M,net.M)/sqrt(net.Nap);[~,index]=max(abs(g(:)).^2);[my,mz]=ind2sub([net.M net.M],index);m=[my-1;mz-1];
end
