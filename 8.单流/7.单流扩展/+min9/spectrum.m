function g=spectrum(net,X)
if strcmp(net.array,'ULA64')
 g=fft(X,net.M,1)/sqrt(net.Nap);
else
 X=reshape(X,32,32,[]);g=fft(fft(X,net.M,1),net.M,2)/sqrt(net.Nap);g=reshape(g,net.M^2,[]);
end
end
