function F=beams(net,m)
m=mod(m,net.M);n=(0:net.Nap-1).';F=exp(1j*2*pi*n*m(:).'/net.M)/sqrt(net.Nap);
end
