function F=beams(net,m)
assert(size(m,1)==2);F=exp(1j*2*pi*(net.ny*double(m(1,:))+net.nz*double(m(2,:)))/net.M)/sqrt(net.Nap);
end
