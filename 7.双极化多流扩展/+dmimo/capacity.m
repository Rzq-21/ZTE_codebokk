function rate=capacity(e,n,rho)
gain=rho*e(1:min(n,numel(e)));gain=gain(gain>max(gain)*1e-14);
for k=numel(gain):-1:1
 level=(1+sum(1./gain(1:k)))/k;
 if level>1/gain(k),power=max(0,level-1./gain(1:k));rate=sum(log2(1+power.*gain(1:k)));return;end
end
rate=0;
end
