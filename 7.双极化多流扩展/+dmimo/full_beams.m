function [F,codes]=full_beams(net,Y,B)
% Joint observed leading subspace, common spatial beam for both polarizations.
[V,~]=dmimo.modes(Y,2);codes=zeros(1+strcmp(net.array,'UPA32x32'),B);
for b=1:B
 ix=(b-1)*2*net.Nap+(1:net.Nap);X=[V(ix,:) V(ix+net.Nap,:)];
 if strcmp(net.array,'ULA64')
  score=sum(abs(fft(X,net.M,1)).^2,2);[~,k]=max(score);codes(b)=k-1;
 else
  X=reshape(X,32,32,[]);Z=fft(fft(X,net.M,1),net.M,2);score=sum(abs(Z).^2,3);
  [~,k]=max(score(:));[y,z]=ind2sub([net.M net.M],k);codes(:,b)=[y-1;z-1];
 end
end
if strcmp(net.array,'ULA64'),F=ula.beams(net,codes);else,F=upa.beams(net,codes);end
end
