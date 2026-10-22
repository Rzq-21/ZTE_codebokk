function E=reduce(net,H,F,dual)
if nargin<4,dual=true;end
B=size(F,2);Nr=size(H,1);P=1+dual;E=zeros(Nr,P*B);
for b=1:B
 for pol=1:P
  ix=(b-1)*P*net.Nap+(pol-1)*net.Nap+(1:net.Nap);
  E(:,P*(b-1)+pol)=H(:,ix)*F(:,b);
 end
end
end
