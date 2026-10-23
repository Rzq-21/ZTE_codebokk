function m=metrics(E,C,rho,e,n)
G=E*C;T=eye(size(G,1))+rho*(G*G');R=chol((T+T')/2);m.joint=2*sum(log2(real(diag(R))));
se=0;sinr=zeros(1,size(C,2));
for j=1:size(C,2)
 other=G;other(:,j)=[];Q=eye(size(G,1))+rho*(other*other');
 sinr(j)=max(0,real(rho*G(:,j)'*(Q\G(:,j))));se=se+log2(1+sinr(j));
end
m.mmse=se;m.sinr=sinr;m.power=sum(abs(G(:)).^2);m.norm=sum(abs(C(:)).^2);
B=size(C,1)/2;bp=zeros(B,1);for b=1:B,bp(b)=norm(C(2*b-1:2*b,:),'fro')^2;end
m.ap_error=max(abs(bp-1/B));m.capacity=dmimo.capacity(e,n,rho);m.ratio=e(2)/e(1);
assert(abs(m.norm-1)<1e-10&&m.ap_error<1e-10);
assert(m.mmse<=m.joint+1e-8&&m.joint<=m.capacity+1e-7);
if n==1,assert(abs(m.mmse-m.joint)<1e-8);end
end
