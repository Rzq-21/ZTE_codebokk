function C=coefficients(E,n,B)
[C,~]=dmimo.modes(E,n);
for b=1:B
 ix=2*b-1:2*b;scale=norm(C(ix,:),'fro');assert(scale>1e-14);
 C(ix,:)=C(ix,:)/(scale*sqrt(B));
end
% Fix one unobservable global phase per stream at public anchor rows 1,2.
for j=1:n,C(:,j)=C(:,j)*exp(-1j*angle(C(j,j)));end
end
