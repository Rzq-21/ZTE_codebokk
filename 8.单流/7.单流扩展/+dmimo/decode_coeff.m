function [C,bits]=decode_coeff(msg)
allowed={'version','B','streams','qA','qP','amplitude_indices','phase_indices'};
assert(isempty(setdiff(fieldnames(msg),allowed)),'dmimo:Boundary','Unexpected coefficient message field');
assert(strcmp(msg.version,'COEFF_V1'));B=msg.B;n=msg.streams;qA=msg.qA;qP=msg.qP;
assert(ismember(n,[1 2])&&B==round(B)&&B>=1&&qA>=1&&qP>=1&&qA==round(qA)&&qP==round(qP));
a=msg.amplitude_indices;ph=msg.phase_indices;
assert(isequal(size(a),[2*B,n])&&numel(ph)==2*B*n-n);
assert(all(a(:)==round(a(:))&a(:)>=0&a(:)<2^qA)&&all(ph(:)==round(ph(:))&ph(:)>=0&ph(:)<2^qP));
mask=true(2*B,n);for j=1:n,mask(j,j)=false;end
P=zeros(2*B,n);P(mask)=ph;C=a/(2^qA-1).*exp(1j*2*pi*P/2^qP);
for b=1:B
 ix=2*b-1:2*b;r=norm(C(ix,:),'fro');assert(r>0,'dmimo:ZeroBlock','All-zero amplitude message');C(ix,:)=C(ix,:)/(r*sqrt(B));
end
bits=numel(a)*qA+numel(ph)*qP;
assert(bits==2*B*n*(qA+qP)-n*qP);
end
