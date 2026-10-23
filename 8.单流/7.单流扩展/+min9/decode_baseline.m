function [F,db,ib]=decode_baseline(net,msg,B)
allowed={'version','kind','codes','amplitude_indices','phase_indices'};
assert(isempty(setdiff(fieldnames(msg),allowed)),'min9:Boundary','Unexpected baseline message field');
assert(strcmp(msg.version,'BASE_DUAL_MIN1'));k=msg.codes;assert(all(k(:)==round(k(:))&k(:)>=0&k(:)<net.M));
d=1+strcmp(net.array,'UPA32x32');F=zeros(net.Nap,B);
if strcmp(msg.kind,'TypeI')
 assert(isequal(size(k),[d B])&&isempty(msg.amplitude_indices)&&isempty(msg.phase_indices));F=min9.beams(net,k);db=B*net.index_bits;ib=0;
else
 assert(strcmp(msg.kind,'TypeII')&&isequal(size(k),[d 2 B]));a=msg.amplitude_indices;p=msg.phase_indices;
 assert(numel(a)==B&&numel(p)==B&&all(a==round(a)&a>=0&a<8)&&all(p==round(p)&p>=0&p<4));
 for b=1:B
  f=min9.beams(net,k(:,:,b))*[1;a(b)/7*exp(1j*2*pi*p(b)/4)];F(:,b)=f/norm(f);
 end
 db=2*B*net.index_bits;ib=5*B;
end
assert(max(abs(sum(abs(F).^2,1)-1))<1e-10);
end
