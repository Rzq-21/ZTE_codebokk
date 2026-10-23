function design=decode(net,message,active,cfg)
% CPU information boundary: public geometry/config + integer coordinate message.
allowed={'method','reference_ap_ids','integer_codes','coordinate_bits','quantizer_id'};
assert(isempty(setdiff(fieldnames(message),allowed)),'posfb:InformationBoundary','Unexpected message field');
assert(strcmp(message.method,'POS')&&strcmp(message.quantizer_id,cfg.id));
k=message.integer_codes;bits=message.coordinate_bits;refs=message.reference_ap_ids;
assert(isequal(size(k),[cfg.dim 1])&&isequal(size(bits),[cfg.dim 1]));
assert(all(k==round(k)&k>=0&k<2.^bits)&&all(bits==round(bits)&bits>=4));
assert(isequal(bits,posfb.allocate(cfg,sum(bits)))&&all(ismember(refs,active))&&active(1)==1);
p=zeros(3,1);p(1:cfg.dim)=cfg.lo+(k+.5).*(cfg.hi-cfg.lo)./2.^bits;
design=posfb.prediction(net,p,refs,active);design.direction_message=message;design.direction_bits=sum(bits);design.feedback_bits=sum(bits);
end
