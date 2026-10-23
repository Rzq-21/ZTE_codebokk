function [message,diagnostic]=encode(position,refs,cfg,total)
% UE encoder: only its continuous estimate, never the true UE location.
validateattributes(position,{'double'},{'size',[3 1],'finite'});
bits=posfb.allocate(cfg,total);count=2.^bits;p=position(1:cfg.dim);
raw=floor((p-cfg.lo)./(cfg.hi-cfg.lo).*count);overflow=any(raw<0|raw>=count);
index=max(0,min(count-1,raw));
message=struct('method','POS','reference_ap_ids',refs(:).','integer_codes',index,'coordinate_bits',bits,'quantizer_id',cfg.id);
diagnostic=struct('overflow',overflow,'unquantized_position',position,'raw_indices',raw);
end
