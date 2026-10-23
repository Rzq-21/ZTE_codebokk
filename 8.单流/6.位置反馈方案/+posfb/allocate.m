function bits=allocate(cfg,total)
% Exhaustive integer allocation; lexicographic tie break, same for CPU/UE/R/B.
persistent memo
if isempty(memo),memo=containers.Map('KeyType','char','ValueType','any');end
key=sprintf('%s_%d_%s',cfg.array,total,sprintf('%.17g,',[cfg.sensitivity_weights;cfg.lo;cfg.hi]));
if isKey(memo,key),bits=memo(key);return;end
best=Inf;bits=[];d=cfg.dim;minimum=cfg.minimum_bits_per_axis;
for bx=minimum:total-minimum*(d-1)
 if d==2,candidates=[bx;total-bx];else
  by=minimum:total-bx-minimum;candidates=[repmat(bx,1,numel(by));by;total-bx-by];
 end
 for j=1:size(candidates,2)
  b=candidates(:,j);cost=sum(cfg.sensitivity_weights.*((cfg.hi-cfg.lo)./2.^b).^2);
  if cost<best,best=cost;bits=b;end
 end
end
assert(numel(bits)==d&&sum(bits)==total);
memo(key)=bits;
end
