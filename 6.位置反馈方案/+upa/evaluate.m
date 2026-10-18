function metric=evaluate(net,sim,Htrue,truth,design,true_position)
% Truth is accessible to the evaluator, not the CPU decoder.
w=design.w;metric.norm_squared=norm(w)^2;assert(abs(metric.norm_squared-1)<1e-10);metric.received_power=norm(Htrue*w)^2;metric.ideal_power=truth.s2;metric.power_gain=metric.received_power/truth.s2;assert(metric.power_gain>=0&&metric.power_gain<=1+1e-8);
metric.SE=log2(1+sim.rho_data*metric.received_power);metric.throughput=sim.bandwidth_Hz/1e6*metric.SE;metric.position_error=NaN;metric.AP1_SSE=NaN;metric.AP1_count=0;metric.pair_SSE=NaN;metric.pair_count=0;
if all(isfinite(design.position)),metric.position_error=norm(design.position-true_position);end
if strcmp(design.method_id,'SVD'),return;end
B=size(design.F,2);power=reshape(abs(w).^2,net.Nap,B);assert(max(abs(sum(power,1)-1/B))<1e-10);a=zeros(B,1);
for slot=1:B,ix=(slot-1)*net.Nap+(1:net.Nap);a(slot)=truth.u'*Htrue(:,ix)*design.F(:,slot);end
% a above is u'*H*f, so its conjugate is f'*H'*u.
a=conj(a);valid=abs(a)>=1e-10*max(abs(a));error=angle(exp(1j*(angle(design.c)-angle(a))));
if valid(1)
 err=angle(exp(1j*(error-error(1))));mask=valid;mask(1)=false;metric.AP1_count=sum(mask);metric.AP1_SSE=sum(err(mask).^2);
end
pair=angle(exp(1j*(error-error.')));mask=triu(valid*valid.',1)>0;metric.pair_count=sum(mask(:));metric.pair_SSE=sum(pair(mask).^2);
end
