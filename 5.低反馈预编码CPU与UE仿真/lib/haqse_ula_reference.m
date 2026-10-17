function f=haqse_ula_reference(S)
% One-dimensional single-tone CF-mQSE sum-of-ratios; four UE antennas are snapshots.
N=size(S,1);nv=(0:N-1)';q=min(1/sqrt(3*N),.3232);c=(pi*cot(pi*q/N)-pi*N*cot(pi*q))/(N*cos(pi*q*(N-1)/N)^2);
[~,ix]=max(sum(abs(fft(S)).^2,2));p=ix-1;if p>N/2,p=p-N;end
f=p/N;
for it=1:4
 plus=exp(-1j*2*pi*nv*(f+q/N)).'*S;minus=exp(-1j*2*pi*nv*(f-q/N)).'*S;
 ratios=(plus-minus)./(plus+minus);ratios(~isfinite(ratios))=0;delta=real(mean(ratios))/c;f=f+delta/N;
end
f=wrap(f);
end

function x=wrap(x),x=mod(x+.5,1)-.5;end
