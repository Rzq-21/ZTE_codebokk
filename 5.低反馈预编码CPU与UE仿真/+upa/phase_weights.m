function [c,index]=phase_weights(a,phase_bits)
theta=angle(exp(1j*(angle(a)-angle(a(1)))));
if isinf(phase_bits),c=exp(1j*theta);index=[];else,L=2^phase_bits;index=mod(round(theta*L/(2*pi)),L);c=exp(1j*2*pi*index/L);end
end
