function msg=encode_coeff(C,qA,qP)
B=size(C,1)/2;n=size(C,2);a=abs(C)*sqrt(B);ph=mod(angle(C),2*pi);
msg=struct('version','COEFF_V1','B',B,'streams',n,'qA',qA,'qP',qP,...
 'amplitude_indices',round(min(1,a)*(2^qA-1)),...
 'phase_indices',mod(round(ph*2^qP/(2*pi)),2^qP));
% Store only transmitted phases. Anchors have fixed zero phase.
mask=true(2*B,n);for j=1:n,mask(j,j)=false;end
msg.phase_indices=msg.phase_indices(mask);
end
