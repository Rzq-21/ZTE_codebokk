function [u,s2,v]=principal(H)
% Only a 4x4 eigendecomposition; avoid a large full SVD.
C=H*H';[U,D]=eig((C+C')/2,'vector');[s2,index]=max(real(D));assert(s2>0&&isfinite(s2));u=U(:,index);if nargout>2,v=H'*u/sqrt(s2);end
end
