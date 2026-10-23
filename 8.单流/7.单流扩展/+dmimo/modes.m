function [V,e]=modes(H,n)
Q=H*H';[U,D]=eig((Q+Q')/2,'vector');[e,ix]=sort(max(0,real(D)),'descend');
n=min(n,numel(e));V=H'*U(:,ix(1:n))./sqrt(max(e(1:n).',eps(max(e))));
end
