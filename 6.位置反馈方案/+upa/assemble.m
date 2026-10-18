function w=assemble(F,c)
B=size(F,2);assert(numel(c)==B);w=reshape(F.*reshape(c,1,B),[],1)/sqrt(B);assert(abs(norm(w)^2-1)<1e-10);
end
