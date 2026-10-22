function codes=codes(net,index)
if strcmp(net.array,'ULA64'),codes=reshape(index-1,1,[]);
else,codes=[reshape(mod(index-1,net.M),1,[]);reshape(floor((index-1)/net.M),1,[])];end
end
