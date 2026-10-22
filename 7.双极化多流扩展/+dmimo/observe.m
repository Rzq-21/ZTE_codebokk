function [H,Y]=observe(blocks,Z,channel,active,snr)
H=reshape(blocks.(channel)(:,:,active),8,[]);Z=reshape(Z(:,:,active),8,[]);
Y=H+sqrt(mean(abs(H(:)).^2)/10^(snr/10))*Z;
end
