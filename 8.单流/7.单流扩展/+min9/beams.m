function F=beams(net,codes)
if strcmp(net.array,'ULA64'),F=ula.beams(net,codes);else,F=upa.beams(net,codes);end
end
