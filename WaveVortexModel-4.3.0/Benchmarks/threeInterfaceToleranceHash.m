function value = threeInterfaceToleranceHash(tolerances)
% Return the portable runtime's 64-bit FNV tolerance-array hash.
low = uint32(hex2dec('739D0383'));
high = uint32(hex2dec('14650FB0'));
mask = uint64(2^32-1);
for tolerance = reshape(double(tolerances),1,[])
    bits = typecast(tolerance,'uint64');
    bits = bitand(bits+uint64(hex2dec('80000')),bitcmp(uint64(hex2dec('FFFFF'))));
    low = bitxor(low,uint32(bitand(bits,mask)));
    high = bitxor(high,uint32(bitshift(bits,-32)));
    lowProduct = double(low)*435;
    carry = floor(lowProduct/2^32);
    shiftedLow = double(bitand(low,uint32(hex2dec('00FFFFFF'))))*256;
    high = uint32(mod(double(high)*435+carry+shiftedLow,2^32));
    low = uint32(mod(lowProduct,2^32));
end
value = string(bitor(bitshift(uint64(high),32),uint64(low)));
end
