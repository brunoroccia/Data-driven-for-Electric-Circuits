function [etilde, stilde] = phiES(e, s, edata, sdata, p)

Se = size(edata);

eones = e*ones(Se);

sones = s*ones(Se);

q = p/(1-p);

dist = 1/p*abs(eones-edata).^p+1/q*abs(sones-sdata).^q;
%dist = 1/p*abs(eones-edata).^p;

[value, index] = min(dist);

etilde = edata(index);

stilde = sdata(index);

end