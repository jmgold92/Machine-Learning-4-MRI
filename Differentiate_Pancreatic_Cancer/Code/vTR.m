function [Signal] = vTR(pars,TR)
Mo=pars(1);
T1=pars(2);

Signal=Mo.*(1-exp(-TR./T1));
end

