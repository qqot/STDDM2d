function [eb] = Eb_func(x,k0)
%incident electric field
eb=exp(-1i*k0*x);
end