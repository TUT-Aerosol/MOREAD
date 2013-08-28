function[out] = lehtinen_factor(GR,CS,Dp_small,Dp_large)
m = -1.7;

% GR in m/s

GR_si = GR.*1e-9./3600;

CoagS_Dpsmall = CS ./ ((0.71e-9./Dp_small).^m);
gamma = (1./(m+1)).*((Dp_large./Dp_small).^(m+1)-1);
factor = exp(-gamma.*Dp_small.*CoagS_Dpsmall./GR_si);

out.gamma = gamma;
out.factor = factor;
out.CO_small = CoagS_Dpsmall;
out.Dp_small = Dp_small;
out.Dp_large = Dp_large;
out.GR = GR;
out.GR_si = GR_si;
out.CS = CS;

end