function[Dp, dN]=disc_conv(Dp_in,N_in)

cs = cumsum(N_in);

newx = logspace(-9,log10(15e-9),100);
Dp = newx(2:end)+(diff(newx)./2);

Ni = interp1(Dp_in,cs,newx,'linear');
dN = diff(Ni);
