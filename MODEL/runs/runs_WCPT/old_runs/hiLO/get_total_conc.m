function[res] = get_total_conc(in,out,tim)

ix = find(out.time >tim,1)

big = out.concs(ix,2);

det = CPCsigmoid(out.wetdiam.*2,2e-9);
Dp3 = find((out.wetdiam.*2)>2e-9,1)


small = sum(out.concs(ix,in.nucsize:end));
N3 = sum(out.concs(ix,Dp3:end).*det(Dp3:end));
Cvap = out.concs(ix,1);


Ntot = big+small;
res.big = big;
res.small = small;
res.Ntot = Ntot;
res.N3 = N3+big;
res.Cvap = Cvap;










