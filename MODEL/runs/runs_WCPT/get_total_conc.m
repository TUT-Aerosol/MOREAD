function[res] = get_total_conc(in,out,tim)

ix = find(out.time >tim,1)

big = out.concs(ix,2);

det = CPCsigmoid(out.wetdiam.*2,2e-9);
det6=CPCsigmoid(out.wetdiam.*2,6e-9);


small = sum(out.concs(ix,in.nucsize:end));
N3 = sum(out.concs(ix,in.nucsize:end).*det(in.nucsize:end));
N6 = sum(out.concs(ix,in.nucsize:end).*det6(in.nucsize:end));



Cvap = out.concs(ix,1);


Ntot = big+small;
res.big = big;
res.small = small;
res.Ntot = Ntot;
res.N3 = N3+big;
res.Cvap = Cvap;
res.N6 = N6+big;










