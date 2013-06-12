function[ou]=plot_disc_auto(in,out)
% make_figs


a = [out.time(:) out.concs];

nus = in.nucsize;
% Dp = ((1:1000).*(6.*8.89e-29.*3./(4.*pi))).^(1./3); 
Dp = out.wetdiam; 
Dpdry = out.drydiam; 

[ro co]  = size(a);

i1 = 1; 
i2 = floor(ro./4);
i3 = floor(ro./2);
i4 = floor(ro.*0.75);
i5 = ro;

sc = 1e6;








% get the nuc rates at different sizes...


dp1 = min(find(Dp>=1e-9));
dp2 = min(find(Dp>=2e-9));
dp3 = min(find(Dp>=3e-9));
dp4 = min(find(Dp>=4e-9));

Ntot = sum(a(:,nus:end),2);
N1 = sum(a(:,dp1:end),2);
N2 = sum(a(:,dp2:end),2);
N3 = sum(a(:,dp3:end),2);
N4 = sum(a(:,dp4:end),2);

dtim = a(2:end,1)-diff(a(:,1))./2;
dNtotdt = diff(Ntot)./diff(a(:,1));
dN1dt = diff(N1)./diff(a(:,1));
dN2dt = diff(N2)./diff(a(:,1));
dN3dt = diff(N3)./diff(a(:,1));
dN4dt = diff(N4)./diff(a(:,1));


ixt = find(dNtotdt>1e-2);
ix1 = find(dN1dt>1e-2);
ix2 = find(dN2dt>1e-2);
ix3 = find(dN3dt>1e-2);

ou.me_dNtot = nanmean(dNtotdt(ixt))./sc;
ou.me_dN1 = nanmean(dN1dt(ix1))./sc;
ou.me_dN2 = nanmean(dN2dt(ix2))./sc;
ou.me_dN3 = nanmean(dN3dt(ix3))./sc;

out.md_dNtot = nanmedian(dNtotdt(ixt))./sc;
ou.md_dN1 = nanmedian(dN1dt(ix1))./sc;
ou.md_dN2 = nanmedian(dN2dt(ix2))./sc;
ou.md_dN3 = nanmedian(dN3dt(ix3))./sc;

ou.ma_dNtot = nanmax(dNtotdt(ixt))./sc;
ou.ma_dN1 = nanmax(dN1dt(ix1))./sc;
ou.ma_dN2 = nanmax(dN2dt(ix2))./sc;
ou.ma_dN3 = nanmax(dN3dt(ix3))./sc;

% scale using Lehtinen et al

ou.gamma1 = (1./(-1.7+1)).*((1e-9./Dp(nus)).^(-1.7)-1);
ou.gamma2 = (1./(-1.7+1)).*((2e-9./Dp(nus)).^(-1.7)-1);
ou.gamma3 = (1./(-1.7+1)).*((3e-9./Dp(nus)).^(-1.7)-1);


v = load('sam20020414.sum');

V = Dlog_to_N(v);

dp = V(1,3:end);
N =  V(67,3:end).*1;

CO = Coags_dR(Dp(nus),dp,N,273); % 1/s
N10 = V(67,3:end).*10;
CO_10 = Coags_dR(Dp(nus),dp,N10,273); %1/s

CO = CO.*3600
CO_10 = CO_10.*3600

ou.factor1 = exp(-ou.gamma1.*Dp(nus).*1e9.*(1.44./(in.cvap_0./1e6./1.4e7)));
ou.factor2 = exp(-ou.gamma2.*Dp(nus).*1e9.*(1.44./(in.cvap_0./1e6./1.4e7)));
ou.factor3 = exp(-ou.gamma3.*Dp(nus).*1e9.*(1.44./(in.cvap_0./1e6./1.4e7)));

ou.factor2_CS10 = exp(-ou.gamma2.*Dp(nus).*1e9.*(CO_10./(in.cvap_0./1e6./1.4e7)));
ou.factor2_CS01 = exp(-ou.gamma2.*Dp(nus).*1e9.*(0.144./(in.cvap_0./1e6./1.4e7)));

ou.kineJ = (in.cvap_0.^2).*(1.0e-20)./sc;
ou.actiJ = (in.cvap_0).*(2.0e-6)./sc;

if isempty(ou.ma_dN1),
    ou.ma_dN1 = NaN;
end
if isempty(ou.ma_dN2),
    ou.ma_dN2 = NaN;
end
if isempty(ou.ma_dN3),
    ou.ma_dN2 = NaN;
end


