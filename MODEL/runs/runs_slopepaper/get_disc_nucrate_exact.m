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


dp1 = find(Dp>=1e-9, 1 );
dp2 = find(Dp>=2e-9, 1 );
dp3 = find(Dp>=3e-9, 1 );
dp4 = find(Dp>=4e-9, 1 );

Ntot = sum(a(:,nus:end),2);
N1 = sum(a(:,dp1:dp2),2); % 1-2
N2 = sum(a(:,dp2:dp3),2); % 2-3
N3 = sum(a(:,dp3:dp4),2); % 3-4
N4 = sum(a(:,dp4:end),2); % 4->


% new way
GR = in.cvap_0./1e6./1.4e7;

v = load('sam20020414.sum');

V = Dlog_to_N(v);

dp = V(1,3:end);
N =  V(67,3:end).*1;

CO = Coags_dR(Dp(nus),dp,N,273); % 1/s
N10 = V(67,3:end).*10;
CO_10 = Coags_dR(Dp(nus),dp,N10,273); %1/s

CO = CO.*3600
CO_10 = CO_10.*3600

CO12 = Coags_dR(Dp(dp1),dp,N,273); % 1/s, overestimate
CO23 = Coags_dR(Dp(dp2),dp,N,273); % 1/s, overestimate
CO34 = Coags_dR(Dp(dp3),dp,N,273); % 1/s, overestimate
%CO45 = Coags_dR(Dp(dp4),dp,N,273); % 1/s, overestimate

Fc1 = N1.*CO12;
Fc2 = N2.*CO23;
Fc3 = N3.*CO34;
%Fc4 = N4.*CO45;

Fg1 = GR.* N1./3600;
Fg2 = GR.* N2./3600;
Fg3 = GR.* N3./3600;
%Fg4 = GR.* N4.*3600;


dtim = a(2:end,1)-diff(a(:,1))./2;
dNtotdt = diff(Ntot)./diff(a(:,1));
dN1dt = diff(N1)./diff(a(:,1));
dN2dt = diff(N2)./diff(a(:,1));
dN3dt = diff(N3)./diff(a(:,1));
dN4dt = diff(N4)./diff(a(:,1));

dN1dt = diff(N1)./diff(a(:,1));
dN2dt = diff(N2)./diff(a(:,1));
dN3dt = diff(N3)./diff(a(:,1));
dN4dt = diff(N4)./diff(a(:,1));

J1 = dN1dt + Fc1(2:end) +Fg1(2:end);
J2 = dN2dt + Fc2(2:end) +Fg2(2:end);
J3 = dN3dt + Fc3(2:end) +Fg3(2:end);
%J4 = dN4dt + Fc4(2:end) +Fg4(2:end);

ixt = find(dNtotdt>1e-2);
ix1 = find(dN1dt>1e-2);
ix2 = find(dN2dt>1e-2);
ix3 = find(dN3dt>1e-2);

figure
plot(a(:,1),Ntot,'m*')
hold on
plot(a(:,1),N1,'k*')
plot(a(:,1),N2,'r*')
plot(a(:,1),N3,'b*')
title('N')

figure
plot(dtim,dNtotdt,'m*')
hold on
plot(dtim,J1,'k*')
plot(dtim,J2,'r*')
plot(dtim,J3,'b*')
%plot(dtim,J4,'g*')
title('J')

figure
plot(dtim,dNtotdt,'m*')
hold on
plot(dtim,dN1dt,'k*')
plot(dtim,dN2dt,'r*')
plot(dtim,dN3dt,'b*')
title('dN')
figure

plot(dtim,Fc1(2:end),'k*')
hold on
plot(dtim,Fc2(2:end),'r*')
plot(dtim,Fc3(2:end),'b*')
title('Coag')
figure
plot(dtim,Fg1(2:end),'k*')
hold on
plot(dtim,Fg2(2:end),'r*')
plot(dtim,Fg3(2:end),'b*')
title('Growth')



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

ou.me_J1 = nanmean(J1(ix1))./sc;
ou.me_J2 = nanmean(J2(ix2))./sc;
ou.me_J3 = nanmean(J3(ix3))./sc;

ou.md_J1 = nanmedian(J1(ix1))./sc;
ou.md_J2 = nanmedian(J2(ix2))./sc;
ou.md_J3 = nanmedian(J3(ix3))./sc;

ou.ma_J1 = nanmax(J1(ix1))./sc;
ou.ma_J2 = nanmax(J2(ix2))./sc;
ou.ma_J3 = nanmax(J3(ix3))./sc;




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
    ou.ma_dN3 = NaN;
end


