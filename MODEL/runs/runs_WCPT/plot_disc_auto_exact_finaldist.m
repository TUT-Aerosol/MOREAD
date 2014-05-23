function plot_disc_auto(in,out)
% make_figs


a = [out.time(:) out.concs];

conc = out.concs;



nus = in.nucsize;
len = in.imax;
% Dp = ((1:1000).*(6.*8.89e-29.*3./(4.*pi))).^(1./3); 
Dp = out.wetdiam.*2; 
Dpdry = out.drydiam.*2; 

[ro co]  = size(a);

i1 = 1; 
i2 = floor(ro./5);
i3 = floor(ro./3);
i4 = floor(ro.*0.5);
%i5 = ro;

sc = 1e6;






l1 = sprintf('T = %4.1f',a(i1,1));
l2 = sprintf('T = %4.1f',a(i2,1));
l3 = sprintf('T = %4.1f',a(i3,1));
l4 = sprintf('T = %4.1f',a(i4,1));
%l5 = sprintf('T = %4.1f',a(i5,1));


% get the nuc rates at different sizes...


dp1 = find(Dp>=1e-9, 1 )
dp2 = find(Dp>=2e-9, 1 )
dp3 = find(Dp>=3e-9, 1 )
dp4 = find(Dp>=4e-9, 1 )

Ntot = sum(conc(:,nus:end),2);
N1 = sum(conc(:,dp1:dp2),2);
N2 = sum(conc(:,dp2:dp3),2);
N3 = sum(conc(:,dp3:dp4),2);
N4 = sum(conc(:,dp4:end),2);

dtim = a(2:end,1)-diff(a(:,1))./2;
dNtotdt = diff(Ntot)./diff(a(:,1));
dN1dt = diff(N1)./diff(a(:,1));
dN2dt = diff(N2)./diff(a(:,1));
dN3dt = diff(N3)./diff(a(:,1));
dN4dt = diff(N4)./diff(a(:,1));


b = a(:,2:end)./sc;


[rr rc] = size(b);
for i = 1:rr,
	[Dpn dN(i,:)] = disc_conv(Dp(nus:end),b(i,nus:end));
end

hold on
plot(Dpn,dN(end,:),'.-')

