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






figure
plot(Dp(nus:len),conc(i1,nus:len)./sc,'m*')
hold on
plot(Dp(nus:len),conc(i2,nus:len)./sc,'r*')
plot(Dp(nus:len),conc(i3,nus:len)./sc,'g*')
plot(Dp(nus:len),conc(i4,nus:len)./sc,'b*')
%plot(Dp(nus:1000),a(i5,nus:1000)./sc,'k*')

l1 = sprintf('T = %4.1f',a(i1,1));
l2 = sprintf('T = %4.1f',a(i2,1));
l3 = sprintf('T = %4.1f',a(i3,1));
l4 = sprintf('T = %4.1f',a(i4,1));
%l5 = sprintf('T = %4.1f',a(i5,1));

%legend({l1 l2 l3 l4 l5})
legend({l1 l2 l3 l4})
%set(gca,'xscale', 'log')
set(gca,'yscale', 'log')
%set(gca,'ylim',[1 1e4])
xlabel('Particle diameter (s)')
ylabel('Concentration (cm^{-3})')

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




figure
subplot(211)
%plot(a(:,1),Ntot./sc,'ms')
hold on
plot(a(:,1),N1./sc,'k-')
plot(a(:,1),N2./sc,'k--')
plot(a(:,1),N3./sc,'k-.')
%plot(a(:,1),N4./sc,'g*')

% legend('Ntot', 'N1', 'N2', 'N3', 'N4')
legend('N_{1.7-2}', 'N_{1.7-2}', 'N_{3-4}')
 ylabel('Concentration (cm^{-3})')
 xlabel('Simulation time (s)')
 
 
subplot(212)
%plot(dtim,dNtotdt./sc,'ms')
hold on
plot(dtim,dN1dt./sc,'k-')
plot(dtim,dN2dt./sc,'k--')
plot(dtim,dN3dt./sc,'k-.')
%plot(dtim,dN4dt./sc,'g*')

%legend('dNtot', 'dN1', 'dN2', 'dN3', 'dN4')
legend( 'dN_{1.7-2}', 'dN_{2-3}', 'dN_{3-4}')
xlabel('Simulation time (s)')
ylabel('dN/dt (cm^{-3}s^{-1})')


b = a(:,2:end)./sc;
figure
pcolor(a(:,1),Dp(nus:end),real(log10(b(:,nus:end)')))
set(gca,'yscale','log')
shading flat
set(gca,'ylim',[1e-9 5e-9])
caxis([0 4])


[rr rc] = size(b);
for i = 1:rr,
	[Dpn dN(i,:)] = disc_conv(Dp(nus:end),b(i,nus:end));
end

figure
pcolor(a(:,1),Dpn,real(log10(dN')))
set(gca,'yscale','log')
shading flat
set(gca,'ylim',[1e-9 5e-9])
caxis([0 4])

