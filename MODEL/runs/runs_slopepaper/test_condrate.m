% compare condensation rate 

load('TEST_MOREAD.mat');
a = load('COND_RATE.TXT');

s = out.wetdiam; 

for i = 1:length(s),
    cr(i) = koag_kernel(s(1),s(i),1.0,273);
    cr17(i) = koag_kernel(s(1),s(i),1.7,273);
    
    cdp(i) = koag_kernel(s(1).*2,s(i).*2,1.0,273);
    cdp17(i) = koag_kernel(s(1).*2,s(i).*2,1.7,273);

    
end

figure
plot(s,a,'r-')
hold on
plot(s,cr,'b--')
plot(s,cr17,'b-')

plot(s,cdp,'k--')
plot(s,cdp17,'k-')

hold off
legend('model','matlab \rho=1.0','matlab \rho=1.7','assuming radius \rho=1.0','assuming radius \rho=1.7')









