 qvap = [0.1 0.5 1 5 10];

cols = 'rgbcmykrgbcmykrgbcmyk'
let = 'ABCDE'
CS = [1e-3 1e-4 1e-2 5e-4 5e-3];

for l = 1:length(let),
    for q = 1:length(qvap)
    clear in out
    fn = sprintf('DR_%s_QVAP_test%02i.mat',let(l),q)
    load(fn)
    figure(1)
    hold on
    plot(CS(l),max(out.concs(:,2)),[cols(l) '*'])
    figure(2)
    hold on
    plot(qvap(q),max(out.concs(:,2)).*1e-6./1800,[cols(l) '*'])
    
end
    
end

figure(1)
set(gca,'xscale','log')
set(gca,'yscale','log')
figure(2)
set(gca,'xscale','log')
set(gca,'yscale','log')