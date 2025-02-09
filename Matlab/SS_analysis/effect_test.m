clear all
% close all
DP_s = 10e-9;

CS = logspace(-6,-1,8);
GR = logspace(-1,1,10);
Q = logspace(4,6,5)
j=1;

for i = 1:length(CS)
        for k = 1:length(Q)
            f = lehtinen_factor((Q(k)./CS(i))./1.4e7,CS(i),1.5e-9,DP_s);
            J10(i,j,k) = 2e-14.*(Q(k).^2./CS(i).^2).*f.factor;
        end
    
end

cols = 'rgbcmykrgbcmyk';
figure(1)
for k = 1:length(Q)
   hold on
   plot(CS,J10(:,j,k),cols(k));
    l{k} = sprintf('Q = %3.1f',Q(k));

set(gca,'yscale','log')
set(gca,'xscale','log')
legend(l)
end

% get the self-regulating 

N = 0;
J = 0;
CSb = 1e-2;
CSn = 0;
CSt = 1e-3;
col = 'b';
tend = 7200

for t = 1:tend,
    
    % get nucleation rate as a function of CS
    J = interp1(CS,J10(:,j,4),CSn+CSb,'linear','extrap');
    J = max([J 0]);
    Jt(t+1) = J;
    N(t+1) = N(t)+J;
    CSn = CS_general(DP_s,N(t+1),288,1.0);
    CSt(t) = CSn;
%     pause
end

figure(2)
hold on
plot(1:tend+1,N,[col '.-'])

figure(3)
hold on
plot(1:tend+1,Jt,[col '.-'])






return

figure
for i = 1:length(CS)
    hold on
    plot(Q,squeeze(J10(i,j,:)),cols(i));
    l{i} = sprintf('CS = %5.3e',CS(i));

set(gca,'yscale','log')
set(gca,'xscale','log')
legend(l)
end



