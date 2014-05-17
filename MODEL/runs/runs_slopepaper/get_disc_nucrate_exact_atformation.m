function[ou]=plot_disc_auto(in,out,varargin)
% make_figs


if nargin>2,
    figs = varargin{1};
else
    figs = 0;
end


a = [out.time(:) out.concs];

nus = in.nucsize;
Dp = out.wetdiam.*2; 
Dpdry = out.drydiam.*2; 

[ro co]  = size(a);

i1 = 1; 
i2 = floor(ro./4);
i3 = floor(ro./2);
i4 = floor(ro.*0.75);
i5 = ro;

sc = 1e6;


% get the nuc rates at different sizes...


dp1 = find(Dp>=1.7e-9, 1 );
dp2 = find(Dp>=2e-9, 1 );
dp3 = find(Dp>=3e-9, 1 );
dp4 = find(Dp>=4e-9, 1 );

Ntot = sum(a(:,nus:end),2);
if Dp(nus)>1e-9,
    dp1 = nus;
end

N1 = sum(a(:,dp1:dp2),2); % 1-2
N2 = sum(a(:,dp2:dp3),2); % 2-3
N3 = sum(a(:,dp3:dp4),2); % 3-4
N4 = sum(a(:,dp4:end),2); % 4->

dtim = a(2:end,1)-diff(a(:,1))./2;
dNtotdt = diff(Ntot)./diff(a(:,1));

dN1dt = diff(N1)./diff(a(:,1));
dN2dt = diff(N2)./diff(a(:,1));
dN3dt = diff(N3)./diff(a(:,1));
dN4dt = diff(N4)./diff(a(:,1));

% getting the midpoint for the growth rate calculations
% algorithm: 
% 1. find max value
% 2. find first value when N reaches 10 % of max
% 3. find first value after this when N reaches 90% of max
% 4. determine the midpoint value of these two values
% 5. find the first value that goes over this

maN1 = max(N1);
maN2 = max(N2);
maN3 = max(N3);
maN4 = max(N4);

v10_N1 = find(N1>(maN1.*0.10),1);
v10_N2 = find(N2>(maN2.*0.10),1);
v10_N3 = find(N3>(maN3.*0.10),1);


v90_N1 = find(N1>(maN1.*0.90),1);
v90_N2 = find(N2>(maN2.*0.90),1);
v90_N3 = find(N3>(maN3.*0.90),1);

mid_N1 = mean([N1(v90_N1) N1(v10_N1)]);
mid_N2 = mean([N2(v90_N2) N2(v10_N2)]);
mid_N3 = mean([N3(v90_N3) N3(v10_N3)]);

mid_ix1 = find(N1>(mid_N1),1);
mid_ix2 = find(N2>(mid_N2),1);
mid_ix3 = find(N3>(mid_N3),1);


if figs

figure
plot(a(:,1),Ntot,'m*')
hold on
plot(a(:,1),N1,'k*')
plot(a(:,1),N2,'r*')
plot(a(:,1),N3,'b*')

% plot the midpoints
plot(a(mid_ix1,1),N1(mid_ix1),'rs')
plot(a(mid_ix2,1),N2(mid_ix2),'rs')
plot(a(mid_ix3,1),N3(mid_ix3),'rs')
title('N')

end




% Growth rates are here
GR_vapor = in.cvap_0./1e6./1.4e7;
delta_t_data  = a(mid_ix3,1)-a(mid_ix2,1); % the time difference in seconds
if isempty(delta_t_data)
    delta_t_data = NaN;
    ou.GR_data = NaN;
    ou.gamma1 = (1./(-1.7+1)).*((2e-9./Dp(nus)).^(-1.7)-1);
    ou.factor1 = NaN;
    ou.J1_avg = NaN;
    return
end

GR_data = 3600./delta_t_data;



% get the condensation and coagulation sinks

CS = CS_general(in.sinkdist(1,:),in.sinkdist(2,:),in.temp,1.0);
CO12 = Coags_dR(Dp(dp1),in.sinkdist(1,:),in.sinkdist(2,:),in.temp);

% and the rate calculations are done here
Fc1 = N1.*CO12;
Fg1 = N1./(0.3.*delta_t_data); % scaled with bin width
J1 = dN1dt + Fc1(2:end) +Fg1(2:end);




% get the time period for which the J2 is above zero

ix_start1 = find(J1>0,1);
numvec = 1:length(J1);
ix_end1 = find((J1<0 & numvec'>5),1); % not at start
if isempty(ix_end1)
    ix_end1 = length(J1);
end



J1_mean = nanmedian(J1(ix_start1:ix_end1));








if figs

figure
plot(a(:,1),N1/1000,'k*')
hold on
plot(dtim,J1,'r-')
plot(dtim,Fc1(2:end),'b-.')
plot(dtim,Fg1(2:end),'k--')
plot(dtim,dN1dt,'g:')
% mark the spot
plot(dtim(ix_start1:ix_end1), J1(ix_start1:ix_end1),'r*')
line([dtim(ix_start1) dtim(ix_end1)],[J1_mean J1_mean])

legend('N1','J1','FCoag','Fgrowth','dNdt')

end
% pause
% 
% 
% % title('J')
% % 
% figure
% plot(dtim,dNtotdt,'m*')
% hold on
% plot(dtim,dN1dt,'k*')
% plot(dtim,dN2dt,'r*')
% plot(dtim,dN3dt,'b*')
% title('dN')
% 
% figure
% plot(dtim,Fc2(2:end),'r*')
% title('Coag')
% 
% figure
% plot(dtim,Fg2(2:end),'r*')
% title('Growth')

ou.GR_data = GR_data;
ou.gamma1 = (1./(-1.7+1)).*((1.7e-9./Dp(nus)).^(-1.7)-1);
ou.factor1 = exp(-ou.gamma1.*Dp(nus).*1e9.*(delta_t_data.*CO12));
ou.J1_avg = J1_mean;

if isempty(J1_mean)
    keyboard
end

    






return

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
% to do later


ou.gamma2 = (1./(-1.7+1)).*((1.7e-9./Dp(nus)).^(-1.7)-1);

ou.factor2 = exp(-ou.gamma2.*Dp(nus).*1e9.*(1./(delta_t_data.*CO23)));




