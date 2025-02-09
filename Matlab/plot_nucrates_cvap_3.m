
for i = 1:40,
   clear in out nu
    load(sprintf('DR_A_KIN_CO_SI_test%02i',i));
   nu= get_disc_nucrate(in,out);
   Cvap(i) = in.cvap_0;
   madNdt1(i)   = nu.ma_dN2;
   medNdt1(i)   = nu.me_dN2;
   mddNdt1(i)   = nu.md_dN2;
   nufact(i)    = nu.factor2;
   kineJ(i)     = nu.kineJ;
end

vapsca = 1.0
% vapsca = 1.4e7
%plot(Cvap./(1e6.*vapsca),madNdt1,'*')
plot(Cvap./(1e6.*vapsca),medNdt1,'r*')
hold on
plot(Cvap./(1e6.*vapsca),medNdt1./nufact,'b*')
% plot(Cvap./(1e6.*vapsca),kineJ,'k*')

%plot(Cvap./(1e6.*vapsca),mddNdt1,'k*')
hold on




for i = 1:40,
   clear in out nu
    load(sprintf('DR_B_KIN_CO_SI_test%02i',i));
   nu= get_disc_nucrate(in,out);
   Cvap_co(i) = in.cvap_0;
   madNdt1_co(i)   = nu.ma_dN2;
   medNdt1_co(i)   = nu.me_dN2;
   mddNdt1_co(i)   = nu.md_dN2;
   nufact_co(i)    = nu.factor2_CS01;
   kineJ_co(i)     = nu.kineJ;
   
end


%plot(Cvap_co./(1e6.*vapsca),madNdt1_co,'o')
hold on
plot(Cvap_co./(1e6.*vapsca),medNdt1_co,'ro')
plot(Cvap_co./(1e6.*vapsca),medNdt1_co./nufact_co,'bo')
% plot(Cvap_co./(1e6.*vapsca),kineJ_co,'ko')
%plot(Cvap_co./(1e6.*vapsca),mddNdt1_co,'ko')

for i = 1:40,
   clear in out nu 
    load(sprintf('DR_C_KIN_CO_SI_test%02i',i));
   nu= get_disc_nucrate(in,out);
   Cvap_co_si(i) = in.cvap_0;
   madNdt1_co_si(i)   = nu.ma_dN2;
   medNdt1_co_si(i)   = nu.me_dN2;
   mddNdt1_co_si(i)   = nu.md_dN2;
   nufact_co_si(i)    = nu.factor2_CS10;
   kineJ_co_si(i)     = nu.kineJ;   
end


%plot(Cvap_co_si./(1e6.*vapsca),madNdt1_co_si,'s')
hold on
plot(Cvap_co_si./(1e6.*vapsca),medNdt1_co_si,'rs')
plot(Cvap_co_si./(1e6.*vapsca),medNdt1_co_si./nufact_co_si,'bs')
plot(Cvap_co_si./(1e6.*vapsca),kineJ_co_si,'ks')
%plot(Cvap_co_si./(1e6.*vapsca),mddNdt1_co_si,'ks')

l1 = 'Mean dN(Dp>2nm)/dt, CS = 4x10^{-4}';
l11= 'Scaled J, CS = 4x10^{-4}';
l2 = 'Mean dN(Dp>2nm)/dt, CS = 4x10^{-5}';
l21= 'Scaled J, CS = 4x10^{-5}';
l3 = 'Mean dN(Dp>2nm)/dt, CS = 4x10^{-3}';
l31= 'Scaled J, CS = 4x10^{-3}';
l4 = 'Actual J'


% l2 = 'Mean dN(Dp>2nm)/dt';
% l3 = 'Median dN(Dp>2nm)/dt';
legend({l1 l11 l2  l21 l3 l31})

xlabel('Vapour concentration (cm^{-3})')
ylabel('dN/dt (cm^{-3})s^{-1}')
    