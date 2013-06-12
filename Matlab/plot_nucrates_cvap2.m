for i = 1:40,
   clear in out nu
    load(sprintf('DR_test%02i',i));
   nu= get_disc_nucrate(in,out);
   Cvap(i) = in.cvap_0;
   madNdt1(i)   = nu.ma_dN2;
   medNdt1(i)   = nu.me_dN2;
   mddNdt1(i)   = nu.md_dN2;
end


plot(Cvap./(1e6.*1.4e7),madNdt1,'*')
hold on
plot(Cvap./(1e6.*1.4e7),medNdt1,'r*')
plot(Cvap./(1e6.*1.4e7),mddNdt1,'k*')

for i = 1:40,
   clear in out nu
    load(sprintf('DR_CO_test%02i',i));
   nu= get_disc_nucrate(in,out);
   Cvap_co(i) = in.cvap_0;
   madNdt1_co(i)   = nu.ma_dN2;
   medNdt1_co(i)   = nu.me_dN2;
   mddNdt1_co(i)   = nu.md_dN2;
end


plot(Cvap_co./(1e6.*1.4e7),madNdt1_co,'o')
hold on
plot(Cvap_co./(1e6.*1.4e7),medNdt1_co,'ro')
plot(Cvap_co./(1e6.*1.4e7),mddNdt1_co,'ko')

for i = 1:40,
   clear in out nu
    load(sprintf('DR_CO_SI_test%02i',i));
   nu= get_disc_nucrate(in,out);
   Cvap_co_si(i) = in.cvap_0;
   madNdt1_co_si(i)   = nu.ma_dN2;
   medNdt1_co_si(i)   = nu.me_dN2;
   mddNdt1_co_si(i)   = nu.md_dN2;
end


plot(Cvap_co_si./(1e6.*1.4e7),madNdt1_co_si,'s')
hold on
plot(Cvap_co_si./(1e6.*1.4e7),medNdt1_co_si,'rs')
plot(Cvap_co_si./(1e6.*1.4e7),mddNdt1_co_si,'ks')

