for i = 1:40,
   clear in out nu
    load(sprintf('DR_test%02i',i));
   nu= get_disc_nucrate(in,out);
   Cvap(i) = in.cvap_0;
   madNdt1(i)   = nu.ma_dN1;
   medNdt1(i)   = nu.me_dN1;
   mddNdt1(i)   = nu.md_dN1;
end


plot(Cvap./(1e6.*1.4e7),madNdt1,'*')
hold on
plot(Cvap./(1e6.*1.4e7),medNdt1,'r*')
plot(Cvap./(1e6.*1.4e7),mddNdt1,'k*')

