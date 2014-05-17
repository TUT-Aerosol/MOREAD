sinks = [1e-5 5e-5 1e-4 5e-4 1e-3 5e-3 1e-2 5e-2];
Cvap  = [1e6 5e6 1e7 5e7 1e8 5e8].*1e6;
% J = [0.01 0.1 1 10].*1e6;

for s = 1:length(sinks),
    for c = 1:length(Cvap),
        
        run_name1 = sprintf('rSLOPE_Jact_CS%i_Cvap%i_Ahigh',s,c);
        run_name2 = sprintf('rSLOPE_Jact_CS%i_Cvap%i_Alow',s,c);
        run_MOREAD('imax',400,'cvap_0',Cvap(c),'condsink_value',sinks(s),'nuc_mech',1,'runfilename',run_name1,'nuc_coeff',6e-7,'nucsize',12)
        run_MOREAD('imax',400,'cvap_0',Cvap(c),'condsink_value',sinks(s),'nuc_mech',1,'runfilename',run_name2,'nuc_coeff',6e-8,'nucsize',12)
        
    end
end












