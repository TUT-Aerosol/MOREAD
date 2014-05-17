% run some runs



sinks = [1e-4 1e-3 1e-2];
Cvap  = [5e6 1e7 5e7 1e8].*1e6;
J = [0.01 0.1 1 10].*1e6;

for s = 1:length(sinks),
    for c = 1:length(Cvap),
        for j = 1:length(J),
            run_name = sprintf('rSLOPE_Jconst_CS%i_Cvap%i_J%i',s,c,j);
            
            run_MOREAD('imax',400,'nucrate',J(j),'cvap_0',Cvap(c),'condsink_value',sinks(s),'nuc_mech',3,'runfilename',run_name,'nucsize',12)
        end
    end
end



sinks = [1e-5 5e-5 1e-4 5e-4 1e-3 5e-3 1e-2 5e-2];
Cvap  = [1e6 5e6 1e7 5e7 1e8 5e8].*1e6;
% J = [0.01 0.1 1 10].*1e6;

for s = 1:length(sinks),
    for c = 1:length(Cvap),
        
        run_name1 = sprintf('rSLOPE_Jkine_CS%i_Cvap%i_Khigh',s,c);
        run_name2 = sprintf('rSLOPE_Jkine_CS%i_Cvap%i_Klow',s,c);
        run_MOREAD('imax',400,'cvap_0',Cvap(c),'condsink_value',sinks(s),'nuc_mech',2,'runfilename',run_name1,'nuc_coeff',1e-20,'nucsize',12)
        run_MOREAD('imax',400,'cvap_0',Cvap(c),'condsink_value',sinks(s),'nuc_mech',2,'runfilename',run_name2,'nuc_coeff',1e-21,'nucsize',12)
        
    end
end


