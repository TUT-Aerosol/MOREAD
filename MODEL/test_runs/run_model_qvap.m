cvap = [0.1];
qvap = [0.1 0.5 1 5 10];
% qvap = 5;
sinksfiles = {'SINKDIST_05.TXT' 'SINKDIST_5.TXT' }
%sinksfiles = {'SINKDIST_0.TXT'}

setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/'); % for MAC

% A = 1, B = 10? C = 0.1; D = 0.5; E = 5

for s= 1:length(sinksfiles),
    delete('SINKDIST.TXT')
    
    str = sprintf('!cp %s SINKDIST.TXT',char(sinksfiles{s}))
    eval(str)
    runns = [67+s '_'];
    
    for i = 1:length(qvap),
        
        
        
        run_name = sprintf('%sQVAP_test%02i',runns,i)
        
        in.imax = 1000;
        in.temp = 273.15;
        in.rh = 0.75;
        in.pres = 1.0;
        in.nucsize = 80;
        in.nucrate = 1000.*1e6;
        in.pulse_length = 1800;
        in.cond_on = 1;
        in.evap_on = 0;
        in.sink_on = 1;
        in.coag_on = 1;
        in.nuc_mech = 1;
        in.nuc_coeff = 1e-14;
        in.nuc_exp = 3;
        in.nuc_coeff_org = 1e-14;
        in.nuc_exp_org = 2;
        in.cvap_0 = cvap(1).*1e7.*1e6;
        in.qvap_0 = qvap(i).*1e5.*1e6;
        
        make_model_setup(in)
        
        delete('DISC_TEST.TXT')
        delete('DISC_DIAMETERS.TXT')
try        
        !./test.bin
        
        a = load('DISC_TEST.TXT');
        b = load('DISC_DIAMETERS.TXT');
        
        out.time  = a(:,1);
        out.concs = a(:,2:end); % note: all under nucsize is rubbish
        out.drydiam = b(1,:);
        out.wetdiam = b(2,:);
        
        sav = sprintf('DR_%s.mat',run_name);
        save(sav,'in','out')
catch
    fprintf('Run %s failed!',run_name)
end
    
    end
end

%  plot_disc_auto_exact(in,out)