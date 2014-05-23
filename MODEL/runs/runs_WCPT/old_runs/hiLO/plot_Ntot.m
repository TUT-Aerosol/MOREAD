% plot total concentrations




mar = '--'
LO = 1;




sc = (98./29)./2.2e19;



sinks = [1e-6 1e-4 5e-2 1e-1 5e-1 ];
Cvap  = [1e8 3e8 5e8 7e8 1e9 5e9 1e10 5e10 1e11].*1e6;
cv = logspace(min(log10(Cvap)),max(log10(Cvap)),200);
lcv = log10(cv./1e6.*sc./LO);



for s = 1:length(sinks),
    for c = 1:length(Cvap),
            clear in out
            run_name = sprintf('wcpt_CS0%i_Cvap%i_rs2_actLO',s,c);
            
            load(run_name) 
            
            r = get_total_conc(in,out,2.6);
            
            Ntot(s,c) = r.Ntot;
            Cvap_d(s,c) = r.Cvap;
            bigP(s,c) = r.big;
            N3(s,c) = r.N3
            
            
            
%             interpolate smoothly
            
            
            
            
      
    end
end

%figure
%plot(log10(Cvap),log10(Ntot(3,:)),'ro')

% 
% 
% plot(log10(Cvap./1e6.*sc),(Ntot(2,:)./1e6),'ks')
% hold on
% plot(log10(Cvap(1,:)./1e6.*sc),(Ntot(1,:)./1e6),'bp')
% plot(log10(Cvap(1,:)./1e6.*sc),(Ntot(3,:)./1e6),'ro')
% plot(log10(Cvap(1,:)./1e6.*sc),(Ntot(4,:)./1e6),'m^')
% 
% % figure
% plot(log10(Cvap./1e6.*sc),(N3(2,:)./1e6),'ks-')
% hold on
% plot(log10(Cvap(1,:)./1e6.*sc),(N3(1,:)./1e6),'bp-')
% plot(log10(Cvap(1,:)./1e6.*sc),(N3(3,:)./1e6),'ro-')
% plot(log10(Cvap(1,:)./1e6.*sc),(N3(4,:)./1e6),'m^-')
% 
% % 


for i = 1:length(sinks),
    lN3 = log10(N3(i,:)./1e6);
    
    lN3_i(i,:) = interp1(log10(Cvap./1e6.*sc./LO),lN3,lcv,'cubic')
end



hold on
plot(log10(Cvap./1e6.*sc./LO),log10(Ntot(2,:)./1e6),'ks')
hold on
plot(log10(Cvap(1,:)./1e6.*sc./LO),log10(Ntot(1,:)./1e6),'bp')
plot(log10(Cvap(1,:)./1e6.*sc./LO),log10(Ntot(3,:)./1e6),'ro')
plot(log10(Cvap(1,:)./1e6.*sc./LO),log10(Ntot(4,:)./1e6),'m^')
plot(log10(Cvap(1,:)./1e6.*sc./LO),log10(Ntot(5,:)./1e6),'g+')


% figure
plot(log10(Cvap./1e6.*sc./LO),log10(N3(2,:)./1e6),['ks'])
hold on
plot(log10(Cvap(1,:)./1e6.*sc./LO),log10(N3(1,:)./1e6),['bp'])
plot(log10(Cvap(1,:)./1e6.*sc./LO),log10(N3(3,:)./1e6),['ro'])
plot(log10(Cvap(1,:)./1e6.*sc./LO),log10(N3(4,:)./1e6),['m^'])
plot(log10(Cvap(1,:)./1e6.*sc./LO),log10(N3(5,:)./1e6),['g+'])

plot((lcv),(lN3_i(2,:)),['k' mar])
hold on
plot((lcv(1,:)),(lN3_i(1,:)),['b' mar])
plot((lcv(1,:)),(lN3_i(3,:)),['r' mar])
plot((lcv(1,:)),(lN3_i(4,:)),['m' mar])
plot((lcv(1,:)),(lN3_i(5,:)),['g' mar])


