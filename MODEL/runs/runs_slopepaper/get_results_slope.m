% make nucleation rate plots for SLOPE runs



% first, constant J plots

d = dir('rSLOPE_Jconst*.*')

for i = 1:length(d)
    clear in out
    load(d(i).name); % in, out
    
    if in.cvap_0>1e12,
        d(i).name
     
        o = get_disc_nucrate_exact(in,out)
        
        CS(i) = in.condsink_value;
        GR(i) = o.GR_data;
        Jgiven(i) = in.nucrate;
        J2_model(i) = o.J2_avg;
        Jcorr(i) = J2_model(i)./o.factor2;
        Cvap(i) = in.cvap_0;
        
        
        
        
        
        
        

            CS(i) = in.condsink_value; 
            GR(i) = NaN;

    end
end
        
        
        
        
    
    
    
    


