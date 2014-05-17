% make nucleation rate plots for SLOPE runs



% first, constant J plots

d = dir('rSLOPE_Jkine*.*')

for i = 1:length(d)
    clear in out
    load(d(i).name); % in, out
    
    
        d(i).name
     
        o = get_disc_nucrate_exact(in,out)
        
        close all
        CS(i) = in.condsink_value;
        GR(i) = o.GR_data;
        Jgiven(i) = in.nucrate;
        J2_model(i) = o.J2_avg;
        Jcorr(i) = J2_model(i)./o.factor2;
        Cvap(i) = in.cvap_0;
        nuc_coeff(i) = in.nuc_coeff;
        
        
        
        
        
        

            CS(i) = in.condsink_value; 
            GR(i) = NaN;

    
end
        
        

% make a figure

six(1,:) = CS==5e-5;
six(2,:) = CS==1e-3;
six(3,:) = CS==5e-2;

hi = nuc_coeff==1e-20;
lo = nuc_coeff==1e-21;


cols = ['rgbcmyk'];

for i = 1:3
ihi = find(hi & six(i,:));
ilo = find(lo & six(i,:));

mat = [Cvap(ihi)' Jcorr(ihi)'];
mat = sortrows(mat,1);
hCvap = mat(:,1);
hJcorr = mat(:,2);

mat2 = [Cvap(ilo)' Jcorr(ilo)'];
mat2 = sortrows(mat2,1);
lCvap = mat2(:,1);
lJcorr = mat2(:,2);

plot(hCvap,hJcorr,[cols(i) 'o-'])
hold on
plot(lCvap./sqrt(10),lJcorr,[cols(i) 's-'])
end
set(gca,'xscale','log')
set(gca,'yscale','log')

cv = logspace(11,15,100);
j = 1e-20.*cv.^2;
plot(cv,j,'k-')
        

        
        
    
    
    
    


