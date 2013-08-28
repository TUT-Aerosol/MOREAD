function[out] = SS_rates(Dp,width,CS,gamma,GR_slope)


wi = width;

GR = 3; % nm/h
GR = 2.*Dp.*1e9-2.8;

GR = GR_slope.*Dp.*1e9+0.2-GR_slope.*1.5;

%GR = GR_slope;


Co = (Dp./0.71e-9).^-1.7.*CS;
%Co = zeros(size(Dp));
G = GR.*1e-9./3600 .*1./(wi); 

%figure
%plot(Dp,Co,'b-')
%hold on
%plot(Dp,G,'r-')


% scale gamma with diffusion coefficient
diff_part = Diff_particle(Dp,293);
x = gamma./diff_part(1);

gamma = x.*diff_part;





N(1) = 1./(Co(1)+G(1)+gamma(1));

for i = 2:length(Dp),
  N(i) = N(i-1).*G(i-1)./(Co(i)+G(i)+gamma(i));    
end

%figure
%plot(Dp,N,'.')
out.Dp = Dp; 
out.Co = Co; 
out.G = G; 
out.N = N; 
out.gamma = gamma; 
out.CS = CS;

