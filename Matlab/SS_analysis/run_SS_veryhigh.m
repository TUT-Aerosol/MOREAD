clear all
%gamma = logspace(-5,-2,10);
gamma = [1e-14 1e-3];
CS = [1 10 20 ];
GR_slope = 0.1:0.2:10;

load MOREAD_Dp.mat;

Dp = MOREAD.Dp,
width = MOREAD.width;

close all


cols = ['rgbcmykrgbcmykrgbcmykrgbcmyk']

for c = 1:length(CS)
 for j =1:length(gamma)
for i =1:length(GR_slope)
    
    
    o(i,j,c) = SS_rates(Dp,width,CS(c),gamma(j),GR_slope(i));
    
J(i,j,c) = o(i,j,c).N(84).*o(i,j,c).G(84);    
    
%     figure(1)
%     hold on
%     plot(o(end,j).Dp,o(end,j).N,cols(j))
%     hold off
%     
%     figure(2)
%     hold on
%     plot(o(i,j).N(21).*(o(i,j).G(21)),o(i,j).N(end).*(o(i,j).G(end)),[cols(j) '*'])
%     hold off
%     
%     figure(3)
%     hold on
%     plot(o(i,j).N(21),o(i,j).N(end),[cols(j) '*'])
%     hold off
% % %     
% %     figure(4)
%     hold on
%     plot(GR_slope(i).*2e7,o(i,j).N(end).*o(i,j).G(end),[cols(j) '*'])
%     plot(GR_slope(i).*2e7,o(i,j).N(21).*o(i,j).G(21),[cols(j) 'o'])
%     
%     
    hold off
    fprintf('%i/%i - %i/%i -%i/%i\n',c,length(CS),j,length(gamma),i,length(GR_slope))
end

 end
end

 
for i= 1:length(CS)

   gamma3 = (1./(-1.7+1)).*((3e-9./1.5e-9).^(-1.7)-1);
   factor3 = exp(-gamma3.*1.5.*CS(i).*(3./1.5).^-1.7./(GR_slope./3600));
    
    
    figure(1)
    hold on
    plot((GR_slope),((GR_slope(:).^2).*J(:,1,i)./factor3(:)),cols(i));
    plot((GR_slope),((GR_slope(:).^2).*J(:,1,i)),[cols(i) ':']);
    
%     plot(GR_slope,(GR_slope.^2).*factor3,[cols(i) ':'])
    plot((GR_slope),((GR_slope(:).^2).*J(:,2,i)./factor3(:)),[cols(i) '--']);
    plot((GR_slope),(GR_slope.^2),'k.-')
    plot((GR_slope),(GR_slope.^3./10),'k.-')
    
    hold off
end

for i= 1:length(CS)

   gamma3 = (1./(-1.7+1)).*((3e-9./1.5e-9).^(-1.7)-1);
   factor3 = exp(-gamma3.*1.5.*CS(i).*(3./1.5).^-1.7./(GR_slope./3600));
    
    
    figure(2)
    hold on
    plot((GR_slope(2:end)),diff(log10((GR_slope(:).^2).*J(:,1,i)./factor3(:)))./diff(log10(GR_slope(:))),cols(i));
%     plot(GR_slope,(GR_slope.^2).*factor3,[cols(i) ':'])
    plot((GR_slope(2:end)),diff(log10((GR_slope(:).^2).*J(:,2,i)./factor3(:)))./diff(log10(GR_slope(:))),[cols(i) '--']);
    plot((GR_slope(2:end)),diff(log10(GR_slope.^2))./diff(log10(GR_slope)),'k.-')
    plot((GR_slope(2:end)),diff(log10(GR_slope.^3./10))./diff(log10(GR_slope)),'k.-')
    
    hold off
end





