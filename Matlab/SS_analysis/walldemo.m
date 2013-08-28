% make wall loss demo figure

load MOREAD_Dp.mat

o1 = SS_rates(MOREAD.Dp,MOREAD.width,1e-3,1e-14,1.7);
o2 = SS_rates(MOREAD.Dp,MOREAD.width,1e-3,0.7e-3,1.7);
o3 = SS_rates_slope(MOREAD.Dp,MOREAD.width,1e-3,1e-14,2);
o4 = SS_rates_slope(MOREAD.Dp,MOREAD.width,1e-3,0.7e-3,2);


figure
plot(o1.Dp,o1.N,'r-')
hold on
plot(o2.Dp,o2.N,'k-')
plot(o3.Dp,o3.N,'b-')
plot(o4.Dp,o4.N,'g-')

figure
plot(o1.Dp,o1.G.*MOREAD.width.*3600.*1e9,'r-')
hold on
plot(o3.Dp,o3.G.*MOREAD.width.*3600.*1e9,'b-')

figure
plot(o2.Dp,o2.G.*MOREAD.width.*3600.*1e9,'k-')
hold on
plot(o4.Dp,o4.G.*MOREAD.width.*3600.*1e9,'g-')
