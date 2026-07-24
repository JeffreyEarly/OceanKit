N0 = 3*2*pi/3600;
L_gm = 1300;
N2 = @(z) N0*N0*exp(2*z/L_gm);
Lz = 4000;

rotationRate = 7.2921E-5;




%%
dx = 10.^(linspace(log10(15),log10(300e3),60));

Lr_scale = (1/(2*rotationRate))*integral(@(z) sqrt(N2(z)),-Lz,0);
Nj_10 = ceil(Lr_scale ./ dx / sind(10));
Nj_30 = ceil(Lr_scale ./ dx / sind(30));
Nj_50 = ceil(Lr_scale ./ dx / sind(50));

figure
plot(dx,Nj_10), hold on
plot(dx,Nj_30)
plot(dx,Nj_50)
xlog, ylog, set(gca, 'XDir', 'reverse')
legend("latitude 10°", "latitude 30°", "latitude 50°",Location="northwest")
xlabel('horizontal resolution (m)')
ylabel('equivalent number of vertical modes')
title("resolution requirements for exponential stratification")
exportgraphics(gca,"resolution-requirements.pdf")