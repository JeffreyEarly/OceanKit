z = linspace(-5000, 0, 500);


imAnalytical = InternalModesExponentialStratification(N0=5.2e-3,b=1300,zIn=[-5000 0],zOut=z,latitude=33,nModes=128);
imAnalytical.normalization = Normalization.kConstant;

profile on
for i=1:10
    [F_analytical,G_analytical,h_analytical] = imAnalytical.ModesAtFrequency( 0.2*imAnalytical.N0 );
end

profile viewer
