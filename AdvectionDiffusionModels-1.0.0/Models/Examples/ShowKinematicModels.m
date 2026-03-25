shouldSaveImages = false;

eddy = TranslatingGaussian();
figure
eddy.plotStreamfunction(), hold on
eddy.plotVelocityField(N=20);

if shouldSaveImages
    exportgraphics('figures/kinematic_model_eddy.png')
end

jet = MeanderingJet();
figure
jet.plotStreamfunction(), hold on
jet.plotVelocityField(N=20);

if shouldSaveImages
    exportgraphics('figures/kinematic_model_jet.png','-dpng')
end

%%

cylinder = CylinderFlow();
figure
cylinder.plotStreamfunction(), hold on
cylinder.plotVelocityField(N=20);

if shouldSaveImages
    exportgraphics('figures/kinematic_model_cylinder.png')
end

%%

strain = LinearVelocityField(1e-6,0,0);
figure
strain.plotStreamfunction(), hold on
strain.plotVelocityField(N=20);

if shouldSaveImages
    exportgraphics('figures/kinematic_model_strain.png')
end

%%
strain = LinearVelocityField(0,0,1e-6);
figure
strain.plotStreamfunction(), hold on
strain.plotVelocityField(N=20);

if shouldSaveImages
    exportgraphics('figures/kinematic_model_vorticity.png')
end