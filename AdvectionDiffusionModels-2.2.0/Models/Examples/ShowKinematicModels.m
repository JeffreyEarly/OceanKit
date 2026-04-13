shouldSaveImages = false;

eddy = TranslatingGaussian();
figureHandle = figure;
eddy.plotStreamfunction();
hold on
eddy.plotVelocityField(numPoints=20);
if shouldSaveImages
    exportgraphics(figureHandle, 'figures/kinematic_model_eddy.png')
end

jet = MeanderingJet();
figureHandle = figure;
jet.plotStreamfunction();
hold on
jet.plotVelocityField(numPoints=20);
if shouldSaveImages
    exportgraphics(figureHandle, 'figures/kinematic_model_jet.png')
end

%%

cylinder = CylinderFlow();
figureHandle = figure;
cylinder.plotStreamfunction();
hold on
cylinder.plotVelocityField(numPoints=20);
if shouldSaveImages
    exportgraphics(figureHandle, 'figures/kinematic_model_cylinder.png')
end

%%

strain = LinearVelocityField(sigma=1e-6);
figureHandle = figure;
strain.plotStreamfunction();
hold on
strain.plotVelocityField(numPoints=20);
if shouldSaveImages
    exportgraphics(figureHandle, 'figures/kinematic_model_strain.png')
end

%%

strain = LinearVelocityField(zeta=1e-6);
figureHandle = figure;
strain.plotStreamfunction();
hold on
strain.plotVelocityField(numPoints=20);
if shouldSaveImages
    exportgraphics(figureHandle, 'figures/kinematic_model_vorticity.png')
end
