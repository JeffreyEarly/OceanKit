% basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";
% basedir = "/Users/Shared/CimRuns_June2025/output/";
% wvd = WVDiagnostics(basedir + replace(getRunParameters(18),"256","512") + ".nc");
% wvt = wvd.wvt;
% timeIndices = 2751:3001;
% timeIndices = 51:251;

%%
% timeIndices = 51:251;
% wvd.createReservoirGroup();

%%
% [transferFlux, forcingFlux, ddt, energy] = wvd.fluxesForReservoirGroup(timeIndices=51:251);

%%
custom_names = configureDictionary("string","cell");
custom_names{"quadratic_bottom_friction"} = "bottom friction";
custom_names{"vertical_diffusivity"} = "diffusivity";
custom_names{"adaptive_damping"} = "damping";
custom_names{"inertial_forcing"} = "NIO";
custom_names{"M2_tidal_forcing"} = "M2 tide";
custom_names{"geostrophic_mean_flow"} = "mean flow";
custom_names{"te_gmda"} = "geostrophic";
custom_names{"damped_geostrophic"} = ["damped", "geostrophic"];
custom_names{"damped_wave"} = ["damped", "wave"];
timeIndices=51:251;

col = configureDictionary("string","cell");
col{"source"} = [191 191 250]/255;
col{"damped_geostrophic"} = [205 253 254]/255;
col{"geostrophic"} = [205 253 254]/255;
col{"wave"} = [205 253 197]/255;
col{"damped_wave"} = [205 253 197]/255;
col{"sink"} = [245 194 193]/255;

order = ["geostrophic", "wave", "damped_geostrophic", "damped_wave"];

[fig, boxDiagram] = wvd.plotSourcesSinksForReservoirGroup(customNames=custom_names,customColors=col,customReservoirOrder=order,shouldShowUnits=true,timeIndices=timeIndices,title="hydrostatic, geostrophic + wave",visible="off");
dampedBoxes = boxDiagram.rows{3};
dampedBoxes(1).Size = [3.5 2.0];
dampedBoxes(2).Size = [3.5 2.0];
dampedBoxes(1).Position = dampedBoxes(1).Position + [4.5 0];
dampedBoxes(2).Position = dampedBoxes(2).Position + [0.5 0];
boxDiagram.arrows(6).LabelOffset = 0.4;
boxDiagram.arrows(9).LabelOffset = 0.4;
boxDiagram.arrows(10).LabelOffset = 0.33;
boxDiagram.arrows(11).LabelOffset = 0.33;

boxDiagram.layoutArrows();
boxDiagram.arrows(4).intermediatePoints = [16, 5.5+.75; 16, -2; 4.5/2, -2];
boxDiagram.arrows(4).LabelPosition = [8 -1.90];

fig = boxDiagram.draw();