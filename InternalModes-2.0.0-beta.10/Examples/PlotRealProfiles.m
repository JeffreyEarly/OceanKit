exampleDir = fileparts(mfilename('fullpath'));
load(fullfile(exampleDir,'SampleLatmixProfiles.mat'));

i=1;
% for i=1:2:length(rhoProfile)
    im = InternalModes(rhoProfile{i},zProfile{i},zProfile{i},latitude,'method','finiteDifference');
    im.upperBoundary = UpperBoundary.rigidLid;
    im.lowerBoundary = LowerBoundary.freeSlip;
    %     im.showLowestModesAtWavenumber(0);
    im.showLowestModesAtFrequency(im.f0);
% end
