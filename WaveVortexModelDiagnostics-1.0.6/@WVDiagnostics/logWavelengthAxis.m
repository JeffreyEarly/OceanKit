function [labels, ticks] = logWavelengthAxis(self,options)
% Produce tick labels and positions for a log-scaled wavelength x-axis.
%
% Returns formatted label strings (pseudo-wavelength in km) and numeric
% tick positions suitable for use with xticks/xticklabels on spectral
% figures. The number of ticks and rounding of the displayed wavelength
% values are controlled via options.
%
% - Topic: Figures — Auxiliary functions
% - Declaration: [labels, ticks] = logWavelengthAxis(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter options.num_ticks: (optional) number of tick labels to generate (default: 6)
% - Parameter options.roundToNearest: (optional) round wavelength labels to this nearest integer (km) (default: 5)
% - Returns labels: cell array of label strings (km)
% - Returns ticks: numeric vector of tick positions (units match the x-axis used in spectral plots)
arguments
    self WVDiagnostics
    options.num_ticks = 6
    options.roundToNearest = 5
end
ticks = logspace(log10(self.kRadial(2)),log10(self.kRadial(end)),options.num_ticks);
ticks = round(2*pi./(1e3.*ticks)/options.roundToNearest)*options.roundToNearest;
labels = cell(length(ticks),1);
for i=1:length(ticks)
    labels{i} = sprintf('%.0f',ticks(i));
end
ticks = 2*pi./(1e3*ticks);
end