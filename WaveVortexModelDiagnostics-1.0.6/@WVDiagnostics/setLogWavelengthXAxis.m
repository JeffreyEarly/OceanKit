function setLogWavelengthXAxis(self,options)
% Configure axis tick labels for a log-scaled wavelength x-axis.
%
% Converts internal radial-wavenumber tick values into human-readable
% pseudo-wavelength labels (kilometers) and applies them to the current
% axes. Useful for spectral plots where the x-axis uses a log-scaled
% wavenumber coordinate but labels should show wavelengths in km.
%
% - Topic: Figures — Auxiliary functions
% - Declaration: setLogWavelengthXAxis(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter options.num_ticks: (optional) Number of tick labels to generate (default: 6)
% - Parameter options.roundToNearest: (optional) Round wavelength labels to this nearest integer (km) (default: 5)
% - Returns: None. Sets xticks and xticklabels on current axes.
arguments
    self WVDiagnostics
    options.num_ticks = 6
    options.roundToNearest = 5
end
[labels_x,ticks_x] = self.logWavelengthAxis(num_ticks=options.num_ticks,roundToNearest=options.roundToNearest);
xscale('log')
xticks(ticks_x)
xticklabels(labels_x)
end