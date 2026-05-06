function overlayFrequencyContours(self,options)
% Overlay contours of nondimensional frequency on the current axes.
%
% Draws contours of the model dispersion frequency (omega/f) converted
% to the plotting coordinates and optionally labels them. Intended for
% spectrum figures to show dispersion curves.
%
% - Topic: Figures — Auxiliary functions
% - Declaration: overlayFrequencyContours(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter options.frequencies: vector of nondimensional frequency values (omega/f) to draw (default: [])
% - Parameter options.lineWidth: contour line width (default: 1)
% - Parameter options.textColor: color for contour labels (default: [0.7 0.7 0.7])
% - Parameter options.labelSpacing: label spacing passed to clabel (default: 1000)
% - Returns: None. Draws contours on the current axes.
arguments
    self
    options.frequencies = [1.01 1.05 1.2 1.5 2 4 8 16]
    options.textColor = [.5,.5,.5]
    options.labelSpacing = 600
    options.lineWidth = 1
end
omegaJK = self.omega_jk;
set(gca,'layer','top'),
hold on
% flipud() and fliplr() help trick clabel into nicer label placement.
% for y-axis, use j+1 so contours line up with pcolor cells.
[C,h] = contour(flipud(2*pi./self.kRadial(2:end)/1000),self.j',fliplr(omegaJK(:,2:end)),options.frequencies,'LineWidth',options.lineWidth,'Color',options.textColor);
clabel(C,h,options.frequencies,'Color',options.textColor,'LabelSpacing',options.labelSpacing)
end