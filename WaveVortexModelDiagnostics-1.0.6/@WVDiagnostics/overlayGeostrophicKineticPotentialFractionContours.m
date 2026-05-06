function overlayGeostrophicKineticPotentialFractionContours(self,options)
% Overlay contour where geostrophic kinetic fraction equals a given value.
%
% Draws the contour line where KE_g / (KE_g + PE_g) equals a chosen
% fractional value (default 0.5). Useful to mark equipartition lines
% on spectral diagrams.
%
% - Topic: Figures — Auxiliary functions
% - Declaration: overlayGeostrophicKineticPotentialFractionContours(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter options.fraction: fractional contour level to draw (default: 0.5)
% - Parameter options.textColor: color for contour labels (default: [0.7 0.7 0.7])
% - Parameter options.labelSpacing: label spacing passed to clabel (default: 1000)
% - Returns: None. Draws contour on the current axes.
arguments
    self
    options.fractions = [.01,.1,.25,.75,.9,.99]
    options.textColor = [.5,.5,.5]
    options.labelSpacing = 600
    options.lineWidth = 1
end
hke = self.geo_hke_jk;
pe = self.geo_pe_jk;
fraction = hke./(hke+pe);
set(gca,'layer','top'),
hold on
% flipud() and fliplr() help trick clabel into nicer label placement.
% for y-axis, use j+1 so contours line up with pcolor cells.
[C,h] = contour(flipud(2*pi./self.kRadial(2:end)/1000),self.j(1:end)'+1,fliplr(fraction(1:end,2:end)),options.fractions,'LineWidth',options.lineWidth,'Color',options.textColor);
clabel(C,h,options.fractions,'Color',options.textColor,'LabelSpacing',options.labelSpacing)
[C,h] = contour(flipud(2*pi./self.kRadial(2:end)/1000),(self.j(1:end)')+1,fliplr(fraction(1:end,2:end)),[.5,.5],'LineWidth',2,'Color',options.textColor);
clabel(C,h,.5,'Color',options.textColor,'LabelSpacing',options.labelSpacing)
end