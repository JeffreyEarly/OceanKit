function overlayGeostrophicKineticPotentialRatioContours(self,options)
% Overlay contours of geostrophic kinetic-to-potential energy ratio.
%
% Adds contours showing the local ratio KE_g/(KE_g+PE_g) on the
% current axes. Useful to indicate balanced vs unbalanced regions on
% energy-spectrum plots.
%
% - Topic: Figures — Auxiliary functions
% - Declaration: overlayGeostrophicKineticPotentialRatioContours(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter options.ratios: vector of ratio levels to contour (default: [.25 .5 .75])
% - Parameter options.lineWidth: contour line width (default: 1)
% - Parameter options.textColor: color for contour labels (default: [0.7 0.7 0.7])
% - Parameter options.labelSpacing: label spacing passed to clabel (default: 1000)
% - Returns: None. Draws contours on the current axes.
arguments
    self
    options.ratios = [-2.5 -2.0 -1.5 -1.0 -0.5 0 0.5 1.0 1.5 2.0 2.5]
    options.textColor = [.5,.5,.5]
    options.labelSpacing = 400
    options.lineWidth = 1
end
hke = self.geo_hke_jk;
pe = self.geo_pe_jk;
ratio = log10(hke./pe);
set(gca,'layer','top'),
hold on
% [C,h] = contour(self.kRadial(2:end),self.j(2:end)',(ratio(2:end,2:end)),options.ratios,'LineWidth',options.lineWidth,'Color',options.textColor);
[C,h] = contour(2*pi./self.kRadial(2:end)/1000,self.j(2:end)',(ratio(2:end,2:end)),options.ratios,'LineWidth',options.lineWidth,'Color',options.textColor);
clabel(C,h,options.ratios,'Color',options.textColor,'LabelSpacing',options.labelSpacing)
[C,h] = contour(2*pi./self.kRadial(2:end)/1000,self.j(1:end)',(ratio(1:end,2:end)),[-3,3],'LineWidth',options.lineWidth,'Color',options.textColor);
clabel(C,h,options.ratios,'Color',options.textColor,'LabelSpacing',options.labelSpacing)
end