function showRossbyRadiusYAxis(self,options)
% Annotate the axes with a Rossby radius y-axis label.
%
% Places a textual annotation (L_r) on the current axes indicating the
% Rossby radius scale in kilometers. Positioning and color can be
% adjusted via options.
%
% - Topic: Figures — Auxiliary functions
% - Declaration: showRossbyRadiusYAxis(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter options.textColor: color for the annotation text (default: 'k')
% - Parameter options.xFrac: fractional x position in axes coordinates (default: 0.02)
% - Parameter options.yFrac: fractional y position in axes coordinates (default: 0.5)
% - Returns: handle to the created text object
arguments
    self
    options.textColor = [.5,.5,.5]
end
set(gca,'Layer','top','TickLength',[0.015 0.015])
% create some nice tick labels to show deformation radius
yticksTemp = yticks;
ticks_y = sqrt(self.wvt.Lr2)./1000;
labels_y = cell(length(yticksTemp),1);
for i=1:length(yticksTemp)
    labels_y{i} = sprintf('%0.1f',ticks_y(yticksTemp(i)+1));
end
text(.7*min(xlim)*ones(size(yticksTemp)),yticksTemp,labels_y,'Color',options.textColor,'HorizontalAlignment','center')
text(.7*min(xlim),1.1*max(ylim),'L_r (km)','Color',options.textColor,'HorizontalAlignment','center')
end