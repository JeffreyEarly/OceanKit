function summarizeDegreesOfFreedom(self)
% Summarize the spatial grid and active spectral degrees of freedom.
%
% The summary lists each primary flow component in lexical `shortName`
% order. Mode counts come from the component's active primary spectral
% masks in the transform's current antialias configuration. Each subtotal
% is the mode count multiplied by the component's degrees of freedom per
% mode; the final line sums those subtotals. The method does not assert an
% equivalence between spatial-grid values and active spectral degrees of
% freedom.
%
% - Topic: Domain attributes — Grid
% - Declaration: summarizeDegreesOfFreedom()
arguments (Input)
    self WVTransform {mustBeNonempty}
end

dimensionNames = string(self.spatialDimensionNames());
gridSize = self.spatialMatrixSize;
gridDescription = strjoin(dimensionNames + "=" + string(gridSize),", ");

components = self.primaryFlowComponents;
[shortNames,componentOrder] = sort(string({components.shortName}));
components = components(componentOrder);
componentNames = string({components.name});
nModesByComponent = arrayfun(@(component) component.nModes,components);
degreesOfFreedomPerMode = arrayfun(@(component) component.degreesOfFreedomPerMode,components);
subtotals = nModesByComponent.*degreesOfFreedomPerMode;

componentWidth = max([strlength("Component") strlength(componentNames)]);
shortNameWidth = max([strlength("Short name") strlength(shortNames)]);

fprintf('Degrees of freedom summary\n');
fprintf('Spatial grid: %s\n\n',gridDescription);
fprintf('Primary flow components (spectral):\n');
fprintf('%-*s  %-*s  %9s  %9s  %9s\n',componentWidth,'Component',shortNameWidth,'Short name','Modes','DOF/mode','Subtotal');
fprintf('%-*s  %-*s  %9s  %9s  %9s\n',componentWidth,repmat('-',1,strlength("Component")),shortNameWidth,repmat('-',1,strlength("Short name")),repmat('-',1,strlength("Modes")),repmat('-',1,strlength("DOF/mode")),repmat('-',1,strlength("Subtotal")));
for iComponent = 1:length(components)
    fprintf('%-*s  %-*s  %9d  %9d  %9d\n',componentWidth,componentNames(iComponent),shortNameWidth,shortNames(iComponent),nModesByComponent(iComponent),degreesOfFreedomPerMode(iComponent),subtotals(iComponent));
end
fprintf('\nTotal active spectral degrees of freedom: %d\n',sum(subtotals));
end
