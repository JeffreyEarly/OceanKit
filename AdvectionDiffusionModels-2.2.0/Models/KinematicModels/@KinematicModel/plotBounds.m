function plotBounds(self, varargin)
% Plot domain boundaries on the current axes.
%
% Finite rectangular domains are drawn as a closed box. Semi-infinite
% domains with finite `ylim` are drawn as vertical boundary lines.
%
% - Topic: Plot model diagnostics
% - Declaration: plotBounds(self,varargin)
% - Parameter varargin: forwarded line or rectangle styling arguments
if all(~isinf(self.xlim)) && all(~isinf(self.ylim)) && ~self.isXPeriodic && ~self.isYPeriodic
    rectangle('Position', [min(self.xlim) min(self.ylim) max(self.xlim) - min(self.xlim) max(self.ylim) - min(self.ylim)] / self.visualScale, varargin{:});
elseif all(~isinf(self.ylim)) && ~self.isYPeriodic
    x = [min(self.xVisualLimits) min(self.xVisualLimits); max(self.xVisualLimits) max(self.xVisualLimits)];
    y = [min(self.yVisualLimits) max(self.yVisualLimits); min(self.yVisualLimits) max(self.yVisualLimits)];
    line(x / self.visualScale, y / self.visualScale, varargin{:});
end
