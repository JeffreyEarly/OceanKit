function fig = plotEnstrophyOverTime(self,options)
% Plot enstrophy over time
%
% Plots the quadratic and APV enstrophy as a function of time.
%
% - Topic: Figures — Potential Enstrophy
% - Declaration: fig = plotEnstrophyOverTime(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter visible: (optional) Figure visibility (default: "on")
% - Parameter timeIndices: (optional) Indices of model times to plot/average (default: Inf -> all times)
% - Returns fig: handle to the generated figure
arguments
    self WVDiagnostics
    options.visible = "on"
    options.timeIndices = Inf;
end
[Z_quadratic, t] = self.quadraticEnstrophyOverTime(timeIndices=options.timeIndices);
[Z_apv, ~] = self.exactEnstrophyOverTime(timeIndices=options.timeIndices);

fig = figure(Visible=options.visible);
plot(t/self.tscale,Z_quadratic/self.zscale,LineWidth=2), hold on
plot(t/self.tscale,Z_apv/self.zscale,LineWidth=2)
legend('quadratic','apv')

xlabel("time (" + self.tscale_units + ")")
ylabel("enstrophy (" + self.zscale_units + ")")
xlim([min(t) max(t)]/self.tscale);
end