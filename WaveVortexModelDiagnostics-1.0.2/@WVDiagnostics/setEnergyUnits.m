function setEnergyUnits(self, units)
% Set the time and energy scaling and units for plotting and output.
%
% Sets tscale, tscale_units, escale, and escale_units based on the specified units.
%
% - Topic: Units
% - Declaration: setEnergyUnits(self, units)
% - Parameter units: must be one of "si", "gm", or "si-yr"
% - Returns: None
arguments
    self
    units {mustBeMember(units, ["si", "gm", "si-yr"])}
end
switch lower(units)
    case "si"
        self.escale = 1;
        self.escale_units = "m^3 s^{-2}";
        self.flux_scale = 1;
        self.flux_scale_units = "m^3 s^{-3}";
    case "gm"
        self.escale = 3.74;
        self.escale_units = "GM";
        self.flux_scale = 3.74/(86400*365);
        self.flux_scale_units = "GM/yr";
    case "si-yr"
        self.escale = 1;
        self.escale_units = "m^3 s^{-2}";
        self.flux_scale = 1/(86400*365);
        self.flux_scale_units = "m^3 s^{-2} yr^{-1}";
    otherwise
        error("Unknown units: %s. Must be 'si', 'gm', or 'si-yr'.", units);
end
end