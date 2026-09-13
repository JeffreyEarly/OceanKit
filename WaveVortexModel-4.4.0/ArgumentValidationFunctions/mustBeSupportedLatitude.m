function mustBeSupportedLatitude(latitude)
% Validate the latitude domain supported by WaveVortexModel 4.2.x.

mustBeGreaterThanOrEqual(abs(latitude),5)
mustBeLessThanOrEqual(abs(latitude),85)
end
