function mustBeDoublyPeriodicFPlane(a)
if ~isa(a,'WVGeometryDoublyPeriodic') || ~isa(a,'WVRotatingFPlane')
    error('mustBeDoublyPeriodicFPlane:invalidClass','This geostrophic component is only valid for doubly periodic geometry on a rotating f-plane.');
end
end
