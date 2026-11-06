function points = referenceBreakpoints(self,zBounds)
% Split reference integration at each native reconstruction knot.
arguments (Input)
    self (1,1) IMBasisSet
    zBounds (1,2) double = self.zDomain
end
arguments (Output)
    points (:,1) double
end
native = self.solver.zNative(:);
points = unique([zBounds(:);native(native > zBounds(1) & native < zBounds(2))]);
end
