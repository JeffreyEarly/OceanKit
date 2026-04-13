function [xCenter, yCenter, q, r] = CenterOfMass(x, y)
if ~isequal(size(x), size(y))
    error('CenterOfMass:SizeMismatch', 'x and y must have the same size.');
end

xCenter = mean(x,2);
yCenter = mean(y,2);
q = x - xCenter;
r = y - yCenter;
end
