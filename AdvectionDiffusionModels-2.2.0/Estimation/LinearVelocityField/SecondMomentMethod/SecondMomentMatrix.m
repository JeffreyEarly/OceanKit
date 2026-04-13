function [varargout] = SecondMomentMatrix(q, r, outputType)
if nargin < 3
    outputType = 'moments';
end

if ~isequal(size(q), size(r))
    error('SecondMomentMatrix:SizeMismatch', 'q and r must have the same size.');
end

Mxx = mean(q.^2,2);
Myy = mean(r.^2,2);
Mxy = mean(q.*r,2);

if strcmp(outputType,'moments')
    varargout = {Mxx, Myy, Mxy};
elseif strcmp(outputType,'eigen')
    traceMoments = Mxx + Myy;
    radiusMoments = sqrt((Mxx - Myy).^2 + 4*Mxy.^2);

    minD = 0.5*(traceMoments - radiusMoments);
    maxD = 0.5*(traceMoments + radiusMoments);
    theta = 0.5*atan2(2*Mxy, Mxx - Myy);

    varargout = {minD, maxD, theta};
else
    error('SecondMomentMatrix:UnknownOutputType', 'Unknown outputType "%s".', outputType);
end
end
