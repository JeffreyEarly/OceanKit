function transformedSpline = affineOutputTransform(self, scale, offset)
% Apply an affine transform to spline outputs without refitting.
%
% - Topic: Transform the spline
% - Declaration: transformedSpline = affineOutputTransform(self,scale,offset)
% - Parameter self: BSpline instance
% - Parameter scale: output scale factor
% - Parameter offset: output offset
% - Returns transformedSpline: BSpline with adjusted output normalization
arguments
    self (1,1) BSpline
    scale (1,1) double
    offset (1,1) double
end

transformedSpline = BSpline(S=self.S, knotPoints=self.knotPoints, xi=self.xi);
transformedSpline.xStd = scale*self.xStd;
transformedSpline.xMean = scale*self.xMean + offset;
