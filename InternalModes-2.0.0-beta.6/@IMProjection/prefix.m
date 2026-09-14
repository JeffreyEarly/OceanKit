function projection = prefix(self,count)
% Construct the projection of a leading set of array columns.
%
% This operation uses array positions. The caller owns whether a leading
% subset is scientifically meaningful; endpoint labels are not mode counts.
%
% - Topic: Create projections
% - Declaration: projection = prefix(self,count)
% - Parameter count: nonnegative number of leading columns
% - Returns projection: projection on the same fixed sample pairing
arguments (Input)
    self (1,1) IMProjection
    count (1,1) double {mustBeInteger,mustBeNonnegative}
end
arguments (Output)
    projection (1,1) IMProjection
end
if count > self.columnCount
    error("IMProjection:InvalidColumnCount","count must not exceed the available column count.");
end
labels = self.columnLabels;
if ~isempty(labels)
    labels = labels(1:count);
end
if self.projectionKind == "prescribedDual"
    projection = IMProjection.fromPrescribedDual(self.samplePairingMatrix(1:count,:),self.gramMatrix(1:count,1:count),self.targetGramMatrix(1:count,1:count),majorantGramMatrix=self.majorantGramMatrix(1:count,1:count),activeColumnMask=self.activeColumnMask(1:count),columnLabels=labels,provenance=self.provenance);
    return
end
projection = IMProjection(self.sampledBasis(:,1:count),self.metricMatrix,self.targetGramMatrix(1:count,1:count),majorantGramMatrix=self.majorantGramMatrix(1:count,1:count),activeColumnMask=self.activeColumnMask(1:count),columnLabels=labels,provenance=self.provenance);
end
