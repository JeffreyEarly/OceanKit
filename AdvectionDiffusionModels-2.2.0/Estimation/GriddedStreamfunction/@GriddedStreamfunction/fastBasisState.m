function fastState = fastBasisState(allT, fastKnotPoints, fastS, shouldComputeDerivative)
fastState = struct();
fastState.B = reshape(BSpline.matrixForDataPoints(allT, knotPoints=fastKnotPoints, S=fastS), numel(allT), []);
fastState.solver = decomposition(fastState.B);

if shouldComputeDerivative
    if fastS > 0
        dBTensor = BSpline.matrixForDataPoints(allT, knotPoints=fastKnotPoints, S=fastS, D=1);
        fastState.dB = reshape(dBTensor(:, :, 2), numel(allT), []);
    else
        fastState.dB = zeros(size(fastState.B));
    end
else
    fastState.dB = [];
end
end
