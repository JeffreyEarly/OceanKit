function psiKnotPoints = defaultPsiKnotPoints(qAll, rAll, tAll, psiS)
qMin = min(qAll);
qMax = max(qAll);
rMin = min(rAll);
rMax = max(rAll);
tMin = min(tAll);
tMax = max(tAll);

qPad = max(50, 0.25 * max(qMax - qMin, eps));
rPad = max(50, 0.25 * max(rMax - rMin, eps));
tPad = max(1, 0.01 * max(tMax - tMin, eps));

lowerBounds = [qMin - qPad, rMin - rPad, tMin - tPad];
upperBounds = [qMax + qPad, rMax + rPad, tMax + tPad];

psiKnotPoints = cell(1, 3);
for iDim = 1:3
    repeatCount = psiS(iDim) + 1;
    psiKnotPoints{iDim} = [ ...
        repmat(lowerBounds(iDim), repeatCount, 1)
        repmat(upperBounds(iDim), repeatCount, 1)
        ];
end
end
