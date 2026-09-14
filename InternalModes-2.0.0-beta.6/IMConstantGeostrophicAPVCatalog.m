classdef (Sealed, Hidden) IMConstantGeostrophicAPVCatalog
    % Exact constant-stratification generalized-energy APV formulas.

    methods (Static)
        function tf = isTargetEVP(evp)
            tf = isa(evp, "IMInternalModes") && string(evp.name) == "geostrophicAPVModes";
        end

        function validateEVP(evp, N0, zDomain)
            if ~IMConstantGeostrophicAPVCatalog.isTargetEVP(evp)
                error("IMConstantStratificationSolution:UnsupportedEVP", "The generalized-energy APV analytical branch requires name=""geostrophicAPVModes"".");
            end
            if evp.formulation ~= "F" || evp.modeFamily ~= "hydrostatic"
                error("IMConstantStratificationSolution:UnsupportedEVP", "Generalized-energy APV constant modes require the hydrostatic F formulation.");
            end
            requiredParameters = ["g0", "gd", "surfaceBoundary"];
            for parameter = requiredParameters
                if ~isfield(evp.parameters, parameter)
                    error("IMConstantStratificationSolution:UnsupportedEVP", "The APV EVP is missing parameters.%s.", parameter);
                end
            end

            g0 = evp.parameters.g0;
            gd = evp.parameters.gd;
            if ~(isscalar(g0) && isreal(g0) && ~isnan(g0) && g0 ~= -Inf)
                error("IMConstantStratificationSolution:UnsupportedEVP", "parameters.g0 must be signed finite, zero, or Inf.");
            end
            if ~(isscalar(gd) && isreal(gd) && ~isnan(gd) && gd ~= -Inf)
                error("IMConstantStratificationSolution:UnsupportedEVP", "parameters.gd must be signed finite, zero, or Inf.");
            end
            surfaceConvention = string(evp.parameters.surfaceBoundary);
            if ~isscalar(surfaceConvention) || ~ismember(surfaceConvention, ["freeSurface", "rigidLid"])
                error("IMConstantStratificationSolution:UnsupportedEVP", "parameters.surfaceBoundary must be ""freeSurface"" or ""rigidLid"".");
            end

            expectedSurface = IMConstantGeostrophicAPVCatalog.surfaceCondition(evp.g, g0, surfaceConvention);
            expectedBottom = IMConstantGeostrophicAPVCatalog.bottomCondition(gd);
            if ~IMConstantGeostrophicAPVCatalog.equivalentBoundary(evp.surfaceBoundary, expectedSurface) || ~IMConstantGeostrophicAPVCatalog.equivalentBoundary(evp.bottomBoundary, expectedBottom)
                error("IMConstantStratificationSolution:UnsupportedBoundary", "The APV boundary coefficients do not match the declared g0, gd, and surface convention.");
            end

            z = linspace(zDomain(1), zDomain(2), 5).';
            expectedN2 = N0*N0*ones(size(z));
            actualN2 = evp.N2(z);
            tolerance = 1e4*eps(max([1; abs(expectedN2(:)); abs(actualN2(:))]));
            if max(abs(actualN2(:) - expectedN2(:))) > tolerance
                error("IMConstantStratificationSolution:UnsupportedStratification", "The APV EVP N2 function does not match this constant analytical solution.");
            end

            context = IMConstantGeostrophicAPVCatalog.context(evp);
            p = IMEigenvalueProblem.evaluateCoefficient(evp.p, z, context);
            q = IMEigenvalueProblem.evaluateCoefficient(evp.q, z, context);
            r = IMEigenvalueProblem.evaluateCoefficient(evp.r, z, context);
            if isscalar(p)
                p = p*ones(size(z));
            end
            if isscalar(q)
                q = q*ones(size(z));
            end
            if isscalar(r)
                r = r*ones(size(z));
            end
            coefficientTolerance = 1e4*eps(max([1; abs(p(:)); abs(q(:)); abs(r(:)); abs(1./expectedN2(:)); 1/evp.g]));
            if max(abs(p(:) - 1./expectedN2(:))) > coefficientTolerance || max(abs(q(:))) > coefficientTolerance || max(abs(r(:) - 1/evp.g)) > coefficientTolerance
                error("IMConstantStratificationSolution:UnsupportedEVP", "The APV EVP coefficients must be p=1/N2, q=0, and r=1/g.");
            end
        end

        function [modeData, h, modeNumber, metadata, diagnostics] = solve(evp, N0, zDomain, nModes)
            IMConstantGeostrophicAPVCatalog.validateEVP(evp, N0, zDomain);
            endpoint = IMConstantGeostrophicAPVCatalog.endpointData(N0, zDomain, evp.g, evp.parameters.g0, evp.parameters.gd, string(evp.parameters.surfaceBoundary));
            [negativeCount, zeroPresent, endpointMatrix] = IMConstantGeostrophicAPVCatalog.endpointInertia(endpoint);

            negativeRoots = IMConstantGeostrophicAPVCatalog.negativeRoots(endpoint, negativeCount, zeroPresent);
            zeroCoefficients = zeros(1,2);
            if zeroPresent
                zeroCoefficients = IMConstantGeostrophicAPVCatalog.zeroCoefficients(endpoint);
            end
            positiveCount = max(0, nModes - negativeCount - double(zeroPresent));
            positiveRoots = IMConstantGeostrophicAPVCatalog.positiveRoots(endpoint, positiveCount, zeroPresent);

            roots = [negativeRoots(:).' NaN(1,double(zeroPresent)) positiveRoots(:).'];
            modeKinds = [repmat("negative",1,length(negativeRoots)) repmat("zero",1,double(zeroPresent)) repmat("positive",1,length(positiveRoots))];
            scale = N0*N0*endpoint.depth*endpoint.depth/evp.g;
            h = [-scale./(negativeRoots(:).'.^2) Inf(1,double(zeroPresent)) scale./(positiveRoots(:).'.^2)];
            coefficients = zeros(length(roots),2);
            for iMode = 1:length(roots)
                if modeKinds(iMode) == "zero"
                    coefficients(iMode,:) = zeroCoefficients;
                else
                    coefficients(iMode,:) = IMConstantGeostrophicAPVCatalog.branchCoefficients(endpoint, roots(iMode), modeKinds(iMode));
                end
            end

            eigenvalues = 1./h;
            [~, order] = sort(eigenvalues, "ascend");
            roots = roots(order);
            modeKinds = modeKinds(order);
            h = h(order);
            coefficients = coefficients(order,:);
            keep = 1:min(nModes,length(h));
            roots = roots(keep);
            modeKinds = modeKinds(keep);
            h = h(keep);
            coefficients = coefficients(keep,:);

            retainedNegativeCount = nnz(modeKinds == "negative");
            modeNumber = zeros(1,length(h));
            modeNumber(1:retainedNegativeCount) = -1:-1:-retainedNegativeCount;
            nextMode = retainedNegativeCount + 1;
            if nextMode <= length(h) && modeKinds(nextMode) == "zero"
                nextMode = nextMode + 1;
            end
            modeNumber(nextMode:end) = 1:(length(h) - nextMode + 1);

            modeData = struct("roots", roots, "modeKinds", modeKinds, "coefficients", coefficients, "endpoint", endpoint, "signFactors", ones(1,length(h)));
            zOrientation = linspace(zDomain(1), zDomain(2), 1025).';
            [FOrientation, ~, ~] = IMConstantGeostrophicAPVCatalog.rawValues(modeData, N0, evp.g, zOrientation, false);
            for iMode = 1:length(h)
                [~, index] = max(abs(FOrientation(:,iMode)));
                if FOrientation(index,iMode) < 0
                    modeData.signFactors(iMode) = -1;
                end
            end

            rootResiduals = NaN(1,length(roots));
            for iMode = 1:length(roots)
                if modeKinds(iMode) ~= "zero"
                    rootResiduals(iMode) = IMConstantGeostrophicAPVCatalog.rootResidual(endpoint, roots(iMode), modeKinds(iMode));
                end
            end
            metadata = struct("roots", roots, "rootResiduals", rootResiduals, "modeKinds", modeKinds, "negativeModeCount", negativeCount, "zeroModePresent", zeroPresent, "endpointMatrix", endpointMatrix, "endpointData", endpoint, "g0", evp.parameters.g0, "gd", evp.parameters.gd, "surfaceBoundary", string(evp.parameters.surfaceBoundary));
            zeroModeStatus = "absent";
            if zeroPresent
                zeroModeStatus = "present";
            end
            diagnostics = struct("negativeModeCount", negativeCount, "zeroModeStatus", zeroModeStatus, "zeroModeCount", double(zeroPresent), "reason", "Exact constant generalized-energy APV endpoint inertia.");
        end

        function values = rawVariable(modeData, N0, g, variable, z)
            [F, G, ~] = IMConstantGeostrophicAPVCatalog.rawValues(modeData, N0, g, z(:), true);
            if string(variable) == "F"
                values = F;
            elseif string(variable) == "G"
                values = G;
            else
                error("IMConstantStratificationSolution:UnsupportedVariable", "Variable must be ""F"" or ""G"".");
            end
        end

        function values = rawFz(modeData, N0, g, z)
            [~, ~, values] = IMConstantGeostrophicAPVCatalog.rawValues(modeData, N0, g, z(:), true);
        end
    end

    methods (Static, Access = private)
        function context = context(evp)
            context.N2 = @(z) evp.N2(z);
            context.g = evp.g;
            context.f0 = evp.f0;
            context.formulation = evp.formulation;
            fields = fieldnames(evp.parameters);
            for iField = 1:length(fields)
                context.(fields{iField}) = evp.parameters.(fields{iField});
            end
        end

        function endpoint = endpointData(N0, zDomain, g, g0, gd, surfaceConvention)
            endpoint.depth = diff(zDomain);
            endpoint.zBottom = zDomain(1);
            endpoint.N0 = N0;
            endpoint.IN = N0*N0*endpoint.depth;
            endpoint.surfaceKind = "active";
            endpoint.bottomKind = "active";
            endpoint.surfaceConvention = surfaceConvention;
            endpoint.g0 = g0;
            endpoint.gd = gd;
            if g0 == 0
                endpoint.surfaceKind = "dirichlet";
                endpoint.betaSurface = NaN;
            elseif surfaceConvention == "freeSurface"
                if isinf(g0)
                    endpoint.betaSurface = 1/g;
                else
                    endpoint.betaSurface = 1/g + 1/g0;
                end
            elseif isinf(g0)
                endpoint.betaSurface = 0;
            else
                endpoint.betaSurface = 1/g0;
            end
            if gd == 0
                endpoint.bottomKind = "dirichlet";
                endpoint.betaBottom = NaN;
            elseif isinf(gd)
                endpoint.betaBottom = 0;
            else
                endpoint.betaBottom = 1/gd;
            end
        end

        function [negativeCount, zeroPresent, matrix] = endpointInertia(endpoint)
            if endpoint.surfaceKind == "active" && endpoint.bottomKind == "active"
                matrix = [1 + endpoint.IN*endpoint.betaSurface, -1; -1, 1 + endpoint.IN*endpoint.betaBottom];
                eigenvalues = eig(matrix);
            elseif endpoint.surfaceKind == "dirichlet" && endpoint.bottomKind == "active"
                matrix = 1 + endpoint.IN*endpoint.betaBottom;
                eigenvalues = matrix;
            elseif endpoint.surfaceKind == "active" && endpoint.bottomKind == "dirichlet"
                matrix = 1 + endpoint.IN*endpoint.betaSurface;
                eigenvalues = matrix;
            else
                matrix = zeros(0,0);
                eigenvalues = zeros(0,1);
            end
            tolerance = 1e3*eps(max([1; abs(eigenvalues(:))]));
            negativeCount = nnz(eigenvalues < -tolerance);
            zeroPresent = any(abs(eigenvalues) <= tolerance);
        end

        function coefficients = zeroCoefficients(endpoint)
            N2 = endpoint.N0*endpoint.N0;
            if endpoint.bottomKind == "dirichlet"
                bottomRow = [1 0];
            else
                bottomRow = [-endpoint.betaBottom 1/N2];
            end
            if endpoint.surfaceKind == "dirichlet"
                surfaceRow = [1 endpoint.depth];
            else
                surfaceRow = [endpoint.betaSurface 1/N2 + endpoint.betaSurface*endpoint.depth];
            end
            residual = det([bottomRow; surfaceRow]);
            tolerance = 1e3*eps(max([1 abs(bottomRow) abs(surfaceRow)]));
            if abs(residual) > tolerance
                error("IMConstantStratificationSolution:ZeroModeMismatch", "The endpoint inertia predicted a zero mode but the exact affine determinant is nonzero.");
            end
            coefficients = [bottomRow(2) -bottomRow(1)];
            coefficients = coefficients/norm(coefficients);
        end

        function roots = negativeRoots(endpoint, count, zeroPresent)
            if count == 0
                roots = zeros(1,0);
                return;
            end
            if count == 2
                strengths = [endpoint.betaSurface endpoint.betaBottom]*endpoint.IN;
                upperBound = max(10, 2*max(abs(strengths)) + 10);
                equalityTolerance = 1e3*eps(max([1 abs(strengths)]));
                if endpoint.surfaceKind == "active" && endpoint.bottomKind == "active" && abs(diff(strengths)) <= equalityTolerance
                    strength = mean(strengths);
                    evenResidual = @(value) value*tanh(value/2) + strength;
                    oddResidual = @(value) value + strength*tanh(value/2);
                    roots = [fzero(evenResidual,[1e-7 upperBound]) fzero(oddResidual,[1e-7 upperBound])];
                    return;
                end
                residual = @(value) IMConstantGeostrophicAPVCatalog.negativeCharacteristic(endpoint, value);
                lowerBound = 1e-7;
                [minimumLocation, minimumValue] = fminbnd(residual, lowerBound, upperBound, optimset("TolX",1e-12));
                while residual(upperBound) <= 0
                    upperBound = upperBound*2;
                end
                if minimumValue < 0 && residual(lowerBound) > 0
                    roots = [fzero(residual,[lowerBound minimumLocation]) fzero(residual,[minimumLocation upperBound])];
                    return;
                end
            end
            upperBound = 10;
            roots = zeros(1,0);
            for iteration = 1:8
                x = logspace(-8, log10(upperBound), 10000).';
                residual = @(value) IMConstantGeostrophicAPVCatalog.branchResidual(endpoint, value, "negative");
                roots = IMConstantGeostrophicAPVCatalog.scanRoots(residual, x);
                if zeroPresent
                    roots = roots(roots > 1e-5);
                end
                if length(roots) >= count
                    roots = roots(end-count+1:end);
                    break;
                end
                upperBound = upperBound*10;
            end
            if length(roots) ~= count
                error("IMConstantStratificationSolution:RootSearchFailed", "Expected %d negative APV roots but found %d.", count, length(roots));
            end
            roots = reshape(roots,1,[]);
        end

        function value = negativeCharacteristic(endpoint, root)
            surfaceStrength = endpoint.betaSurface*endpoint.IN;
            bottomStrength = endpoint.betaBottom*endpoint.IN;
            if endpoint.surfaceKind == "active" && endpoint.bottomKind == "active"
                value = root*(surfaceStrength + bottomStrength) + (root*root + surfaceStrength*bottomStrength)*tanh(root);
            elseif endpoint.surfaceKind == "dirichlet"
                value = root + bottomStrength*tanh(root);
            elseif endpoint.bottomKind == "dirichlet"
                value = root + surfaceStrength*tanh(root);
            else
                value = sinh(root);
            end
        end

        function roots = positiveRoots(endpoint, count, zeroPresent)
            if count == 0
                roots = zeros(1,0);
                return;
            end
            upperBound = max(20, (count + 4)*pi);
            roots = zeros(1,0);
            for iteration = 1:8
                x = unique([logspace(-8,0,1024) linspace(1,upperBound,max(8192,512*(count + 4)))]).';
                residual = @(value) IMConstantGeostrophicAPVCatalog.branchResidual(endpoint, value, "positive");
                roots = IMConstantGeostrophicAPVCatalog.scanRoots(residual, x);
                if zeroPresent
                    roots = roots(roots > 1e-5);
                end
                if length(roots) >= count
                    roots = roots(1:count);
                    break;
                end
                upperBound = upperBound*1.75;
            end
            if length(roots) ~= count
                error("IMConstantStratificationSolution:RootSearchFailed", "Expected %d positive APV roots but found %d.", count, length(roots));
            end
            roots = reshape(roots,1,[]);
        end

        function value = branchResidual(endpoint, root, branch)
            [bottomRow, surfaceRow] = IMConstantGeostrophicAPVCatalog.branchRows(endpoint, root, branch);
            bottomNorm = norm(bottomRow);
            surfaceNorm = norm(surfaceRow);
            if bottomNorm == 0 || surfaceNorm == 0 || ~isfinite(bottomNorm) || ~isfinite(surfaceNorm)
                value = NaN;
                return;
            end
            value = det([bottomRow/bottomNorm; surfaceRow/surfaceNorm]);
        end

        function value = rootResidual(endpoint, root, branch)
            if branch == "negative" && endpoint.surfaceKind == "active" && endpoint.bottomKind == "active"
                strengths = [endpoint.betaSurface endpoint.betaBottom]*endpoint.IN;
                equalityTolerance = 1e3*eps(max([1 abs(strengths)]));
                if abs(diff(strengths)) <= equalityTolerance
                    strength = mean(strengths);
                    if root >= -strength
                        value = (root*tanh(root/2) + strength)/max([1 abs(root) abs(strength)]);
                    else
                        value = (root + strength*tanh(root/2))/max([1 abs(root) abs(strength)]);
                    end
                    return;
                end
            end
            value = IMConstantGeostrophicAPVCatalog.branchResidual(endpoint, root, branch);
        end

        function coefficients = branchCoefficients(endpoint, root, branch)
            if branch == "negative" && endpoint.surfaceKind == "active" && endpoint.bottomKind == "active"
                strengths = [endpoint.betaSurface endpoint.betaBottom]*endpoint.IN;
                equalityTolerance = 1e3*eps(max([1 abs(strengths)]));
                if abs(diff(strengths)) <= equalityTolerance
                    if root >= -mean(strengths)
                        coefficients = [1 0];
                    else
                        coefficients = [0 1];
                    end
                    return;
                end
            end
            [bottomRow, ~] = IMConstantGeostrophicAPVCatalog.branchRows(endpoint, root, branch);
            coefficients = [bottomRow(2) -bottomRow(1)];
            coefficients = coefficients/norm(coefficients);
        end

        function [bottomRow, surfaceRow] = branchRows(endpoint, root, branch)
            N2 = endpoint.N0*endpoint.N0;
            if branch == "positive"
                valueBottom = [1 0];
                fluxBottom = root/(endpoint.depth*N2)*[0 1];
                valueSurface = [cos(root) sin(root)];
                fluxSurface = root/(endpoint.depth*N2)*[-sin(root) cos(root)];
            else
                scaledCosh = (1 + exp(-root))/2;
                scaledSinh = (1 - exp(-root))/2;
                valueBottom = [scaledCosh -scaledSinh];
                fluxBottom = root/(endpoint.depth*N2)*[-scaledSinh scaledCosh];
                valueSurface = [scaledCosh scaledSinh];
                fluxSurface = root/(endpoint.depth*N2)*[scaledSinh scaledCosh];
            end
            if endpoint.bottomKind == "dirichlet"
                bottomRow = valueBottom;
            else
                bottomRow = fluxBottom - endpoint.betaBottom*valueBottom;
            end
            if endpoint.surfaceKind == "dirichlet"
                surfaceRow = valueSurface;
            else
                surfaceRow = fluxSurface + endpoint.betaSurface*valueSurface;
            end
        end

        function roots = scanRoots(residual, x)
            y = arrayfun(residual, x);
            roots = zeros(length(x)-1,1);
            rootCount = 0;
            for iPoint = 1:(length(x)-1)
                if ~isfinite(y(iPoint)) || ~isfinite(y(iPoint+1))
                    continue;
                end
                if y(iPoint) == 0
                    rootCount = rootCount + 1;
                    roots(rootCount) = x(iPoint);
                elseif sign(y(iPoint)) ~= sign(y(iPoint+1))
                    try
                        root = fzero(residual, [x(iPoint) x(iPoint+1)]);
                        if isfinite(root) && abs(residual(root)) <= 1e-8
                            rootCount = rootCount + 1;
                            roots(rootCount) = root;
                        end
                    catch
                    end
                end
            end
            roots = IMConstantGeostrophicAPVCatalog.deduplicateRoots(roots(1:rootCount));
        end

        function roots = deduplicateRoots(roots)
            roots = sort(real(roots(:)));
            roots = roots(isfinite(roots) & roots > 0);
            if isempty(roots)
                return;
            end
            tolerance = 1e-8*max(1,max(abs(roots)));
            roots = roots([true; abs(diff(roots)) > tolerance]);
        end

        function [F, G, Fz] = rawValues(modeData, N0, g, z, applySigns)
            z = z(:);
            x = z - modeData.endpoint.zBottom;
            xi = x/modeData.endpoint.depth;
            centeredXi = xi - 0.5;
            F = zeros(length(z),length(modeData.roots));
            Fz = zeros(size(F));
            for iMode = 1:length(modeData.roots)
                coefficients = modeData.coefficients(iMode,:);
                if modeData.modeKinds(iMode) == "zero"
                    F(:,iMode) = coefficients(1) + coefficients(2)*x;
                    Fz(:,iMode) = coefficients(2);
                elseif modeData.modeKinds(iMode) == "positive"
                    phase = modeData.roots(iMode)*xi;
                    F(:,iMode) = coefficients(1)*cos(phase) + coefficients(2)*sin(phase);
                    Fz(:,iMode) = modeData.roots(iMode)/modeData.endpoint.depth*(-coefficients(1)*sin(phase) + coefficients(2)*cos(phase));
                else
                    root = modeData.roots(iMode);
                    positiveExponential = exp(root*(centeredXi - 0.5));
                    negativeExponential = exp(-root*(centeredXi + 0.5));
                    scaledCosh = (positiveExponential + negativeExponential)/2;
                    scaledSinh = (positiveExponential - negativeExponential)/2;
                    F(:,iMode) = coefficients(1)*scaledCosh + coefficients(2)*scaledSinh;
                    Fz(:,iMode) = root/modeData.endpoint.depth*(coefficients(1)*scaledSinh + coefficients(2)*scaledCosh);
                end
            end
            G = -(g/(N0*N0))*Fz;
            if applySigns
                F = F.*modeData.signFactors;
                G = G.*modeData.signFactors;
                Fz = Fz.*modeData.signFactors;
            end
        end

        function boundary = surfaceCondition(g, g0, surfaceConvention)
            if g0 == 0
                boundary = IMBoundaryCondition.dirichlet();
            elseif surfaceConvention == "freeSurface"
                if isinf(g0)
                    boundary = IMBoundaryCondition(a=-1/g, b=1);
                else
                    boundary = IMBoundaryCondition(a=-(1/g + 1/g0), b=1);
                end
            elseif isinf(g0)
                boundary = IMBoundaryCondition.neumann();
            else
                boundary = IMBoundaryCondition(a=-1/g0, b=1);
            end
        end

        function boundary = bottomCondition(gd)
            if gd == 0
                boundary = IMBoundaryCondition.dirichlet();
            elseif isinf(gd)
                boundary = IMBoundaryCondition.neumann();
            else
                boundary = IMBoundaryCondition(a=1/gd, b=1);
            end
        end

        function tf = equivalentBoundary(actual, expected)
            actualCoefficients = [actual.a actual.b actual.c actual.d];
            expectedCoefficients = [expected.a expected.b expected.c expected.d];
            active = find(abs(expectedCoefficients) > 0, 1, "first");
            if isempty(active)
                tf = false;
                return;
            end
            scale = actualCoefficients(active)/expectedCoefficients(active);
            tolerance = 1e3*eps(max([1 abs(actualCoefficients) abs(scale*expectedCoefficients)]));
            tf = isfinite(scale) && scale ~= 0 && max(abs(actualCoefficients - scale*expectedCoefficients)) <= tolerance;
        end
    end
end
