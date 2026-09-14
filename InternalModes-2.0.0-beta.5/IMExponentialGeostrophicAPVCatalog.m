classdef (Sealed, Hidden) IMExponentialGeostrophicAPVCatalog
    % Exact exponential-stratification generalized-energy APV formulas.

    methods (Static)
        function tf = isTargetEVP(evp)
            tf = isa(evp, "IMInternalModes") && string(evp.name) == "geostrophicAPVModes";
        end

        function validateEVP(evp, N0, b, zDomain)
            if ~IMExponentialGeostrophicAPVCatalog.isTargetEVP(evp)
                error("IMExponentialStratificationSolution:UnsupportedEVP", ...
                    "The generalized-energy APV analytical branch requires name=""geostrophicAPVModes"".");
            end
            if evp.formulation ~= "F" || evp.modeFamily ~= "hydrostatic"
                error("IMExponentialStratificationSolution:UnsupportedEVP", ...
                    "Generalized-energy APV exponential modes require the hydrostatic F formulation.");
            end
            requiredParameters = ["g0", "gd", "surfaceBoundary"];
            for parameter = requiredParameters
                if ~isfield(evp.parameters, parameter)
                    error("IMExponentialStratificationSolution:UnsupportedEVP", ...
                        "The APV EVP is missing parameters.%s.", parameter);
                end
            end

            g0 = evp.parameters.g0;
            gd = evp.parameters.gd;
            if ~(isscalar(g0) && isreal(g0) && ~isnan(g0) && g0 ~= -Inf)
                error("IMExponentialStratificationSolution:UnsupportedEVP", ...
                    "parameters.g0 must be signed finite, zero, or Inf.");
            end
            if ~(isscalar(gd) && isreal(gd) && ~isnan(gd) && gd ~= -Inf)
                error("IMExponentialStratificationSolution:UnsupportedEVP", ...
                    "parameters.gd must be signed finite, zero, or Inf.");
            end
            surfaceConvention = string(evp.parameters.surfaceBoundary);
            if ~isscalar(surfaceConvention) || ~ismember(surfaceConvention, ["freeSurface", "rigidLid"])
                error("IMExponentialStratificationSolution:UnsupportedEVP", ...
                    "parameters.surfaceBoundary must be ""freeSurface"" or ""rigidLid"".");
            end

            expectedSurface = IMExponentialGeostrophicAPVCatalog.surfaceCondition( ...
                evp.g, g0, surfaceConvention);
            expectedBottom = IMExponentialGeostrophicAPVCatalog.bottomCondition(gd);
            if ~IMExponentialGeostrophicAPVCatalog.equivalentBoundary(evp.surfaceBoundary, expectedSurface) ...
                    || ~IMExponentialGeostrophicAPVCatalog.equivalentBoundary(evp.bottomBoundary, expectedBottom)
                error("IMExponentialStratificationSolution:UnsupportedBoundary", ...
                    "The APV boundary coefficients do not match the declared g0, gd, and surface convention.");
            end

            z = linspace(zDomain(1), zDomain(2), 5).';
            expectedN2 = N0*N0*exp(2*z/b);
            actualN2 = evp.N2(z);
            tolerance = 1e4*eps(max([1; abs(expectedN2(:)); abs(actualN2(:))]));
            if max(abs(actualN2(:) - expectedN2(:))) > tolerance
                error("IMExponentialStratificationSolution:UnsupportedStratification", ...
                    "The APV EVP N2 function does not match this exponential analytical solution.");
            end

            context = IMExponentialGeostrophicAPVCatalog.context(evp);
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
            if max(abs(p(:) - 1./expectedN2(:))) > coefficientTolerance ...
                    || max(abs(q(:))) > coefficientTolerance ...
                    || max(abs(r(:) - 1/evp.g)) > coefficientTolerance
                error("IMExponentialStratificationSolution:UnsupportedEVP", ...
                    "The APV EVP coefficients must be p=1/N2, q=0, and r=1/g.");
            end
        end

        function [modeData, h, modeNumber, metadata, diagnostics] = solve(evp, N0, b, zDomain, nModes)
            IMExponentialGeostrophicAPVCatalog.validateEVP(evp, N0, b, zDomain);
            endpoint = IMExponentialGeostrophicAPVCatalog.endpointData( ...
                N0, b, zDomain, evp.g, evp.parameters.g0, evp.parameters.gd, ...
                string(evp.parameters.surfaceBoundary));
            [negativeCount, zeroPresent, endpointMatrix] = ...
                IMExponentialGeostrophicAPVCatalog.endpointInertia(endpoint);

            negativeRoots = IMExponentialGeostrophicAPVCatalog.negativeRoots(endpoint, negativeCount);
            zeroCoefficients = zeros(1,2);
            if zeroPresent
                zeroCoefficients = IMExponentialGeostrophicAPVCatalog.zeroCoefficients(endpoint);
            end
            positiveCount = max(0, nModes - negativeCount - double(zeroPresent));
            positiveRoots = IMExponentialGeostrophicAPVCatalog.positiveRoots( ...
                endpoint, positiveCount, zeroPresent);

            roots = [negativeRoots(:).' NaN(1,double(zeroPresent)) positiveRoots(:).'];
            modeKinds = [repmat("negative",1,length(negativeRoots)) ...
                repmat("zero",1,double(zeroPresent)) repmat("positive",1,length(positiveRoots))];
            h = [-(b*N0./negativeRoots(:).').^2/evp.g ...
                Inf(1,double(zeroPresent)) (b*N0./positiveRoots(:).').^2/evp.g];
            coefficients = zeros(length(roots),2);
            sReference = NaN(1,length(roots));
            for iMode = 1:length(roots)
                if modeKinds(iMode) == "zero"
                    coefficients(iMode,:) = zeroCoefficients;
                else
                    [coefficients(iMode,:), sReference(iMode)] = ...
                        IMExponentialGeostrophicAPVCatalog.branchCoefficients( ...
                        endpoint, roots(iMode), modeKinds(iMode));
                end
            end

            eigenvalues = 1 ./ h;
            [~, order] = sort(eigenvalues, "ascend");
            roots = roots(order);
            modeKinds = modeKinds(order);
            h = h(order);
            coefficients = coefficients(order,:);
            sReference = sReference(order);
            keep = 1:min(nModes,length(h));
            roots = roots(keep);
            modeKinds = modeKinds(keep);
            h = h(keep);
            coefficients = coefficients(keep,:);
            sReference = sReference(keep);

            retainedNegativeCount = nnz(modeKinds == "negative");
            modeNumber = zeros(1,length(h));
            modeNumber(1:retainedNegativeCount) = -1:-1:-retainedNegativeCount;
            nextMode = retainedNegativeCount + 1;
            if nextMode <= length(h) && modeKinds(nextMode) == "zero"
                nextMode = nextMode + 1;
            end
            modeNumber(nextMode:end) = 1:(length(h) - nextMode + 1);

            modeData = struct();
            modeData.roots = roots;
            modeData.modeKinds = modeKinds;
            modeData.coefficients = coefficients;
            modeData.sReference = sReference;
            modeData.endpoint = endpoint;
            modeData.signFactors = ones(1,length(h));
            zOrientation = linspace(zDomain(1), zDomain(2), 1025).';
            [FOrientation, ~, ~] = IMExponentialGeostrophicAPVCatalog.rawValues( ...
                modeData, N0, b, evp.g, zOrientation, false);
            for iMode = 1:length(h)
                [~, index] = max(abs(FOrientation(:,iMode)));
                if FOrientation(index,iMode) < 0
                    modeData.signFactors(iMode) = -1;
                end
            end

            rootResiduals = NaN(1,length(roots));
            for iMode = 1:length(roots)
                if modeKinds(iMode) ~= "zero"
                    rootResiduals(iMode) = IMExponentialGeostrophicAPVCatalog.branchResidual( ...
                        endpoint, roots(iMode), modeKinds(iMode));
                end
            end
            metadata = struct( ...
                "roots", roots, ...
                "rootResiduals", rootResiduals, ...
                "modeKinds", modeKinds, ...
                "negativeModeCount", negativeCount, ...
                "zeroModePresent", zeroPresent, ...
                "endpointMatrix", endpointMatrix, ...
                "endpointData", endpoint, ...
                "g0", evp.parameters.g0, ...
                "gd", evp.parameters.gd, ...
                "surfaceBoundary", string(evp.parameters.surfaceBoundary));
            zeroModeStatus = "absent";
            if zeroPresent
                zeroModeStatus = "present";
            end
            diagnostics = struct( ...
                "negativeModeCount", negativeCount, ...
                "zeroModeStatus", zeroModeStatus, ...
                "zeroModeCount", double(zeroPresent), ...
                "reason", "Exact exponential generalized-energy APV endpoint inertia.");
        end

        function values = rawVariable(modeData, N0, b, g, variable, z)
            [F, G, ~] = IMExponentialGeostrophicAPVCatalog.rawValues( ...
                modeData, N0, b, g, z(:), true);
            if string(variable) == "F"
                values = F;
            elseif string(variable) == "G"
                values = G;
            else
                error("IMExponentialStratificationSolution:UnsupportedVariable", ...
                    "Variable must be ""F"" or ""G"".");
            end
        end

        function values = rawFz(modeData, N0, b, g, z)
            [~, ~, values] = IMExponentialGeostrophicAPVCatalog.rawValues( ...
                modeData, N0, b, g, z(:), true);
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

        function endpoint = endpointData(N0, b, zDomain, g, g0, gd, surfaceConvention)
            depth = diff(zDomain);
            endpoint.depth = depth;
            endpoint.N0 = N0;
            endpoint.b = b;
            endpoint.expMinusDOverB = exp(-depth/b);
            endpoint.IN = N0*N0*b*(1 - exp(-2*depth/b))/2;
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
                matrix = [1 + endpoint.IN*endpoint.betaSurface, -1; ...
                    -1, 1 + endpoint.IN*endpoint.betaBottom];
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
            if endpoint.bottomKind == "dirichlet"
                bottomRow = [1 0];
            else
                bottomRow = [-endpoint.betaBottom 1];
            end
            if endpoint.surfaceKind == "dirichlet"
                surfaceRow = [1 endpoint.IN];
            else
                surfaceRow = [endpoint.betaSurface 1 + endpoint.betaSurface*endpoint.IN];
            end
            residual = det([bottomRow; surfaceRow]);
            tolerance = 1e3*eps(max([1 abs(bottomRow) abs(surfaceRow)]));
            if abs(residual) > tolerance
                error("IMExponentialStratificationSolution:ZeroModeMismatch", ...
                    "The endpoint inertia predicted a zero mode but the exact zero determinant is nonzero.");
            end
            coefficients = [bottomRow(2) -bottomRow(1)];
            coefficients = coefficients/norm(coefficients);
        end

        function roots = negativeRoots(endpoint, count)
            if count == 0
                roots = zeros(1,0);
                return;
            end
            upperBound = 10;
            roots = zeros(1,0);
            for iteration = 1:6
                x = logspace(-8, log10(upperBound), 12000).';
                residual = @(value) IMExponentialGeostrophicAPVCatalog.branchResidual( ...
                    endpoint, value, "negative");
                roots = IMExponentialGeostrophicAPVCatalog.scanRoots(residual, x);
                if length(roots) >= count
                    break;
                end
                upperBound = upperBound*10;
            end
            if length(roots) ~= count
                error("IMExponentialStratificationSolution:RootSearchFailed", ...
                    "Expected %d negative APV roots but found %d.", count, length(roots));
            end
            roots = reshape(roots,1,[]);
        end

        function roots = positiveRoots(endpoint, count, zeroPresent)
            if count == 0
                roots = zeros(1,0);
                return;
            end
            phaseSpan = max(1 - endpoint.expMinusDOverB, 1e-3);
            upperBound = max(20, (count + 4)*pi/phaseSpan);
            roots = zeros(1,0);
            for iteration = 1:8
                nLinear = max(8192, 512*(count + 4));
                x = unique([logspace(-8,0,1024) linspace(1,upperBound,nLinear)]).';
                residual = @(value) IMExponentialGeostrophicAPVCatalog.branchResidual( ...
                    endpoint, value, "positive");
                roots = IMExponentialGeostrophicAPVCatalog.scanRoots(residual, x);
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
                error("IMExponentialStratificationSolution:RootSearchFailed", ...
                    "Expected %d positive APV roots but found %d.", count, length(roots));
            end
            roots = reshape(roots,1,[]);
        end

        function value = branchResidual(endpoint, s0, branch)
            if branch == "negative"
                % Cancel the known exponential factors before forming the
                % characteristic determinant. Row normalization otherwise
                % turns a localized root into an unresolvable sign jump.
                sd=endpoint.expMinusDOverB*s0;
                x=[sd;s0]; N2=endpoint.N0^2*exp(2*[-endpoint.depth;0]/endpoint.b);
                values=[besseli(1,x,1),besselk(1,x,1)];
                flux=x./(endpoint.b*N2).*[besseli(0,x,1),-besselk(0,x,1)];
                beta=[-endpoint.betaBottom;endpoint.betaSurface];
                rows=zeros(2); rowScale=zeros(2,1);
                kinds=[string(endpoint.bottomKind),string(endpoint.surfaceKind)];
                for row=1:2
                    if kinds(row)=="dirichlet"
                        rows(row,:)=values(row,:); rowScale(row)=norm(values(row,:));
                    else
                        rows(row,:)=flux(row,:)+beta(row)*values(row,:);
                        rowScale(row)=norm(abs(flux(row,:))+abs(beta(row)*values(row,:)));
                    end
                end
                rows=rows./rowScale;
                value=exp(-2*(s0-sd))*rows(1,1)*rows(2,2)-rows(1,2)*rows(2,1);
                return;
            end
            [bottomRow, surfaceRow] = IMExponentialGeostrophicAPVCatalog.branchRows( ...
                endpoint, s0, branch);
            bottomNorm = norm(bottomRow);
            surfaceNorm = norm(surfaceRow);
            if bottomNorm == 0 || surfaceNorm == 0 || ~isfinite(bottomNorm) || ~isfinite(surfaceNorm)
                value = NaN;
                return;
            end
            value = det([bottomRow/bottomNorm; surfaceRow/surfaceNorm]);
        end

        function [coefficients, sReference] = branchCoefficients(endpoint, s0, branch)
            [bottomRow, ~, sReference] = IMExponentialGeostrophicAPVCatalog.branchRows( ...
                endpoint, s0, branch);
            coefficients = [bottomRow(2) -bottomRow(1)];
            coefficients = coefficients/norm(coefficients);
        end

        function [bottomRow, surfaceRow, sReference] = branchRows(endpoint, s0, branch)
            sd = endpoint.expMinusDOverB*s0;
            sReference = NaN;
            if branch == "positive"
                [bottomValue, bottomFlux] = IMExponentialGeostrophicAPVCatalog.positiveRows( ...
                    sd, -endpoint.depth, endpoint);
                [surfaceValue, surfaceFlux] = IMExponentialGeostrophicAPVCatalog.positiveRows( ...
                    s0, 0, endpoint);
            else
                sReference = (s0 + sd)/2;
                [bottomValue, bottomFlux] = IMExponentialGeostrophicAPVCatalog.negativeRows( ...
                    sd, -endpoint.depth, endpoint, sReference);
                [surfaceValue, surfaceFlux] = IMExponentialGeostrophicAPVCatalog.negativeRows( ...
                    s0, 0, endpoint, sReference);
            end
            if endpoint.bottomKind == "dirichlet"
                bottomRow = bottomValue;
            else
                bottomRow = bottomFlux - endpoint.betaBottom*bottomValue;
            end
            if endpoint.surfaceKind == "dirichlet"
                surfaceRow = surfaceValue;
            else
                surfaceRow = surfaceFlux + endpoint.betaSurface*surfaceValue;
            end
        end

        function [valueRow, fluxRow] = positiveRows(s, z, endpoint)
            N2 = endpoint.N0*endpoint.N0*exp(2*z/endpoint.b);
            valueRow = s*[besselj(1,s) bessely(1,s)];
            fluxRow = (s*s/(endpoint.b*N2))*[besselj(0,s) bessely(0,s)];
        end

        function [valueRow, fluxRow] = negativeRows(s, z, endpoint, sReference)
            N2 = endpoint.N0*endpoint.N0*exp(2*z/endpoint.b);
            I1 = exp(s - sReference)*besseli(1,s,1);
            K1 = exp(sReference - s)*besselk(1,s,1);
            I0 = exp(s - sReference)*besseli(0,s,1);
            K0 = exp(sReference - s)*besselk(0,s,1);
            valueRow = s*[I1 K1];
            fluxRow = (s*s/(endpoint.b*N2))*[I0 -K0];
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
            roots = roots(1:rootCount);
            roots = IMExponentialGeostrophicAPVCatalog.deduplicateRoots(roots);
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

        function [F, G, Fz] = rawValues(modeData, N0, b, g, z, applySigns)
            z = z(:);
            N2 = N0*N0*exp(2*z/b);
            depth = modeData.endpoint.depth;
            Q = N0*N0*b*(exp(2*z/b) - exp(-2*depth/b))/2;
            F = zeros(length(z),length(modeData.roots));
            G = zeros(size(F));
            Fz = zeros(size(F));
            for iMode = 1:length(modeData.roots)
                coefficients = modeData.coefficients(iMode,:);
                if modeData.modeKinds(iMode) == "zero"
                    F(:,iMode) = coefficients(1) + coefficients(2)*Q;
                    Fz(:,iMode) = coefficients(2)*N2;
                else
                    s = modeData.roots(iMode)*exp(z/b);
                    if modeData.modeKinds(iMode) == "positive"
                        F(:,iMode) = s.*(coefficients(1)*besselj(1,s) + coefficients(2)*bessely(1,s));
                        Fz(:,iMode) = (s.*s/b).*(coefficients(1)*besselj(0,s) + coefficients(2)*bessely(0,s));
                    else
                        sReference = modeData.sReference(iMode);
                        I1 = exp(s - sReference).*besseli(1,s,1);
                        K1 = exp(sReference - s).*besselk(1,s,1);
                        I0 = exp(s - sReference).*besseli(0,s,1);
                        K0 = exp(sReference - s).*besselk(0,s,1);
                        F(:,iMode) = s.*(coefficients(1)*I1 + coefficients(2)*K1);
                        Fz(:,iMode) = (s.*s/b).*(coefficients(1)*I0 - coefficients(2)*K0);
                    end
                end
                G(:,iMode) = -(g./N2).*Fz(:,iMode);
            end
            if applySigns
                F = F .* modeData.signFactors;
                G = G .* modeData.signFactors;
                Fz = Fz .* modeData.signFactors;
            end
        end

        function boundary = surfaceCondition(g, g0, surfaceConvention)
            if g0 == 0
                boundary = IMBoundaryCondition.dirichlet();
            elseif surfaceConvention == "freeSurface"
                if isinf(g0)
                    beta = 1/g;
                else
                    beta = 1/g + 1/g0;
                end
                boundary = IMBoundaryCondition(a=-beta, b=1);
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
            tf = isfinite(scale) && scale ~= 0 ...
                && max(abs(actualCoefficients - scale*expectedCoefficients)) <= tolerance;
        end
    end
end
