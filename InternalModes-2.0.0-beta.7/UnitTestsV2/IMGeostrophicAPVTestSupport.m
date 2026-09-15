classdef IMGeostrophicAPVTestSupport
    % Shared reference construction for generalized-energy APV tests.

    methods (Static)
        function evp = manualEVP(options)
            arguments
                options.N2 (1,1) function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.g0 (1,1) double {mustBeReal}
                options.gd (1,1) double {mustBeReal}
                options.surfaceBoundary {mustBeTextScalar, mustBeMember(options.surfaceBoundary, ["freeSurface", "rigidLid"])} = "freeSurface"
            end

            surfaceBoundary = IMGeostrophicAPVTestSupport.surfaceCondition( ...
                options.g, options.g0, string(options.surfaceBoundary));
            bottomBoundary = IMGeostrophicAPVTestSupport.bottomCondition(options.gd);
            parameters = struct( ...
                "g0", options.g0, ...
                "gd", options.gd, ...
                "surfaceBoundary", string(options.surfaceBoundary));
            evp = IMInternalModes( ...
                name="geostrophicAPVModes", ...
                formulation="F", ...
                modeFamily="hydrostatic", ...
                N2=options.N2, ...
                zDomain=options.zDomain, ...
                p=@(z,ctx) 1 ./ ctx.N2(z), ...
                q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) ones(size(z))/ctx.g, ...
                surfaceBoundary=surfaceBoundary, ...
                bottomBoundary=bottomBoundary, ...
                g=options.g, ...
                parameters=parameters);
        end

        function data = endpointData(options)
            arguments
                options.N0 (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
                options.b (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.g0 (1,1) double {mustBeReal}
                options.gd (1,1) double {mustBeReal}
                options.surfaceBoundary {mustBeTextScalar, mustBeMember(options.surfaceBoundary, ["freeSurface", "rigidLid"])} = "freeSurface"
            end

            depth = diff(sort(options.zDomain));
            data.IN = options.N0*options.N0*options.b*(1 - exp(-2*depth/options.b))/2;
            data.surfaceKind = "active";
            data.bottomKind = "active";
            if options.g0 == 0
                data.surfaceKind = "dirichlet";
                data.betaSurface = NaN;
            elseif string(options.surfaceBoundary) == "freeSurface"
                if isinf(options.g0)
                    data.betaSurface = 1/options.g;
                else
                    data.betaSurface = 1/options.g + 1/options.g0;
                end
            elseif isinf(options.g0)
                data.betaSurface = 0;
            else
                data.betaSurface = 1/options.g0;
            end

            if options.gd == 0
                data.bottomKind = "dirichlet";
                data.betaBottom = NaN;
            elseif isinf(options.gd)
                data.betaBottom = 0;
            else
                data.betaBottom = 1/options.gd;
            end
        end

        function parameters = parametersFromStrengths(options)
            arguments
                options.N0 (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
                options.b (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.surfaceStrength (1,1) double {mustBeReal, mustBeFinite}
                options.bottomStrength (1,1) double {mustBeReal, mustBeFinite}
                options.surfaceBoundary {mustBeTextScalar, mustBeMember(options.surfaceBoundary, ["freeSurface", "rigidLid"])} = "freeSurface"
            end

            depth = diff(sort(options.zDomain));
            IN = options.N0*options.N0*options.b*(1 - exp(-2*depth/options.b))/2;
            betaSurface = options.surfaceStrength/IN;
            betaBottom = options.bottomStrength/IN;
            if string(options.surfaceBoundary) == "freeSurface"
                denominator = betaSurface - 1/options.g;
            else
                denominator = betaSurface;
            end
            if denominator == 0
                g0 = Inf;
            else
                g0 = 1/denominator;
            end
            if betaBottom == 0
                gd = Inf;
            else
                gd = 1/betaBottom;
            end
            parameters = struct( ...
                "g0", g0, ...
                "gd", gd, ...
                "surfaceBoundary", string(options.surfaceBoundary), ...
                "IN", IN, ...
                "betaSurface", betaSurface, ...
                "betaBottom", betaBottom);
        end

        function [F, G, factors] = depthNormalize(z, F, G)
            depth = max(z) - min(z);
            factors = sqrt(trapz(z, F.*F, 1)/depth);
            F = F ./ factors;
            G = G ./ factors;
            [F, G] = IMGeostrophicAPVTestSupport.orientColumns(F, G);
        end

        function [F, G] = orientColumns(F, G)
            for iMode = 1:size(F,2)
                [~, index] = max(abs(F(:,iMode)));
                if F(index,iMode) < 0
                    F(:,iMode) = -F(:,iMode);
                    G(:,iMode) = -G(:,iMode);
                end
            end
        end

        function errorValue = subspaceError(z, leftModes, rightModes)
            z = z(:);
            weights = zeros(size(z));
            weights(1) = (z(2) - z(1))/2;
            weights(end) = (z(end) - z(end-1))/2;
            weights(2:end-1) = (z(3:end) - z(1:end-2))/2;
            weightedLeft = sqrt(weights).*leftModes;
            weightedRight = sqrt(weights).*rightModes;
            [leftQ, ~] = qr(weightedLeft, 0);
            [rightQ, ~] = qr(weightedRight, 0);
            singularValues = svd(leftQ.'*rightQ);
            errorValue = sqrt(max(0, 1 - min(singularValues)^2));
        end
    end

    methods (Static, Access = private)
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
    end
end
