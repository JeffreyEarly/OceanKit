classdef (Hidden, Sealed) IMGeostrophicZeroAPVFormTools
    % Share canonical endpoint forms and their scientific coordinate changes.
    methods (Static)
        function [responseMetric,physicalEnergy,surfaceBuoyancy,bottomBuoyancy] = formMatrices(FSurface,FBottom,surfaceResponse,bottomResponse,selectedResponse,endpoints,k,f0,g)
            nEndpoints = numel(endpoints);
            nK = numel(k);
            responseMetric = zeros(nEndpoints,nEndpoints,nK);
            physicalEnergy = zeros(nEndpoints,nEndpoints,nK);
            surfaceBuoyancy = zeros(nEndpoints,nEndpoints,nK);
            bottomBuoyancy = zeros(nEndpoints,nEndpoints,nK);

            for iK = 1:nK
                response = selectedResponse(:,:,iK);
                responseMetric(:,:,iK) = IMGeostrophicZeroAPVFormTools.symmetrize(response.'*response);
                energyScale = f0^2/(2*g*k(iK)^4);
                boundaryForm = -FSurface(:,:,iK).'*surfaceResponse(:,:,iK) + FBottom(:,:,iK).'*bottomResponse(:,:,iK);
                physicalEnergy(:,:,iK) = energyScale*IMGeostrophicZeroAPVFormTools.symmetrize(boundaryForm);
                buoyancyScale = f0^2/(2*g^2*k(iK)^4);
                if ismember("surface",endpoints)
                    surfaceBuoyancy(:,:,iK) = buoyancyScale*IMGeostrophicZeroAPVFormTools.symmetrize(surfaceResponse(:,:,iK).'*surfaceResponse(:,:,iK));
                end
                if ismember("bottom",endpoints)
                    bottomBuoyancy(:,:,iK) = buoyancyScale*IMGeostrophicZeroAPVFormTools.symmetrize(bottomResponse(:,:,iK).'*bottomResponse(:,:,iK));
                end
            end
        end
        function matrix = generalizedEnergy(energy,surface,bottom,endpoints,g0,gd,errorPrefix)
            for item = {g0,gd}
                coefficient = item{1};
                if ~isscalar(coefficient) || ~isreal(coefficient) || isnan(coefficient) || coefficient == -Inf
                    error(errorPrefix+":InvalidEnergyCoefficient","Use finite signed coefficients; positive infinity is allowed only for an absent endpoint.");
                end
            end
            if (isinf(g0) && ismember("surface",endpoints)) || (isinf(gd) && ismember("bottom",endpoints))
                error(errorPrefix+":InactiveEndpointRequired","An infinite coefficient cannot remove a retained endpoint. Choose endpoints explicitly in the original problem.");
            end
            matrix = energy;
            if ismember("surface",endpoints), matrix = matrix+g0*surface; end
            if ismember("bottom",endpoints), matrix = matrix+gd*bottom; end
            if any(~isfinite(matrix),"all")
                error(errorPrefix+":UnrepresentableQuadraticForm","Generalized energy must remain finite in these coefficient coordinates.");
            end
            matrix = IMGeostrophicZeroAPVFormTools.symmetrizePages(matrix);
        end

        function [map,depth,evidence] = orthogonalizeEnergy(matrix,k,errorPrefix)
            n = size(matrix,1); pages = numel(k);
            map = zeros(n,n,pages); depth = zeros(n,pages);
            residual = zeros(n,pages); orthogonality = zeros(1,pages);
            for page = 1:pages
                if n == 0, continue; end
                H = matrix(:,:,page);
                scale = max(abs(H),[],"all");
                if scale == 0
                    C = eye(n); lambda = zeros(n,1);
                else
                    unit = H/scale;
                    if any(H~=0 & unit==0,"all")
                        error(errorPrefix+":UnrepresentableGeneralizedEnergy","Scaling generalized energy lost a nonzero matrix entry.");
                    end
                    if norm(unit-unit.',"fro") > 1e-12*norm(unit,"fro")
                        error(errorPrefix+":NonSymmetricGeneralizedEnergy","The canonical generalized-energy form must be symmetric.");
                    end
                    [C,lambda] = eig(IMGeostrophicZeroAPVFormTools.symmetrize(unit),"vector");
                    [lambda,order] = sort(lambda,"descend"); C = C(:,order);
                end
                for j = 1:n
                    [~,index] = max(abs(C(:,j)));
                    if C(index,j)<0, C(:,j) = -C(:,j); end
                    if scale>0
                        difference = (H/scale)*C(:,j)-lambda(j)*C(:,j);
                        denominator = (norm(H/scale,"fro")+abs(lambda(j)))*norm(C(:,j));
                        if denominator>0, residual(j,page) = norm(difference)/denominator; end
                    end
                end
                orthogonality(page) = norm(C.'*C-eye(n),"fro")/max(1,sqrt(n));
                if max(residual(:,page))>1e-10 || orthogonality(page)>1e-10
                    error(errorPrefix+":UnqualifiedOrthogonality","Energy residuals and orthonormal endpoint geometry must pass the relative qualification checks.");
                end
                map(:,:,page) = C;
                if scale>0
                    physical = IMGeostrophicZeroAPVFormTools.restoreScale(lambda,scale);
                    depth(:,page) = (2*k(page)^2)*physical;
                    if any(~isfinite(depth(:,page))) || any(physical~=0 & depth(:,page)==0)
                        error(errorPrefix+":UnrepresentableBoundaryDepth","Signed boundary depths must be representable in the supplied units.");
                    end
                end
            end
            evidence = struct(relativeEnergyResidual=residual,endpointOrthogonalityError=orthogonality,method="real symmetric canonical generalized-energy eigenproblem");
        end

        function result = transformForm(canonical,map)
            n = size(map,2); pages = size(map,3);
            result = zeros(n,n,pages);
            for page = 1:pages
                result(:,:,page) = IMGeostrophicZeroAPVFormTools.symmetrize(map(:,:,page).'*canonical(:,:,page)*map(:,:,page));
            end
        end

        function result = applyMap(canonical,map)
            result = zeros(size(canonical,1),size(map,2),size(map,3));
            for page = 1:size(map,3), result(:,:,page) = canonical(:,:,page)*map(:,:,page); end
        end

        function factors = normalizationFactors(basis,name,page)
            n = basis.nColumns;
            switch name
                case "endpointResponse"
                    directions = basis.endpointUnscaledResponseMatrix();
                    factors = reshape(sqrt(sum(abs(directions(:,:,page)).^2,1)),n,1);
                    return

            end
            [z,weights] = basis.endpointNormalizationQuadrature();
            if endsWith(name,"F"), variable="F"; else, variable="G"; end
            values = basis.endpointUnscaledProfiles(variable,z,page);
            if startsWith(name,"rms")
                factors = sqrt(sum(weights.*abs(values).^2,1)/diff(basis.zDomain)).';
            else
                points = [z;basis.zDomain.'];
                values = basis.endpointUnscaledProfiles(variable,points,page);
                evaluator = @(zz) basis.endpointUnscaledProfiles(variable,zz,page);
                factors = IMInternalModesNormalizationTools.maximumAmplitudes(evaluator,points,values).';
            end
        end

        function matrix = symmetrize(matrix)
            if isequal(matrix,matrix.'), return; end
            diagonal = diag(matrix);
            matrix = .5*matrix+.5*matrix.';
            matrix(1:size(matrix,1)+1:end) = diagonal;
        end

        function pages = symmetrizePages(pages)
            for page = 1:size(pages,3), pages(:,:,page) = IMGeostrophicZeroAPVFormTools.symmetrize(pages(:,:,page)); end
        end
    end
    methods (Static, Access = private)
        function values = restoreScale(values,scale)
            [fraction,exponent] = log2(scale);
            [valueFraction,valueExponent] = log2(abs(values));
            restored = sign(values).*pow2(valueFraction*fraction,valueExponent+exponent);
            if any(~isfinite(restored)) || any(values~=0 & restored==0)
                error("IMGeostrophicZeroAPVFormTools:UnrepresentableEnergySpectrum","Generalized-energy eigenvalues must be representable in the supplied units.");
            end
            values = restored;
        end
    end
end
