classdef (Hidden, Sealed) IMGeostrophicZeroAPVFormTools
    % Share representation-independent zero-APV form and rotation operations.

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

        function pages = expandMatrixPages(matrix,nEndpoints,nK,name,errorPrefix)
            if size(matrix,1) ~= nEndpoints || size(matrix,2) ~= nEndpoints
                error(errorPrefix + ":InvalidPencilShape", "%s must have nEndpoints rows and columns.",name);
            end
            if size(matrix,3) == 1
                pages = repmat(matrix,1,1,nK);
            elseif size(matrix,3) == nK
                pages = matrix;
            else
                error(errorPrefix + ":InvalidPencilShape", "%s must have one page or nK pages.",name);
            end
            for iK = 1:nK
                scale = max(1,norm(pages(:,:,iK),"fro"));
                if norm(pages(:,:,iK)-pages(:,:,iK).',"fro") > 1e-12*scale
                    error(errorPrefix + ":NonSymmetricPencil", "%s must be symmetric on every wavenumber page.",name);
                end
                pages(:,:,iK) = IMGeostrophicZeroAPVFormTools.symmetrize(pages(:,:,iK));
            end
        end

        function [vectors,eigenvalues,rightSignatures,residuals] = solveSymmetricPencil(leftMatrix,rightMatrix,nEndpoints,nK,ordering,errorPrefix)
            leftMatrix = IMGeostrophicZeroAPVFormTools.expandMatrixPages(leftMatrix,nEndpoints,nK,"leftMatrix",errorPrefix);
            rightMatrix = IMGeostrophicZeroAPVFormTools.expandMatrixPages(rightMatrix,nEndpoints,nK,"rightMatrix",errorPrefix);
            vectors = zeros(nEndpoints,nEndpoints,nK);
            eigenvalues = zeros(nEndpoints,nK);
            rightSignatures = zeros(nEndpoints,nK);
            relativeResidual = zeros(nEndpoints,nK);
            rightOrthogonalityError = zeros(1,nK);
            rotationConditionNumber = zeros(1,nK);

            for iK = 1:nK
                left = leftMatrix(:,:,iK);
                right = rightMatrix(:,:,iK);
                [C,Lambda] = eig(left,right);
                lambda = diag(Lambda);
                imaginaryScale = max(1,max(abs(real(lambda))));
                if any(~isfinite(lambda)) || any(abs(imag(lambda)) > 1e-10*imaginaryScale)
                    error(errorPrefix + ":InvalidPencilSpectrum", "The matrix pencil must have a complete finite real spectrum on every wavenumber page.");
                end
                C = real(C);
                lambda = real(lambda);
                for iMode = 1:nEndpoints
                    rightNorm = C(:,iMode).'*right*C(:,iMode);
                    normScale = max(1,norm(right,"fro")*norm(C(:,iMode))^2);
                    if abs(rightNorm) <= 1e-12*normScale
                        error(errorPrefix + ":NullRightMetric", "A pencil eigenvector has a numerically zero right-metric norm.");
                    end
                    C(:,iMode) = C(:,iMode)/sqrt(abs(rightNorm));
                end

                switch string(ordering)
                    case "nonzeroFirst"
                        zeroTolerance = 1e-11*max(1,max(abs(lambda)));
                        isZero = abs(lambda) <= zeroTolerance;
                        [~,sortIndex] = sortrows([isZero(:),-abs(lambda(:)),-lambda(:)],[1 2 3]);
                    otherwise
                        [~,sortIndex] = sort(lambda,"descend");
                end
                lambda = lambda(sortIndex);
                C = C(:,sortIndex);

                for iMode = 1:nEndpoints
                    [~,referenceIndex] = max(abs(C(:,iMode)));
                    if C(referenceIndex,iMode) < 0
                        C(:,iMode) = -C(:,iMode);
                    end
                    rightSignatures(iMode,iK) = sign(C(:,iMode).'*right*C(:,iMode));
                    residual = left*C(:,iMode)-lambda(iMode)*right*C(:,iMode);
                    residualScale = (norm(left,"fro")+abs(lambda(iMode))*norm(right,"fro"))*norm(C(:,iMode));
                    relativeResidual(iMode,iK) = norm(residual)/max(residualScale,eps);
                end

                targetRightGram = diag(rightSignatures(:,iK));
                rightOrthogonalityError(iK) = norm(C.'*right*C-targetRightGram,"fro")/max(1,norm(targetRightGram,"fro"));
                rotationConditionNumber(iK) = cond(C);
                vectors(:,:,iK) = C;
                eigenvalues(:,iK) = lambda;
            end
            residuals = struct(relativePencilResidual=relativeResidual,rightOrthogonalityError=rightOrthogonalityError,rotationConditionNumber=rotationConditionNumber);
        end

        function tf = pagesEqual(left,right)
            scale = max([1,norm(left(:)),norm(right(:))]);
            tf = isequal(size(left),size(right)) && norm(left(:)-right(:)) <= 1e-12*scale;
        end

        function matrix = symmetrize(matrix)
            matrix = 0.5*(matrix+matrix.');
        end

        function pages = symmetrizePages(pages)
            for iPage = 1:size(pages,3)
                pages(:,:,iPage) = IMGeostrophicZeroAPVFormTools.symmetrize(pages(:,:,iPage));
            end
        end

        function validateRotationPages(rotation,nEndpoints,nK,errorPrefix)
            if ~isequal(size(rotation),[nEndpoints nEndpoints nK]) && ~(nK == 1 && isequal(size(rotation),[nEndpoints nEndpoints]))
                error(errorPrefix + ":InvalidRotationShape", "A rotation must have nEndpoints rows and columns and one page per wavenumber.");
            end
            for iK = 1:nK
                if rcond(rotation(:,:,iK)) <= sqrt(eps)
                    error(errorPrefix + ":SingularCoordinateRotation", "The coordinate rotation is singular or ill-conditioned on page %d.",iK);
                end
            end
        end
    end
end
