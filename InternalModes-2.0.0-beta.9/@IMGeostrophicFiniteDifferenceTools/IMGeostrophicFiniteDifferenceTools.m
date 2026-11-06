classdef IMGeostrophicFiniteDifferenceTools
    % Assemble conservative staggered geostrophic finite-difference pairs.
    % F uses the supplied primal grid. G uses its cell midpoints and both
    % endpoints, retaining the discrete integration-by-parts identity.

    methods (Static)
        function [A,B,samples] = assemble(evp,solver)
            data = IMGeostrophicFiniteDifferenceTools.configuration(evp,solver);
            [A,B] = IMGeostrophicFiniteDifferenceTools.inverseDepthPencil(evp,data);
            samples = struct("N2",IMGeostrophicFiniteDifferenceTools.stratification(evp.N2,data.zF));
            if evp.name == "geostrophicGeneralizedPotentialEnstrophyModes"
                rows = data.n+(1:data.n);
                A(rows,:) = (evp.f0^2/evp.g)*A(rows,:)+evp.parameters.k^2*B(rows,:);
                A(data.n,:) = 0;
                A(data.n,1) = -1;
                A(data.n,data.n+1) = 1;
                B(data.n,:) = 0;
                B(data.n,1) = -evp.g/evp.parameters.alpha0;
                A(end,:) = 0;
                A(end,end) = 1;
                B(end,:) = 0;
                B(end,data.n) = evp.g/evp.parameters.alphaD;
            end
        end

        function [A,B,map] = preparePencil(evp,solver)
            data = IMGeostrophicFiniteDifferenceTools.configuration(evp,solver);
            if data.IN == 0
                [r,V] = IMGeostrophicFiniteDifferenceTools.fullyNeutralMode(evp,data);
                rScale = 1/diff(evp.zDomain);
                A = r/rScale;
                B = 1;
                sizes = 1;
                columnScale = 1;
                expand = @(coefficients) V*coefficients;
            else
                [A,B] = IMGeostrophicFiniteDifferenceTools.inverseDepthPencil(evp,data);
                rScale = evp.g/(data.IN*diff(evp.zDomain));
                B = rScale*B;
                sizes = [data.n data.n+1];
                columnScale = [data.IN/evp.g*ones(1,data.n),ones(1,data.n+1)];
                expand = @(coefficients) coefficients;
            end
            map = struct("nativeFieldSizes",sizes,"hasInverseDepth",true, ...
                "nativeColumnScale",columnScale,"restoreEigenvalues",@(lambda) rScale*lambda, ...
                "expandEigenvectors",expand,"publicEigenvalues",@(V,r) IMGeostrophicFiniteDifferenceTools.publicEigenvalues(evp,solver,V,r));
        end

        function [V,lambda] = qualify(evp,solver,A,B,V,lambda)
            data = IMGeostrophicFiniteDifferenceTools.configuration(evp,solver);
            finite = isfinite(real(lambda)) & isfinite(imag(lambda));
            V = V(:,finite);
            lambda = lambda(finite);
            if isempty(lambda), return; end
            realPair = abs(imag(lambda))<=1e-10*max(abs(real(lambda)),1);
            for j = 1:numel(lambda)
                [~,index] = max(abs(V(:,j)));
                V(:,j) = V(:,j)*exp(-1i*angle(V(index,j)));
            end
            realPair = realPair & (vecnorm(imag(V),2,1).'<=1e-8*vecnorm(real(V),2,1).');
            V = real(V(:,realPair));
            lambda = real(lambda(realPair));
            if isempty(lambda), return; end
            scale = ones(size(V,1),1);
            if data.IN>0, scale(1:data.n) = data.IN/evp.g; end
            balancedV = V./scale;
            metric = vecnorm(B*V,2,1)./max(norm(B.*scale.',1)*vecnorm(balancedV,2,1),realmin);
            valid = metric>1e3*eps;
            if data.IN>0
                valid = valid & vecnorm(V(1:data.n,:),2,1)/(data.IN/evp.g)>1e3*eps*vecnorm(balancedV,2,1);
            end
            V = V(:,valid);
            lambda = lambda(valid);
            diagnostics = evp.modeSelectionDiagnostics(solver,[]);
            if diagnostics.zeroInverseDepthCount>0 && ~isempty(lambda)
                [~,j] = min(abs(lambda));
                zeroResidual = norm(A*V(:,j))/max(norm(A.*scale.',1)*norm(V(:,j)./scale),realmin);
                if zeroResidual<1e-9 && abs(lambda(j))<1e-7, lambda(j) = 0; end
            end
            if data.IN>0
                [V,lambda] = IMGeostrophicFiniteDifferenceTools.refineSmallModes(evp,data,V,lambda);
            end
        end

        function mapped = publicEigenvalues(evp,solver,V,r)
            data = IMGeostrophicFiniteDifferenceTools.configuration(evp,solver);
            valid = IMGeostrophicFiniteDifferenceTools.qualifiedPairs(evp,data,V,r);
            if evp.name == "geostrophicAPVModes"
                mapped = struct("values",r,"qualified",valid);
                return;
            end
            F = V(1:data.n,:);
            G = V(data.n+1:end,:);
            volumeF = sum(data.WF.*F.^2,1);
            volumeG = sum(data.cellIntegral.*G(2:end-1,:).^2,1);
            numerator = evp.parameters.k^2*volumeF+(evp.f0^2/evp.g^2)*volumeG+(evp.f0^2/evp.g)*F(1,:).^2;
            denominator = volumeF+(evp.f0^2/evp.parameters.alpha0)*F(1,:).^2+(evp.f0^2/evp.parameters.alphaD)*F(end,:).^2;
            values = (numerator./denominator).';
            values(r==0) = evp.parameters.k^2;
            affine = evp.parameters.k^2+(evp.f0^2/evp.g)*r;
            affineScale = evp.parameters.k^2+(evp.f0^2/evp.g)*abs(r);
            [A,B] = IMGeostrophicFiniteDifferenceTools.assemble(evp,solver);
            residual = vecnorm(A*V-(B*V).*values.',2,1) ...
                ./max((norm(A,1)+abs(values.')*norm(B,1)).*vecnorm(V,2,1),realmin);
            valid = valid & isfinite(values) & values>0 & residual.'<=1e-7 & abs(values-affine)<=1e-7*affineScale;
            values(~valid) = max(affine(~valid),realmin);
            mapped = struct("values",values,"qualified",valid);
        end

        function basis = makeBasis(evp,solver,V,lambda,labels,diagnostics)
            data = IMGeostrophicFiniteDifferenceTools.configuration(evp,solver);
            if isfield(diagnostics,"geostrophicInverseDepth")
                r = diagnostics.geostrophicInverseDepth;
                diagnostics = rmfield(diagnostics,"geostrophicInverseDepth");
            elseif evp.name == "geostrophicAPVModes"
                r = lambda;
            else
                r = evp.g*(lambda-evp.parameters.k^2)/evp.f0^2;
            end
            basis = IMGeostrophicModesBasis(solver=solver,evp=evp,nativeModes=V(1:data.n,:), ...
                nativeGModes=V(data.n+1:end,:),companionSolver=IMSolverFiniteDifference(z=data.zG), ...
                eigenvalues=lambda,inverseDepth=reshape(r,1,[]),modeNumber=labels, ...
                modeSelectionDiagnostics=diagnostics,metadata=evp.parameters);
        end

        function basis = solveZeroAPV(problem,solver,solveSystems)
            data = IMGeostrophicFiniteDifferenceTools.configuration(problem,solver);
            nEndpoints = numel(problem.endpoints);
            rhs = zeros(2*data.n+1,nEndpoints);
            for j = 1:nEndpoints
                if problem.endpoints(j)=="surface", rhs(data.n,j) = 1; else, rhs(end,j) = 1; end
            end
            F = zeros(data.n,nEndpoints,numel(problem.k));
            G = zeros(data.n+1,nEndpoints,numel(problem.k));
            solvedK = zeros(1,numel(problem.k));
            solved = cell(1,numel(problem.k));
            count = 0;
            for iK = 1:numel(problem.k)
                index = find(solvedK(1:count)==problem.k(iK),1);
                if isempty(index)
                    count = count+1;
                    solvedK(count) = problem.k(iK);
                    a = problem.g*problem.k(iK)^2/problem.f0^2;
                    sigma = double(problem.surfaceBoundary=="freeSurface");
                    if data.IN==0
                        small = [-(a*diff(problem.zDomain)+sigma),1;0,1];
                        response = solveSystems(small,rhs([data.n,end],:));
                        nativeF = ones(data.n,1)*response(1,:);
                        nativeG = ones(data.n+1,1)*response(2,:)-a*(data.zG-problem.zDomain(1))*response(1,:);
                        solved{count} = [nativeF;nativeG];
                    else
                        [A,B] = IMGeostrophicFiniteDifferenceTools.interiorPencil(data,problem.g);
                        A = A+a*B;
                        A(data.n,:) = 0;
                        A(data.n,1) = -sigma;
                        A(data.n,data.n+1) = 1;
                        A(end,:) = 0;
                        A(end,end) = 1;
                        fieldScale = sqrt(max(data.cellIntegral./data.h)/(problem.g*a));
                        columns = [fieldScale*ones(1,data.n),ones(1,data.n+1)];
                        A = A.*columns;
                        rowScale = max(abs(A),[],2);
                        rowScale(rowScale==0) = 1;
                        solved{count} = columns.'.*solveSystems(A./rowScale,rhs./rowScale);
                    end
                    index = count;
                end
                F(:,:,iK) = solved{index}(1:data.n,:);
                G(:,:,iK) = solved{index}(data.n+1:end,:);
            end
            metadata = problem.metadata;
            metadata.solutionKind = "geostrophicZeroAPVModes";
            metadata.k = problem.k;
            metadata.endpoints = problem.endpoints;
            metadata.surfaceBoundary = problem.surfaceBoundary;
            metadata.modesPerWavenumber = nEndpoints;
            metadata.factorizations = count;
            metadata.discretization = "staggeredFiniteVolume";
            basis = IMGeostrophicZeroAPVModesBasis(problem=problem,solver=solver,nativeModes=F,nativeGModes=G, ...
                companionSolver=IMSolverFiniteDifference(z=data.zG),metadata=metadata);
        end
    end

    methods (Static, Access = private)
        function data = configuration(problem,solver)
            z = solver.zNative(:);
            if numel(z)<3 || any(~isfinite(z)) || any(diff(z)>=0)
                error("IMGeostrophicFiniteDifferenceTools:InvalidGrid","The primal grid must contain at least three distinct finite depths.");
            end
            h = -diff(z);
            WF = [h(1)/2;(h(1:end-1)+h(2:end))/2;h(end)/2];
            cellIntegral = zeros(numel(h),1);
            IMGeostrophicFiniteDifferenceTools.stratification(problem.N2,z);
            for i = 1:numel(h)
                cellIntegral(i) = integral(@(zz) reshape(IMGeostrophicFiniteDifferenceTools.stratification(problem.N2,zz),size(zz)), ...
                    z(i+1),z(i),AbsTol=0,RelTol=1e-12);
            end
            data = struct("zF",z,"zG",[z(1);(z(1:end-1)+z(2:end))/2;z(end)], ...
                "n",numel(z),"h",h,"WF",WF,"cellIntegral",cellIntegral,"IN",sum(cellIntegral));
        end

        function values = stratification(N2,z)
            values = N2(z);
            if numel(values)~=numel(z) || ~isreal(values) || any(~isfinite(values(:))) || any(values(:)<0)
                error("IMGeostrophicModes:InvalidStratification","N2 must return finite real nonnegative values at every requested depth.");
            end
            values = values(:);
        end

        function [A,B] = interiorPencil(data,g)
            n = data.n;
            A = zeros(2*n+1);
            B = zeros(2*n+1);
            for j = 1:n-1
                A(j,j) = 1/data.h(j);
                A(j,j+1) = -1/data.h(j);
                A(j,n+j+1) = data.cellIntegral(j)/(data.h(j)*g);
            end
            for j = 1:n
                A(n+j,n+j) = 1/data.WF(j);
                A(n+j,n+j+1) = -1/data.WF(j);
                B(n+j,j) = 1;
            end
        end

        function [A,B] = inverseDepthPencil(evp,data)
            [A,B] = IMGeostrophicFiniteDifferenceTools.interiorPencil(data,evp.g);
            n = data.n;
            if evp.name=="geostrophicAPVModes"
                sigma = double(evp.parameters.surfaceBoundary=="freeSurface");
                g0 = evp.parameters.g0;
                gd = evp.parameters.gd;
                if isinf(g0)
                    A(n,1) = -sigma;
                    A(n,n+1) = 1;
                else
                    A(n,1) = -(evp.g+sigma*g0);
                    A(n,n+1) = g0;
                end
                if isinf(gd)
                    A(end,end) = 1;
                else
                    A(end,n) = evp.g;
                    A(end,end) = gd;
                end
            else
                A(n,1) = evp.g*evp.parameters.k^2/evp.parameters.alpha0-1;
                A(n,n+1) = 1;
                B(n,1) = -evp.f0^2/evp.parameters.alpha0;
                A(end,n) = -evp.g*evp.parameters.k^2/evp.parameters.alphaD;
                A(end,end) = 1;
                B(end,n) = evp.f0^2/evp.parameters.alphaD;
            end
        end

        function [r,V] = fullyNeutralMode(evp,data)
            F = ones(data.n,1);
            if evp.name=="geostrophicAPVModes"
                if evp.parameters.g0==0 || evp.parameters.gd==0
                    error("IMSolver:NoValidEigenvalues","The fully neutral APV problem has no F-normalizable mode when an endpoint imposes F=0.");
                end
                betaS = double(evp.parameters.surfaceBoundary=="freeSurface")/evp.g+1/evp.parameters.g0;
                betaB = 1/evp.parameters.gd;
                r = evp.g*(betaS+betaB)/diff(evp.zDomain);
                G = -evp.g*betaB+r*(data.zG-evp.zDomain(1));
            else
                a = 1/evp.parameters.alpha0+1/evp.parameters.alphaD;
                cS = 1-evp.g*evp.parameters.k^2/evp.parameters.alpha0;
                if isinf(evp.parameters.alphaD) && (evp.parameters.alpha0==evp.g*evp.parameters.k^2 ...
                        || evp.parameters.alpha0==(evp.g*evp.parameters.k)*evp.parameters.k)
                    cS = 0;
                end
                cB = evp.g*evp.parameters.k^2/evp.parameters.alphaD;
                r = (cS-cB)/(diff(evp.zDomain)+evp.f0^2*a);
                G = (evp.g*evp.parameters.k^2+evp.f0^2*r)/evp.parameters.alphaD+r*(data.zG-evp.zDomain(1));
            end
            V = [F;G];
        end

        function [V,lambda] = refineSmallModes(evp,data,V,lambda)
            small = find(abs(lambda)<1e-6);
            if isempty(small), return; end
            F = V(1:data.n,small);
            G = V(data.n+1:end,small);
            volumeF = sum(data.WF.*F.^2,1);
            volumeG = sum(data.cellIntegral.*G(2:end-1,:).^2,1)/evp.g;
            if evp.name=="geostrophicAPVModes"
                sigma = double(evp.parameters.surfaceBoundary=="freeSurface");
                surface = zeros(1,numel(small));
                bottom = surface;
                if evp.parameters.g0~=0, surface = (sigma+evp.g/evp.parameters.g0)*F(1,:).^2; end
                if evp.parameters.gd~=0, bottom = evp.g/evp.parameters.gd*F(end,:).^2; end
                r = (volumeG+surface+bottom)./volumeF;
            else
                cS = 1-evp.g*evp.parameters.k^2/evp.parameters.alpha0;
                cB = evp.g*evp.parameters.k^2/evp.parameters.alphaD;
                r = (volumeG+cS*F(1,:).^2-cB*F(end,:).^2) ...
                    ./(volumeF+evp.f0^2/evp.parameters.alpha0*F(1,:).^2+evp.f0^2/evp.parameters.alphaD*F(end,:).^2);
            end
            r(lambda(small)==0) = 0;
            primitive = [zeros(1,numel(small));cumsum(data.WF.*F,1)];
            if evp.name=="geostrophicAPVModes" && evp.parameters.gd==0
                Gs = (sigma+evp.g/evp.parameters.g0)*F(1,:);
                G = Gs-r.*primitive;
            else
                if evp.name=="geostrophicAPVModes"
                    Gb = -evp.g/evp.parameters.gd*F(end,:);
                else
                    Gb = (evp.g*evp.parameters.k^2+evp.f0^2*r)/evp.parameters.alphaD.*F(end,:);
                end
                G = Gb+r.*(primitive(end,:)-primitive);
            end
            V(data.n+1:end,small) = G;
            lambda(small) = (r*(data.IN*diff(evp.zDomain)/evp.g)).';
        end

        function valid = qualifiedPairs(evp,data,V,r)
            [A,B] = IMGeostrophicFiniteDifferenceTools.inverseDepthPencil(evp,data);
            if data.IN>0
                fieldScale = [data.IN/evp.g*ones(data.n,1);ones(data.n+1,1)];
                rScale = evp.g/(data.IN*diff(evp.zDomain));
            else
                fieldScale = ones(2*data.n+1,1);
                rScale = 1/diff(evp.zDomain);
            end
            B = rScale*B;
            rowScale = max(max(abs(A.*fieldScale.'),[],2),max(abs(B.*fieldScale.'),[],2));
            rowScale(rowScale==0) = 1;
            A = A./rowScale;
            B = B./rowScale;
            balancedV = V./fieldScale;
            lambda = r/rScale;
            residual = vecnorm(A*V-(B*V).*lambda.',2,1) ...
                ./max((norm(A.*fieldScale.',1)+abs(lambda.')*norm(B.*fieldScale.',1)).*vecnorm(balancedV,2,1),realmin);
            valid = isfinite(r) & residual.'<=1e-8;
        end
    end
end
