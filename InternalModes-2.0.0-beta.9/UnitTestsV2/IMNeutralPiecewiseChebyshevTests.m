classdef IMNeutralPiecewiseChebyshevTests < matlab.unittest.TestCase
    % Qualify the study-only interval solver against independent references.

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot,fullfile(repoRoot,"UnitTestsV2","TestSupport"));
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function neutralAPVModesRetainBothEndpointChoicesAndSignedSectors(testCase)
            p = neutralGeostrophicProfile();
            endpoints = [0 -.01;-p.integratedN2 p.integratedN2];
            for i = 1:size(endpoints,1)
                nIndexedModes = 6;
                if i == 1
                    nIndexedModes = 40;
                end
                problem = IMInternalModes.geostrophicAPVModes(N2=p.N2,zDomain=p.zDomain, ...
                    latitude=45,g0=endpoints(i,1),gd=endpoints(i,2),surfaceBoundary="freeSurface");
                reference = neutralGeostrophicReference(problem,nModes=nIndexedModes,includeSignedSectors=true,tolerance=1e-13);
                candidate = neutralGeostrophicPiecewiseBasis(problem,384,nModes=numel(reference.h),joins=p.joins);
                testCase.verifyIndependentAgreement(reference,candidate);
                testCase.verifyMathematicalIdentities(candidate);
                testCase.verifyNeutralCap(candidate);
                testCase.verifyTrue(any(candidate.inverseDepth < 0))
                testCase.verifyEqual(nnz(candidate.modeNumber > 0),nIndexedModes)
                testCase.verifyBudgetAndInterfaces(candidate,384)
            end
        end

        function neutralMDAModesRetainTheConstantAndSignedSectors(testCase)
            p = neutralGeostrophicProfile();
            endpoints = [0 -.01;-p.integratedN2 p.integratedN2];
            for i = 1:size(endpoints,1)
                problem = IMInternalModes.meanDensityAnomalyModes(N2=p.N2,zDomain=p.zDomain, ...
                    latitude=45,g0=endpoints(i,1),gd=endpoints(i,2));
                reference = neutralGeostrophicReference(problem,nModes=6,includeSignedSectors=true,tolerance=1e-13);
                candidate = neutralGeostrophicPiecewiseBasis(problem,384,nModes=numel(reference.h),joins=p.joins);
                testCase.verifyIndependentAgreement(reference,candidate);
                testCase.verifyMathematicalIdentities(candidate);
                testCase.verifyNeutralCap(candidate);
                testCase.verifyBudgetAndInterfaces(candidate,384)
                values = candidate.evaluate([-1000;-700;-400;-200;-100;0]);
                testCase.verifyEqual(values.F(end,:),zeros(1,numel(candidate.h)),AbsTol=1e-10)
                zero = candidate.modeNumber == 0;
                testCase.verifyEqual(nnz(zero),1)
                testCase.verifyEqual(candidate.inverseDepth(zero),0,AbsTol=0)
                testCase.verifyEqual(candidate.h(zero),Inf,AbsTol=0)
                testCase.verifyEqual(values.G(:,zero),repmat(values.G(1,zero),6,1),RelTol=1e-10,AbsTol=1e-10)
                testCase.verifyEqual(nnz(candidate.modeNumber > 0),6)
                testCase.verifyTrue(any(candidate.inverseDepth < 0))
            end
        end

        function neutralEnstrophyModesRetainBothWavelengthsAndAllEndpointWeights(testCase)
            p = neutralGeostrophicProfile();
            f0 = 2*7.2921e-5*sind(45);
            endpoints = {struct(),struct("alpha0",2*f0^2/1000,"alphaD",3*f0^2/1000), ...
                struct("alpha0",Inf,"alphaD",Inf)};
            for wavelength = [10e3 1e3]
                for i = 1:numel(endpoints)
                    options = endpoints{i};
                    options.N2 = p.N2;
                    options.zDomain = p.zDomain;
                    options.k = 2*pi/wavelength;
                    options.latitude = 45;
                    args = namedargs2cell(options);
                    problem = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes(args{:});
                    nModes = 6;
                    if wavelength == 1e3 && i == 2
                        nModes = 40;
                    end
                    reference = neutralGeostrophicReference(problem,nModes=nModes,tolerance=1e-13);
                    candidate = neutralGeostrophicPiecewiseBasis(problem,384,nModes=nModes,joins=p.joins);
                    testCase.verifyIndependentAgreement(reference,candidate);
                    testCase.verifyMathematicalIdentities(candidate);
                    testCase.verifyNeutralCap(candidate);
                    testCase.verifyBudgetAndInterfaces(candidate,384)
                    testCase.verifyGreaterThan(min(candidate.eigenvalues),0)
                end
            end
        end

        function zeroAPVHasCanonicalResponsesAtBothWavelengths(testCase)
            p = neutralGeostrophicProfile();
            for wavelength = [10e3 1e3]
                problem = IMGeostrophicZeroAPVModes(N2=p.N2,zDomain=p.zDomain,k=2*pi/wavelength,latitude=45);
                reference = neutralGeostrophicReference(problem,tolerance=1e-13);
                candidate = neutralGeostrophicPiecewiseBasis(problem,384,joins=p.joins);
                testCase.verifyIndependentAgreement(reference,candidate);
                testCase.verifyMathematicalIdentities(candidate);
                testCase.verifyNeutralCap(candidate);
                testCase.verifyBudgetAndInterfaces(candidate,384)
                values = candidate.evaluate(p.zDomain.');
                testCase.verifyEqual([values.G(end,:)-values.F(end,:);values.G(1,:)],eye(2),AbsTol=1e-10)
            end
        end

        function constantStratificationAgreesWithAnalyticalModesWithoutAddingPoints(testCase)
            N0 = 5e-3;
            problem = IMInternalModes.geostrophicAPVModes(N2=@(z) N0^2*ones(size(z)), ...
                zDomain=[-1000 0],latitude=45,g0=0,gd=-.01,surfaceBoundary="freeSurface");
            candidate = neutralGeostrophicPiecewiseBasis(problem,128,nModes=6);
            analytical = IMConstantStratificationSolution(N0=N0,zDomain=[-1000 0],latitude=45).internalModes(problem,nModes=6);
            [z,weights] = neutralGeostrophicQuadrature([-1000 -400 -200 0],257);
            values = candidate.evaluate(z);
            exactF = analytical.F(z);
            exactG = analytical.G(z);
            amplitudes = sum(weights.*values.F.*exactF,1)./sum(weights.*values.F.^2,1);
            testCase.verifyEqual(candidate.modeNumber,analytical.modeNumber)
            testCase.verifyEqual(candidate.h,analytical.h,RelTol=1e-7)
            testCase.verifyLessThan(max(vecnorm(values.F.*amplitudes-exactF)./vecnorm(exactF)),1e-7)
            testCase.verifyLessThan(max(vecnorm(values.G.*amplitudes-exactG)./vecnorm(exactG)),1e-7)
            testCase.verifyBudgetAndInterfaces(candidate,128)
            testCase.verifyMathematicalIdentities(candidate);
        end

        function exactEnstrophyZeroInverseDepthSurvivesIntervalMatching(testCase)
            p = neutralGeostrophicProfile();
            k = 2*pi/500e3;
            g = 9.81;
            problem = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                N2=p.N2,zDomain=p.zDomain,k=k,g=g,latitude=45,alpha0=g*k*k,alphaD=Inf);
            candidate = neutralGeostrophicPiecewiseBasis(problem,128,nModes=6,joins=p.joins);
            zero = candidate.inverseDepth == 0;
            testCase.verifyEqual(nnz(zero),1)
            testCase.verifyEqual(candidate.h(zero),Inf,AbsTol=0)
            testCase.verifyEqual(candidate.eigenvalues(zero),k*k,AbsTol=0)
            values = candidate.evaluate(linspace(-1000,0,101).');
            testCase.verifyEqual(values.F(:,zero),repmat(values.F(1,zero),101,1),RelTol=1e-10,AbsTol=1e-10)
            testCase.verifyEqual(values.G(:,zero),zeros(101,1),AbsTol=1e-10)
            testCase.verifyMathematicalIdentities(candidate);
        end
    end

    methods (Access = private)
        function verifyIndependentAgreement(testCase,reference,candidate)
            [z,weights] = neutralGeostrophicQuadrature(reference.joins,513);
            metrics = neutralGeostrophicComparison(reference,candidate,z,weights);
            testCase.verifyTrue(metrics.labelsMatch)
            testCase.verifyLessThanOrEqual(metrics.maximumFH1Error,1e-7)
            testCase.verifyLessThanOrEqual(metrics.maximumGH1Error,1e-7)
            testCase.verifyLessThanOrEqual(metrics.maximumRelativeHError,1e-7)
        end

        function verifyBudgetAndInterfaces(testCase,candidate,nEVP)
            testCase.verifyEqual(sum(candidate.metadata.nPerInterval),nEVP)
            testCase.verifyNumElements(candidate.metadata.nPerInterval,3)
            testCase.verifyLessThanOrEqual(max(abs(candidate.metadata.interfaceResiduals),[],"all"),1e-7)
        end

        function verifyNeutralCap(testCase,candidate)
            z = (-200:20:0).';
            values = candidate.evaluate(z);
            [zq,~] = neutralGeostrophicQuadrature(candidate.configuration.joins,257);
            fullDomain = candidate.evaluate(zq);
            % A trapped mode can have an exponentially small cap trace;
            % scale its interval identity by the physical mode amplitude.
            FScale = max(max(abs(fullDomain.F),[],1),1);
            GScale = max(max(abs(fullDomain.G),[],1),1);
            testCase.verifyLessThanOrEqual(max(abs(values.F-values.F(1,:))./FScale,[],"all"),1e-10)
            if candidate.family == "mda"
                expectedG = values.G(1,:)+(z-z(1)).*values.Gz(1,:);
            else
                expectedG = values.G(1,:)+(z-z(1)).*(candidate.inverseDepth.*values.F(1,:));
            end
            testCase.verifyLessThanOrEqual(max(abs(values.G-expectedG)./GScale,[],"all"),1e-10)
        end

        function verifyMathematicalIdentities(testCase,candidate)
            c = candidate.configuration;
            [z,weights] = neutralGeostrophicQuadrature([-1000 -400 -200 0],513);
            values = candidate.evaluate(z);
            r = candidate.inverseDepth;
            D = diff(c.zDomain);
            FScale = max(abs(values.F),[],1);
            GScale = max(abs(values.G),[],1);
            first = values.Fz+(c.N2(z)/c.g).*values.G;
            firstScale = max(FScale/D,max(abs((c.N2(z)/c.g).*values.G),[],1));
            if candidate.family == "mda"
                second = values.Gzz+r.*(c.N2(z)/c.g).*values.G;
                secondScale = max(GScale/D^2,max(abs(r.*(c.N2(z)/c.g).*values.G),[],1));
            else
                second = values.Gz-r.*values.F;
                secondScale = max(GScale/D,max(abs(r.*values.F),[],1));
            end
            testCase.verifyLessThanOrEqual(max(abs(first)./max(firstScale,realmin),[],"all"),1e-7)
            testCase.verifyLessThanOrEqual(max(abs(second)./max(secondScale,realmin),[],"all"),1e-7)
            testCase.verifyEndpointConditions(candidate,values,FScale,GScale,D);
            FGram = values.F.'*(weights.*values.F);
            GGram = values.G.'*(weights.*(c.N2(z)/c.g).*values.G);
            Fs = values.F(end,:);
            Fb = values.F(1,:);
            Gs = values.G(end,:);
            Gb = values.G(1,:);
            if candidate.family == "apv"
                if c.g0 ~= 0
                    GGram = GGram+(c.sigma+c.g/c.g0)*(Fs.'*Fs);
                end
                if c.gd ~= 0
                    GGram = GGram+(c.g/c.gd)*(Fb.'*Fb);
                end
                testCase.verifyLessThanOrEqual(norm(FGram/D-eye(numel(r)),2),2e-7)
                scale = sqrt(max(abs(r),1/D));
                testCase.verifyLessThanOrEqual(norm((GGram/D-diag(r))./(scale.'*scale),2),2e-7)
            elseif candidate.family == "gpe"
                FGram = FGram+(c.f0^2/c.alpha0)*(Fs.'*Fs)+(c.f0^2/c.alphaD)*(Fb.'*Fb);
                GGram = GGram+(1-c.g*c.k^2/c.alpha0)*(Fs.'*Fs)-(c.g*c.k^2/c.alphaD)*(Fb.'*Fb);
                testCase.verifyLessThanOrEqual(norm(FGram/D-eye(numel(r)),2),2e-7)
                scale = sqrt(max(abs(r),1/D));
                testCase.verifyLessThanOrEqual(norm((GGram/D-diag(r))./(scale.'*scale),2),2e-7)
            elseif candidate.family == "mda"
                if isfinite(c.g0)
                    GGram = GGram+(c.g0/c.g)*(Gs.'*Gs);
                end
                if isfinite(c.gd)
                    GGram = GGram+(c.gd/c.g)*(Gb.'*Gb);
                end
                testCase.verifyLessThanOrEqual(norm(GGram-diag(candidate.signatures),2),2e-7)
            else
                testCase.verifyEqual([Gs-c.sigma*Fs;Gb],eye(2),AbsTol=1e-10)
            end
        end

        function verifyEndpointConditions(testCase,candidate,values,FScale,GScale,D)
            c = candidate.configuration;
            Fs = values.F(end,:);
            Fb = values.F(1,:);
            Gs = values.G(end,:);
            Gb = values.G(1,:);
            r = candidate.inverseDepth;
            if candidate.family == "apv"
                if isinf(c.g0)
                    surface = Gs-c.sigma*Fs;
                    surfaceScale = GScale+c.sigma*FScale;
                else
                    surface = c.g0*Gs-(c.g+c.sigma*c.g0)*Fs;
                    surfaceScale = abs(c.g0)*GScale+abs(c.g+c.sigma*c.g0)*FScale;
                end
                if isinf(c.gd)
                    bottom = Gb;
                    bottomScale = GScale;
                else
                    bottom = c.g*Fb+c.gd*Gb;
                    bottomScale = c.g*FScale+abs(c.gd)*GScale;
                end
            elseif candidate.family == "gpe"
                if isinf(c.alpha0)
                    surface = Gs-Fs;
                    surfaceScale = GScale+FScale;
                else
                    surface = c.alpha0*(Gs-Fs)+c.g*candidate.eigenvalues.*Fs;
                    surfaceScale = c.alpha0*(GScale+FScale)+c.g*abs(candidate.eigenvalues).*FScale;
                end
                if isinf(c.alphaD)
                    bottom = Gb;
                    bottomScale = GScale;
                else
                    bottom = c.alphaD*Gb-c.g*candidate.eigenvalues.*Fb;
                    bottomScale = c.alphaD*GScale+c.g*abs(candidate.eigenvalues).*FScale;
                end
            elseif candidate.family == "mda"
                surface = values.Gz(end,:)-(c.g0/c.g)*r.*Gs;
                bottom = values.Gz(1,:)+(c.gd/c.g)*r.*Gb;
                surfaceScale = GScale/D+abs(c.g0*r/c.g).*GScale;
                bottomScale = GScale/D+abs(c.gd*r/c.g).*GScale;
            else
                surface = Gs-c.sigma*Fs-[1 0];
                bottom = Gb-[0 1];
                surfaceScale = max(1,GScale+c.sigma*FScale);
                bottomScale = max(1,GScale);
            end
            errors = [abs(surface)./max(surfaceScale,realmin) abs(bottom)./max(bottomScale,realmin)];
            testCase.verifyLessThanOrEqual(max(errors),1e-7)
        end
    end
end
