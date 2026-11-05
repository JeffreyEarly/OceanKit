classdef (Hidden, Sealed) IMBasisContracts
    % Describe continuous scientific bases without evaluating profiles.
    methods (Static)
        function value = capabilities(basis)
            variables = "u"; orders = {0:1};
            interval = isa(basis,"IMNeutralLabelModesBasis") || isa(basis,"IMZeroDepthWaveModesBasis");
            endpoint = isa(basis,"IMGeostrophicZeroAPVModesBasis") || isa(basis,"IMAnalyticalGeostrophicZeroAPVModesBasis");
            if isa(basis,"IMNeutralLabelModesBasis")
                variables = ["G","APV"]; orders = {0:2,0};
            elseif endpoint || isa(basis,"IMZeroDepthWaveModesBasis")
                variables = ["F","G"]; orders = {0:1,0:2};
            elseif isa(basis,"IMInternalModesBasis") || isa(basis,"IMAnalyticalInternalModesBasis")
                variables = ["F","G"]; orders = {0,0};
                evp = basis.evp;
                solved = find(variables == string(evp.formulation));
                orders{solved} = 0:1;
                if isa(basis,"IMInternalModesBasis"), orders{solved} = 0:2; end
                if IMInternalModesNormalizationTools.isWave(evp) || IMInternalModesNormalizationTools.hasHydrostaticDerivativeRelations(evp) || isa(basis,"IMMeanDensityAnomalyModesBasis") || isa(basis,"IMGeostrophicModesBasis")
                    orders = {0:1,0:2};
                end
            end
            derivativeConvention = "classicalPhysical"; traceConvention = "physicalDomainEvaluation";
            if interval, derivativeConvention = "classicalIntervalInterior"; traceConvention = "insideIntervalEndpoints"; end
            spectra = strings(1,0);
            for name = ["h","eigenvalues","frequency"]
                if isprop(basis,name), spectra(end+1) = name; end %#ok<AGROW>
            end
            horizontal = IMBasisContracts.horizontalIdentity(basis);
            if horizontal.kind == "perColumn", spectra = [spectra,"kSquared"]; end
            if endpoint && isprop(basis,"metadata") && isfield(basis.metadata,"basisChoice") && basis.metadata.basisChoice == "generalizedEnergyOrthogonal"
                spectra = [spectra,"boundaryDepth"];
            end
            operations = struct("variable",{},"kind",{},"availability",{},"requirements",{});
            if isa(basis,"IMNeutralLabelModesBasis")
                operations = struct(variable="APV",kind="prescribedContinuousDual",availability="conditional",requirements=["fixed APV observations and dz quadrature","nonzero f0 and finite nonzero retained APV norms","independent particle-label information for a complete label state"]);
            elseif isprop(basis,"evp") && isa(basis.evp,"IMInternalModes") && ~IMInternalModesNormalizationTools.isWave(basis.evp) && basis.evp.name ~= "geostrophicGeneralizedPotentialEnstrophyModes"
                evp = basis.evp;
                for variable = ["F","G"]
                    if evp.modeFamily == "meanDensityAnomaly" && variable == "F", continue; end
                    spec = evp.innerProduct(variable);
                    if ~spec.hasInnerProduct, continue; end
                    operations(end+1) = struct(variable=variable,kind="familyPhysicalDual",availability="conditional",requirements=["independent observations required by the physical pairing","fixed physical samples and dz quadrature","separate explicit numerical recovery policy"]); %#ok<AGROW>
                end
            end
            recovery = struct(synthesisVariables=variables,preparation="IMDiscreteTransform.fromObservationMatrices",availability="requiresPreparedSystem",requirements=["explicit observation matrix and scientific column activity","explicit fitting and observation-noise meanings","successful prepared-system rank qualification"],noiseEvidence="requiresDeclaredObservationFactor");
            value = struct(variables=variables,derivativeOrders={orders},derivativeConvention=derivativeConvention,traceConvention=traceConvention,traceOrders={orders},spectrumNames=spectra,normalizationNames=reshape(string(basis.normalizationNames()),1,[]),normalizationRequirements="Every retained column needs a finite positive divisor; signed scientific norms retain their signs separately.",projectionOperations=operations,sampledRecovery=recovery,supportsEmptyColumns=true);
            if endpoint
                value.quadraticFormNames = ["energyMatrix","generalizedEnergyMatrix","surfaceBuoyancyMatrix","bottomBuoyancyMatrix"];
                value.basisOperations = "orthogonalizeGeneralizedEnergy";
            end
        end

        function value = horizontalIdentity(basis)
            value = struct(kind="unspecified",k=zeros(1,0),kSquared=zeros(1,0));
            if isa(basis,"IMZeroDepthWaveModesBasis") || isa(basis,"IMGeostrophicZeroAPVModesBasis") || isa(basis,"IMAnalyticalGeostrophicZeroAPVModesBasis")
                value.kind = "fixedWavenumber"; value.k = reshape(basis.k,1,[]);
            elseif isa(basis,"IMNeutralLabelModesBasis")
                value.kind = "independent";
            elseif isprop(basis,"evp") && isa(basis.evp,"IMInternalModes")
                evp = basis.evp;
                if evp.name == "waveModesAtFrequency"
                    value.kind = "perColumn";
                    omega = evp.parameters.omega;
                    value.kSquared = ((omega-abs(evp.f0))*(omega+abs(evp.f0))/evp.g)*reshape(basis.eigenvalues,1,[]);
                elseif isfield(evp.parameters,"k")
                    value.kind = "fixedWavenumber"; value.k = reshape(evp.parameters.k,1,[]);
                else
                    value.kind = "independent";
                end
            end
        end

        function state = initialState(basis)
            n = basis.nColumns;
            state = struct(originalColumnCount=n,initialColumnIndices=1:n,activeSourceColumns=1:n,sourceEigenpairIndices=zeros(1,0),sourceModeNumbers=zeros(1,0),requestedScientific=struct(),requestedSolver=struct(),request=struct(),selectionDiagnostics=struct(),sourceClass=string(class(basis)),constructionSkipped=false);
            if isprop(basis,"modeNumber"), state.sourceModeNumbers = reshape(basis.modeNumber,1,[]); end
            if isprop(basis,"modeSelectionDiagnostics"), state.selectionDiagnostics = basis.modeSelectionDiagnostics; end
            if isprop(basis,"solver"), state.requestedSolver = IMBasisContracts.solverRecord(basis.solver,true); end
            if isprop(basis,"evp")
                evp = basis.evp;
                if isprop(evp,"requestedSpecification") && ~isempty(fieldnames(evp.requestedSpecification)), state.requestedScientific = evp.requestedSpecification; end
            elseif isprop(basis,"problem")
                problem = basis.problem;
                state.requestedScientific = struct(k=problem.k,endpoints=problem.endpoints,surfaceBoundary=problem.surfaceBoundary,zDomain=problem.zDomain,N2=problem.N2);
                state.request = struct(endpoints=problem.endpoints);
                if ~isempty(problem.background), state.physicalBackground = IMBasisContracts.physicalBackground(problem.background); end
            end
        end

        function state = selectState(state,columns)
            if ~isempty(state.initialColumnIndices), state.initialColumnIndices = state.initialColumnIndices(columns); end
            if ~isempty(state.activeSourceColumns), state.activeSourceColumns = state.activeSourceColumns(columns); end
            if isfield(state,"sourceDirections"), state.sourceDirections = state.sourceDirections(:,columns,:); end
        end

        function value = scientificSpecification(basis)
            horizontal = IMBasisContracts.horizontalIdentity(basis);
            physical = struct(zDomain=basis.zDomain);
            evpRecord = struct(); pairings = struct(); energy = struct(); role = "notSpecified";
            profileRelation = "sourceDefined"; intervals = zeros(0,2); intervalKind = "none";
            if isprop(basis,"evp")
                evp = basis.evp;
                evpRecord = struct(className=string(class(evp)),name=evp.name,parameters=evp.parameters,p=evp.p,q=evp.q,r=evp.r,surfaceBoundary=IMBasisContracts.boundaryRecord(evp.surfaceBoundary),bottomBoundary=IMBasisContracts.boundaryRecord(evp.bottomBoundary));
                if isa(evp,"IMInternalModes")
                    physical = IMBasisContracts.backgroundRecord(evp.background);
                    evpRecord.formulation = evp.formulation; evpRecord.modeFamily = evp.modeFamily;
                    pairings.F = evp.innerProduct("F"); pairings.G = evp.innerProduct("G");
                    if IMInternalModesNormalizationTools.isWave(evp), profileRelation = "F=D*Gz"; else, profileRelation = "familyDefinedCoupledFields"; end
                    if isfield(evp.parameters,"g0"), energy = struct(g0=evp.parameters.g0,gd=evp.parameters.gd); role = "definesPairingAndEndpointConditions"; end
                end
            elseif isprop(basis,"background")
                physical = IMBasisContracts.backgroundRecord(basis.background);
            else
                for name = ["N2","rotationRate","latitude","f0","g"]
                    if isprop(basis,name), physical.(name) = basis.(name); end
                end
            end
            if isa(basis,"IMNeutralLabelModesBasis")
                intervals = basis.neutralIntervals; intervalKind = "neutral"; profileRelation = "APV=-f0*Gz";
                energy = struct(g0=basis.g0,gd=basis.gd); role = "provenanceOnly";
            elseif isa(basis,"IMZeroDepthWaveModesBasis")
                intervals = basis.criticalIntervals; intervalKind = "critical"; profileRelation = "F=D*Gz";
            end
            state = basis.basisContractState;
            if isfield(state,"physicalBackground"), physical = state.physicalBackground; end
            normalization = struct(name=string(basis.normalization),scalesStatus="lazy; request normalizationFactors explicitly",signaturesStatus="lazy; request normSignatures explicitly");
            resolved = struct(evp=evpRecord,pairings=pairings,energyBoundary=energy,energyBoundaryRole=role,intervals=intervals,intervalKind=intervalKind,intervalBoundaryTypes=repmat("Dirichlet",size(intervals)),endpointDerivativeTrace="notApplicable",profileRelation=profileRelation,horizontal=horizontal,columnLabels=basis.columnLabels,modeInfo=basis.modeInfo,sourceColumns=state.activeSourceColumns,normalization=normalization,includesAnalyticalBarotropic=false);
            if ~isempty(intervals), resolved.endpointDerivativeTrace = "insideDeclaredInterval"; end
            if isprop(basis,"modeInfo") && ismember("branch",string(basis.modeInfo.Properties.VariableNames)), resolved.includesAnalyticalBarotropic = any(basis.modeInfo.branch == "barotropic"); end
            if isprop(basis,"canonicalEndpointLabels")
                resolved.endpoints = basis.canonicalEndpointLabels;
                if isprop(basis,"metadata")
                    if isfield(basis.metadata,"basisChoice"), resolved.basisChoice = basis.metadata.basisChoice; end
                    if isfield(basis.metadata,"generalizedEnergyCoefficients")
                        resolved.energyBoundary = basis.metadata.generalizedEnergyCoefficients;
                        resolved.energyBoundaryRole = "definesGeneralizedEnergyOrthogonality";
                    end
                    if isfield(basis.metadata,"orthogonality"), resolved.orthogonality = basis.metadata.orthogonality; end
                    if isfield(basis.metadata,"basisChoice") && basis.metadata.basisChoice == "generalizedEnergyOrthogonal", resolved.boundaryDepth = basis.boundaryDepth; end
                end
                if isfield(state,"sourceDirections"), resolved.sourceDirections = state.sourceDirections; end
                resolved.profileRelation = "Gz=-g*k^2*F/f0^2; Fz=-N2*G/g";
                resolved.boundaryResponseMatrixStatus = "available on demand in current normalization";
            end
            value = struct(schemaVersion=1,family=basis.family,background=physical,requested=state.requestedScientific,resolved=resolved);
        end

        function value = constructionRecord(basis)
            state = basis.basisContractState;
            backend = "intrinsic"; actual = struct();
            if isprop(basis,"solver")
                actual = IMBasisContracts.solverRecord(basis.solver,true); backend = "numerical";
            elseif isa(basis,"IMAnalyticalInternalModesBasis") || isa(basis,"IMAnalyticalGeostrophicZeroAPVModesBasis")
                backend = "analytical"; actual = struct(className="IMSolverAnalytical",coordinateKind="analytical",nEVP=[],nativeRepresentationAvailable=false);
            elseif isa(basis,"IMZeroDepthWaveModesBasis")
                backend = "analyticalDDKernel";
            end
            selection = struct(request=state.request,sourceEigenpairIndices=state.sourceEigenpairIndices,sourceModeNumbers=state.sourceModeNumbers,diagnostics=state.selectionDiagnostics,activeSourceColumns=state.activeSourceColumns,initialColumnIndices=state.initialColumnIndices,columnLabels=basis.columnLabels,originalColumnCount=state.originalColumnCount);
            if isfield(state,"sourceDirections"), selection.sourceDirections = state.sourceDirections; end
            evidence = struct(solveAccuracy=struct(status="unverified"),referenceConvergence=struct(status="unverified"));
            value = struct(schemaVersion=1,requested=struct(solver=state.requestedSolver),resolved=struct(backendKind=backend,sourceClass=state.sourceClass,solver=actual,constructionSkipped=state.constructionSkipped),selection=selection,evidence=evidence);
        end

        function value = normalizationSignatures(basis)
            % Signs belong to the active normalization, not to an inferred state energy.
            value = ones(1,basis.nColumns);
            if basis.nColumns == 0, return; end
            if string(basis.normalization) == "meanDensityAnomalyEnergy"
                gram = basis.gramMatrix(variable="G");
                value = sign(real(diag(gram))).';
            end
        end

        function value = requestedArguments(args)
            value = struct(alpha0=[],alphaD=[]);
            for index = 1:2:numel(args)
                name = string(args{index});
                value.(name) = args{index+1};
            end
            if isfield(value,"g0") || isfield(value,"gd")
                energy = struct();
                if isfield(value,"g0"), energy.g0 = value.g0; end
                if isfield(value,"gd"), energy.gd = value.gd; end
                value.energyBoundary = energy;
            end
        end

        function value = physicalBackground(background)
            value = IMBasisContracts.backgroundRecord(background);
        end

        function validateEvaluation(basis,variable,order)
            description = basis.capabilities;
            index = find(description.variables == variable,1);
            if isempty(index)
                error(string(class(basis))+":UnsupportedVariable","The basis does not provide channel %s.",variable);
            end
            if ~ismember(order,description.derivativeOrders{index})
                error(string(class(basis))+":UnsupportedDerivative","Derivative order %d of %s is unavailable.",order,variable);
            end
        end

        function validateNormalization(basis,name)
            name = IMInternalModesNormalizationTools.name(name,string(class(basis)));
            if ~ismember(name,string(basis.normalizationNames()))
                error(string(class(basis))+":UnsupportedNormalization","The basis does not define normalization %s.",string(name));
            end
        end

        function value = solverRecord(solver,hasNative)
            value = struct(className=string(class(solver)),nEVP=[],coordinateKind="notApplicable",zNative=[],xNative=[],nativeRepresentationAvailable=false);
            if isprop(solver,"nEVP"), value.nEVP = solver.nEVP; end
            if isprop(solver,"coordinateKind"), value.coordinateKind = solver.coordinateKind; end
            if isa(solver,"IMSolverFiniteDifference"), value.coordinateKind = "z"; hasNative = true; end
            if isa(solver,"IMSolverAnalytical"), value.coordinateKind = "analytical"; end
            if hasNative
                if isprop(solver,"zNative"), value.zNative = solver.zNative; end
                if isprop(solver,"xNative"), value.xNative = solver.xNative; end
                value.nativeRepresentationAvailable = ~isempty(value.zNative);
            end
        end
    end
    methods (Static, Access = private)
        function value = backgroundRecord(background)
            value = struct();
            for name = ["N2","zDomain","rotationRate","latitude","f0","g","surfaceCondition","profileKind","profileParameters"]
                value.(name) = background.(name);
            end
        end
        function value = boundaryRecord(boundary)
            value = struct(a=boundary.a,b=boundary.b,c=boundary.c,d=boundary.d);
        end
    end
end
