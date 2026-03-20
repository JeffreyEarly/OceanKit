basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";
basedir = "/Users/Shared/CimRuns_June2025/output/";
wvd = WVDiagnostics(basedir + replace(getRunParameters(18),"256","512") + ".nc");
wvt = wvd.wvt;
% timeIndices = 2751:3001;


%%
diagfile = wvd.diagfile;
groupName = "reservoir-damped-wave-geo";

%%
% we want to save name and 
svv = wvt.forcingWithName("adaptive damping");
NoDampMask = (wvt.Kh < (svv.k_damp+svv.k_no_damp)/2) & (wvt.J < (svv.j_damp + svv.j_no_damp)/2);

i = 1;
flowComponents(i) = wvt.flowComponentWithName("wave") + wvt.flowComponentWithName("inertial");
flowComponents(i).name = "wave";
flowComponents(i).maskAp = flowComponents(i).maskAp & NoDampMask;
flowComponents(i).maskAm = flowComponents(i).maskAm & NoDampMask;
componentEnergy{i} =  @(Ep,Em,E0) sum(flowComponents(i).maskAp(:).*Ep(:) + flowComponents(i).maskAm(:).*Em(:));

i = i+1;
flowComponents(i) = wvt.flowComponentWithName("geostrophic") + wvt.flowComponentWithName("mda");
flowComponents(i).name = "geostrophic";
flowComponents(i).maskA0 = flowComponents(i).maskA0 & NoDampMask;
componentEnergy{i} =  @(Ep,Em,E0) sum(flowComponents(i).maskA0(:).*E0(:));

i = i+1;
flowComponents(i) = wvt.flowComponentWithName("wave") + wvt.flowComponentWithName("inertial");
flowComponents(i).name = "damped wave";
flowComponents(i).maskAp = flowComponents(i).maskAp & ~NoDampMask;
flowComponents(i).maskAm = flowComponents(i).maskAm & ~NoDampMask;
componentEnergy{i} =  @(Ep,Em,E0) sum(flowComponents(i).maskAp(:).*Ep(:) + flowComponents(i).maskAm(:).*Em(:));

i = i+1;
flowComponents(i) = wvt.flowComponentWithName("geostrophic") + wvt.flowComponentWithName("mda");
flowComponents(i).name = "damped geostrophic";
flowComponents(i).maskA0 = flowComponents(i).maskA0 & ~NoDampMask;
componentEnergy{i} =  @(Ep,Em,E0) sum(flowComponents(i).maskA0(:).*E0(:));

% componentEnergy{i} =  @(Ep,Em,E0) sum(flowComponents(i).maskAp(:).*Ep(:) + flowComponents(i).maskAm(:).*Em(:) + flowComponents(i).maskA0(:).*E0(:));

%%

iTriad = 0;
for i=1:length(flowComponents)
    for j=1:i
        for k=1:length(flowComponents)
            iTriad = iTriad + 1;
            triadName(iTriad) = "T_" + i + "_" + j + "_" + k;
        end
    end
end
%%


triadVar = configureDictionary("string","cell");
if diagfile.hasGroupWithName(groupName)
    group = diagfile.groupWithName(groupName);
    for iTriad=1:length(triadName)
        triadVar{triadName(iTriad)} = group.variableWithName(triadName(iTriad));
    end
else
    group = diagfile.addGroup(groupName);
    group.addAttribute('flow-components',reshape([flowComponents.name],[],1));
    for iTriad=1:length(triadName)
        triadVar{triadName(iTriad)} = group.addVariable(triadName(iTriad),"t",type="double",isComplex=false);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Loop over the the requested time indices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

integrationLastInformWallTime = datetime('now');
loopStartTime = integrationLastInformWallTime;
integrationLastInformLoopNumber = 1;
integrationInformTime = 10;
fprintf("Starting loop to compute reservoir fluxes for %d time indices.\n",length(timeIndices));
for timeIndex = 1:length(timeIndices)
    deltaWallTime = datetime('now')-integrationLastInformWallTime;
    if ( seconds(deltaWallTime) > integrationInformTime)
        wallTimePerLoopTime = deltaWallTime / (timeIndex - integrationLastInformLoopNumber);
        wallTimeRemaining = wallTimePerLoopTime*(length(timeIndices) - timeIndex + 1);
        fprintf('Time index %d of %d. Estimated time to finish is %s (%s)\n', timeIndex, length(timeIndices), wallTimeRemaining, datetime(datetime('now')+wallTimeRemaining,TimeZone='local',Format='d-MMM-y HH:mm:ss Z')) ;
        integrationLastInformWallTime = datetime('now');
        integrationLastInformLoopNumber = timeIndex;
    end

    outputIndex = timeIndices(timeIndex);
    wvd.iTime = timeIndices(timeIndex);

    for i=1:length(flowComponents)
        for j=1:i
            if i==j
                [Fp,Fm,F0] = wvt.nonlinearFluxForFlowComponents(flowComponents(i),flowComponents(j));
            else
                [Fp_ij,Fm_ij,F0_ij] = wvt.nonlinearFluxForFlowComponents(flowComponents(i),flowComponents(j));
                [Fp_ji,Fm_ji,F0_ji] = wvt.nonlinearFluxForFlowComponents(flowComponents(j),flowComponents(i));
                Fp = Fp_ji + Fp_ij;
                Fm = Fm_ji + Fm_ij;
                F0 = F0_ji + F0_ij;
            end
            [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0);
            for k=1:length(flowComponents)
                var = triadVar{"T_" + i + "_" + j + "_" + k};
                dE = componentEnergy{k}(Ep,Em,E0);
                var.setValueAlongDimensionAtIndex(dE,'t',outputIndex);
            end
        end
    end

end

%%
if ~diagfile.hasGroupWithName(groupName)
    error("Unable to find the group "+groupName+".");
end
group = diagfile.groupWithName(groupName);
flowComponentNames = reshape(group.attributes('flow-components'),[],1);

% loop over T_i_j_k
transfers = zeros(length(flowComponentNames),length(flowComponentNames));
for k=1:length(flowComponentNames)
    % Flux *into* the k-th reservoir
    for i=1:length(flowComponentNames)
        for j=1:i
            
            if i==j
                % has the form T_i_i_k OR T_k_k_k
                E = group.readVariables("T_" + i + "_" + j + "_" + k);
                E = mean(E(timeIndices));
                transfers(i,k) = transfers(i,k) + E;
            elseif i==k % our for-loop is such that i will reach k, but never let j get that high
                % has the form T_k_j_k, so we need -T_k_k_j
                E = group.readVariables("T_" + k + "_" + k + "_" + j);
                E = mean(E(timeIndices));
                transfers(j,k) = transfers(j,k) - E;
            elseif j==k % our for-loop is such that i will reach k, but never let j get that high
                % has the form T_i_k_k, so we need -T_k_k_i
                E = group.readVariables("T_" + k + "_" + k + "_" + i);
                E = mean(E(timeIndices));
                transfers(i,k) = transfers(i,k) - E;
            else
                % now we have a true mixed triad
                E_i_j_k = group.readVariables("T_" + i + "_" + j + "_" + k);
                E_j_k_i = group.readVariables("T_" + max(j,k) + "_" + min(j,k) + "_" + i);
                E_k_i_j = group.readVariables("T_" + max(i,k) + "_" + min(i,k) + "_" + j);

                E_i_j_k = mean(E_i_j_k(timeIndices));
                E_j_k_i = mean(E_j_k_i(timeIndices));
                E_k_i_j = mean(E_k_i_j(timeIndices));

                d = (E_j_k_i + E_k_i_j + E_i_j_k)/wvd.flux_scale; % / max([abs(E_j_k_i),abs(E_k_i_j),abs(E_i_j_k)]);
                disp( "T_" + i + "_" + j + "_" + k + " with relative error " + string(d))

                % a -> i, b -> j, c -> k
                [Tji, Tki, Tkj] = transfersFromFluxes(E_j_k_i,E_k_i_j,E_i_j_k);
                transfers(i,k) = transfers(i,k) - Tki;
                transfers(j,k) = transfers(j,k) - Tkj;
            end
        end
    end
end

function [Tba, Tca, Tcb] = transfersFromFluxes(dA,dB,dC)
if abs(dC) >= abs(dA) && abs(dC) >= abs(dB)
    Tba = 0;
    Tca = dA;
    Tcb = dB;
elseif abs(dB) >= abs(dA) && abs(dB) >= abs(dC)
    Tba = dA;
    Tca = 0;
    Tcb = -dC;
else
    Tba = -dB;
    Tca = -dC;
    Tcb = 0;
end
end

% function [Tba, Tca, Tcb] = transfersFromFluxes(dA,dB,dC)
% if abs(dA) < abs(dB) && abs(dA) < abs(dC)
%     alpha = abs(dB/(dC-dB));
%     Tba = alpha*dA;
%     Tca = (1-alpha)*dA;
%     Tcb = dB + alpha*dA;
% elseif abs(dB) < abs(dA) && abs(dB) < abs(dC)
%     alpha = abs(dA/(dA-dC));
%     Tba =- alpha*dB;
%     Tca = dA + alpha*dB;
%     Tcb = (1-alpha)*dB;
% else
%     alpha = abs(dA/(dA-dB));
%     Tba = dA + alpha*dC;
%     Tca = -alpha*dC;
%     Tcb = -(1-alpha)*dC;
% end
% end