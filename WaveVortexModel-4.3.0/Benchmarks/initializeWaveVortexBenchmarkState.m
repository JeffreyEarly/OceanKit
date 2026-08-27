function state = initializeWaveVortexBenchmarkState(wvt,seed)
% Initialize and capture the deterministic canonical benchmark state.
arguments
    wvt WVTransform
    seed (1,1) double
end
rng(seed,"twister");
wvt.initWithRandomFlow(uvMax=0.01);
state = struct("t0",wvt.t,"Ap",[],"Am",[],"A0",[]);
if ~isempty(wvt.Ap)
    state.Ap = wvt.Ap;
end
if ~isempty(wvt.Am)
    state.Am = wvt.Am;
end
if ~isempty(wvt.A0)
    state.A0 = wvt.A0;
end
end
