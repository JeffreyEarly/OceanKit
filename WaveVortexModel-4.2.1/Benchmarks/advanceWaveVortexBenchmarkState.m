function advanceWaveVortexBenchmarkState(wvt,state,iState)
% Advance through public state setters without explicitly clearing caches.
arguments
    wvt WVTransform
    state (1,1) struct
    iState (1,1) double {mustBeInteger,mustBeNonnegative}
end
theta = 0.017*iState;
wvt.t = state.t0 + 30*iState;
if ~isempty(state.Ap)
    wvt.Ap = state.Ap.*exp(1i*theta);
end
if ~isempty(state.Am)
    wvt.Am = state.Am.*exp(-1i*theta);
end
if ~isempty(state.A0)
    wvt.A0 = state.A0.*exp(1i*theta/3);
end
end
