function failure = wvCompiledBackendFailure(stage,exception)
stack = strings(numel(exception.stack),1);
for iFrame = 1:numel(exception.stack)
    frame = exception.stack(iFrame);
    stack(iFrame) = string(frame.name)+":"+frame.line;
end
failure = struct("stage",string(stage),"identifier",string(exception.identifier),"message",string(exception.message),"stack",stack);
end
