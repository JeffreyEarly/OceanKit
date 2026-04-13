function spectrum = coherenceSpectrumFromStorage(frequency, coherence)
frequency = reshape(frequency, [], 1);
coherence = reshape(coherence, [], 1);

if isempty(frequency) && isempty(coherence)
    spectrum = struct("frequency", zeros(0, 1), "coherence", zeros(0, 1));
    return
end

if numel(frequency) == 1 && numel(coherence) == 1 && isnan(frequency) && isnan(coherence)
    spectrum = struct("frequency", zeros(0, 1), "coherence", zeros(0, 1));
    return
end

spectrum = struct("frequency", frequency, "coherence", coherence);
end
