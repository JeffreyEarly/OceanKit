function storageSpectrum = coherenceSpectrumStorage(spectrum)
frequency = reshape(spectrum.frequency, [], 1);
coherence = reshape(spectrum.coherence, [], 1);

if isempty(frequency) && isempty(coherence)
    storageSpectrum = struct("frequency", NaN, "coherence", NaN);
    return
end

storageSpectrum = struct("frequency", frequency, "coherence", coherence);
end
