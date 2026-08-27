function u = transformToSpatialDomainWithFourier(self,u_bar)
layout = self.fourierStorageLayout;
self.complexBufferRows(layout.fourierRowsForDirectWVIndices,:) = u_bar(:,layout.directWVIndices).';
if ~isempty(layout.hermitianCompletionRows)
    self.complexBufferRows(layout.hermitianCompletionRows,:) = conj(u_bar(:,layout.hermitianSourceWVIndices).');
end
if ~isempty(layout.selfConjugateFourierRows)
    self.complexBufferRows(layout.selfConjugateFourierRows,:) = real(self.complexBufferRows(layout.selfConjugateFourierRows,:));
end
complexBuffer = layout.reshapeFourierRowsToStorage(self.complexBufferRows);
if self.wvg.conjugateDimension == 1
    u = ifft(ifft(complexBuffer,self.wvg.Ny,2),self.wvg.Nx,1,'symmetric')*(self.wvg.Nx*self.wvg.Ny);
else
    u = ifft(ifft(complexBuffer,self.wvg.Nx,1),self.wvg.Ny,2,'symmetric')*(self.wvg.Nx*self.wvg.Ny);
end
end
