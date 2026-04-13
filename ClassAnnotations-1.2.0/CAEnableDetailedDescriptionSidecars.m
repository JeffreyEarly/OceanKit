function cleanup = CAEnableDetailedDescriptionSidecars()
% Enable markdown sidecar loading during annotation construction.
%
% Use this helper around documentation-building code that constructs
% `CAPropertyAnnotation` or `CAMethodAnnotation` instances and wants them
% to load same-named markdown sidecars relative to their defining call
% sites.
%
% ```matlab
% cleanup = CAEnableDetailedDescriptionSidecars();
% annotations = MyClass.classDefinedPropertyAnnotations();
% clear cleanup
% ```

stateKey = 'CAEnableDetailedDescriptionSidecars_isEnabled';
hadPreviousState = isappdata(0, stateKey);
if hadPreviousState
    previousValue = getappdata(0, stateKey);
else
    previousValue = false;
end

setappdata(0, stateKey, true);
cleanup = onCleanup(@() restoreDetailedDescriptionSidecarsState(stateKey, hadPreviousState, previousValue));
end

function restoreDetailedDescriptionSidecarsState(stateKey, hadPreviousState, previousValue)
if hadPreviousState
    setappdata(0, stateKey, previousValue);
elseif isappdata(0, stateKey)
    rmappdata(0, stateKey);
end
end
