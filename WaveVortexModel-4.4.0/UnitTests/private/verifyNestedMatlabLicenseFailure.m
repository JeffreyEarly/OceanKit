function verifyNestedMatlabLicenseFailure(testCase,failure,expectedIdentifier)
% Verify the expected nested-MATLAB licensing failure on GitHub Actions.
arguments
    testCase (1,1) matlab.unittest.TestCase
    failure (1,1) struct
    expectedIdentifier (1,1) string
end

diagnostic = "Nested MATLAB worker failed: " + string(failure.identifier) + ": " + string(failure.message);
testCase.verifyEqual(string(failure.identifier),expectedIdentifier,diagnostic);
licenseMessages = ["License checkout failed" "Unable to find a license for MATLAB"];
testCase.verifyTrue(any(contains(string(failure.message),licenseMessages)),diagnostic);
end
