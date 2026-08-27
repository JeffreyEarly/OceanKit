function seedRandomNumberGenerator(testCase,seed)
% Seed the global random stream and restore its previous state after the test.

arguments
    testCase (1,1) matlab.unittest.TestCase
    seed (1,1) double {mustBeInteger,mustBeNonnegative}
end

previousRandomState = rng;
testCase.addTeardown(@()rng(previousRandomState));
rng(seed,"twister")
end
