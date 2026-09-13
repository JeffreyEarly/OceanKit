function writePortableRunRequest(path,modelFiles,options)
% Write a portable-runtime request for a MATLAB-authored NetCDF bundle.
%
% The referenced NetCDF files remain the authoritative scientific model,
% including transform configuration, state, forcing, observers, schedules,
% and restart progress. This method writes only execution choices and file
% routing for `wave-vortex-run --request`.
% Run-request v2 is the default. Omitted execution controls select MATLAB's
% standard `ode78` configuration: relative tolerance `1e-3`, absolute
% tolerance scale `1e-6`, an initial step from CFL `0.5` after state
% restoration, and a maximum step equal to one tenth of the continuation
% interval. The standalone runtime selects native FFTW with an automatically
% bounded thread count. These defaults do not recover custom MATLAB-session
% settings that were never persisted.
%
% The metadata-only writer accepts supported
% `WVTransformConstantStratification` bundles with complete `Ap`, `Am`, and
% `A0` restart state and `WVTransformBarotropicQG` bundles with one compact
% `A0` stream declared by `WVCoefficients`. An Eulerian field named `A0`
% does not by itself make a Barotropic QG file restart-capable.
%
% Relative model, output, and report paths are interpreted relative to the
% request document. Output destinations are keyed by the stable identifiers
% returned in validation errors or stored as `portableFileIdentifier` in a
% portable-runtime-authored file.
%
% - Topic: Write model output
% - Declaration: WVModel.writePortableRunRequest(path,modelFiles,options)
% - Parameter path: destination JSON request path
% - Parameter modelFiles: ordered complete set of source NetCDF paths
% - Parameter options.schemaVersion: exact request schema version; default `2`
% - Parameter options.method: `fixed-rk4`, `adaptive-rk23`, `adaptive-rk45`, or `adaptive-rk78`; v2 default `adaptive-rk78`
% - Parameter options.finalTime: requested final integration time
% - Parameter options.initialStep: explicit RK4 step or adaptive initial step
% - Parameter options.cfl: CFL number for schema-v2 CFL-selected RK4
% - Parameter options.timeStepConstraint: `advective`, `oscillatory`, or `min`
% - Parameter options.maximumStep: adaptive maximum step
% - Parameter options.relativeTolerance: adaptive relative tolerance
% - Parameter options.absoluteToleranceScale: adaptive absolute-tolerance scale
% - Parameter options.outputPolicy: `create`, `replace`, or `append`; default `append`
% - Parameter options.destinations: string-to-string dictionary keyed by output-file identifier
% - Parameter options.fftProvider: `native-fftw` or `reference`; v2 default `native-fftw`
% - Parameter options.threads: positive execution thread count; v2 default automatically bounded hardware concurrency
% - Parameter options.reportPath: report path; default `<request-name>-report.json`

arguments
    path (1,1) string {mustBeNonzeroLengthText}
    modelFiles string {mustBeNonempty}
    options.schemaVersion (1,1) double {mustBeMember(options.schemaVersion,[1 2])} = 2
    options.method (1,1) string {mustBeMember(options.method,["","fixed-rk4","adaptive-rk23","adaptive-rk45","adaptive-rk78"])} = ""
    options.finalTime (1,1) double {mustBeFinite}
    options.initialStep (1,1) double = NaN
    options.cfl (1,1) double = NaN
    options.timeStepConstraint (1,1) string = ""
    options.maximumStep (1,1) double = NaN
    options.relativeTolerance (1,1) double = NaN
    options.absoluteToleranceScale (1,1) double = NaN
    options.outputPolicy (1,1) string {mustBeMember(options.outputPolicy,["create","replace","append"])} = "append"
    options.destinations (1,1) dictionary = configureDictionary("string","string")
    options.fftProvider (1,1) string {mustBeMember(options.fftProvider,["","native-fftw","reference"])} = ""
    options.threads (1,1) double = NaN
    options.reportPath (1,1) string = ""
end

configuration = struct(options);
WVInternal.writePortableRunRequest(path,modelFiles,configuration,"");
end
