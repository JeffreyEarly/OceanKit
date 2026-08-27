cmake_minimum_required(VERSION 3.20)

if(NOT DEFINED WV_REPOSITORY_ROOT)
    message(FATAL_ERROR "WV_REPOSITORY_ROOT is required.")
endif()

function(require_token source token message_text)
    string(FIND "${source}" "${token}" position)
    if(position EQUAL -1)
        message(FATAL_ERROR "${message_text}")
    endif()
endfunction()

function(forbid_token source token message_text)
    string(FIND "${source}" "${token}" position)
    if(NOT position EQUAL -1)
        message(FATAL_ERROR "${message_text}: ${token}")
    endif()
endfunction()

file(READ
    "${WV_REPOSITORY_ROOT}/PortableRuntime/include/WaveVortexRuntime/WVOutputSchedule.hpp"
    schedule_header)
foreach(token
        "WVMaximumOutputSchedulePayloadBytes = 4096"
        "WVOutputSchedulePayloadSchema"
        "std::array<std::uint8_t, WVMaximumOutputSchedulePayloadBytes>"
        "valueFingerprint() const noexcept"
        "virtual const WVOutputSchedulePayloadSchema &payloadSchema()")
    require_token("${schedule_header}" "${token}"
        "The schedule occurrence payload lost a required compact resolved contract")
endforeach()
string(FIND "${schedule_header}" "enum class WVOutputSchedulePayloadType" payload_type_begin)
string(FIND "${schedule_header}" "struct WVOutputSchedulePayloadField" payload_type_end)
if(payload_type_begin EQUAL -1 OR payload_type_end EQUAL -1 OR
   payload_type_end LESS payload_type_begin)
    message(FATAL_ERROR "The schedule occurrence payload type declaration is malformed.")
endif()
math(EXPR payload_type_length "${payload_type_end} - ${payload_type_begin}")
string(SUBSTRING "${schedule_header}" ${payload_type_begin}
       ${payload_type_length} payload_types)
forbid_token("${payload_types}" "text"
    "Schedule occurrence payloads must remain numeric, integer, or Boolean")

file(READ
    "${WV_REPOSITORY_ROOT}/PortableRuntime/include/WaveVortexRuntime/WVObserverOutputProvider.hpp"
    provider_header)
foreach(token
        "WVObserverOccurrencePreparationContext"
        "WVObserverOccurrenceWorkspace"
        "occurrencePayloadSchema"
        "occurrenceStateBlocks"
        "occurrencePositionSets"
        "occurrenceValues"
        "resolvedAdditionalStateBlockIndex"
        "resolvedXValueSlot"
        "resolvedYValueSlot"
        "value(std::size_t resolvedValueSlot"
        "scheduleOrdinal() const noexcept")
    require_token("${provider_header}" "${token}"
        "The observer occurrence provider lost a construction-resolved API")
endforeach()
forbid_token("${provider_header}" "value(const std::string &variableIdentifier"
    "The event evaluation context restored name-based value lookup")
forbid_token("${provider_header}" "eventOrdinal() const noexcept"
    "Observer providers must use the schedule occurrence ordinal, not the driver-global event ordinal")

file(READ
    "${WV_REPOSITORY_ROOT}/PortableRuntime/include/WaveVortexRuntime/WVFieldEvaluationService.hpp"
    field_header)
foreach(token
        "WVEventFieldEvaluationPlan"
        "WVEventFieldEvaluationBatchEntry"
        "WVPreparedFieldGeometry"
        "createEventPlan("
        "prepareEventGeometry("
        "evaluateEvent("
        "evaluateEventBatch(")
    require_token("${field_header}" "${token}"
        "The central field service lost an event-geometry stage")
endforeach()
string(FIND "${field_header}" "struct WVEventFieldRequest" event_field_begin)
string(FIND "${field_header}" "struct WVFieldEvaluationMetrics" event_field_end)
if(event_field_begin EQUAL -1 OR event_field_end EQUAL -1 OR
   event_field_end LESS event_field_begin)
    message(FATAL_ERROR "The event field contract declaration is malformed.")
endif()
math(EXPR event_field_length "${event_field_end} - ${event_field_begin}")
string(SUBSTRING "${field_header}" ${event_field_begin}
       ${event_field_length} event_field_contract)
forbid_token("${event_field_contract}" "virtual"
    "Event geometry and field operations must use resolved concrete plans")

file(READ
    "${WV_REPOSITORY_ROOT}/PortableRuntime/include/WaveVortexRuntime/WVModelOutputNetCDF.hpp"
    sample_source_header)
foreach(token
        "WVObservationOccurrenceIdentity"
        "preparationOwner"
        "preparationGeneration"
        "preparedOccurrenceSlot"
        "resolvedObserverRecord"
        "logicalScheduleRecord"
        "schedulePayloadSchema"
        "proposedScheduleCursor"
        "resolvedSchedulePayload"
        "samePreparedObservationOccurrenceIdentity("
        "preparedOccurrenceIdentity("
        "const WVObservationOccurrenceIdentity &identity")
    require_token("${sample_source_header}" "${token}"
        "The sink/source boundary lost destination-independent occurrence identity")
endforeach()
forbid_token("${sample_source_header}"
    "observationBatch(const WVObserverRecord &observer"
    "The event batch API restored observer-only cache identity")

file(READ
    "${WV_REPOSITORY_ROOT}/PortableRuntime/src/WVObserverSampleSource.cpp"
    sample_source_body)
string(FIND "${sample_source_body}"
    "bool sameObservationOccurrenceIdentity("
    semantic_identity_begin)
string(FIND "${sample_source_body}"
    "bool samePreparedObservationOccurrenceIdentity("
    semantic_identity_end)
if(semantic_identity_begin EQUAL -1 OR semantic_identity_end EQUAL -1 OR
   semantic_identity_end LESS semantic_identity_begin)
    message(FATAL_ERROR "The exact semantic occurrence comparison body is unavailable for policy inspection.")
endif()
math(EXPR semantic_identity_length
    "${semantic_identity_end} - ${semantic_identity_begin}")
string(SUBSTRING "${sample_source_body}" ${semantic_identity_begin}
       ${semantic_identity_length} semantic_identity)
foreach(token
        "sameOutputObserverSemanticIdentity("
        "sameLogicalOutputScheduleIdentity("
        "sameOutputSchedulePayloadSchema("
        "samePortableTypedRecordValue("
        "resolvedSchedulePayload->sameValue(")
    require_token("${semantic_identity}" "${token}"
        "Semantic occurrence identity lost an exact construction or event value comparison")
endforeach()
foreach(token "observerOrdinal" "semanticScheduleOrdinal"
              "scheduleCursorIdentity" "payloadFingerprint"
              "geometryFingerprint" "fieldPlanFingerprint"
              "preparationOwner" "preparationGeneration"
              "preparedOccurrenceSlot"
              "std::find_if" ".find(" "std::map" "std::unordered_map")
    forbid_token("${semantic_identity}" "${token}"
        "Exact semantic occurrence comparison restored a plan-local, fingerprint-only, or named identity")
endforeach()

file(READ
    "${WV_REPOSITORY_ROOT}/PortableRuntime/src/WVObserverOutputEvaluationService.cpp"
    evaluator_source)
string(FIND "${evaluator_source}"
    "WVKernelStatus WVObserverOutputEvaluationService::prepare("
    evaluator_prepare_begin)
string(FIND "${evaluator_source}"
    "void WVObserverOutputEvaluationService::complete("
    evaluator_prepare_end)
if(evaluator_prepare_begin EQUAL -1 OR evaluator_prepare_end EQUAL -1 OR
   evaluator_prepare_end LESS evaluator_prepare_begin)
    message(FATAL_ERROR "The observer occurrence preparation body is unavailable for policy inspection.")
endif()
math(EXPR evaluator_prepare_length
    "${evaluator_prepare_end} - ${evaluator_prepare_begin}")
string(SUBSTRING "${evaluator_source}" ${evaluator_prepare_begin}
       ${evaluator_prepare_length} evaluator_prepare)
foreach(token "std::find_if" ".find(" "std::map" "std::unordered_map")
    forbid_token("${evaluator_prepare}" "${token}"
        "Observer occurrence preparation performs event-loop name lookup")
endforeach()
foreach(token "sameOutputSchedulePayloadSchema(" "validateObservationBatch("
              "preparedOccurrenceOrdinal")
    forbid_token("${evaluator_prepare}" "${token}"
        "Observer occurrence preparation restored named/deep event-path work or unstable identity")
endforeach()
foreach(token
        "sameTime(event.state.waveVortex.t, event.scheduledTime)"
        "stateBlock.resolvedAdditionalStateBlockIndex"
        "occurrence.observerStateViews.data()")
    require_token("${evaluator_prepare}" "${token}"
        "Observer occurrence preparation lost a resolved trigger-state invariant")
endforeach()
foreach(token "samePortableTypedRecordValue(" "candidate.payload.sameValue(")
    require_token("${evaluator_source}" "${token}"
        "Prepared occurrence lookup lost exact cursor or payload comparison")
endforeach()
foreach(token
        "occurrence.identity.resolvedObserverRecord = binding.record"
        "occurrence.identity.logicalScheduleRecord ="
        "occurrence.identity.schedulePayloadSchema ="
        "occurrence.identity.proposedScheduleCursor ="
        "occurrence.identity.resolvedSchedulePayload =")
    require_token("${evaluator_source}" "${token}"
        "Prepared occurrences lost an exact borrowed semantic identity view")
endforeach()
require_token("${evaluator_source}"
    "if (identity == nullptr) {\n    status = validateObservationBatch"
    "The legacy named-schema validator is no longer isolated to initial delivery")

file(READ
    "${WV_REPOSITORY_ROOT}/PortableRuntime/src/WVObserverContracts.cpp"
    observer_contract_source)
string(FIND "${observer_contract_source}"
    "WVKernelStatus WVObservingSystem::observationBatch("
    observer_batch_begin)
string(FIND "${observer_contract_source}"
    "WVKernelStatus WVResolvedObserver::outputPlan("
    observer_batch_end)
if(observer_batch_begin EQUAL -1 OR observer_batch_end EQUAL -1 OR
   observer_batch_end LESS observer_batch_begin)
    message(FATAL_ERROR "The default observer batch body is unavailable for policy inspection.")
endif()
math(EXPR observer_batch_length
    "${observer_batch_end} - ${observer_batch_begin}")
string(SUBSTRING "${observer_contract_source}" ${observer_batch_begin}
       ${observer_batch_length} observer_batch)
foreach(token "std::find_if" ".find(" "std::map" "std::unordered_map"
              "observationVariable(")
    forbid_token("${observer_batch}" "${token}"
        "Default occurrence batch assembly performs event-loop name lookup")
endforeach()
require_token("${observer_batch}"
    "if (kind == WVObservationBatchKind::initial)"
    "The legacy named batch validator is no longer isolated from event delivery")

file(READ
    "${WV_REPOSITORY_ROOT}/PortableRuntime/src/WVModelOutputNetCDFWriter.cpp"
    sink_source)
string(FIND "${sink_source}"
    "validateResolvedObservationBatch("
    sink_validate_begin)
string(FIND "${sink_source}"
    "WVCheckpointStatus prepareObservationBatches("
    sink_validate_end)
if(sink_validate_begin EQUAL -1 OR sink_validate_end EQUAL -1 OR
   sink_validate_end LESS sink_validate_begin)
    message(FATAL_ERROR "The resolved sink batch validator is unavailable for policy inspection.")
endif()
math(EXPR sink_validate_length
    "${sink_validate_end} - ${sink_validate_begin}")
string(SUBSTRING "${sink_source}" ${sink_validate_begin}
       ${sink_validate_length} sink_validate)
foreach(token "std::find_if" ".find(" "std::map" "std::unordered_map"
              "variableIdentifier" "schemaIdentifier")
    forbid_token("${sink_validate}" "${token}"
        "Resolved NetCDF event validation performs event-loop name lookup")
endforeach()
require_token("${sink_validate}" "value.resolvedVariableIndex"
    "Resolved NetCDF event validation lost numeric variable-slot identity")

string(FIND "${sink_source}"
    "WVCheckpointStatus prepareObservationBatches("
    sink_prepare_begin)
string(FIND "${sink_source}" "void completePreparedRoute("
    sink_prepare_end)
if(sink_prepare_begin EQUAL -1 OR sink_prepare_end EQUAL -1 OR
   sink_prepare_end LESS sink_prepare_begin)
    message(FATAL_ERROR "The sink occurrence preparation body is unavailable for policy inspection.")
endif()
math(EXPR sink_prepare_length "${sink_prepare_end} - ${sink_prepare_begin}")
string(SUBSTRING "${sink_source}" ${sink_prepare_begin}
       ${sink_prepare_length} sink_prepare)

string(FIND "${sink_source}" "WVCheckpointStatus writeRoute("
    sink_write_begin)
string(FIND "${sink_source}" "WVCheckpointStatus close() noexcept"
    sink_write_end)
if(sink_write_begin EQUAL -1 OR sink_write_end EQUAL -1 OR
   sink_write_end LESS sink_write_begin)
    message(FATAL_ERROR "The sink event-write body is unavailable for policy inspection.")
endif()
math(EXPR sink_write_length "${sink_write_end} - ${sink_write_begin}")
string(SUBSTRING "${sink_source}" ${sink_write_begin}
       ${sink_write_length} sink_write)
foreach(event_path IN ITEMS sink_prepare sink_write)
    foreach(token "std::find_if" ".find(" "std::map" "std::unordered_map")
        forbid_token("${${event_path}}" "${token}"
            "NetCDF occurrence persistence performs event-loop name lookup")
    endforeach()
endforeach()
require_token("${sink_prepare}"
    "samePreparedObservationOccurrenceIdentity("
    "The NetCDF event cache is not keyed by the exact prepared-occurrence token")
forbid_token("${sink_prepare}"
    "sameObservationOccurrenceIdentity("
    "The NetCDF event cache must use authoritative prepared tokens, not cross-source semantic equality")

file(READ
    "${WV_REPOSITORY_ROOT}/PortableRuntime/src/WVModelOutputConfiguration.cpp"
    output_configuration_source)
string(FIND "${output_configuration_source}"
    "WVCheckpointStatus WVModelOutputConfiguration::openNetCDFSink("
    output_open_begin)
string(FIND "${output_configuration_source}"
    "WVModelOutputConfiguration::descriptor() const noexcept"
    output_open_end)
if(output_open_begin EQUAL -1 OR output_open_end EQUAL -1 OR
   output_open_end LESS output_open_begin)
    message(FATAL_ERROR "The NetCDF output-construction body is unavailable for policy inspection.")
endif()
math(EXPR output_open_length "${output_open_end} - ${output_open_begin}")
string(SUBSTRING "${output_configuration_source}" ${output_open_begin}
       ${output_open_length} output_open)
forbid_token("${output_open}" "sampleSource->preflight("
    "The configuration facade must delegate source preflight to every public NetCDF factory")
string(REGEX MATCHALL "impl_->descriptor,[ \t\r\n]*impl_->plan"
    output_plan_forwarding "${output_open}")
list(LENGTH output_plan_forwarding output_plan_forwarding_count)
if(NOT output_plan_forwarding_count EQUAL 3)
    message(FATAL_ERROR
        "The configuration facade must pass its compiled plan to create, replace, and append factories.")
endif()

string(FIND "${sink_source}"
    "WVCheckpointStatus WVModelOutputNetCDFSink::createNew("
    sink_create_begin)
string(FIND "${sink_source}"
    "WVCheckpointStatus WVModelOutputNetCDFSink::replaceExisting("
    sink_replace_begin)
string(FIND "${sink_source}"
    "WVCheckpointStatus WVModelOutputNetCDFSink::openAppend("
    sink_append_begin)
string(FIND "${sink_source}"
    "WVKernelStatus WVModelOutputNetCDFSink::preflight("
    sink_preflight_begin)
if(sink_create_begin EQUAL -1 OR sink_replace_begin EQUAL -1 OR
   sink_append_begin EQUAL -1 OR sink_preflight_begin EQUAL -1 OR
   sink_replace_begin LESS sink_create_begin OR
   sink_append_begin LESS sink_replace_begin OR
   sink_preflight_begin LESS sink_append_begin)
    message(FATAL_ERROR "The public NetCDF factory bodies are unavailable for policy inspection.")
endif()
math(EXPR sink_create_length "${sink_replace_begin} - ${sink_create_begin}")
math(EXPR sink_replace_length "${sink_append_begin} - ${sink_replace_begin}")
math(EXPR sink_append_length "${sink_preflight_begin} - ${sink_append_begin}")
string(SUBSTRING "${sink_source}" ${sink_create_begin}
       ${sink_create_length} sink_create)
string(SUBSTRING "${sink_source}" ${sink_replace_begin}
       ${sink_replace_length} sink_replace)
string(SUBSTRING "${sink_source}" ${sink_append_begin}
       ${sink_append_length} sink_append)
foreach(factory_body IN ITEMS sink_create sink_replace sink_append)
    string(FIND "${${factory_body}}" "sampleSource->preflight(plan)"
        factory_source_preflight_position)
    string(FIND "${${factory_body}}" "candidate->validateConfiguration()"
        factory_schema_discovery_position)
    string(FIND "${${factory_body}}"
        "candidate->validateAndCompilePlan(plan)"
        factory_plan_compilation_position)
    string(FIND "${${factory_body}}" "candidate->validateDestinations("
        factory_destination_validation_position)
    if(factory_source_preflight_position EQUAL -1 OR
       factory_schema_discovery_position EQUAL -1 OR
       factory_source_preflight_position GREATER
           factory_schema_discovery_position)
        message(FATAL_ERROR
            "Every public NetCDF factory must preflight the observer source before schema discovery or destination mutation.")
    endif()
    if(factory_plan_compilation_position EQUAL -1 OR
       factory_destination_validation_position EQUAL -1 OR
       factory_plan_compilation_position GREATER
           factory_destination_validation_position)
        message(FATAL_ERROR
            "Every public NetCDF factory must validate and numerically compile its plan before destination inspection or mutation.")
    endif()
endforeach()
foreach(token
        "WVKernelStatus validateAndCompilePlan(const WVOutputPlan &plan)"
        "plan.metrics().fileCount != outputFiles.size()"
        "plan.groupCount() != expectedGroupCount"
        "compileEventBindings(transient, route)"
        "files[route.fileOrdinal].groups[route.groupOrdinal], route)")
    require_token("${sink_source}" "${token}"
        "The NetCDF factory lost full graph validation or numeric event-binding compilation")
endforeach()

set(observer_implementation_paths
    "PortableRuntime/include/WaveVortexRuntime/WVObservingSystem.hpp"
    "PortableRuntime/include/WaveVortexRuntime/WVObserverOutputProvider.hpp"
    "PortableRuntime/src/WVObserverAdapter.cpp"
    "PortableRuntime/src/WVObserverContracts.cpp")
foreach(relative_path IN LISTS observer_implementation_paths)
    file(READ "${WV_REPOSITORY_ROOT}/${relative_path}" source)
    foreach(token
            "netcdf.h" "nc_open(" "nc_create(" "nc_def_" "nc_put_" "nc_get_"
            "createEventPlan(" "prepareEventGeometry(" "evaluateEvent("
            "evaluateMoving(")
        forbid_token("${source}" "${token}"
            "Observer implementations must not own NetCDF or private field evaluation in ${relative_path}")
    endforeach()
endforeach()

file(GLOB_RECURSE runtime_production_sources LIST_DIRECTORIES false
    "${WV_REPOSITORY_ROOT}/PortableRuntime/include/*.hpp"
    "${WV_REPOSITORY_ROOT}/PortableRuntime/src/*.cpp"
    "${WV_REPOSITORY_ROOT}/PortableRuntime/src/*.hpp"
    "${WV_REPOSITORY_ROOT}/PortableRuntime/app/*.cpp"
    "${WV_REPOSITORY_ROOT}/PortableRuntime/app/*.hpp")
foreach(source_path IN LISTS runtime_production_sources)
    file(READ "${source_path}" source)
    string(TOLOWER "${source}" source_lower)
    foreach(token
            "alongtrack" "along-track" "along_track" "adcp" "glider"
            "argo" "satellite" "shipboard")
        forbid_token("${source_lower}" "${token}"
            "Portable runtime contains deferred scientific-observer policy in ${source_path}")
    endforeach()
endforeach()

message(STATUS "Observation occurrence source policy passed.")
