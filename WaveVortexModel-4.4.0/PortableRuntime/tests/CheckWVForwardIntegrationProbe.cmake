if(NOT DEFINED PROBE OR NOT DEFINED FIXTURE OR NOT DEFINED WORK)
    message(FATAL_ERROR "Probe contract requires executable, fixture and work paths")
endif()
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")

foreach(invalid_count IN ITEMS 0 -1 1x 999999999999999999999999999999999999999999)
    execute_process(COMMAND "${PROBE}" --stop-boundary accepted-step
        --stop-after "${invalid_count}" --request missing.json
        RESULT_VARIABLE status OUTPUT_QUIET ERROR_QUIET)
    if(NOT status EQUAL 2)
        message(FATAL_ERROR "Invalid stop count ${invalid_count} was not rejected")
    endif()
endforeach()
execute_process(COMMAND "${PROBE}" --stop-boundary rhs --stop-after 1
    --request missing.json RESULT_VARIABLE status OUTPUT_QUIET ERROR_QUIET)
if(NOT status EQUAL 2)
    message(FATAL_ERROR "Unknown stop boundary was not rejected")
endif()

foreach(boundary IN ITEMS accepted-step output-occurrence)
    set(output "${WORK}/${boundary}.nc")
    set(report "${WORK}/${boundary}.json")
    execute_process(COMMAND "${PROBE}" --stop-boundary "${boundary}"
        --stop-after 1 "${FIXTURE}" "${output}" --restart-mode coefficients
        --output-policy create --steps 4 --delta-t 1e-5 --fft-provider reference
        --benchmark-dense-outputs-per-step 1 --report "${report}"
        RESULT_VARIABLE status OUTPUT_QUIET ERROR_VARIABLE errors)
    if(NOT status EQUAL 0 OR NOT EXISTS "${output}" OR NOT EXISTS "${report}")
        message(FATAL_ERROR "Source-linked probe failed ${boundary}: ${errors}")
    endif()
    file(READ "${report}" json)
    string(JSON run_status GET "${json}" status)
    string(JSON requested_boundary GET "${json}" termination requestedBoundary)
    string(JSON steps GET "${json}" state stepCount)
    string(JSON final_time GET "${json}" state finalTime)
    string(JSON accepted_time GET "${json}" termination finalAcceptedTime)
    string(JSON request_time GET "${json}" termination requestedAtAcceptedTime)
    if(NOT run_status STREQUAL "stopped" OR
       NOT requested_boundary STREQUAL boundary OR NOT steps EQUAL 1 OR
       NOT final_time EQUAL accepted_time OR NOT final_time EQUAL request_time)
        message(FATAL_ERROR "Probe did not retain the selected accepted boundary")
    endif()
endforeach()

# Positive counts beyond the run's available boundaries must complete normally;
# the probe is a one-shot callback, not a change to runner integration policy.
execute_process(COMMAND "${PROBE}" --stop-boundary accepted-step --stop-after 9
    "${FIXTURE}" "${WORK}/complete.nc" --restart-mode coefficients
    --output-policy create --steps 2 --delta-t 1e-5 --fft-provider reference
    --report "${WORK}/complete.json" RESULT_VARIABLE status
    OUTPUT_QUIET ERROR_VARIABLE errors)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "Continuing source-linked probe failed: ${errors}")
endif()
file(READ "${WORK}/complete.json" json)
string(JSON run_status GET "${json}" status)
string(JSON steps GET "${json}" state stepCount)
if(NOT run_status STREQUAL "complete" OR NOT steps EQUAL 2)
    message(FATAL_ERROR "Unreached probe stop must preserve complete integration")
endif()
message(STATUS "Forward integration probe source/CLI contract passed")
