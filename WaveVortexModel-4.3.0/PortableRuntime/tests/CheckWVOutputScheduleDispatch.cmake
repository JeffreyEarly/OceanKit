file(GLOB runtime_sources "${WV_REPOSITORY_ROOT}/PortableRuntime/src/*.cpp")
foreach(source IN LISTS runtime_sources)
  file(READ "${source}" contents)
  if(NOT source MATCHES "WVOutputOrchestration.cpp$" AND
     contents MATCHES "\\.eventCount\\(|\\.event\\(")
    message(FATAL_ERROR "Production source enumerates output-plan events: ${source}")
  endif()
  if(NOT source MATCHES "WVOutputSchedule.cpp$|WVExtensionCatalog.cpp$" AND
     contents MATCHES "WVEvenlySpacedOutputSchedule|WVStateTriggeredOutputSchedule|WVTestQuadraticOutputSchedule")
    message(FATAL_ERROR "Production source dispatches on a schedule provider identity: ${source}")
  endif()
  if(contents MATCHES "typeIdentifier[ \t\r\n]*[!=]=[ \t\r\n]*(WVEvenlySpacedOutputScheduleType|WVStateTriggeredOutputScheduleType)" OR
     contents MATCHES "(WVEvenlySpacedOutputScheduleType|WVStateTriggeredOutputScheduleType)[ \t\r\n]*[!=]=[ \t\r\n]*[^;\n]*typeIdentifier")
    message(FATAL_ERROR "Production source compares a schedule provider identity: ${source}")
  endif()
endforeach()
