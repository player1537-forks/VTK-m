##============================================================================
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##============================================================================

#-----------------------------------------------------------------------------
# Adds a performance benchmark test
#
# Usage:
#   add_benchmark_test(benchmark
#     [ NAME <name> ]
#     [ ARGS <args...> ]
#     [ REGEX <benchmark_regex...> ]
#     )
#
# benchmark:    Target of an executable that uses Google Benchmark.
#
# NAME:         The name given to the CMake tests. The benchmark target name is used
#               if NAME is not specified.
#
# ARGS:         Extra arguments passed to the benchmark executable when run.
#
# REGEX:        Regular expressions that select the specific benchmarks within the binary
#               to be used. It populates the Google Benchmark
#               --benchmark_filter parameter. When multiple regexes are passed
#               as independent positional arguments, they are joined using the "|"
#               regex operator before populating the  `--benchmark_filter` parameter.
#
function(add_benchmark_test benchmark)

  # We need JSON support among other things for this to work
  if (CMAKE_VERSION VERSION_LESS 3.19)
    message(FATAL_ERROR "Performance regression testing needs CMAKE >= 3.19")
    return()
  endif()

  ###TEST VARIABLES############################################################

  set(options)
  set(one_value_keywords NAME)
  set(multi_value_keywords ARGS REGEX)
  cmake_parse_arguments(PARSE_ARGV 1 VTKm_PERF "${options}" "${one_value_keywords}" "${multi_value_keywords}")
  if (VTKm_PERF_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Bad arguments to add_benchmark_test (${VTKm_PERF_UNPARSED_ARGUMENTS}).")
  endif()

  if (NOT VTKm_PERF_NAME)
    set(VTKm_PERF_NAME ${benchmark})
  endif()

  if (VTKm_PERF_REGEX)
    string(REPLACE ";" "|" VTKm_PERF_REGEX "${VTKm_PERF_REGEX}")
  else()
    set(VTKm_PERF_REGEX ".*")
  endif()

  set(VTKm_PERF_REMOTE_URL "https://gitlab.kitware.com/vbolea/vtk-m-benchmark-records.git")

  # Parameters for the null hypothesis test
  set(VTKm_PERF_ALPHA  0.05)
  set(VTKm_PERF_REPETITIONS 10)
  set(VTKm_PERF_MIN_TIME 1)
  set(VTKm_PERF_DIST "normal")

  set(VTKm_PERF_REPO           "${CMAKE_BINARY_DIR}/vtk-m-benchmark-records")
  set(VTKm_PERF_COMPARE_JSON   "${CMAKE_BINARY_DIR}/nocommit_${VTKm_PERF_NAME}.json")
  set(VTKm_PERF_STDOUT         "${CMAKE_BINARY_DIR}/benchmark_${VTKm_PERF_NAME}.stdout")
  set(VTKm_PERF_COMPARE_STDOUT "${CMAKE_BINARY_DIR}/compare_${VTKm_PERF_NAME}.stdout")

  if (DEFINED ENV{CI_COMMIT_SHA})
    set(VTKm_PERF_COMPARE_JSON "${CMAKE_BINARY_DIR}/$ENV{CI_COMMIT_SHA}_${VTKm_PERF_NAME}.json")
  endif()

  set(test_name "PerformanceTest${VTKm_PERF_NAME}")

  ###TEST INVOKATIONS##########################################################
  add_test(NAME "${test_name}Run"
    COMMAND ${CMAKE_COMMAND}
    "-DVTKm_PERF_BENCH_DEVICE=Any"
    "-DVTKm_PERF_BENCH_PATH=${CMAKE_BINARY_DIR}/bin/${benchmark}"
    "-DVTKm_PERF_ARGS=${VTKm_PERF_ARGS}"
    "-DVTKm_PERF_REGEX=${VTKm_PERF_REGEX}"
    "-DVTKm_PERF_REPETITIONS=${VTKm_PERF_REPETITIONS}"
    "-DVTKm_PERF_MIN_TIME=${VTKm_PERF_MIN_TIME}"
    "-DVTKm_PERF_COMPARE_JSON=${VTKm_PERF_COMPARE_JSON}"
    "-DVTKm_PERF_STDOUT=${VTKm_PERF_STDOUT}"
    "-DVTKm_SOURCE_DIR=${VTKm_SOURCE_DIR}"
    -P "${VTKm_SOURCE_DIR}/CMake/testing/VTKmPerformanceTestRun.cmake"
    )

  add_test(NAME "${test_name}Fetch"
    COMMAND ${CMAKE_COMMAND}
    "-DVTKm_PERF_REPO=${VTKm_PERF_REPO}"
    "-DVTKm_SOURCE_DIR=${VTKm_SOURCE_DIR}"
    "-DVTKm_PERF_REMOTE_URL=${VTKm_PERF_REMOTE_URL}"
    -P "${VTKm_SOURCE_DIR}/CMake/testing/VTKmPerformanceTestFetch.cmake"
    )

  add_test(NAME "${test_name}Upload"
    COMMAND ${CMAKE_COMMAND}
    "-DVTKm_PERF_REPO=${VTKm_PERF_REPO}"
    "-DVTKm_PERF_COMPARE_JSON=${VTKm_PERF_COMPARE_JSON}"
    "-DVTKm_SOURCE_DIR=${VTKm_SOURCE_DIR}"
    -P "${VTKm_SOURCE_DIR}/CMake/testing/VTKmPerformanceTestUpload.cmake"
    )

  add_test(NAME "${test_name}Report"
    COMMAND ${CMAKE_COMMAND}
    "-DBENCHMARK_NAME=${benchmark}"
    "-DVTKm_PERF_ALPHA=${VTKm_PERF_ALPHA}"
    "-DVTKm_PERF_DIST=${VTKm_PERF_DIST}"
    "-DVTKm_PERF_REPO=${VTKm_PERF_REPO}"
    "-DVTKm_PERF_COMPARE_JSON=${VTKm_PERF_COMPARE_JSON}"
    "-DVTKm_SOURCE_DIR=${VTKm_SOURCE_DIR}"
    "-DVTKm_BINARY_DIR=${VTKm_BINARY_DIR}"
    "-DVTKm_PERF_COMPARE_STDOUT=${VTKm_PERF_COMPARE_STDOUT}"
    -P "${VTKm_SOURCE_DIR}/CMake/testing/VTKmPerformanceTestReport.cmake"
    )

  add_test(NAME "${test_name}CleanUp"
    COMMAND ${CMAKE_COMMAND} -E rm -rf "${VTKm_PERF_REPO}"
    )

  ###TEST PROPERTIES###########################################################
  set_tests_properties("${test_name}Report" "${test_name}Upload"
    PROPERTIES
    FIXTURE_REQUIRED "${test_name}Run;${test_name}Fetch"
    FIXTURE_CLEANUP  "${test_name}CleanUp"
    REQUIRED_FILES "${VTKm_PERF_COMPARE_JSON}")

  set_tests_properties("${test_name}Run"
                        "${test_name}Report"
                        "${test_name}Upload"
                        "${test_name}Fetch"
                        "${test_name}CleanUp"
                        PROPERTIES RUN_SERIAL ON)

  set_tests_properties(${test_name}Run PROPERTIES TIMEOUT 1800)

  # Only upload when we are inside a CI build
  if (NOT DEFINED ENV{CI_COMMIT_SHA} OR NOT DEFINED ENV{VTKM_BENCH_RECORDS_TOKEN})
    set_tests_properties(${test_name}Upload PROPERTIES DISABLED TRUE)
  endif()
endfunction()
