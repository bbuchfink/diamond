SET(CMD "./diamond ${ARGS} -o ${NAME}.out")
#separate_arguments (SEP NATIVE_COMMAND PROGRAM SEPARATE_ARGS ${CMD})
separate_arguments (SEP NATIVE_COMMAND ${CMD})
execute_process(COMMAND ${SEP} RESULT_VARIABLE CMD_RESULT)
if(CMAKE_VERSION VERSION_LESS 3.14)
  # No compare_files --ignore-eol yet: fall back to an external diff.
  if(CMAKE_HOST_WIN32)
    execute_process(COMMAND busybox diff ${TEST_DIR}/${NAME}.out ${NAME}.out RESULT_VARIABLE DIFF_RESULT)
  else()
    execute_process(COMMAND diff ${TEST_DIR}/${NAME}.out ${NAME}.out RESULT_VARIABLE DIFF_RESULT)
  endif()
else()
  # Needs no external diff tool, and ignores CRLF/LF differences between the
  # checked-out reference and the freshly written output.
  execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol
    ${TEST_DIR}/${NAME}.out ${NAME}.out RESULT_VARIABLE DIFF_RESULT)
endif()
if(NOT ${DIFF_RESULT} EQUAL 0)
  message(FATAL_ERROR "${NAME} failed. Compare ${TEST_DIR}/${NAME}.out with ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.out")
endif()
