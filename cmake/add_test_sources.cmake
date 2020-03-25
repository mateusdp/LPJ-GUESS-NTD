# By including this file all files in ${headers} and ${source} is
# added to test_sources.
foreach(file ${headers} ${source})
  list(APPEND these_sources ${CMAKE_CURRENT_SOURCE_DIR}/${file})
endforeach(file)

set(test_sources ${test_sources} ${these_sources} PARENT_SCOPE)
