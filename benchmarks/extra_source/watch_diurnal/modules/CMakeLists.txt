# Include all .h and .cpp files. Using wildcards like this goes
# against cmake best practices (you should list each file individually).
# However, doing it like this means we don't need to update this file
# everytime a new file is added in modules.

file(GLOB headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.h")

file(GLOB source RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")

include(add_guess_sources)
