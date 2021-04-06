cmake_minimum_required(VERSION 3.6.3)
project(testMetagraphInterface VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_VERBOSE_MAKEFILE ON)

# set warnings for everything
if(MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -Weffc++ -Woverloaded-virtual
        -Wuninitialized -Wmaybe-uninitialized -pedantic)# -Werror)
endif()

# tell cmake where to put the compiled binaries (usually in build/ dir)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# find libraries
# create interface library targets from other libs
add_library(catch2 INTERFACE)
target_include_directories(catch2 SYSTEM INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib/Catch2/single_include")
#add_library(cpp-json-outstream INTERFACE)
#target_include_directories(cpp-json-outstream SYSTEM INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib/cpp-json-outstream/include")
#add_library(cpp-progress-bar INTERFACE)
#target_include_directories(cpp-progress-bar SYSTEM INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib/cpp-progress-bar/include")
#add_library(cpp-timestep INTERFACE)
#target_include_directories(cpp-timestep SYSTEM INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib/cpp-timestep/include")
#add_library(cxx-prettyprint INTERFACE)
#target_include_directories(cxx-prettyprint SYSTEM INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib/cxx-prettyprint")
#add_library(hopscotch-map INTERFACE)
#target_include_directories(hopscotch-map SYSTEM INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib/hopscotch-map/include")
add_library(json INTERFACE)
target_include_directories(json SYSTEM INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib/json/single_include")

# make targets in subdirs visible
add_subdirectory(lib/mabl3)

# the main program
add_executable(testMetagraphInterface src/main.cpp
                                      src/testMetagraphInterface.cpp
                                      src/loadExpectedTestdata.h)

target_include_directories(testMetagraphInterface PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")


# set testfile location (sets -DTESTDATAPATH=\"/home/.../testdata\" that is then used in loadExpectedTestdata.h)
target_compile_definitions(testMetagraphInterface PRIVATE TESTDATAPATH="${CMAKE_CURRENT_SOURCE_DIR}/testdata")

# link sources
target_link_libraries(testMetagraphInterface PRIVATE metagraphInterface)

# link third party libraries
target_link_libraries(testMetagraphInterface PRIVATE catch2)
target_link_libraries(testMetagraphInterface PRIVATE json)

# create test files from python
set(PYSCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/testdata/generateTestdata.py")
set(PYARG1 "--metagraph")
set(PYARG2 "${METAGRAPH_DIR_PARENT}/build/metagraph")
message("[INFO] -- Running python3 ${PYSCRIPT} ${PYARG1} ${PYARG2} in working directory ${CMAKE_CURRENT_SOURCE_DIR}/testdata/")
execute_process(COMMAND "python3" ${PYSCRIPT} ${PYARG1} ${PYARG2}
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/testdata/"
                OUTPUT_VARIABLE CMDOUT
                ERROR_VARIABLE CMDOUT)
message(${CMDOUT})
