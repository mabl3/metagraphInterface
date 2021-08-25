project(metagraphInterface)

# find metagraph
set(METAGRAPH_PROJECT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../projects2014-metagenome" CACHE PATH "Metagraph project root directory")
message("[INFO] -- METAGRAPH_PROJECT_ROOT: ${METAGRAPH_PROJECT_ROOT}")
set(METAGRAPH_DIR "${METAGRAPH_PROJECT_ROOT}/metagraph")
set(METAGRAPH_DIR_PARENT "${METAGRAPH_PROJECT_ROOT}/metagraph" PARENT_SCOPE)
message("[INFO] -- METAGRAPH_DIR: ${METAGRAPH_DIR}")
message("[INFO] -- Searching Metagraph in ${METAGRAPH_DIR}")
if (${METAGRAPH_PROJECT_ROOT} STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}/../projects2014-metagenome")
    message("          If this is not successful, point to the right project directory with -DMETAGRAPH_PROJECT_ROOT")
endif()

# find metagraph library file
set(METAGRAPH_LIB METAGRAPH_LIB-NOTFOUND)
find_library(METAGRAPH_LIB
             NAMES metagraph-core
             PATHS "${METAGRAPH_DIR}/build"
             NO_DEFAULT_PATH)
if (METAGRAPH_LIB)
    message("[INFO] -- Found metagraph: ${METAGRAPH_LIB}")
else()
    message("[INFO] -- Metagraph not found")
endif()

# get graph type
file(READ "${METAGRAPH_DIR}/build/compile_commands.json" METAGRAPH_COMPILE_COMMAND_JSON)
string(REGEX MATCH "-D_[A-Z_]+_GRAPH" METAGRAPH_GRAPH_TYPE_FLAG ${METAGRAPH_COMPILE_COMMAND_JSON})
message("[INFO] -- METAGRAPH_GRAPH_TYPE_FLAG: ${METAGRAPH_GRAPH_TYPE_FLAG}")

find_library(LIBFBVECTOR fbvector PATHS "${METAGRAPH_DIR}/build" NO_DEFAULT_PATH)
if (NOT LIBFBVECTOR)
    message("[INFO] -- Library fbvector not found.")
else()
    message("[INFO] -- Found fbvector: ${LIBFBVECTOR}")
endif()
find_library(LIBMERSENNE mersenne_twister PATHS "${METAGRAPH_DIR}/build" NO_DEFAULT_PATH)
if (NOT LIBMERSENNE)
    message("[INFO] -- Library mersenne-twister not found.")
else()
    message("[INFO] -- Found mersenne-twister: ${LIBMERSENNE}")
endif()
find_library(LIBOMP OpenMP omp omp5 HINTS "/usr/lib/llvm-10/lib" "~/.linuxbrew/lib" "~/.linuxbrew/opt/libomp/lib")
if (NOT LIBOMP)
    message("[INFO] -- Library OpenMP not found.")
else()
    message("[INFO] -- Found OpenMP: ${LIBOMP}")
endif()
find_library(LIBPROGRESSBAR progress_bar PATHS "${METAGRAPH_DIR}/build" NO_DEFAULT_PATH)
if (NOT LIBPROGRESSBAR)
    message("[INFO] -- Library progress-bar not found.")
else()
    message("[INFO] -- Found progress-bar: ${LIBPROGRESSBAR}")
endif()
find_library(LIBSDSL sdsl PATHS "${METAGRAPH_DIR}/external-libraries/sdsl-lite/build/lib" NO_DEFAULT_PATH)
if (NOT LIBSDSL)
    message("[INFO] -- Library sdsl-lite not found.")
else()
    message("[INFO] -- Found sdsl-lite: ${LIBSDSL}")
endif()
find_library(LIBSPDLOG spdlog PATHS "${METAGRAPH_DIR}/build/external-libraries/spdlog" NO_DEFAULT_PATH)
if (NOT LIBSPDLOG)
    message("[INFO] -- Library spdlog not found.")
else()
    message("[INFO] -- Found spdlog: ${LIBSPDLOG}")
endif()

# g++ needs the flag -latomic but CMake does not find atomic for whatever reason.
# Thus, query the GCC major version to construct the correct hint path
string(REGEX MATCH "^[0-9]+" CXXMAJOR ${CMAKE_CXX_COMPILER_VERSION}) # CMAKE_CXX_COMPILER_VERSION is e.g. "10.3.0", we only need the "10"
find_library(ATOMIC atomic HINTS "/usr/lib/gcc/x86_64-linux-gnu/${CXXMAJOR}")
if (NOT ATOMIC)
    message("[INFO] -- Library atomic not found.")
else()
    message("[INFO] -- Found atomic: ${ATOMIC}")
endif()

# assemble metagraph target
add_library(metagraph INTERFACE)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(metagraph INTERFACE tbb)
target_link_libraries(metagraph INTERFACE Threads::Threads)

# jemalloc required for folly/fbvector
list(APPEND CMAKE_MODULE_PATH "${METAGRAPH_DIR}") # metagraph contains custom FindJemalloc.cmake
find_package(Jemalloc REQUIRED)
target_link_libraries(metagraph INTERFACE ${JEMALLOC_LIBRARIES})

target_include_directories(metagraph SYSTEM INTERFACE "${METAGRAPH_DIR}/src")
target_include_directories(metagraph SYSTEM INTERFACE "${METAGRAPH_DIR}/external-libraries")
target_include_directories(metagraph SYSTEM INTERFACE "${METAGRAPH_DIR}/external-libraries/DYNAMIC/include")
target_include_directories(metagraph SYSTEM INTERFACE "${METAGRAPH_DIR}/external-libraries/hopscotch-map/include")
target_include_directories(metagraph SYSTEM INTERFACE "${METAGRAPH_DIR}/external-libraries/ips4o")
target_include_directories(metagraph SYSTEM INTERFACE "${METAGRAPH_DIR}/external-libraries/rollinghashcpp")
target_include_directories(metagraph SYSTEM INTERFACE "${METAGRAPH_DIR}/external-libraries/spdlog/include")
target_include_directories(metagraph SYSTEM INTERFACE "${METAGRAPH_DIR}/external-libraries/sdsl-lite/include")

target_link_libraries(metagraph INTERFACE ${METAGRAPH_LIB})
target_link_libraries(metagraph INTERFACE ${LIBFBVECTOR})
target_link_libraries(metagraph INTERFACE ${LIBMERSENNE})
target_link_libraries(metagraph INTERFACE ${LIBOMP})
target_link_libraries(metagraph INTERFACE ${LIBPROGRESSBAR})
target_link_libraries(metagraph INTERFACE ${LIBSDSL})
target_link_libraries(metagraph INTERFACE ${LIBSPDLOG})
target_link_libraries(metagraph INTERFACE ${ATOMIC})

target_compile_definitions(metagraph INTERFACE ${METAGRAPH_GRAPH_TYPE_FLAG})

# add interface target
add_library(metagraphInterface STATIC MetagraphInterface.cpp)
target_link_libraries(metagraphInterface PRIVATE tbb)
target_include_directories(metagraphInterface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# set C++ standard
target_compile_features(metagraphInterface PUBLIC cxx_std_17)

# link metagraph stuff
target_link_libraries(metagraphInterface PUBLIC metagraph)
