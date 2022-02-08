cmake_minimum_required(VERSION 3.14)
project(project CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O1 -g -Wall -Wextra -mavx")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -Wall -Wextra -mavx -fsanitize=address -fno-omit-frame-pointer")

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -g -Wall -Wextra -mavx")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -g -Wall -Wextra -mavx")

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

message("Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

add_executable(serialbh
    utility.cpp serial.cpp interface.cpp main.cpp)
target_compile_definitions(serialbh
    PUBLIC BUILDSERIAL)

add_executable(ompbh
        utility.cpp visualizer.cpp serial.cpp interface.cpp omp_interface.cpp ompbh.cpp main.cpp)
target_compile_definitions(ompbh
        PUBLIC BUILDOMP)
target_link_libraries(ompbh PRIVATE OpenMP::OpenMP_CXX)

add_executable(mpibh
        utility.cpp visualizer.cpp serial.cpp interface.cpp mpi_interface.cpp mpibh.cpp main.cpp)
target_compile_definitions(mpibh
        PUBLIC BUILDMPI)
target_link_libraries(mpibh PRIVATE MPI::MPI_CXX)

add_executable(visualbh
        utility.cpp visualizer.cpp serial.cpp interface.cpp main.cpp)
target_compile_definitions(visualbh
        PUBLIC BUILDSERIAL VISUALIZATION CONVERT_PPM)

add_executable(test_omp test_omp.cpp ompbh.cpp utility.cpp serial.cpp)
target_link_libraries(test_omp gtest_main OpenMP::OpenMP_CXX)

add_executable(test_mpi test_mpi.cpp mpibh.cpp utility.cpp serial.cpp)
target_link_libraries(test_mpi gtest_main gmock MPI::MPI_CXX)

include(GoogleTest)
gtest_discover_tests(test_omp)
gtest_discover_tests(test_mpi)
