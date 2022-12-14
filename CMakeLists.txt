cmake_minimum_required(VERSION 3.21)
project(test)

SET(CMAKE_CXX_FLAGS "-O0")
SET(CMAKE_C_FLAGS "-O0")

set(CMAKE_CXX_STANDARD 20)

find_package (Python REQUIRED
        COMPONENTS Interpreter Development)

        add_executable(test test.cxx)
        #add_executable(test test_class.cxx)

target_include_directories(${PROJECT_NAME} PRIVATE
        ${Python_INCLUDE_DIRS})

target_link_directories(${PROJECT_NAME} PRIVATE
        ${Python_LIBRARY_DIRS})

target_link_libraries(${PROJECT_NAME} PRIVATE
        ${Python_LIBRARIES})