cmake_minimum_required(VERSION 3.10)
project(Extreme_Env_Nav_Suite LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

add_executable(nav_suite_run src/main.cpp)
target_include_directories(nav_suite_run PRIVATE include)
target_link_libraries(nav_suite_run PRIVATE Eigen3::Eigen)
