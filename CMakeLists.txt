cmake_minimum_required(VERSION 3.22)
set(PROJECT_NAME sat-nav)
project(${PROJECT_NAME} VERSION 1.0.0 LANGUAGES CXX)

find_package(Eigen3 REQUIRED NO_MODULE)

add_library(${PROJECT_NAME}_options INTERFACE)
target_compile_options(${PROJECT_NAME}_options INTERFACE 
    -Werror -Wall -Wextra -pedantic -Wconversion -Wsign-conversion
)
target_compile_features(${PROJECT_NAME}_options INTERFACE cxx_std_23)

file(GLOB_RECURSE LIB_SOURCES "src/*.cpp")
add_library(${PROJECT_NAME}_core STATIC ${LIB_SOURCES})
target_link_libraries(${PROJECT_NAME}_core PUBLIC 
    ${PROJECT_NAME}_options
    Eigen3::Eigen
)
target_include_directories(${PROJECT_NAME}_core PUBLIC src)

file(GLOB APP_SOURCES "apps/*.cpp")

foreach(app_source ${APP_SOURCES})
    get_filename_component(target_name ${app_source} NAME_WE)
    add_executable(${target_name} ${app_source})
    target_link_libraries(${target_name} PRIVATE 
        ${PROJECT_NAME}_core
        ${PROJECT_NAME}_options
    )
endforeach()