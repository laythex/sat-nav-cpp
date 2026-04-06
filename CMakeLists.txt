cmake_minimum_required(VERSION 3.22)
set(PROJECT_NAME sat-nav)
project(${PROJECT_NAME})

file(GLOB SOURCE src/*.cpp)
add_executable(${PROJECT_NAME} ${SOURCE})
target_compile_options(${PROJECT_NAME} PRIVATE -Werror -Wall -Wextra -pedantic -Wconversion -Wsign-conversion) # -Wshadow опущен, вернуть
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_23)