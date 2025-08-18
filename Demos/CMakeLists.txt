cmake_minimum_required(VERSION 3.0)

project(Demos LANGUAGES CXX)

set( easycppogl_DIR "../Externals/easycppogl/Install/EasyCppOGL/lib/cmake/easycppogl")
find_package(easycppogl REQUIRED)

set(CMAKE_DEBUG_POSTFIX "_d")

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/Install/Demos" CACHE PATH "Install path" FORCE)

macro(compile)
  add_executable(${ARGV0} ${ARGV0}.cpp)
  target_link_libraries(${ARGV0} easycppogl)
  
  #install(TARGETS ${ARGV0} DESTINATION bin)

  message(STATUS "CMAKE_CURRENT_SOURCE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}")

  set(RAW_DATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../Data")
  message(STATUS "RAW_DATA_PATH = ${RAW_DATA_PATH}")

  get_filename_component(DATA_PATH_ABS "${CMAKE_CURRENT_SOURCE_DIR}/../Data" ABSOLUTE)
  message(STATUS "Resolved absolute path to data: ${DATA_PATH_ABS}")

#  target_compile_definitions(${ARGV0} PRIVATE 
#    DEF_SHADERS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/"
#    DEF_DATA_PATH="${DATA_PATH_ABS}/"
#  )

  target_compile_definitions(${ARGV0} PRIVATE 
    MYAPP_SHADERS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/"
    MYAPP_DATA_PATH=${DATA_PATH_ABS}
  )

endmacro()

# ----------------------------------------------------------------------------------------

macro(compile_demo)
  add_executable(${ARGV0} ${ARGV0}.cpp WaveNoise.cpp WaveNoise.h hvPair.h)
  target_link_libraries(${ARGV0} easycppogl)
  
  #install(TARGETS ${ARGV0} DESTINATION bin)

  message(STATUS "CMAKE_CURRENT_SOURCE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}")

  set(RAW_DATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../Data")
  message(STATUS "RAW_DATA_PATH = ${RAW_DATA_PATH}")

  get_filename_component(DATA_PATH_ABS "${CMAKE_CURRENT_SOURCE_DIR}/../Data" ABSOLUTE)
  message(STATUS "Resolved absolute path to data: ${DATA_PATH_ABS}")

#  target_compile_definitions(${ARGV0} PRIVATE 
#    DEF_SHADERS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/"
#    DEF_DATA_PATH="${DATA_PATH_ABS}/"
#  )

  target_compile_definitions(${ARGV0} PRIVATE 
    MYAPP_SHADERS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/"
    MYAPP_DATA_PATH=${DATA_PATH_ABS}
  )

endmacro()

# ----------------------------------------------------------------------------------------

#add_definitions("-DDEF_SHADERS_PATH=${CMAKE_CURRENT_SOURCE_DIR}/")
#add_definitions("-DDEF_DATA_PATH=${CMAKE_CURRENT_SOURCE_DIR}/../Data")

set( demos_list 
  demo_gl_01
  wave_noise_01
)

foreach(example ${demos_list})
  compile(${example})
endforeach()

compile_demo( pg_wave_noise_demo )

# nice vs solution
set_target_properties(${demos_list} PROPERTIES FOLDER Demos)

set_target_properties(pg_wave_noise_demo PROPERTIES FOLDER Demos)

