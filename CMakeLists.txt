cmake_minimum_required(VERSION 3.10)
project(genetic)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${GLFW_INCLUDE_DIR})

#通过file(GLOB命令匹配文件，并存入变量SRC_FILES中
file(GLOB SRC_FILES 
  "${PROJECT_SOURCE_DIR}/src/*.h"
  "${PROJECT_SOURCE_DIR}/src/*.cpp" 
  "${PROJECT_SOURCE_DIR}/src/*.c"
  "${PROJECT_SOURCE_DIR}/src/*.cc")


#CMAKE_PROJECT_NAME宏替换前面的project(Blackhole)
add_executable(${CMAKE_PROJECT_NAME} ${SRC_FILES})

# Copy assets files after build.
add_custom_command(
  TARGET ${CMAKE_PROJECT_NAME}
  #编译后的操作
  POST_BUILD
  #将根目录下的assets复制到输出路径下
  COMMAND ${CMAKE_COMMAND} -E copy_directory "${PROJECT_SOURCE_DIR}/assets"
          "$<TARGET_FILE_DIR:${CMAKE_PROJECT_NAME}>/assets")