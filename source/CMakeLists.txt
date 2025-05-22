cmake_minimum_required(VERSION 3.22)
project(chroma_ex LANGUAGES C CXX)

find_package(Chroma REQUIRED)
find_package(CURL   REQUIRED)

message(STATUS "Use Chroma from: ${Chroma_CONFIG}")

# ----------- 源文件分组 -------------
set(SOURCE_CHROMA ${CMAKE_CURRENT_SOURCE_DIR}/source/chroma.cc)
file(GLOB_RECURSE SOURCE CONFIGURE_DEPENDS source/*.cc)
list(REMOVE_ITEM SOURCE ${SOURCE_CHROMA})

add_library(chroma_ex_libs OBJECT ${SOURCE})
target_link_libraries(chroma_ex_libs PUBLIC Chroma::chromalib CURL::libcurl)

add_executable(chroma ${SOURCE_CHROMA})
target_link_libraries(chroma PRIVATE chroma_ex_libs)

set_property(TARGET chroma chroma_ex_libs PROPERTY CXX_STANDARD 17)
set_target_properties(chroma PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")

install(TARGETS chroma
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)
