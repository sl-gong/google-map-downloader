cmake_minimum_required(VERSION 3.15)

project(GoogleMapDownloader)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 导入 GDAL 和 libcurl
find_package(GDAL REQUIRED)
find_package(CURL REQUIRED)

# 添加源文件
set(SOURCES
    downloader.cpp
    )

# 创建可执行文件
add_executable(${PROJECT_NAME} ${SOURCES})

# 链接 GDAL 和 libcurl
target_include_directories(${PROJECT_NAME} PRIVATE ${GDAL_INCLUDE_DIRS} ${CURL_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${GDAL_LIBRARIES} ${CURL_LIBRARIES})
