cmake_minimum_required(VERSION 3.20)
project(satellitePRJ)

set(CMAKE_CXX_STANDARD 20)

# 添加平台特定的编译定义
if(WIN32)
    # Windows 特定的定义
    add_compile_definitions(_WIN32_WINNT=0x0A00)
    add_compile_definitions(NOMINMAX)
    add_compile_definitions(WIN32_LEAN_AND_MEAN)
else()
    # Linux 特定的定义（如果有）
    add_compile_definitions(LINUX_PLATFORM)
endif()

add_executable(satellitePRJ main.cpp
        src/visibility_calculator.h
        src/visibility_calculator.cpp
        src/Funs.cpp
        src/visibility_service.cpp
        src/delay_service.cpp
        src/link_margin_service.cpp
        src/communication_quality_service.cpp
        src/interference_assessment_service.cpp
        src/link_assessment_service.cpp
        src/simplified_setup_time_service.cpp)

# 添加平台特定的链接库
if(WIN32)
    # Windows 网络库
    target_link_libraries(satellitePRJ ws2_32 wsock32)
else()
    # Linux 通常不需要显式链接网络库
    # 但如果需要，可以添加其他库，例如：
    # target_link_libraries(satellitePRJ pthread)
endif()