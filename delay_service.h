#ifndef DELAY_SERVICE_H
#define DELAY_SERVICE_H

#include "nlohmann/json.hpp"
#include "Funs.h"  // 包含已有的公共函数头文件

namespace DelayService {

// 常量定义
const double SPEED_OF_LIGHT = 299792.458; // 光速 (km/s)
const double QUEUEING_DELAY = 0.010;     // 排队时延10ms

// 计算单个地面站与单个卫星的通信时延
nlohmann::json calculateTerminalSatelliteDelay(
    const nlohmann::json& terminal, 
    const nlohmann::json& satellite);

// 计算多个地面站与多个卫星的通信时延
nlohmann::json calculateMultiDelay(
    const nlohmann::json& terminals, 
    const nlohmann::json& satellites);

} // namespace DelayService

#endif // DELAY_SERVICE_H