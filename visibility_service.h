#ifndef VISIBILITY_SERVICE_H
#define VISIBILITY_SERVICE_H

#include "nlohmann/json.hpp"
#include "Funs.h"  // 包含已有的公共函数头文件

namespace VisibilityService {

// 常量定义
const double SPEED_OF_LIGHT = 299792.458; // 光速 (km/s)

// 计算单个地面站与单个卫星的可见性
nlohmann::json calculateTerminalSatelliteVisibility(
    const nlohmann::json& terminal, 
    const nlohmann::json& satellite);

// 计算多个地面站与多个卫星的可见性
nlohmann::json calculateMultiVisibility(
    const nlohmann::json& terminals, 
    const nlohmann::json& satellites);

} // namespace VisibilityService

#endif // VISIBILITY_SERVICE_H