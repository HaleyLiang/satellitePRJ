#ifndef VISIBILITY_CALCULATOR_H
#define VISIBILITY_CALCULATOR_H

#include "../nlohmann/json.hpp"
#include "../src/Funs.h"

using json = nlohmann::json;

// 常量定义
const double SPEED_OF_LIGHT = 299792.458; // 光速 (km/s)

// 计算卫星与地面终端的可见性
json calculateVisibility(const json& terminal, const json& satellites);

#endif // VISIBILITY_CALCULATOR_H