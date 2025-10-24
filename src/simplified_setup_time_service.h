#ifndef SIMPLIFIED_SETUP_TIME_SERVICE_H
#define SIMPLIFIED_SETUP_TIME_SERVICE_H

#include "../nlohmann/json.hpp"
#include <cmath>
#include <string>
#include <vector>
#include <map>

namespace SimplifiedSetupTimeService {

// 常量定义
const double SPEED_OF_LIGHT = 299792.458; // km/s

// 计算建链时长
double calculateSetupTime(const nlohmann::json& link);

// 计算所有链路的建链时长
nlohmann::json calculateAllSetupTimes(const nlohmann::json& links);

} // namespace SimplifiedSetupTimeService

#endif // SIMPLIFIED_SETUP_TIME_SERVICE_H