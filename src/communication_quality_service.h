#ifndef COMMUNICATION_QUALITY_SERVICE_H
#define COMMUNICATION_QUALITY_SERVICE_H

#include "../nlohmann/json.hpp"
#include <cmath>
#include <string>

namespace CommunicationQualityService {

// 计算通信质量
nlohmann::json calculateCommunicationQuality(double info_rate, double delay);

// 计算重要信息可靠性
nlohmann::json calculateInformationReliability(double ber);

// 计算受干扰情况
nlohmann::json calculateInterferenceLevel(double snr, double sinr);

// 计算总体质量评估
std::string calculateOverallQuality(
    const nlohmann::json& comm_quality,
    const nlohmann::json& info_reliability,
    const nlohmann::json& interference_level);

// 计算所有链路的通信质量
nlohmann::json calculateAllLinkQualities(const nlohmann::json& links);

} // namespace CommunicationQualityService

#endif // COMMUNICATION_QUALITY_SERVICE_H