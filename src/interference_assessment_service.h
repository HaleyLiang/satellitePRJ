#ifndef INTERFERENCE_ASSESSMENT_SERVICE_H
#define INTERFERENCE_ASSESSMENT_SERVICE_H

#include "../nlohmann/json.hpp"
#include <cmath>
#include <string>
#include <vector>

namespace InterferenceAssessmentService {

// 计算通信容量（香农公式）
double calculateCapacity(double snr);

// 计算通信容量下降比例
double calculateCapacityReduction(double snr, double sinr);

// 评估受干扰情况
nlohmann::json assessInterferenceLevel(double snr, double sinr);

// 评估所有链路的受干扰情况
nlohmann::json assessAllLinksInterference(const nlohmann::json& links);

} // namespace InterferenceAssessmentService

#endif // INTERFERENCE_ASSESSMENT_SERVICE_H