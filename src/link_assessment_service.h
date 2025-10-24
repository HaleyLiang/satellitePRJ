#ifndef LINK_ASSESSMENT_SERVICE_H
#define LINK_ASSESSMENT_SERVICE_H

#include "../nlohmann/json.hpp"
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

namespace LinkAssessmentService {

// 计算链路稳定性
double calculateStabilityScore(const nlohmann::json& link);

// 计算聚合指标
nlohmann::json calculateAggregateMetrics(const nlohmann::json& link);

// 计算总体评分
double calculateOverallScore(const nlohmann::json& link);

// 生成链路评估
nlohmann::json assessLink(const nlohmann::json& link);

// 评估所有链路并推荐优选链路
nlohmann::json assessAllLinks(const nlohmann::json& links);

} // namespace LinkAssessmentService

#endif // LINK_ASSESSMENT_SERVICE_H