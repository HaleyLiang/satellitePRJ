#include "link_assessment_service.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace LinkAssessmentService {

// 计算链路稳定性
double calculateStabilityScore(const nlohmann::json& link) {
    double stability = 0.0;
    
    // 1. 中继窗口数量稳定性 (权重0.4)
    if (link["type"] == "relay_link") {
        int relayWindows = link["relay_windows_count"].get<int>();
        // 假设理想中继窗口数为200，最大为300
        double windowStability = std::min(1.0, relayWindows / 200.0);
        stability += windowStability * 0.4;
    } else {
        // 直接链路稳定性基础值
        stability += 0.8 * 0.4;
    }
    
    // 2. 仰角稳定性 (权重0.3)
    double minElevation = link["min_elevation"].get<double>();
    double elevationStability = std::min(1.0, minElevation / 45.0); // 45度以上为理想
    stability += elevationStability * 0.3;
    
    // 3. 持续时间稳定性 (权重0.2)
    double duration = link["duration_minutes"].get<double>();
    double durationStability = std::min(1.0, duration / 1440.0); // 24小时为理想
    stability += durationStability * 0.2;
    
    // 4. SNR稳定性 (权重0.1)
    double snr = link["avg_snr"].get<double>();
    double snrStability = std::min(1.0, snr / 50.0); // 50dB为理想
    stability += snrStability * 0.1;
    
    return std::min(1.0, stability);
}

// 计算聚合指标
nlohmann::json calculateAggregateMetrics(const nlohmann::json& link) {
    nlohmann::json metrics;
    
    // 1. 信号质量 (SNR和链路余量)
    double snr = link["avg_snr"].get<double>();
    double margin = link["avg_margin"].get<double>();
    double signalQuality = 0.6 * std::min(1.0, snr / 50.0) + 
                          0.4 * std::min(1.0, margin / 40.0);
    metrics["signal_quality"] = signalQuality;
    
    // 2. 可靠性 (仰角和持续时间)
    double minElevation = link["min_elevation"].get<double>();
    double duration = link["duration_minutes"].get<double>();
    double reliability = 0.7 * std::min(1.0, minElevation / 45.0) + 
                         0.3 * std::min(1.0, duration / 1440.0);
    metrics["reliability"] = reliability;
    
    // 3. 可用性 (中继窗口或直接链路特性)
    double availability = 0.0;
    if (link["type"] == "relay_link") {
        int relayWindows = link["relay_windows_count"].get<int>();
        availability = std::min(1.0, relayWindows / 200.0);
    } else {
        // 直接链路可用性基于仰角
        availability = std::min(1.0, minElevation / 45.0);
    }
    metrics["availability"] = availability;
    
    // 4. 性能 (综合指标)
    double performance = 0.4 * signalQuality + 
                         0.3 * reliability + 
                         0.3 * availability;
    metrics["performance"] = performance;
    
    return metrics;
}

// 计算总体评分
double calculateOverallScore(const nlohmann::json& link) {
    double stability = calculateStabilityScore(link);
    auto metrics = calculateAggregateMetrics(link);
    
    // 权重分配：稳定性40%，信号质量30%，可靠性20%，可用性10%
    double overall = 0.4 * stability +
                     0.3 * metrics["signal_quality"].get<double>() +
                     0.2 * metrics["reliability"].get<double>() +
                     0.1 * metrics["availability"].get<double>();
    
    return std::min(1.0, overall);
}

// 生成链路评估
nlohmann::json assessLink(const nlohmann::json& link) {
    nlohmann::json result;
    
    result["link_id"] = link["link_id"];
    result["type"] = link["type"];
    
    // 计算稳定性评分
    double stability = calculateStabilityScore(link);
    result["stability_score"] = stability;
    
    // 计算聚合指标
    result["aggregate_metrics"] = calculateAggregateMetrics(link);
    
    // 计算总体评分
    double overallScore = calculateOverallScore(link);
    result["overall_score"] = overallScore;
    
    // 生成推荐等级
    if (overallScore >= 0.9) {
        result["recommendation"] = "优选";
    } else if (overallScore >= 0.8) {
        result["recommendation"] = "次选";
    } else if (overallScore >= 0.7) {
        result["recommendation"] = "备选";
    } else {
        result["recommendation"] = "不推荐";
    }
    
    // 分析优缺点
    std::vector<std::string> strengths;
    std::vector<std::string> weaknesses;
    
    // 信号质量分析
    double signalQuality = result["aggregate_metrics"]["signal_quality"].get<double>();
    if (signalQuality >= 0.9) {
        strengths.push_back("极高信号质量");
    } else if (signalQuality >= 0.8) {
        strengths.push_back("高信号质量");
    } else if (signalQuality < 0.7) {
        weaknesses.push_back("信号质量一般");
    }
    
    // 可靠性分析
    double reliability = result["aggregate_metrics"]["reliability"].get<double>();
    if (reliability >= 0.9) {
        strengths.push_back("高可靠性");
    } else if (reliability < 0.8) {
        weaknesses.push_back("可靠性有待提升");
    }
    
    // 可用性分析
    double availability = result["aggregate_metrics"]["availability"].get<double>();
    if (availability >= 0.9) {
        strengths.push_back("高可用性");
    } else if (availability < 0.8) {
        weaknesses.push_back("可用性略低");
    }
    
    // 链路类型特定分析
    if (link["type"] == "relay_link") {
        strengths.push_back("多跳传输稳定性");
        weaknesses.push_back("多跳传输时延较高");
    } else {
        strengths.push_back("低时延直接传输");
        weaknesses.push_back("单点故障风险");
    }
    
    result["strengths"] = strengths;
    result["weaknesses"] = weaknesses;
    
    return result;
}

// 评估所有链路并推荐优选链路
nlohmann::json assessAllLinks(const nlohmann::json& links) {
    nlohmann::json results = nlohmann::json::array();
    nlohmann::json summary;
    
    int totalLinks = 0;
    double totalStability = 0.0;
    double totalSignalQuality = 0.0;
    double maxScore = 0.0;
    nlohmann::json preferredLink;
    
    for (const auto& link : links) {
        auto result = assessLink(link);
        results.push_back(result);
        
        // 更新统计信息
        totalLinks++;
        totalStability += result["stability_score"].get<double>();
        totalSignalQuality += result["aggregate_metrics"]["signal_quality"].get<double>();
        
        // 检查是否是最优选
        double overallScore = result["overall_score"].get<double>();
        if (overallScore > maxScore) {
            maxScore = overallScore;
            preferredLink = result;
        }
    }
    
    // 创建汇总信息
    summary["total_links"] = totalLinks;
    summary["avg_stability"] = totalStability / totalLinks;
    summary["avg_signal_quality"] = totalSignalQuality / totalLinks;
    
    // 创建优选链路信息
    nlohmann::json preferred;
    preferred["link_id"] = preferredLink["link_id"];
    preferred["type"] = preferredLink["type"];
    preferred["overall_score"] = preferredLink["overall_score"];
    
    // 生成推荐理由
    std::string reason = "综合评分最高";
    if (preferredLink["type"] == "relay_link") {
        reason += "，中继链路提供更高稳定性";
    } else {
        reason += "，直接链路提供更低时延";
    }
    preferred["reason"] = reason;
    
    // 生成总体建议
    std::string overallRecommendation = "建议使用" + preferredLink["link_id"].get<std::string>() + 
                                       "作为主链路";
    
    // 如果有次优选链路，建议作为备份
    for (const auto& result : results) {
        if (result["link_id"] != preferredLink["link_id"] && 
            result["overall_score"].get<double>() >= 0.8) {
            overallRecommendation += "，" + result["link_id"].get<std::string>() + "作为备份";
            break;
        }
    }
    summary["recommendation"] = overallRecommendation;
    
    return {
        {"results", results},
        {"preferred_link", preferred},
        {"summary", summary}
    };
}

} // namespace LinkAssessmentService