#include "communication_quality_service.h"
#include <cmath>
#include <algorithm>

namespace CommunicationQualityService {

// 计算通信质量
nlohmann::json calculateCommunicationQuality(double info_rate, double delay) {
    nlohmann::json result;
    
    if (delay < 200 && info_rate >= 10.0) {
        result["level"] = "高";
        result["description"] = "可稳定传输视频";
    } else if (delay < 300 && info_rate >= 1.0) {
        result["level"] = "中";
        result["description"] = "可断续传输语音或稳定传输语音";
    } else if (delay < 500 && info_rate < 1.0) {
        result["level"] = "低";
        result["description"] = "可断续传输语音或稳定传输短信";
    } else {
        result["level"] = "极低";
        result["description"] = "通信质量差，可能无法正常通信";
    }
    
    result["info_rate"] = info_rate;
    result["delay"] = delay;
    
    return result;
}

// 计算重要信息可靠性
nlohmann::json calculateInformationReliability(double ber) {
    nlohmann::json result;
    
    if (ber <= 1e-5) {
        result["level"] = "高";
        result["description"] = "BER ≤ 10⁻⁵";
    } else if (ber > 1e-5 && ber <= 1e-4) {
        result["level"] = "中";
        result["description"] = "10⁻⁵ < BER ≤ 10⁻⁴";
    } else if (ber > 1e-4 && ber <= 1e-3) {
        result["level"] = "低";
        result["description"] = "10⁻⁴ < BER ≤ 10⁻³";
    } else {
        result["level"] = "极低";
        result["description"] = "BER > 10⁻³，通信不可靠";
    }
    
    result["ber"] = ber;
    return result;
}

// 计算受干扰情况
nlohmann::json calculateInterferenceLevel(double snr, double sinr) {
    nlohmann::json result;
    
    // 将dB值转换为线性值
    double snr_linear = std::pow(10.0, snr / 10.0);
    double sinr_linear = std::pow(10.0, sinr / 10.0);
    
    // 计算通信容量下降比例
    double capacity_snr = std::log2(1 + snr_linear);
    double capacity_sinr = std::log2(1 + sinr_linear);
    double reduction = (capacity_snr - capacity_sinr) / capacity_snr;
    
    result["capacity_reduction"] = reduction;
    
    if (reduction <= 0.1) {
        result["level"] = "无";
        result["description"] = "通信速率基本无下降";
    } else if (reduction > 0.1 && reduction <= 0.3) {
        result["level"] = "弱";
        result["description"] = "通信速率下降10%-30%";
    } else if (reduction > 0.3 && reduction <= 0.5) {
        result["level"] = "中";
        result["description"] = "通信速率下降30%-50%";
    } else {
        result["level"] = "强";
        result["description"] = "通信速率下降超过50%";
    }
    
    return result;
}

// 计算总体质量评估
std::string calculateOverallQuality(
    const nlohmann::json& comm_quality,
    const nlohmann::json& info_reliability,
    const nlohmann::json& interference_level) {
    
    std::string comm_level = comm_quality["level"].get<std::string>();
    std::string info_level = info_reliability["level"].get<std::string>();
    std::string interf_level = interference_level["level"].get<std::string>();
    
    // 所有指标均为高/无干扰
    if (comm_level == "高" && info_level == "高" && interf_level == "无") {
        return "优秀";
    }
    
    // 通信质量和信息可靠性至少为中等，干扰不超过弱
    if ((comm_level == "高" || comm_level == "中") &&
        (info_level == "高" || info_level == "中") &&
        (interf_level == "无" || interf_level == "弱")) {
        return "良好";
    }
    
    // 通信质量至少为低，信息可靠性至少为低，干扰不超过中
    if (comm_level != "极低" && info_level != "极低" && interf_level != "强") {
        return "一般";
    }
    
    // 其他情况
    return "较差";
}

// 计算所有链路的通信质量
nlohmann::json calculateAllLinkQualities(const nlohmann::json& links) {
    nlohmann::json results = nlohmann::json::array();
    
    // 统计信息
    int total_links = 0;
    int excellent_count = 0;
    int good_count = 0;
    int fair_count = 0;
    int poor_count = 0;
    double total_info_rate = 0.0;
    double total_delay = 0.0;
    double total_ber = 0.0;
    
    for (const auto& link : links) {
        std::string link_id = link["link_id"].get<std::string>();
        double info_rate = link["info_rate"].get<double>();
        double delay = link["delay"].get<double>();
        double ber = link["ber"].get<double>();
        double snr = link["snr"].get<double>();
        double sinr = link["sinr"].get<double>();
        
        // 计算各项指标
        auto comm_quality = calculateCommunicationQuality(info_rate, delay);
        auto info_reliability = calculateInformationReliability(ber);
        auto interference_level = calculateInterferenceLevel(snr, sinr);
        
        // 计算总体质量
        std::string overall_quality = calculateOverallQuality(
            comm_quality, info_reliability, interference_level);
        
        // 更新统计信息
        total_links++;
        total_info_rate += info_rate;
        total_delay += delay;
        total_ber += ber;
        
        if (overall_quality == "优秀") excellent_count++;
        else if (overall_quality == "良好") good_count++;
        else if (overall_quality == "一般") fair_count++;
        else if (overall_quality == "较差") poor_count++;
        
        // 构建结果对象
        nlohmann::json result = {
            {"link_id", link_id},
            {"communication_quality", comm_quality},
            {"information_reliability", info_reliability},
            {"interference_level", interference_level},
            {"overall_quality", overall_quality}
        };
        
        results.push_back(result);
    }
    
    // 创建汇总信息
    nlohmann::json summary = {
        {"total_links", total_links},
        {"excellent_quality", excellent_count},
        {"good_quality", good_count},
        {"fair_quality", fair_count},
        {"poor_quality", poor_count},
        {"avg_info_rate", total_info_rate / total_links},
        {"avg_delay", total_delay / total_links},
        {"avg_ber", total_ber / total_links}
    };
    
    return {
        {"results", results},
        {"summary", summary}
    };
}

} // namespace CommunicationQualityService