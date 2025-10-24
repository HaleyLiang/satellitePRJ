#include "interference_assessment_service.h"
#include <cmath>
#include <algorithm>

namespace InterferenceAssessmentService {

// 计算通信容量（香农公式）
double calculateCapacity(double snr) {
    // 将dB值转换为线性值
    double snr_linear = std::pow(10.0, snr / 10.0);
    // 香农公式：C = log2(1 + SNR)
    return std::log2(1.0 + snr_linear);
}

// 计算通信容量下降比例
double calculateCapacityReduction(double snr, double sinr) {
    double capacity_snr = calculateCapacity(snr);
    double capacity_sinr = calculateCapacity(sinr);
    
    if (capacity_snr <= 0.0) {
        return 1.0; // 完全下降
    }
    
    return (capacity_snr - capacity_sinr) / capacity_snr;
}

// 评估受干扰情况
nlohmann::json assessInterferenceLevel(double snr, double sinr) {
    nlohmann::json result;
    
    // 计算通信容量
    double capacity_snr = calculateCapacity(snr);
    double capacity_sinr = calculateCapacity(sinr);
    double reduction = calculateCapacityReduction(snr, sinr);
    
    result["capacity_snr"] = capacity_snr;
    result["capacity_sinr"] = capacity_sinr;
    result["capacity_reduction"] = reduction;
    
    // 根据标准判断干扰等级
    if (capacity_sinr > (1.0 - 0.1) * capacity_snr) {
        // log(1+SINR) > (1-10%)*log(1+SNR)
        result["interference_level"] = "无干扰";
        result["description"] = "通信速率基本无下降";
        result["impact_assessment"] = "干扰影响可忽略";
    }
    else if (capacity_sinr > (1.0 - 0.3) * capacity_snr && 
             capacity_sinr <= (1.0 - 0.1) * capacity_snr) {
        // (1-30%)*log(1+SNR) < log(1+SINR) ≤ (1-10%)*log(1+SNR)
        result["interference_level"] = "弱干扰";
        result["description"] = "通信速率下降10%-30%";
        result["impact_assessment"] = "影响轻微，通信质量基本不受影响";
    }
    else if (capacity_sinr > (1.0 - 0.5) * capacity_snr && 
             capacity_sinr <= (1.0 - 0.3) * capacity_snr) {
        // (1-50%)*log(1+SNR) < log(1+SINR) ≤ (1-30%)*log(1+SNR)
        result["interference_level"] = "中干扰";
        result["description"] = "通信速率下降30%-50%";
        result["impact_assessment"] = "影响中等，建议优化通信参数";
    }
    else {
        // log(1+SINR) ≤ (1-50%)*log(1+SNR)
        result["interference_level"] = "强干扰";
        result["description"] = "通信速率下降超过50%";
        result["impact_assessment"] = "影响严重，需要采取抗干扰措施";
    }
    
    return result;
}

// 评估所有链路的受干扰情况
nlohmann::json assessAllLinksInterference(const nlohmann::json& links) {
    nlohmann::json results = nlohmann::json::array();
    
    // 统计信息
    int total_links = 0;
    int weak_interference = 0;
    int medium_interference = 0;
    int strong_interference = 0;
    int no_interference = 0;
    double total_reduction = 0.0;
    
    for (const auto& link : links) {
        std::string link_id = link["link_id"].get<std::string>();
        double snr = link["snr"].get<double>();
        double sinr = link["sinr"].get<double>();
        
        // 评估受干扰情况
        auto interference_result = assessInterferenceLevel(snr, sinr);
        
        // 添加链路基本信息
        nlohmann::json result;
        result["link_id"] = link_id;
        result["snr"] = snr;
        result["sinr"] = sinr;
        result.merge_patch(interference_result);
        
        results.push_back(result);
        
        // 更新统计信息
        total_links++;
        total_reduction += interference_result["capacity_reduction"].get<double>();
        
        std::string level = interference_result["interference_level"].get<std::string>();
        if (level == "无干扰") no_interference++;
        else if (level == "弱干扰") weak_interference++;
        else if (level == "中干扰") medium_interference++;
        else if (level == "强干扰") strong_interference++;
    }
    
    // 创建汇总信息
    nlohmann::json summary;
    summary["total_links"] = total_links;
    summary["no_interference"] = no_interference;
    summary["weak_interference"] = weak_interference;
    summary["medium_interference"] = medium_interference;
    summary["strong_interference"] = strong_interference;
    summary["avg_capacity_reduction"] = total_reduction / total_links;
    
    return {
        {"results", results},
        {"summary", summary}
    };
}

} // namespace InterferenceAssessmentService