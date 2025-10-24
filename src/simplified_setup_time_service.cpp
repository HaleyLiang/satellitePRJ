#include "simplified_setup_time_service.h"
#include <cmath>
#include <algorithm>
#include <vector>

namespace SimplifiedSetupTimeService {

// 计算传播时延 (ms)
double calculatePropagationDelay(double distance) {
    return (distance / SPEED_OF_LIGHT) * 1000; // 转换为毫秒
}

// 计算处理时延 (ms)
double calculateProcessingDelay(const nlohmann::json& link) {
    std::string modulation = link["modulation"].get<std::string>();
    double codingRate = link["coding_rate"].get<double>();
    std::string terminalType = link["terminal_type"].get<std::string>();
    
    // 基础处理时延
    double baseDelay = 10.0; // ms
    
    // 调制方式影响
    std::map<std::string, double> modulationFactors = {
        {"BPSK", 0.8},
        {"QPSK", 1.0},
        {"8PSK", 1.2},
        {"16QAM", 1.5},
        {"64QAM", 2.0},
        {"256QAM", 2.5}
    };
    double modFactor = modulationFactors.count(modulation) ? 
        modulationFactors[modulation] : 1.0;
    
    // 编码率影响
    double codingFactor = 1.0 + (1.0 - codingRate) * 0.5;
    
    // 终端类型影响
    double terminalFactor = 1.0;
    if (terminalType == "mobile") terminalFactor = 1.2;
    else if (terminalType == "portable") terminalFactor = 1.5;
    
    return baseDelay * modFactor * codingFactor * terminalFactor;
}

// 计算同步时延 (ms)
double calculateSynchronizationDelay(const nlohmann::json& link) {
    std::string antennaType = link["antenna_type"].get<std::string>();
    double minElevation = link["min_elevation"].get<double>();
    double frequency = link["frequency"].get<double>();
    
    // 基础同步时延
    double baseDelay = 5.0; // ms
    
    // 天线类型影响
    std::map<std::string, double> antennaFactors = {
        {"parabolic", 0.8},
        {"phased_array", 1.0},
        {"omnidirectional", 1.5}
    };
    double antennaFactor = antennaFactors.count(antennaType) ? 
        antennaFactors[antennaType] : 1.0;
    
    // 仰角影响
    double elevationFactor = 1.0 + (30.0 - std::min(30.0, minElevation)) / 30.0 * 0.5;
    
    // 频率影响
    double frequencyFactor = 1.0 - (frequency / 10e9) * 0.2;
    frequencyFactor = std::max(0.7, frequencyFactor);
    
    return baseDelay * antennaFactor * elevationFactor * frequencyFactor;
}

// 计算大气时延 (ms)
double calculateAtmosphericDelay(const nlohmann::json& link) {
    std::string conditions = link["atmospheric_conditions"].get<std::string>();
    double frequency = link["frequency"].get<double>();
    
    // 基础大气时延
    double baseDelay = 0.5; // ms
    
    // 大气条件影响
    std::map<std::string, double> conditionFactors = {
        {"clear", 1.0},
        {"cloudy", 1.5},
        {"rain", 2.5},
        {"storm", 4.0}
    };
    double conditionFactor = conditionFactors.count(conditions) ? 
        conditionFactors[conditions] : 1.0;
    
    // 频率影响
    double frequencyFactor = 1.0 + (frequency / 10e9) * 0.5;
    
    return baseDelay * conditionFactor * frequencyFactor;
}

// 计算建链时长
double calculateSetupTime(const nlohmann::json& link) {
    // 计算各组成部分时延
    double propagationDelay = calculatePropagationDelay(link["distance"].get<double>());
    double processingDelay = calculateProcessingDelay(link);
    double syncDelay = calculateSynchronizationDelay(link);
    double atmosDelay = calculateAtmosphericDelay(link);
    
    // 计算总时延
    return propagationDelay + processingDelay + syncDelay + atmosDelay;
}

// 计算所有链路的建链时长
nlohmann::json calculateAllSetupTimes(const nlohmann::json& links) {
    nlohmann::json results = nlohmann::json::array();
    
    // 统计信息
    int totalLinks = 0;
    double totalDelay = 0.0;
    double minDelay = 1000.0;
    double maxDelay = 0.0;
    
    for (const auto& link : links) {
        double setupTime = calculateSetupTime(link);
        
        // 创建结果对象
        nlohmann::json result;
        result["link_id"] = link["link_id"];
        result["setup_time_ms"] = setupTime;
        results.push_back(result);
        
        // 更新统计信息
        totalLinks++;
        totalDelay += setupTime;
        minDelay = std::min(minDelay, setupTime);
        maxDelay = std::max(maxDelay, setupTime);
    }
    
    // 创建汇总信息
    nlohmann::json summary;
    summary["total_links"] = totalLinks;
    summary["avg_setup_time_ms"] = totalDelay / totalLinks;
    summary["min_setup_time_ms"] = minDelay;
    summary["max_setup_time_ms"] = maxDelay;
    
    return {
        {"results", results},
        {"summary", summary}
    };
}

} // namespace SimplifiedSetupTimeService