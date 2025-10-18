#include "link_margin_service.h"
#include <cmath>
#include <algorithm>

namespace LinkMarginService {

// 获取调制方式对应的Eb/N0要求
double getRequiredEbNo(const std::string& modulation, double ber) {
    // 查找调制方式
    auto modIt = MODULATION_EBN0_MAP.find(modulation);
    if (modIt == MODULATION_EBN0_MAP.end()) {
        // 默认使用QPSK
        modIt = MODULATION_EBN0_MAP.find("QPSK");
        if (modIt == MODULATION_EBN0_MAP.end()) {
            return 10.5; // 默认值
        }
    }
    
    // 查找最接近的BER要求
    const auto& berMap = modIt->second;
    auto lower = berMap.lower_bound(ber);
    
    if (lower == berMap.begin()) {
        return lower->second;
    }
    
    if (lower == berMap.end()) {
        return berMap.rbegin()->second;
    }
    
    auto prev = std::prev(lower);
    
    // 线性插值
    double ber1 = prev->first;
    double ebno1 = prev->second;
    double ber2 = lower->first;
    double ebno2 = lower->second;
    
    double weight = (ber - ber1) / (ber2 - ber1);
    return ebno1 + weight * (ebno2 - ebno1);
}

// 计算链路余量
nlohmann::json calculateLinkMargin(const nlohmann::json& link) {
    // 获取输入参数
    std::string linkId = link["link_id"].get<std::string>();
    double signalPower = link["signal_power"].get<double>(); // dBm
    double noisePower = link["noise_power"].get<double>(); // dBm
    double interferencePower = link["interference_power"].get<double>(); // dBm
    double bandwidth = link["bandwidth"].get<double>(); // MHz
    double bitRate = link["bit_rate"].get<double>(); // Mbps
    double berRequirement = link["ber_requirement"].get<double>();
    std::string modulation = link["modulation"].get<std::string>();
    
    // 计算载噪比 (C/N)
    double cnLinear = pow(10, (signalPower - noisePower) / 10.0);
    
    // 计算载干噪比 (C/(N+I))
    double interferenceLinear = pow(10, interferencePower / 10.0);
    double noiseLinear = pow(10, noisePower / 10.0);
    double totalInterferenceNoise = noiseLinear + interferenceLinear;
    double cnirLinear = pow(10, signalPower / 10.0) / totalInterferenceNoise;
    double cnir = 10 * log10(cnirLinear); // dB
    
    // 计算等效Eb/(N0+NI)
    double equivalentEbNo = cnir + 10 * log10(bandwidth * 1e6 / (bitRate * 1e6)); // dB
    
    // 获取所需Eb/N0
    double requiredEbNo = getRequiredEbNo(modulation, berRequirement);
    
    // 计算链路余量
    double linkMargin = equivalentEbNo - requiredEbNo;
    
    // 判断链路余量是否足够
    bool isSufficient = (linkMargin >= 2.0);
    std::string marginStatus = isSufficient ? "足够" : "不足";
    
    // 返回结果
    return {
        {"link_id", linkId},
        {"equivalent_ebno", equivalentEbNo},
        {"required_ebno", requiredEbNo},
        {"link_margin", linkMargin},
        {"is_sufficient", isSufficient},
        {"margin_status", marginStatus}
    };
}

// 计算多个链路的链路余量
nlohmann::json calculateAllLinkMargins(const nlohmann::json& links) {
    nlohmann::json results = nlohmann::json::array();
    double minMargin = 100.0;
    double maxMargin = -100.0;
    int sufficientCount = 0;
    int insufficientCount = 0;
    
    // 遍历所有链路
    for (const auto& link : links) {
        nlohmann::json result = calculateLinkMargin(link);
        results.push_back(result);
        
        // 更新统计信息
        double margin = result["link_margin"].get<double>();
        if (margin < minMargin) minMargin = margin;
        if (margin > maxMargin) maxMargin = margin;
        
        if (result["is_sufficient"].get<bool>()) {
            sufficientCount++;
        } else {
            insufficientCount++;
        }
    }
    
    // 创建汇总信息
    nlohmann::json summary = {
        {"total_links", links.size()},
        {"sufficient_links", sufficientCount},
        {"insufficient_links", insufficientCount},
        {"min_margin", minMargin},
        {"max_margin", maxMargin}
    };
    
    // 返回完整结果
    return {
        {"results", results},
        {"summary", summary}
    };
}

} // namespace LinkMarginService