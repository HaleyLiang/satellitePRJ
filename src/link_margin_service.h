#ifndef LINK_MARGIN_SERVICE_H
#define LINK_MARGIN_SERVICE_H

#include "../nlohmann/json.hpp"
#include <string>
#include <map>

namespace LinkMarginService {

// 调制方式与Eb/N0要求映射
const std::map<std::string, std::map<double, double>> MODULATION_EBN0_MAP = {
    {"BPSK", {
        {1e-3, 6.8},
        {1e-4, 8.4},
        {1e-5, 9.6},
        {1e-6, 10.5}
    }},
    {"QPSK", {
        {1e-3, 6.8},
        {1e-4, 8.4},
        {1e-5, 9.6},
        {1e-6, 10.5}
    }},
    {"8PSK", {
        {1e-3, 9.8},
        {1e-4, 11.6},
        {1e-5, 12.8},
        {1e-6, 13.5}
    }},
    {"16QAM", {
        {1e-3, 11.2},
        {1e-4, 13.2},
        {1e-5, 14.5},
        {1e-6, 15.2}
    }},
    {"64QAM", {
        {1e-3, 14.8},
        {1e-4, 17.0},
        {1e-5, 18.5},
        {1e-6, 19.5}
    }}
};

// 计算链路余量
nlohmann::json calculateLinkMargin(const nlohmann::json& link);

// 计算多个链路的链路余量
nlohmann::json calculateAllLinkMargins(const nlohmann::json& links);

} // namespace LinkMarginService

#endif // LINK_MARGIN_SERVICE_H