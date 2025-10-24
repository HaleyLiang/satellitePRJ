#include "delay_service.h"
#include <cmath>
#include <vector>

namespace DelayService {

// 计算单个地面站与单个卫星的通信时延
nlohmann::json calculateTerminalSatelliteDelay(
    const nlohmann::json& terminal, 
    const nlohmann::json& satellite) {
    
    double grdLon = terminal["longitude"].get<double>();
    double grdLat = terminal["latitude"].get<double>();
    double grdHeight = terminal.contains("height") ? terminal["height"].get<double>() : 0.0;
    
    int satId = satellite["satellite_id"].get<int>();
    double satLon = satellite["longitude"].get<double>();
    double satLat = satellite["latitude"].get<double>();
    double satAlt = satellite["altitude"].get<double>(); // 千米
    
    // 计算卫星到地面终端的直线距离
    double distance = CalculateDistanceEarthToAir(satLon, satLat, grdLon, grdLat, satAlt);
    
    // 计算传播时延（秒）
    double propagationDelay = CalTimeDealy(distance);
    
    // 计算排队时延（固定10ms）
    double queueingDelay = QUEUEING_DELAY;
    
    // 计算总时延
    double totalDelay = propagationDelay + queueingDelay;
    
    // 计算仰角
    double elevation = calElevation(satAlt, satLon, satLat, grdLon, grdLat);
    
    // 返回结果
    return {
        {"satellite_id", satId},
        {"distance", distance},
        {"elevation", elevation},
        {"propagation_delay", propagationDelay},
        {"queueing_delay", queueingDelay},
        {"total_delay", totalDelay}
    };
}

// 计算多个地面站与多个卫星的通信时延
nlohmann::json calculateMultiDelay(
    const nlohmann::json& terminals, 
    const nlohmann::json& satellites) {
    
    nlohmann::json results = nlohmann::json::array();
    
    // 遍历所有地面站
    for (const auto& terminal : terminals) {
        int terminalId = terminal["id"].get<int>();
        nlohmann::json terminalResult;
        terminalResult["terminal_id"] = terminalId;
        
        nlohmann::json satelliteDelays = nlohmann::json::array();
        
        // 遍历所有卫星
        for (const auto& satellite : satellites) {
            nlohmann::json delayResult = calculateTerminalSatelliteDelay(terminal, satellite);
            satelliteDelays.push_back(delayResult);
        }
        
        terminalResult["satellite_delays"] = satelliteDelays;
        results.push_back(terminalResult);
    }
    
    return results;
}

} // namespace DelayService