#include "visibility_service.h"
#include <cmath>
#include <vector>
#include <algorithm>

namespace VisibilityService {

// 计算单个地面站与单个卫星的可见性
nlohmann::json calculateTerminalSatelliteVisibility(
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
    
    // 计算仰角
    double elevation = calElevation(satAlt, satLon, satLat, grdLon, grdLat);
    
    // 存储波束距离信息
    std::vector<std::pair<int, double>> beamDistances;
    bool visible = false;
    int visibleBeamId = -1;
    
    // 遍历卫星的所有波束
    for (const auto& beam : satellite["beams"]) {
        int beamId = beam["beam_id"].get<int>();
        double beamLon = beam["center_longitude"].get<double>();
        double beamLat = beam["center_latitude"].get<double>();
        
        // 计算地面终端到波束中心的距离
        double beamDistance = CalculateDistanceEarthToEarth(grdLon, grdLat, beamLon, beamLat);
        beamDistances.push_back({beamId, beamDistance});
        
        // 计算卫星到波束中心的距离
        double satToBeamDist = CalculateDistanceEarthToAir(satLon, satLat, beamLon, beamLat, satAlt);
        
        // 计算波束覆盖半径
        double beamAngle = beam["beam_angle_3db"].get<double>();
        double beamRadius = satToBeamDist * tan(deg2Rad(beamAngle / 2));
        
        // 检查可见性
        if (beamDistance <= beamRadius) {
            visible = true;
            visibleBeamId = beamId;
            
            // 计算通信参数
            double frequency = satellite["frequency"].get<double>();
            double EIRP = satellite["EIRP"].get<double>();
            double recAntGain = terminal["receiver_gain"].get<double>();
            
            // 计算自由空间损耗
            double FSPL = CalculateFSPL(distance, frequency);
            
            // 计算气候衰减
            double climaticLoss = 0.0;
            if (satellite.contains("climatic_params")) {
                const auto& clim = satellite["climatic_params"];
                climaticLoss = CalClimaticLoss(
                    frequency, satAlt, satLon, satLat, grdLon, grdLat, grdHeight,
                    clim["rain"].get<double>(), clim["temperature"].get<double>(),
                    clim["air_pressure"].get<double>(), clim["liquid_water_density"].get<double>(),
                    clim["cloud_vertical_thickness"].get<double>(), clim["rho"].get<double>()
                );
            }
            
            // 计算接收信号功率
            double recSignalPower = CalRecSignalPower(EIRP, recAntGain, FSPL, climaticLoss);
            
            // 计算噪声功率
            double bandwidth = terminal["bandwidth"].get<double>();
            double noisePower = CalNoisePower(bandwidth);
            
            // 计算信噪比
            double SNR = CalSNR(recSignalPower, noisePower);
            
            // 计算信息速率
            double infoRate = CalInfoRate(bandwidth, SNR);
            
            // 计算时延
            double timeDelay = CalTimeDealy(distance);
            
            // 计算通信质量等级
            int commQuality = calLinkCommLevel(timeDelay, infoRate);
            
            // 返回可见结果
            return {
                {"satellite_id", satId},
                {"visible", true},
                {"beam_id", visibleBeamId},
                {"distance", distance},
                {"elevation", elevation},
                {"signal_power", recSignalPower},
                {"noise_power", noisePower},
                {"SNR", SNR},
                {"info_rate", infoRate},
                {"time_delay", timeDelay},
                {"communication_quality", commQuality}
            };
        }
    }
    
    // 如果没有可见波束，返回不可见结果
    return {
        {"satellite_id", satId},
        {"visible", false},
        {"distance", distance},
        {"elevation", elevation}
    };
}

// 计算多个地面站与多个卫星的可见性
nlohmann::json calculateMultiVisibility(
    const nlohmann::json& terminals, 
    const nlohmann::json& satellites) {
    
    nlohmann::json results = nlohmann::json::array();
    
    // 遍历所有地面站
    for (const auto& terminal : terminals) {
        int terminalId = terminal["id"].get<int>();
        nlohmann::json terminalResult;
        terminalResult["terminal_id"] = terminalId;
        
        nlohmann::json satelliteResults = nlohmann::json::array();
        
        // 遍历所有卫星
        for (const auto& satellite : satellites) {
            nlohmann::json satResult = calculateTerminalSatelliteVisibility(terminal, satellite);
            satelliteResults.push_back(satResult);
        }
        
        terminalResult["satellite_results"] = satelliteResults;
        results.push_back(terminalResult);
    }
    
    return results;
}

} // namespace VisibilityService