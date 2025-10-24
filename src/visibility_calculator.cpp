#include "visibility_calculator.h"
#include <algorithm>
#include <vector>
#include <utility>

// 计算卫星与地面终端的可见性
json calculateVisibility(const json& terminal, const json& satellites) {
    json results = json::array();
    double grdLon = terminal["longitude"].get<double>();
    double grdLat = terminal["latitude"].get<double>();
    double grdHeight = terminal.contains("height") ? terminal["height"].get<double>() : 0.0;

    // 遍历所有卫星
    for (const auto& sat : satellites) {
        int satId = sat["satellite_id"].get<int>();
        double satLon = sat["longitude"].get<double>();
        double satLat = sat["latitude"].get<double>();
        double satAlt = sat["altitude"].get<double>(); // 千米
        
        // 计算卫星到地面终端的直线距离
        double distance = CalculateDistanceEarthToAir(satLon, satLat, grdLon, grdLat, satAlt);
        
        // 计算仰角
        double elevation = calElevation(satAlt, satLon, satLat, grdLon, grdLat);
        
        // 存储波束距离信息
        std::vector<std::pair<int, double>> beamDistances;
        bool visible = false;
        int visibleBeamId = -1;
        double bestSNR = -1000.0; // 初始化为极低值
        
        // 遍历卫星的所有波束
        for (const auto& beam : sat["beams"]) {
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
                
                // 计算通信参数
                double frequency = sat["frequency"].get<double>();
                double EIRP = sat["EIRP"].get<double>();
                double recAntGain = terminal["receiver_gain"].get<double>();
                
                // 计算自由空间损耗
                double FSPL = CalculateFSPL(distance, frequency);
                
                // 计算气候衰减
                double climaticLoss = 0.0;
                if (sat.contains("climatic_params")) {
                    const auto& clim = sat["climatic_params"];
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
                
                // 记录最佳波束
                if (SNR > bestSNR) {
                    bestSNR = SNR;
                    visibleBeamId = beamId;
                    
                    // 添加到结果
                    json satResult = {
                        {"satellite_id", satId},
                        {"visible", true},
                        {"beam_id", beamId},
                        {"distance", distance},
                        {"elevation", elevation},
                        {"signal_power", recSignalPower},
                        {"noise_power", noisePower},
                        {"SNR", SNR},
                        {"info_rate", infoRate},
                        {"time_delay", timeDelay},
                        {"communication_quality", commQuality}
                    };
                    
                    // 如果已有结果，替换为更好的波束
                    bool found = false;
                    for (auto& result : results) {
                        if (result["satellite_id"] == satId) {
                            result = satResult;
                            found = true;
                            break;
                        }
                    }
                    
                    if (!found) {
                        results.push_back(satResult);
                    }
                }
            }
        }
        
        // 如果没有可见波束，添加不可见结果
        if (!visible) {
            results.push_back({
                {"satellite_id", satId},
                {"visible", false},
                {"distance", distance},
                {"elevation", elevation}
            });
        }
    }
    
    return results;
}