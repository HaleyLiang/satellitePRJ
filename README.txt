卫星通信系统接口文档
1. 多个卫星单个地面站通信可见性接口
**URL:** /check_visibility
请求参数说明
{
  "terminal": {
    "id": 1,                       // 地面终端唯一标识符
    "longitude": 120.0,            // 地面终端经度（度）
    "latitude": 30.0,              // 地面终端纬度（度）
    "height": 0.0,                 // 地面终端海拔高度（米）
    "receiver_gain": 10.0,         // 接收机增益（dBi）
    "bandwidth": 10.0              // 通信带宽（MHz）
  },
  "satellites": [
    {
      "satellite_id": 101,         // 卫星唯一标识符
      "longitude": 100.0,          // 卫星经度（度）
      "latitude": 40.0,            // 卫星纬度（度）
      "altitude": 36000.0,         // 卫星海拔高度（千米）
      "frequency": 2.4e9,          // 通信频率（Hz）
      "EIRP": 50.0,               // 等效全向辐射功率（dBW）
      "climatic_params": {         // 气候参数
        "rain": 10.0,              // 降雨率（mm/h）
        "temperature": 25.0,       // 温度（摄氏度）
        "air_pressure": 1013.25,   // 大气压（hPa）
        "liquid_water_density": 5.0, // 液态水密度（g/m³）
        "cloud_vertical_thickness": 2.0, // 云层垂直厚度（km）
        "rho": 7.5                 // 水汽密度参数
      },
      "beams": [                   // 卫星波束配置
        {
          "beam_id": 1,            // 波束唯一标识符
          "center_longitude": 100.5, // 波束中心经度（度）
          "center_latitude": 40.5,   // 波束中心纬度（度）
          "beam_angle_3db": 10.0     // 波束3dB角（度）
        }
      ]
    }
  ]
}
响应参数说明
[
  {
    "satellite_id": 101,           // 卫星标识符
    "visible": true,               // 可见性标志（布尔值）
    "beam_id": 2,                  // 最佳波束ID
    "distance": 3678.45,           // 卫星到终端距离（千米）
    "elevation": 45.2,             // 卫星仰角（度）
    "signal_power": -85.3,         // 接收信号功率（dBm）
    "noise_power": -100.2,         // 噪声功率（dBm）
    "SNR": 14.9,                   // 信噪比（dB）
    "info_rate": 25.6,             // 信息速率（Mbps）
    "time_delay": 0.0123,          // 传输时延（秒）
    "communication_quality": 1     // 通信质量等级（0:高, 1:中, 2:低）
  }
]


2. 多卫星多地面站通信可见性接口
**URL:** /check_multi_visibility
请求参数说明
{
  "terminals": [                   // 多个地面终端
    {
      "id": 1,                     // 终端标识符
      "longitude": 120.0,          // 经度
      "latitude": 30.0,            // 纬度
      "height": 0.0,               // 海拔
      "receiver_gain": 10.0,       // 接收增益
      "bandwidth": 10.0            // 带宽
    }
  ],
  "satellites": [                  // 多个卫星
    {
      "satellite_id": 101,         // 卫星标识符
      // ... 其他参数同单卫星接口
    }
  ]
}
响应参数说明
[
  {
    "terminal_id": 1,              // 终端标识符
    "satellite_results": [         // 该终端对所有卫星的结果
      {
        // 卫星通信结果参数（同单卫星接口）
      }
    ]
  }
]


3. 时延计算接口
**URL:** /calculate_delays
请求参数说明
{
  "terminals": [
    {
      "id": 1,                     // 终端标识符
      "longitude": 120.0,          // 经度
      "latitude": 30.0,            // 纬度
      "height": 0.0                // 海拔
    }
  ],
  "satellites": [
    {
      "satellite_id": 101,         // 卫星标识符
      "longitude": 100.0,          // 经度
      "latitude": 40.0,            // 纬度
      "altitude": 36000.0          // 海拔高度
    }
  ]
}
响应参数说明
{
   "results": [  // 链路评估结果数组，包含所有被评估链路的详细信息
     {
       "communication_quality": {  // 通信质量评估
       "delay": 120.0,        // 传输时延，单位：毫秒(ms)
       "description": "可稳定传输视频",  // 通信质量描述
       "info_rate": 25.6,     // 信息传输速率，单位：Mbps
       "level": "高"          // 通信质量等级：高/中/低
     },
       "information_reliability": {  // 信息可靠性评估
       "ber": 1e-06,          // 误码率(Bit Error Rate)，百万分之一
       "description": "BER ≤ 10⁻⁵",  // 可靠性标准描述
       "level": "高"          // 可靠性等级：高/中/低
     },
       "interference_level": {    // 干扰水平评估
       "capacity_reduction": 0.08967262268017845,  // 容量下降比例，约8.97%
       "description": "通信速率基本无下降",  // 干扰影响描述
       "level": "无"          // 干扰等级：无/弱/中/强
     },
       "link_id": "link1",        // 链路唯一标识符
       "overall_quality": "优秀"  // 链路整体质量评级：优秀/良好/一般/较差
     }
   ],
   "summary": {  // 总体统计摘要
     "avg_ber": 1e-06,              // 平均误码率
     "avg_delay": 120.0,            // 平均传输时延，单位：毫秒(ms)
     "avg_info_rate": 25.6,         // 平均信息传输速率，单位：Mbps
     "excellent_quality": 1,        // 质量为"优秀"的链路数量
     "fair_quality": 0,             // 质量为"一般"的链路数量
     "good_quality": 0,             // 质量为"良好"的链路数量
     "poor_quality": 0,             // 质量为"较差"的链路数量
     "total_links": 1               // 评估的链路总数
   }
}


4. 链路余量计算接口
**URL:** /calculate_link_margins
请求参数说明
{
  "links": [
    {
      "link_id": "link1",           // 链路标识符
      "signal_power": -85.3,       // 信号功率（dBm）
      "noise_power": -100.2,       // 噪声功率（dBm）
      "interference_power": -105.0, // 干扰功率（dBm）
      "bandwidth": 10.0,            // 带宽（MHz）
      "bit_rate": 5.0,             // 比特率（Mbps）
      "ber_requirement": 1e-6,     // 误码率要求
      "modulation": "QPSK"         // 调制方式
    }
  ]
}
响应参数说明
{
  "results": [  // 链路余量评估结果数组
    {
      "equivalent_ebno": 16.66809158422522,  // 等效信噪比（Eb/N0），单位：dB
      "is_sufficient": true,                 // 链路余量是否足够（布尔值）
      "link_id": "link1",                    // 链路唯一标识符
      "link_margin": 6.16809158422522,       // 链路余量，单位：dB
      "margin_status": "足够",               // 余量状态描述
      "required_ebno": 10.5                 // 系统要求的Eb/N0门限，单位：dB
    }
  ],
  "summary": {  // 总体统计摘要
    "insufficient_links": 0,      // 余量不足的链路数量
    "max_margin": 6.16809158422522,  // 最大链路余量，单位：dB
    "min_margin": 6.16809158422522,  // 最小链路余量，单位：dB
    "sufficient_links": 1,        // 余量足够的链路数量
    "total_links": 1              // 评估的链路总数
  }
}


5. 链路质量评估接口
**URL:** /calculate_communication_quality
请求参数说明
{
  "links": [
    {
      "link_id": "link1",           // 链路标识符
      "info_rate": 25.6,           // 信息速率（Mbps）
      "delay": 120,                // 时延（ms）
      "ber": 1e-6,                 // 误码率
      "snr": 14.9,                 // 信噪比（dB）
      "sinr": 13.5                 // 信干噪比（dB）
    }
  ]
}
响应参数说明
{
  "results": [  // 链路质量评估结果数组，包含所有被评估链路的详细信息
    {
      "communication_quality": {  // 通信质量评估部分
        "delay": 120.0,        // 传输时延，单位：毫秒(ms)，值越小越好
        "description": "可稳定传输视频",  // 通信质量描述，说明该质量等级支持的应用场景
        "info_rate": 25.6,     // 信息传输速率，单位：Mbps（兆比特每秒），值越大越好
        "level": "高"          // 通信质量等级：高/中/低，基于时延和速率综合评估
      },
      "information_reliability": {  // 信息可靠性评估部分
        "ber": 1e-06,          // 误码率(Bit Error Rate)，表示传输错误的比特比例，值越小越可靠
        "description": "BER ≤ 10⁻⁵",  // 可靠性标准描述，说明满足的误码率要求
        "level": "高"          // 可靠性等级：高/中/低，基于误码率评估
      },
        "interference_level": {    // 干扰水平评估部分
        "capacity_reduction": 0.08967262268017845,  // 容量下降比例（0-1之间的小数），8.97%的下降
        "description": "通信速率基本无下降",  // 干扰影响描述
        "level": "无"          // 干扰等级：无/弱/中/强，基于容量下降比例分类
      },
      "link_id": "link1",        // 链路唯一标识符，用于区分不同的通信链路
      "overall_quality": "优秀"  // 链路整体质量评级：优秀/良好/一般/较差，综合所有指标的评价
    }
  ],
    "summary": {  // 总体统计摘要信息
      "avg_ber": 1e-06,              // 平均误码率，所有链路误码率的平均值
      "avg_delay": 120.0,            // 平均传输时延，单位：毫秒(ms)
      "avg_info_rate": 25.6,         // 平均信息传输速率，单位：Mbps
      "excellent_quality": 1,        // 质量为"优秀"的链路数量
      "fair_quality": 0,             // 质量为"一般"的链路数量
      "good_quality": 0,             // 质量为"良好"的链路数量
      "poor_quality": 0,             // 质量为"较差"的链路数量
      "total_links": 1               // 评估的链路总数
   }
}

6. 受干扰情况计算接口
**URL:** /assess_interference
请求参数说明
{
  "links": [
    {
      "link_id": "SAT001-GS001",   // 链路标识符（卫星-地面站）
      "snr": 46.67,               // 信噪比（dB）
      "sinr": 40.67               // 信干噪比（dB）
    }
  ]
}
响应参数说明
{
  "results": [
    {
      "link_id": "SAT001-GS001",   // 链路标识符
      "snr": 46.67,               // 信噪比
      "sinr": 40.67,              // 信干噪比
      "capacity_snr": 15.23,       // 基于SNR的容量（bps/Hz）
      "capacity_sinr": 13.45,      // 基于SINR的容量（bps/Hz）
      "capacity_reduction": 0.117, // 容量下降比例
      "interference_level": "弱干扰", // 干扰等级
      "description": "通信速率下降10%-30%", // 干扰描述
      "impact_assessment": "影响轻微，通信质量基本不受影响" // 影响评估
    }
  ],
  "summary": {                    // 汇总信息
    "total_links": 2,             // 总链路数
    "weak_interference": 1,       // 弱干扰链路数
    "medium_interference": 1,     // 中干扰链路数
    "strong_interference": 0,     // 强干扰链路数
    "avg_capacity_reduction": 0.157 // 平均容量下降比例
  }
}
7. 优选链路评估接口
**URL:** /assess_links
请求参数说明
{
  "links": [
    {
      "type": "relay_link",               // 链路类型（中继/直连）
      "link_id": "SAT001-RELAY001-GS001", // 链路标识符
      "satellite": "SAT001",              // 卫星标识
      "relay": "RELAY001",               // 中继站标识（仅中继链路）
      "ground_station": "GS001",          // 地面站标识
      "relay_windows_count": 145,        // 中继窗口数量
      "start_time": "2025-09-09T00:00:00", // 通信开始时间
      "end_time": "2025-03-09T12:00:00",  // 通信结束时间
      "duration_minutes": 720.0,         // 通信持续时间（分钟）
      "min_elevation": 53.67,            // 最小仰角（度）
      "max_elevation": 53.67,            // 最大仰角（度）
      "avg_snr": 46.67,                  // 平均信噪比（dB）
      "avg_margin": 36.67                // 平均链路余量（dB）
    }
  ]
}
响应参数说明
{
  "results": [
    {
      "link_id": "SAT001-RELAY001-GS001", // 链路标识符
      "type": "relay_link",               // 链路类型
      "stability_score": 0.92,           // 稳定性评分（0-1）
      "aggregate_metrics": {              // 聚合指标
        "signal_quality": 0.95,          // 信号质量评分
        "reliability": 0.98,             // 可靠性评分
        "availability": 0.93,            // 可用性评分
        "performance": 0.96              // 性能评分
      },
      "overall_score": 0.94,              // 综合评分
      "recommendation": "优选",           // 推荐等级
      "strengths": ["高稳定性", "高可靠性"], // 优势
      "weaknesses": ["多跳传输时延较高"]   // 劣势
    }
  ],
  "preferred_link": {                    // 优选链路
    "link_id": "SAT001-RELAY001-GS001",
    "type": "relay_link",
    "overall_score": 0.94,
    "reason": "综合稳定性和可靠性最优"
  },
  "summary": {                          // 汇总信息
    "total_links": 2,                   // 总链路数
    "avg_stability": 0.90,              // 平均稳定性
    "avg_signal_quality": 0.955,        // 平均信号质量
    "recommendation": "建议使用中继链路..." // 总体建议
  }
}
8. 建链时长计算接口
**URL:** /calculate_setup_times
请求参数说明
{
  "links": [
    {
      "link_id": "SAT001-GS001",       // 链路标识符
      "satellite_type": "LEO",         // 卫星类型（LEO/MEO/GEO）
      "distance": 3678.45,             // 通信距离（km）
      "min_elevation": 53.67,          // 最小仰角（度）
      "frequency": 2.4e9,              // 频率（Hz）
      "modulation": "QPSK",            // 调制方式
      "coding_rate": 0.75,             // 编码率
      "terminal_type": "fixed",        // 终端类型（固定/移动）
      "antenna_type": "parabolic",     // 天线类型
      "atmospheric_conditions": "clear" // 大气条件
    }
  ]
}
响应参数说明
{
  "results": [
    {
      "link_id": "SAT001-GS001",       // 链路标识符
      "setup_time_ms": 36.96           // 建链时间（毫秒）
    }
  ],
  "summary": {                        // 汇总信息
    "total_links": 2,                 // 总链路数
    "avg_setup_time_ms": 96.23,       // 平均建链时间
    "min_setup_time_ms": 36.96,       // 最小建链时间
    "max_setup_time_ms": 155.5         // 最大建链时间
  }
}
技术说明
干扰等级判定标准：
**弱干扰**: (1-30%)log(1+SNR) < log(1+SINR) ≤ (1-10%)log(1+SNR)
**中干扰**: (1-50%)log(1+SNR) < log(1+SINR) ≤ (1-30%)log(1+SNR)
**强干扰**: log(1+SINR) ≤ (1-50%)*log(1+SNR)
通信质量等级：
**0**: 高质量通信
**1**: 中等质量通信
**2**: 低质量通信
卫星类型：
**LEO**: 低地球轨道卫星
**MEO**: 中地球轨道卫星
**GEO**: 地球同步轨道卫星