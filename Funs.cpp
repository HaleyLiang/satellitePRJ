
#include "Funs.h"
#include <chrono>

using namespace std;

// 角度转弧度
double deg2Rad(double deg) {
	return deg * M_PI / 180.0;
}

// 错误信息输出
//void logError(const std::string& message)
//{
//	// 以追加模式打开日志文件,保证后续的错误不覆盖前面的错误
//	std::ofstream error_log("error.log", std::ios::app);       
//	if (error_log.is_open()) {
//
//		// 获取当前时间并格式化
//		time_t raw_time;
//		struct tm time_info;
//		char buffer[26];
//
//		time(&raw_time);
//		localtime_s(&time_info, &raw_time);
//		asctime_s(buffer, sizeof(buffer), &time_info);
//
//		// 写入时间和错误信息
//		error_log << "[" << buffer << "] " << message << std::endl;
//		error_log.close();
//	}
//	else {
//		std::cerr << "Unable to open log file" << std::endl;
//	}
//}

// 经纬度转换为笛卡尔坐标系
void computeCartesian(double lat, double lon, double height, double& x, double& y, double& z)
{
	lat = deg2Rad(lat);
	lon = deg2Rad(lon);
	double radius = EARTH_RADIUS + height;
	x = radius * cos(lat) * cos(lon);
	y = radius * cos(lat) * sin(lon);
	z = radius * sin(lat);
};

//以笛卡尔坐标系 计算两点之间距离
double computeDistance(double x1, double y1, double z1, double x2, double y2, double z2)
{
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
};

// 使用经纬度 计算空中目标和地面目标两点之间的距离
double CalculateDistanceEarthToAir(double air_lon, double air_lat, double grd_lon, double grd_lat, double airHeight)
{														//	飞行器经度			飞行器纬度		地面终端经度		地面终端纬度		飞行器高度
	// 计算空中飞行物与地面终端两点的距离（km）
	float RE = 6371.0;			//地球半径 km

	// 由度数转换为弧度
	double lon_A = air_lon * M_PI / 180;		
	double lat_A = air_lat * M_PI / 180;
	double lon_G = grd_lon * M_PI / 180;
	double lat_G = grd_lat * M_PI / 180;

	double dlon = lon_A - lon_G;
	double dlat = lat_A - lat_G;

	// double a = ((sin(dlat / 2)) * (sin(dlat / 2))) + cos(air_lat) * cos(grd_lat) * (sin(dlon / 2) * sin(dlon / 2));		
	// 输入的应该是弧度，而不是度数，这一句中有错误
	double a = ((sin(dlat / 2)) * (sin(dlat / 2))) + cos(lat_A) * cos(lat_G) * (sin(dlon / 2) * sin(dlon / 2));
	if (a < -1) {
		a = -1;
	}
	else if (a > 1) {
		a = 1;
	}
	// 20240813 -nan(ind)出现这个结果的原因就是因为a会是一个负数，而sqrt没办法计算负数
	double C = 2 * asin(sqrt(a));			
	double disEarthToEarth = RE * C;

	double dis = sqrt(disEarthToEarth * disEarthToEarth + airHeight * airHeight);

	return dis;
}

// 使用经纬度 计算地面目标和地面目标两点之间的距离
// 20240815 新增，通过经纬度计算地面两点之间的距离
double CalculateDistanceEarthToEarth(double air_lon, double air_lat, double grd_lon, double grd_lat)		
{																//	飞行器经度		飞行器纬度		地面终端经度	地面终端纬度
	//计算空中飞行物与地面终端两点的距离（km）
	float RE = 6371.0;			//地球半径 km

	double lon_A = air_lon * M_PI / 180;		//由度数转换为弧度
	double lat_A = air_lat * M_PI / 180;
	double lon_G = grd_lon * M_PI / 180;
	double lat_G = grd_lat * M_PI / 180;

	double dlon = lon_A - lon_G;
	double dlat = lat_A - lat_G;

	// double a = ((sin(dlat / 2)) * (sin(dlat / 2))) + cos(air_lat) * cos(grd_lat) * (sin(dlon / 2) * sin(dlon / 2));			
	// 输入的应该是弧度，而不是度数，这一句中有错误
	double a = ((sin(dlat / 2)) * (sin(dlat / 2))) + cos(lat_A) * cos(lat_G) * (sin(dlon / 2) * sin(dlon / 2));
	if (a < -1) {
		a = -1;
	}
	else if (a > 1) {
		a = 1;
	}
	// 20240813 -nan(ind)出现这个结果的原因就是因为a会是一个负数，而sqrt没办法计算负数
	double C = 2 * asin(sqrt(a));			
	double disEarthToEarth = RE * C;

	return disEarthToEarth;
}

// 自由空间传播损耗				
//20240814 MCS修改：增加非法参数判断
double CalculateFSPL(double Distance, double fre)
{
	// 判断参数是否为非法参数
	if (Distance <= 0 || fre <= 0)
	{
		//非法的传入参数
		return 0;
	}
	else
	{
		//计算电磁波在自由空间内的路径损耗（dB）（距离单位km，频率单位MHz）
		double FSPL = 20.0 * log10(Distance) + 20.0 * log10(fre * 1000) + 32.45;
		return FSPL; // dB
	}
}

// 计算噪声功率(3.4)				// 20240814 MCS修改：将功率单位转换为mW再编写函数
double CalNoisePower(double BW)
{	// BW带宽：MHz
	//计算噪声功率 sysNoiseTem：系统噪声温度 ， BW：系统带宽 噪声功率 = 玻尔兹曼常数 * 系统噪声温度 * 通信信号带宽
	double Bolz = 1.23e-23;	//玻尔兹曼常数约为1.23*10^-23J/K
	double sysNoiseTem = 25 + 273.15;		//系统噪声温度 一般是室温25°，转换为开尔文温度
	double noisePower = Bolz * sysNoiseTem * BW * 1e6;
	noisePower = noisePower * 1000;			// 由于接口中单位改为dBm，噪声功率转换为mW
	return noisePower;  //(mW)
}

// 计算接收端信号功率 （3.5）		// 20240814 MCS修改：将功率单位转换为dBm再编写函数
double CalRecSignalPower(double EIRP, double recAntGain, double FSPL, double ClimaticLoss)
{	//                    EIRP(dBm)	  接收天线增益（dB） 自由传播路径损耗(dB)	气候损耗(dB)
	// 接收端信号功率 = EIRP * 接收天线增益 / 自由空间损耗 * 气候损耗（卫星端EIRP=发射功率+发射端天线增益）
	double recSignalPower = pow(10, EIRP / 10) * pow(10, recAntGain / 10) / pow(10, FSPL / 10) * pow(10, ClimaticLoss / 10);
	return recSignalPower;		// 接收端信号功率(mW)(修改后单位由W--->mW)
}

// 计算接收端干扰功率 （3.6）		// 20240814 MCS修改：将功率单位转换为dBm再编写函数
double CalRecDisturbPower(double jamSendPower, double recAntGain, double jam2RecFSPL)
{	//					   干扰端发射功率dBm	接收端天线增益dB	干扰信号的FSPL(dB)
	//接收端干扰功率 = 发射功率 * 接收端天线增益 / 干扰机到接收端传播损耗
	double recDisturbPower = pow(10, jamSendPower / 10) * pow(10, recAntGain / 10) / pow(10, jam2RecFSPL / 10);
	return recDisturbPower;		// 接收端干扰功率(m)
}

// 计算干信比：干扰信号功率/噪声信号功率(3.7)			//20240814 MCS修改：增加非法参数判断
double CalJSR(double signalPower, double jamPower)
{	//			信号功率W			噪声功率W
	double Ratio = jamPower / signalPower;
	if (Ratio <= 0)
	{
		// 传入的信号功率或噪声功率有非法数值
		return 0;
	}
	else
	{
		double JSR = 10 * log(jamPower / signalPower);
		return JSR;					//干信比 dB
	}
}

// 计算扩频增益						//20240814 MCS修改：增加非法参数判断
double CalDSSSGain(int DSRatio)
{//计算扩频增益（dB）	若输入的扩频比为1，则不需扩频，为其他值即可计算扩频增益
	double DSSSGain;
	if (DSRatio <= 0)
	{
		// 非法的参数
		return 0;
	}
	else if (DSRatio == 1)
	{
		DSSSGain = 0;
	}
	else {
		DSSSGain = 10 * log10(DSRatio);
	}
	return DSSSGain;
}

// 计算跳频增益						//20240814 MCS修改：增加非法参数判断，删除带宽传入参数
double CalFHGain(int FHRatio)
//double CalFHGain(int FHRatio, double BW)
{		//计算跳频增益 若跳频比为1，则不需跳频，为其他值则计算跳频
	double FHGain;
	if (FHRatio <= 0)
	{
		// 非法的传入参数
		return 0;
	}
	else if (FHRatio == 1)
	{
		return FHGain = 0;
	}
	else {
		return FHGain = 10 * log10(FHRatio);
	}

}

//计算信噪比						//20240814 MCS修改：增加非法参数判断
double CalSNR(double signalPower, double noisePower)
{
	double SNR = 0;
	if (isinf(signalPower) || noisePower == 0)
	{
		// 非法的传入参数
		return 0;
	}
	else
	{
		double SNR = signalPower / noisePower;			//线性值,非dB值
		return SNR;
	}

}

// 计算信息速率（3.13）
double CalInfoRate(double BW, double SNR)
{//					通信带宽MHz	  信噪比（线性值）
	double InfoRate = BW * 1e6 * log2(1 + SNR);			//bps
	return InfoRate;
}

// 计算时延
double CalTimeDealy(double disPointToPoint)
{//                  (两点之间的距离，km)
	//计算时延（传播时延）（毫秒）
	double c = 3e8;
	double timeDelay = (disPointToPoint * 1000 / c) * 1000;
	return timeDelay;
}

// SNR转换为Eb/N0					//20240814 MCS修改：增加非法参数判断
//double CalBER(double SNR, double BW, double Rb)
double CalEbN0(double SNR, double BW, double Rb)
{//             SNR信噪比 BW为通信带宽Mbps，Rb为信息速率bps
	// 将SNR转换为Eb/N0(单位dB)	 BW为通信带宽，Rb为信息速率
	if (Rb == 0)
	{
		// 非法的传入参数
		return 0;
	}
	else
	{
		double EbN0 = SNR * BW * 1e6 / Rb;
		return EbN0;
	}
}

// 计算干扰之后对信号影响			// 20240814 此函数估计用不太上
double CalDisturb(double BW, double JamPower, double JamBW, double EbN0, double recGain, double jamFSPL)
//实时输入中给的是干扰机的上下频率，未给出带宽
{
	//				通信带宽MHz		干扰功率dBW    干扰带宽MHz		Eb/N0	   接收端天线增益	干扰信号路径传播损耗

	// 噪声功率
	double k = 1.381e-23;	//玻尔兹曼常数(Boltzmann constant),约为 1.38*10^-23 J/K
	double Ta = 298.15;		//系统噪声温度 (Antenna temperature),单位是开尔文(K) 采用室温25度，即298.15K
	double N = k * Ta * (BW * 1e6);	//噪声功率W

	// 噪声功率谱密度N0和Eb
	double N0 = N / BW * 1e6;
	double Eb = EbN0 * N0;

	// 接收端干扰信号的功率谱密度
	double Pj = JamPower * pow(10, recGain / 10) / pow(10, jamFSPL / 10);		//接收端干扰功率W
	double Nj = Pj / JamBW * 1e6;				//干扰信号的功率谱密度

	// 受到干扰后的EbN0
	double EbN0j = Eb / (N0 + Nj);
	return EbN0j;
}

// 通信质量判断
int calLinkCommLevel(double delay, double infoRate)
{
	int commLevel = 1;
	// infoRate = infoRate / 1e6;	// 20240813 MCS修改，将传入进来的infoRate转换为Mbps才能进行函数中的判断  // 20240815 修改，不需要进行转换了，直接传入的参数单位就是Mbps
	if (delay < 200 && infoRate > 10) {
		commLevel = 0;						//高：时延<200ms,信息速率>10Mbps；
	}
	else if (delay < 300 && infoRate > 1) {
		commLevel = 1;						//中：时延<300ms,信息速率>1Mbps；
	}
	else if (delay < 500 && infoRate < 1) {
		commLevel = 2;						//低：时延<500ms,信息速率<1Mbps
	}
	else {
		commLevel = 1;						//默认值，表示没有匹配到任何条件
	}
	return commLevel;
}

//根据BER计算重要信息反击ZD(0高，1中，2低)
int calSatReliability(double ber)
{
	int reliability;
	if (ber < 1e-5) {
		reliability = 0;						//高：BER<10^-5；
	}
	else if (1e-5 < ber && ber < 1e-4)
	{
		reliability = 1;						//中：BER<10^-4；
	}
	else if (1e-4 < ber && ber < 1e-3)
	{
		reliability = 2;						//低：BER<10^-3
	}
	else
	{
		reliability = 1;						//默认值，表示没有匹配到任何条件
	}
	return reliability;
}

// 干扰机J2000坐标转换为经纬度

// 将弧度转换为度数 pi-->180°
double radToDeg(double rad) {
	return rad * 180.0 / M_PI;
}
// 将J2000坐标转换为经纬度
vector<double> j2000ToWGS84(double J2000Pos[3])
{	//对这个函数在原来的版本上进行了修改，将传入参数改为了一个double[3]数组，再从中取出使用
	// 简化的J2000转换旋转矩阵
	double j2000Matrix[3][3] = {
		{0.999999999999, -0.000000440360, -0.000000190919},
		{0.000000440360, 0.999999999998, -0.000000344478},
		{0.000000190919, 0.000000344478, 0.999999999999}
	};

	// 计算地心笛卡尔坐标
	double x = j2000Matrix[0][0] * J2000Pos[0] + j2000Matrix[0][1] * J2000Pos[1] + j2000Matrix[0][2] * J2000Pos[2];
	double y = j2000Matrix[1][0] * J2000Pos[0] + j2000Matrix[1][1] * J2000Pos[1] + j2000Matrix[1][2] * J2000Pos[2];
	double z = j2000Matrix[2][0] * J2000Pos[0] + j2000Matrix[2][1] * J2000Pos[1] + j2000Matrix[2][2] * J2000Pos[2];

	double p = sqrt(x * x + y * y);
	double theta = atan2(z * EARTH_RADIUS, p * EARTH_RADIUS * (1 - FLATTENING));
	double sinTheta = sin(theta);
	double cosTheta = cos(theta);

	double latitude = atan2(z + ECCENTRICITY_SQUARED * EARTH_RADIUS * sinTheta * sinTheta * sinTheta, p - ECCENTRICITY_SQUARED * EARTH_RADIUS * cosTheta * cosTheta * cosTheta);
	double longitude = atan2(y, x);
	double N = EARTH_RADIUS / sqrt(1 - ECCENTRICITY_SQUARED * sin(latitude) * sin(latitude));
	double height = p / cos(latitude) - N;

	// 将高度调整为正值，表示相对于地球表面的海拔高度
	height = fabs(height);

	// 将纬度和经度从弧度转换为度数
	latitude = radToDeg(latitude);
	longitude = radToDeg(longitude);

	vector<double> WGSPos;
	WGSPos.push_back(longitude);
	WGSPos.push_back(latitude);
	WGSPos.push_back(height);
	return WGSPos;
}

// 计算地面站和卫星之间的仰角
double calElevation(double satHeight, double satLon, double satLat, double grdLon, double grdLat)
{//                  卫星高度(km)       卫星经度         卫星纬度        地面站经度      地面站纬度
	double RE = 6371.0;             //地球半径km
	double rs = RE + satHeight;

	// 角度转换为弧度
	double satLonDeg = satLon * M_PI / 180;
	double satLatDeg = satLat * M_PI / 180;
	double grdLonDeg = grdLon * M_PI / 180;
	double grdLatDeg = grdLat * M_PI / 180;

	// 计算仰角
	double gamma1 = sin(satLatDeg) * sin(grdLatDeg) + cos(grdLatDeg) * cos(satLatDeg) * cos(satLonDeg - grdLonDeg);
	if (gamma1 < -1)
	{
		gamma1 = -1;
	}
	else if (gamma1 > 1)
	{
		gamma1 = 1;
	}
	double gamma = acos(gamma1);
	double eleangle = acos(sin(gamma) / sqrt(1 + pow((RE / rs), 2) - 2 * (RE / rs) * cos(gamma)));
	double elevation = radToDeg(eleangle);              //由弧度转换为度数
	return elevation;
}

//*************************************************环境衰减*************************************************//
//1.计算雨衰包含以下函数**********
//
//struct TableRow {     // 定义表格数据结构     
//	double frequency; // 频率
//	double kH;        // 水平极化参数 k
//	double alphaH;    // 水平极化参数 alpha
//	double kV;        // 垂直极化参数 k
//	double alphaV;    // 垂直极化参数 alpha
//};
//
//// 从文件中读取数据，并存储为二维向量
//vector<vector<double>> readFile(const string& filename) {
//	std::ifstream file(filename); // 打开文件
//	std::vector<std::vector<double>> data; // 存储数据
//	std::string line; // 用于读取文件行
//	while (std::getline(file, line)) { // 逐行读取文件
//		std::stringstream ss(line); // 将行转换为字符串流
//		std::vector<double> row; // 存储行数据
//		double value; // 存储每个值
//		while (ss >> value) { // 读取每个值
//			row.push_back(value); // 添加到行数据
//		}
//		data.push_back(row); // 添加到数据
//	}
//	return data; // 返回读取的数据
//}
//
//// 双线性插值函数
//double bilinearInterpolation(double x, double y, double x1, double y1, double x2, double y2, double Q11, double Q12, double Q21, double Q22) {
//	double R1 = ((x2 - x) / (x2 - x1)) * Q11 + ((x - x1) / (x2 - x1)) * Q21; // 计算 R1
//	double R2 = ((x2 - x) / (x2 - x1)) * Q12 + ((x - x1) / (x2 - x1)) * Q22; // 计算 R2
//	double P = ((y2 - y) / (y2 - y1)) * R1 + ((y - y1) / (y2 - y1)) * R2; // 计算 P
//	return P; // 返回插值结果
//}
//
//// 查找最近的网格点
//void findNearestGridPoints(const vector<vector<double>>& latitudes, const vector<std::vector<double>>& longitudes, double inputLat, double inputLon,
//	int& x1, int& y1, int& x2, int& y2) {
//	double minDist = std::numeric_limits<double>::max(); // 最小距离初始化为最大值
//	for (size_t i = 0; i < latitudes.size() - 1; ++i) { // 遍历纬度
//		for (size_t j = 0; j < latitudes[i].size() - 1; ++j) { // 遍历经度
//			double lat1 = latitudes[i][j]; // 获取纬度点1
//			double lat2 = latitudes[i + 1][j]; // 获取纬度点2
//			double lon1 = longitudes[i][j]; // 获取经度点1
//			double lon2 = longitudes[i][j + 1]; // 获取经度点2
//			if (lat1 <= inputLat && lat2 >= inputLat && lon1 <= inputLon && lon2 >= inputLon) { // 如果找到最近点
//				x1 = i; y1 = j; x2 = i + 1; y2 = j + 1; // 更新最近点
//				return; // 退出函数
//			}
//			double dist = std::hypot(inputLat - lat1, inputLon - lon1); // 计算距离
//			if (dist < minDist) { // 如果距离更小
//				x1 = i; y1 = j; // 更新点1
//				x2 = i + 1; y2 = j + 1; // 更新点2
//				minDist = dist; // 更新最小距离
//			}
//		}
//	}
//}
//
//// 读取文本数据
//vector<TableRow> readTxtData(const string& filePath) {
//	std::vector<TableRow> tableData; // 存储表格数据
//	std::ifstream file(filePath); // 打开文件
//	std::string line; // 用于读取文件行
//
//	if (!file.is_open()) { // 如果文件无法打开
//		std::cerr << "无法打开文件: " << filePath << std::endl; // 输出错误信息
//		return tableData; // 返回空数据
//	}
//
//	while (std::getline(file, line)) { // 逐行读取文件
//		if (line[0] == '#') continue; // 跳过注释行
//		std::istringstream iss(line); // 将行转换为字符串流
//		TableRow rowData = { 0 }; // 存储行数据
//		if (iss >> rowData.frequency >> rowData.kH >> rowData.alphaH >> rowData.kV >> rowData.alphaV) { // 读取行数据
//			tableData.push_back(rowData); // 添加到表格数据
//		}
//	}
//
//	file.close(); // 关闭文件
//	return tableData; // 返回表格数据
//}
//
//// 查找最接近的行数据
//TableRow findClosestRow(const vector<TableRow>& tableData, double frequency) {
//	auto it = std::min_element(tableData.begin(), tableData.end(),
//		[frequency](const TableRow& a, const TableRow& b) {
//			return std::abs(a.frequency - frequency) < std::abs(b.frequency - frequency); // 查找最接近的频率
//		});
//	return *it; // 返回最接近的行数据
//}

//// 计算 k
//double calculate_k(double kH, double kV) {
//	return (kH + kV) / 2.0;
//}

// 计算 alpha
//double calculate_a(double kH, double kV, double aH, double aV, double k) {
//	return (kH * aH + kV * aV) / (2.0 * k);
//}
//
// 计算 gamma_R
double calculate_gamma_R(double k, double R, double alpha) {
	return k * pow(R, alpha);
}

// 计算 Ls
double calculateLs(double hR, double hS, double theta, double Re) {
	if (theta >= 5) {
		return (hR - hS) / sin(theta * M_PI / 180.0);
	}
	else {
		double sinTheta = sin(theta * M_PI / 180.0);
		double numerator = 2 * (hR - hS);
		double denominator = sqrt(pow(sinTheta, 2) + 2 * (hR - hS) / Re) + sinTheta;
		return numerator / denominator;
	}
}

// 计算 Lg
double calculateLg(double Ls, double theta) {
	return Ls * cos(theta * M_PI / 180.0);
}

// 计算 r001
double calculateR001(double Lg, double gammaR, double f) {
	double exponent = -2 * Lg;
	double numerator = 1;
	double denominator = 1 + 0.78 * sqrt(Lg * gammaR / f) - 0.38 * (1 - exp(exponent));
	return numerator / denominator;
}

// 计算 zeta
double calculateZeta(double hR, double hS, double Lg, double r001) {
	return atan(((hR - hS) / (Lg * r001)) * M_PI / 180.0);
}

// 计算 Lr
double calculateLr(double hR, double hS, double Lg, double r001, double theta, double zeta) {
	if (zeta > 0) {
		return Lg * r001 / cos(theta * M_PI / 180.0);
	}
	else {
		return (hR - hS) / sin(theta * M_PI / 180.0);
	}
}

// 计算 v001
double calculateV001(double theta, double phi, double Lr, double gammaR, double f) {
	double chi = fabs(phi) < 36 ? 36 - fabs(phi) : 0;
	double a = theta * M_PI / 180.0;
	double sinTheta = sin(a);
	double x = chi * M_PI / 180.0;
	double b = -(theta / (1 + chi));
	double numerator = 1;
	double denominator = 1 + (pow(sinTheta, 0.5) * (31 * (1 - exp(b)) * pow(Lr * gammaR, 0.5) / (f * f) - 0.45));
	return numerator / denominator;
}

// 计算 Le
double calculateLe(double Lr, double v001) {
	return Lr * v001;
}

// 计算 A001
double calculateA001(double gammaR, double Le) {
	return gammaR * Le;
}

// 计算 Ap
double calculateAp(double A001, double phi, double theta) {
	double p = 0.01;
	double beta;
	if (p >= 1 || fabs(phi) >= 36) {
		beta = 0;
	}
	else if (p < 1 && fabs(phi) < 36 && theta >= 25) {
		beta = -0.005 * (fabs(phi) - 36);
	}
	else {
		beta = -0.005 * (fabs(phi) - 36) + 1.8 - 4.25 * sin(theta * M_PI / 180.0);
	}
	return A001 * pow(p / 0.01, -(0.655 + 0.033 * log(p) - 0.045 * log(A001) - beta * (1 - p) * sin(theta * M_PI / 180.0)));
}

// 1.雨衰
double calculateRainAttenuation(double inputLatitude, double inputLongitude, double Rain, double hs, double theta, double f) {
	//                              地面站纬度             地面站经度     年均单点降雨量(mm/h)  地面站高度     仰角   信号频率(GHz)  
	double Re = 8500; // 地球半径
	double Ap = 0.0; // 初始衰减
	double hR = 3.30; //年平均降雨高度
	double k = 0.03; //频率相关参数
	double a = 1.10; // 频率相关参数

	if (hR - hs <= 0 || Rain == 0) { // 如果高度差小于等于0或无降雨量
		return 0; // 返回0
	}
	else {
		double Ls = calculateLs(hR, hs, theta, Re); // 计算 Ls
		double Lg = calculateLg(Ls, theta); // 计算 Lg

		double R = Rain; // 定义降雨量
		double frequency = f; // 定义频率

		double gammaR = calculate_gamma_R(k, R, a); // 计算 gammaR

		double r001 = calculateR001(Lg, gammaR, f); // 计算 r001
		double zeta = calculateZeta(hR, hs, Lg, r001); // 计算 zeta
		double Lr = calculateLr(hR, hs, Lg, r001, theta, zeta); // 计算 Lr
		double v001 = calculateV001(theta, inputLatitude, Lr, gammaR, f); // 计算 v001

		double Le = calculateLe(Lr, v001); // 计算 Le
		double A001 = calculateA001(gammaR, Le); // 计算 A001
		Ap = calculateAp(A001, inputLatitude, theta); // 计算 Ap
	}

	return Ap; // 返回 Ap
}


// 2.云衰**********
double calCloudAttenuation(double frequency, double temperature, double airPressure, double liquidWaterDensity, double elevationAngle, double cloudVerticalThickness)
{   //                       信号频率(GHz)      温度(摄氏度)        气压(百帕hPa)      云雾的液态水密度（g/m³）        仰角（度）               云层垂直厚度(km)
	//介电常数
	// 根据ITU-R P.840建议书计算介电常数的实部和虚部
	double epsilon_s = 77.66 + 103.3 * temperature / (temperature + 273.15); // 静态介电常数
	double epsilon_inf = 5.48; // 高频极限介电常数
	double alpha = 0.00041 * (temperature + 273.15); // 温度相关系数
	double beta = 87.91 - 0.404 * (temperature + 273.15) + 0.00087 * (temperature + 273.15) * (temperature + 273.15); // 频率相关系数
	double delta = 0.0001 * airPressure; // 压力相关系数

	double epsilon_real = epsilon_s - (epsilon_s - epsilon_inf) / (1 + pow(frequency / beta, 2)) + delta;           //介电常数实部
	double epsilon_imag = (epsilon_s - epsilon_inf) * (frequency / beta) / (1 + pow(frequency / beta, 2)) + alpha * frequency;      //介电常数虚部

	//云雾衰减系数
	double gamma_m = 0.819 * frequency * liquidWaterDensity * (epsilon_real - 1) / (pow(epsilon_real - 1, 2) + pow(epsilon_imag, 2));

	//路径长度
	double pathLength = cloudVerticalThickness / sin(elevationAngle * M_PI / 180.0);

	//云雾衰减（B）d
	double cloudAttenuation = gamma_m * pathLength;
	return cloudAttenuation;
}

// 3.气衰**********
double calGasAttenuation(double f, double P, double T, double rho) {
	//                频率(GHz)  气压(hpa百帕) 温度(摄氏度)   水汽压(g/m^3)
	double T_kelvin = T + 273.15;
	double theta = 300.0 / T_kelvin;
	double P_d = P * 100.0; // 转换为Pa
	double e = rho * T / 216.7; // 水蒸气分压

	// 氧气的吸收线参数
	struct oxygen_line {
		double f_center; // 中心频率
		double a1; // 参数a1
		double a2; // 参数a2
		double a3; // 参数a3
		double a4; // 参数a4
		double a5; // 参数a5
		double a6; // 参数a6
	};

	struct oxygen_line oxygen_lines[] = {
		{60.4348, 2438.000, 0.386, 13.390, 0.0, 6.342, -2.825},
		{118.7503, 940.300, 0.010, 16.640, 0.0, -0.439, 0.079},
		{368.4982, 67.400, 0.048, 16.400, 0.0, 0.000, 0.000}
	};

	// 水汽的吸收线参数
	struct water_vapor_line {
		double f_center; // 中心频率
		double b1; // 参数b1
		double b2; // 参数b2
		double b3; // 参数b3
		double b4; // 参数b4
		double b5; // 参数b5
		double b6; // 参数b6
	};

	struct water_vapor_line water_vapor_lines[] = {
		{22.2351, 0.1079, 2.144, 26.38, 0.76, 5.087, 1.00},
		{183.3101, 2.273, 0.668, 29.06, 0.77, 5.022, 0.85},
		{325.153, 1.514, 1.541, 28.23, 0.64, 4.893, 0.74}
	};

	double oxygen_sum = 0.0;
	double water_sum = 0.0;

	// 计算氧气引起的衰减
	for (int i = 0; i < sizeof(oxygen_lines) / sizeof(oxygen_lines[0]); i++) {
		struct oxygen_line line = oxygen_lines[i];
		double df = 0.002 * P_d * pow(theta, 0.8);
		double S = 1e-7 * P_d * line.a1 * pow(theta, 3) * exp(line.a2 * (1 - theta));
		double delta = line.a5 * P_d * pow(theta, line.a6);
		double gamma = line.a3 * P_d * pow(theta, line.a4);
		double f_i = line.f_center;

		double term1 = S / ((f_i * f_i) - (f * f) + delta * delta);
		double term2 = gamma * (f * f) / ((f_i * f_i) - (f * f) + gamma * gamma + delta * delta);
		oxygen_sum += term1 + term2;
	}

	// 计算水汽引起的衰减
	for (int i = 0; i < sizeof(water_vapor_lines) / sizeof(water_vapor_lines[0]); i++) {
		struct water_vapor_line line = water_vapor_lines[i];
		double df = 0.003 * P_d * pow(theta, 0.8);
		double S = 1e-1 * e * line.b1 * pow(theta, 3.5) * exp(line.b2 * (1 - theta));
		double delta = line.b5 * P_d * pow(theta, line.b6);
		double gamma = line.b3 * P_d * pow(theta, line.b4);
		double f_i = line.f_center;

		double term1 = S / ((f_i * f_i) - (f * f) + delta * delta);
		double term2 = gamma * (f * f) / ((f_i * f_i) - (f * f) + gamma * gamma + delta * delta);
		water_sum += term1 + term2;
	}

	// 计算氧气的衰减
	double N_pp_oxygen = 6.14e-5 * P_d * pow(theta, 2.5);
	double A_o = 0.0118 * N_pp_oxygen * (1 - 0.023 * log10(N_pp_oxygen));
	double A_o_f = A_o * (1 + 0.012 * (f / 50.0));
	double oxygen_attenuation = 0.1820 * f * oxygen_sum + A_o_f;

	// 计算水汽的衰减
	double N_pp_water = 1.51e-4 * rho * pow(theta, 4.5);
	double A_w = 0.0011 * N_pp_water * (1 - 0.016 * log10(N_pp_water));
	double A_w_f = A_w * (1 + 0.007 * (f / 50.0));
	double water_attenuation = 0.1820 * f * water_sum + A_w_f;

	// 返回总衰减
	return oxygen_attenuation + water_attenuation;
}

// 4.闪烁与多径衰减*******
double getSignalAttenuation() {
	bool hasCloud = 1;       // 是否有云的布尔值
	int season = 2;
	double attenuation = 0.0;

	// 如果有云，增加衰减值
	if (hasCloud) {
		attenuation += 2.0;
	}

	// 根据不同的季节增加不同的衰减值
	if (season == 0) {           //0 春季
		attenuation += 0.2;
	}
	else if (season == 1) {      //1 夏季
		attenuation += 1.0;
	}
	else if (season == 2) {      //2 秋季
		attenuation += 0.3;
	}
	else if (season == 3) {      //3 冬季
		attenuation += 0.2;
	}
	else {
		std::cerr << "无效的季节输入！" << std::endl;
		return -1.0; // 返回错误值
	}

	return attenuation;
}

// 总的气候损耗
double CalClimaticLoss(double frequency, double satHeight, double satLon, double satLat, double grdLon, double grdLat, double grdHeight,
	double Rain, double temperature, double airPressure, double liquidWaterDensity, double cloudVerticalThickness, double rho)
{
	//                  通信频率GHz         卫星高度(km)      卫星经度      卫星纬度       地面站经度     地面站纬度      地面站高度  
 //年均单点降雨量(mm/h) 温度(摄氏度)     气压(百帕hPa)    云雾的液态水密度（g/m³）        云层垂直厚度(km)        水汽压(g/m^3)
 //liquidWaterDensity :云雾的液态水密度通常在0.1到1.0 g / m³之间。 
 // rho:在标准条件下（0°C, 1 atm），水汽压大约是 6.11 mbar（611 Pa）。计算出的质量浓度约为 7.5 g/m³。
	// 计算仰角
	double el = calElevation(satHeight, satLon, satLat, grdLon, grdLat);
	// 雨衰
	double RainAttenuation = calculateRainAttenuation(grdLat, grdLon, Rain, grdHeight, el, frequency);
	// 云衰
	double CloudAttenuation = calCloudAttenuation(frequency, temperature, airPressure, liquidWaterDensity, el, cloudVerticalThickness);
	// 气衰
	double GasAttenuation = calGasAttenuation(frequency, airPressure, temperature, rho);
	// 闪烁与多径衰减
	double flickerAttenuation = getSignalAttenuation();

	//总衰减
	double ClimaticLoss = RainAttenuation + CloudAttenuation + GasAttenuation + flickerAttenuation;
	return ClimaticLoss;
}

/********************************************气候衰减*********************************************************/

// Q函数
double Q(double x)
{
	return 0.5 * erfc(x / sqrt(2));
}

// 根据调制方式不同计算BER
double CalBEROfDiffMod(int modType, double EbN0)
{
	double BER;
	switch (modType)
	{
	case 0:	//BPSK
		BER = Q(sqrt(2 * EbN0));
		break;
	case 1:	//DPSK
		BER = 0.5 * exp(-EbN0);
		break;
	case 2:	//QPSK
		BER = Q(sqrt(EbN0));
		break;
	case 3:	//16QAM
		BER = 3 * Q(sqrt((4.0 / 5) * EbN0));
		break;
	case 4:	//2FSK
		BER = Q(sqrt(EbN0));
		break;
	default:
		BER = Q(sqrt(EbN0));
		break;
	}
	return BER;
}

// 波束功率初始化
vector<double> iniBeamsPower(double totalPower, int numBeams)
{
	totalPower = pow(10, totalPower / 10);      //dBW--->W
	vector<double> beamsPower;
	beamsPower.resize(numBeams, totalPower / numBeams);
	return beamsPower;
}

// 波束功率固定增强/减弱10%
vector<double> enhanceBeam(int index, int operation, double totalPower, vector<double>& beamsPower)
{
	// index：波束索引(要改变的波束)     operation：操作索引（0：不变，1：增强 ， 2：减弱）     totalPower：总功率       totalPower：存放各个波束功率
	if (operation == 0) {
		// 功率不变，直接返回
		return beamsPower;
	}
	else if (operation == 1) {
		// 增强指定波束的功率
		beamsPower[index] += beamsPower[index] * 0.1;
		totalPower += beamsPower[index] * 0.1; // 更新总功率
	}
	else if (operation == 2) {
		// 减小指定波束的功率
		beamsPower[index] -= beamsPower[index] * 0.1;
		totalPower -= beamsPower[index] * 0.1; // 更新总功率
	}
	else {
		return beamsPower;              // 未匹配到任何数值
	}

	// 重新分配剩余波束的功率
	double remainingPower = totalPower - beamsPower[index]; // 剩余的总功率
	int remainingBeams = beamsPower.size() - 1; // 剩余的波束数量

	for (size_t i = 0; i < beamsPower.size(); ++i) {
		if (i != static_cast<size_t>(index)) {
			beamsPower[i] = remainingPower / remainingBeams;
		}
	}
	return beamsPower; // 返回修改后的 vector
}

// 随机的故障概率函数，用卫星终端参数中，模仿自然情况下故障概率
double breakdownFun(double number) {
	// 参数             载荷故障概率（范围(0,1)，可定典型值，模拟由空间环境因素引起的故障）
	// 返回值若为1，则说明产生了自然故障，若为0，则载荷正常
	   // 判断一个小数后面有几位
	ostringstream oss;
	oss << number;
	string str = oss.str();

	size_t pos = str.find('.');
	if (pos == string::npos) return 0;  // No decimal point

	double decimalNumes = str.length() - pos - 1;

	// 将一个小于1的小数部分转换为整数
	double decimal_part_as_int = static_cast<int>(round(number * pow(10, decimalNumes)));

	// 转换为的整数部分，一次性生成多少个随机数
	vector<int> numbers1;
	vector<int> numbers2;

	// 创建随机数生成器
	random_device rd;
	mt19937 gen(rd());

	// 定义随机数分布范围
	double max = pow(10, decimalNumes);
	uniform_int_distribution<> dis(1, max);

	// 生成第一个随机数数组
	for (int i = 0; i < decimal_part_as_int; ++i) {
		numbers1.push_back(dis(gen));
	}

	// 生成第二个随机数数组
	for (int i = 0; i < decimal_part_as_int; ++i) {
		numbers2.push_back(dis(gen));
	}

	// 使用 set 存储第一个随机数数组
	set<int> numberSet(numbers1.begin(), numbers1.end());

	// 检查第二组随机数是否存在于 set 中
	for (int num : numbers2) {
		if (numberSet.find(num) != numberSet.end()) {
			// 找到相同的数字
			return 1;  // 表示找到重复
		}
	}

	// 没有找到相同的数字
	return 0;  // 表示没有重复

}

// 20240814 MCS新增 编码增益函数
double getCodeGain(int codeType, double codeRate) {
	double codeGain = 0;

	switch (codeType) {      // 根据编码方式不同
	case 0:         // 无编码
		codeGain = 0;
		break;
	case 1:         // 卷积码
		if (codeRate == 1.0 / 2) {
			codeGain = 3.5;
		}
		else if (codeRate == 1.0 / 3) {
			codeGain = 4.9;
		}
		else if (codeRate == 2.0 / 3) {
			codeGain = 2.8;
		}
		break;
	case 2:         // Turbo码
		if (codeRate == 1.0 / 2) {
			codeGain = 2.1;
		}
		else if (codeRate == 1.0 / 3) {
			codeGain = 3.8;
		}
		else if (codeRate == 2.0 / 3) {
			codeGain = 0.9;
		}
		break;
	case 3:         // LDPC码
		if (codeRate == 1.0 / 2) {
			codeGain = 3.1;
		}
		else if (codeRate == 1.0 / 3) {
			codeGain = 5.8;
		}
		else if (codeRate == 2.0 / 3) {
			codeGain = 0.8;
		}
		break;
	default:
		codeGain = 2.1;
		break;
	}

	return codeGain;
}
