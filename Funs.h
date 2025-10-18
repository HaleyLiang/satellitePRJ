#ifndef  __FUNS_H
#define __FUNS_H

#include<cmath>
#include<vector>
#include<string>
#include<iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>
#include <set>

using namespace std;

// 定义 M_PI
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// 定义地球的平均半径（单位：千米）
const double EARTH_RADIUS = 6371.0;
const double FLATTENING = 1.0 / 298.257223563;
const double ECCENTRICITY_SQUARED = 2 * FLATTENING - FLATTENING * FLATTENING;
// 包含公共函数的头文件

// 角度转弧度
double deg2Rad(double deg);

// 经纬度转换为笛卡尔坐标系
void computeCartesian(double lat, double lon, double height, double& x, double& y, double& z);

// 计算两点之间距离
double computeDistance(double x1, double y1, double z1, double x2, double y2, double z2);

//计算卫星到地面站的距离（经纬度方式）
double CalculateDistanceEarthToAir(double air_lon, double air_lat, double grd_lon, double grd_lat, double airHeight);

// 使用经纬度 计算地面目标和地面目标两点之间的距离
double CalculateDistanceEarthToEarth(double air_lon, double air_lat, double grd_lon, double grd_lat);

// 自由空间传播损耗
double CalculateFSPL(double DistanceAirToGrd, double fre);

// 计算噪声功率
double CalNoisePower(double BW);

// 计算接收端信号功率 （3.5）
double CalRecSignalPower(double EIRP, double recAntGain, double FSPL, double ClimaticLoss);

// 计算接收端干扰功率 （3.6）
double CalRecDisturbPower(double jamSendPower, double recAntGain, double jam2RecFSPL);

// 计算干信比：干扰信号功率/噪声信号功率(3.7)
double CalJSR(double signalPower, double noisePower);

// 计算扩频增益
double CalDSSSGain(int DSRatio);

// 计算跳频增益;
double CalFHGain(int FHRatio);

//计算信噪比 
double CalSNR(double signalPower, double noisePower);

// 计算信息速率（3.13）
double CalInfoRate(double BW, double SNR);

// 计算时延
double CalTimeDealy(double disPointToPoint);

// SNR转换为Eb/N0
double CalEbN0(double SNR, double BW, double Rb);

//根据时延和信息速率计算评价通信质量
int calLinkCommLevel(double delay, double infoRate);

//根据BER计算重要信息反击ZD(0高，1中，2低)
int calSatReliability(double ber);

// 计算干扰对信号影响
double CalDisturb(double BW, double JamPower, double JamBW, double EbN0, double recGain, double jamFSPL);

// 将弧度转换为度数
double radToDeg(double rad);

// 将J2000坐标转换为经纬度
vector <double> j2000ToWGS84(double J2000Pos[3]);

// 求解仰角
double calElevation(double satHeight, double satLon, double satLat, double grdLon, double grdLat);

//*************************************************气候衰减*************************************************//
// 1.雨衰
double calculateRainAttenuation(double inputLatitude, double inputLongitude, double Rain, double hs, double theta, double f);

// 2.云衰**********
double calCloudAttenuation(double frequency, double temperature, double airPressure, double liquidWaterDensity, double elevationAngle, double cloudVerticalThickness);

// 3.气衰**********
double calGasAttenuation(double f, double P, double T, double rho);

// 4.闪烁多径衰减
double getSignalAttenuation();

// 总的气候损耗
double CalClimaticLoss(double frequency, double satHeight, double satLon, double satLat, double grdLon, double grdLat, double grdHeight, double Rain, double temperature, double airPressure, double liquidWaterDensity, double cloudVerticalThickness, double rho);
//*************************************************气候衰减*************************************************//

// Q函数
double Q(double x);

// 根据调制方式不同计算BER
double CalBEROfDiffMod(int modType, double EbN0);

// 波束功率初始化
vector<double> iniBeamsPower(double totalPower, int numBeams);

// 波束功率固定增强/减弱10%
vector<double> enhanceBeam(int index, int operation, double totalPower, vector<double>& beamsPower);

// 随机的故障概率函数，用卫星终端参数中，模仿自然情况下故障概率
double breakdownFun(double number);

// 20240814 MCS新增 编码增益函数
double getCodeGain(int codeType, double codeRate);

#endif      // ifndef __FUNS_H