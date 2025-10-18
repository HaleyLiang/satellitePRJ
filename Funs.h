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

// ���� M_PI
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ��������ƽ���뾶����λ��ǧ�ף�
const double EARTH_RADIUS = 6371.0;
const double FLATTENING = 1.0 / 298.257223563;
const double ECCENTRICITY_SQUARED = 2 * FLATTENING - FLATTENING * FLATTENING;
// ��������������ͷ�ļ�

// �Ƕ�ת����
double deg2Rad(double deg);

// ��γ��ת��Ϊ�ѿ�������ϵ
void computeCartesian(double lat, double lon, double height, double& x, double& y, double& z);

// ��������֮�����
double computeDistance(double x1, double y1, double z1, double x2, double y2, double z2);

//�������ǵ�����վ�ľ��루��γ�ȷ�ʽ��
double CalculateDistanceEarthToAir(double air_lon, double air_lat, double grd_lon, double grd_lat, double airHeight);

// ʹ�þ�γ�� �������Ŀ��͵���Ŀ������֮��ľ���
double CalculateDistanceEarthToEarth(double air_lon, double air_lat, double grd_lon, double grd_lat);

// ���ɿռ䴫�����
double CalculateFSPL(double DistanceAirToGrd, double fre);

// ������������
double CalNoisePower(double BW);

// ������ն��źŹ��� ��3.5��
double CalRecSignalPower(double EIRP, double recAntGain, double FSPL, double ClimaticLoss);

// ������ն˸��Ź��� ��3.6��
double CalRecDisturbPower(double jamSendPower, double recAntGain, double jam2RecFSPL);

// ������űȣ������źŹ���/�����źŹ���(3.7)
double CalJSR(double signalPower, double noisePower);

// ������Ƶ����
double CalDSSSGain(int DSRatio);

// ������Ƶ����;
double CalFHGain(int FHRatio);

//��������� 
double CalSNR(double signalPower, double noisePower);

// ������Ϣ���ʣ�3.13��
double CalInfoRate(double BW, double SNR);

// ����ʱ��
double CalTimeDealy(double disPointToPoint);

// SNRת��ΪEb/N0
double CalEbN0(double SNR, double BW, double Rb);

//����ʱ�Ӻ���Ϣ���ʼ�������ͨ������
int calLinkCommLevel(double delay, double infoRate);

//����BER������Ҫ��Ϣ����ZD(0�ߣ�1�У�2��)
int calSatReliability(double ber);

// ������Ŷ��ź�Ӱ��
double CalDisturb(double BW, double JamPower, double JamBW, double EbN0, double recGain, double jamFSPL);

// ������ת��Ϊ����
double radToDeg(double rad);

// ��J2000����ת��Ϊ��γ��
vector <double> j2000ToWGS84(double J2000Pos[3]);

// �������
double calElevation(double satHeight, double satLon, double satLat, double grdLon, double grdLat);

//*************************************************����˥��*************************************************//
// 1.��˥
double calculateRainAttenuation(double inputLatitude, double inputLongitude, double Rain, double hs, double theta, double f);

// 2.��˥**********
double calCloudAttenuation(double frequency, double temperature, double airPressure, double liquidWaterDensity, double elevationAngle, double cloudVerticalThickness);

// 3.��˥**********
double calGasAttenuation(double f, double P, double T, double rho);

// 4.��˸�ྶ˥��
double getSignalAttenuation();

// �ܵ��������
double CalClimaticLoss(double frequency, double satHeight, double satLon, double satLat, double grdLon, double grdLat, double grdHeight, double Rain, double temperature, double airPressure, double liquidWaterDensity, double cloudVerticalThickness, double rho);
//*************************************************����˥��*************************************************//

// Q����
double Q(double x);

// ���ݵ��Ʒ�ʽ��ͬ����BER
double CalBEROfDiffMod(int modType, double EbN0);

// �������ʳ�ʼ��
vector<double> iniBeamsPower(double totalPower, int numBeams);

// �������ʹ̶���ǿ/����10%
vector<double> enhanceBeam(int index, int operation, double totalPower, vector<double>& beamsPower);

// ����Ĺ��ϸ��ʺ������������ն˲����У�ģ����Ȼ����¹��ϸ���
double breakdownFun(double number);

// 20240814 MCS���� �������溯��
double getCodeGain(int codeType, double codeRate);

#endif      // ifndef __FUNS_H