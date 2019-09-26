#include "Block4D.h"
#include <iostream>
#include <fstream>

using namespace std;

class Statistics 
{
public:
	ofstream yResiduesFile, cbResiduesFile, crResiduesFile, refPlaneCoefficients, allCoefficients, DC_coeff, AC_coeff;
	int predictionType, maxRefPlane, minRefPlane, maxOtherPlanes, minOtherPlanes;
	int yCounterRP, cbCounterRP, crCounterRP, yCounterOP, cbCounterOP, crCounterOP, counterOPCoeffs, counterRPCoeffs, counterAllCoeffs;
	int *yRefPlaneSamples;
	int *cbRefPlaneSamples;
	int *crRefPlaneSamples;
	int *yOtherPlanesSamples;
	int *cbOtherPlanesSamples;
	int *crOtherPlanesSamples;
	int *refPlaneCoeff;
	int *otherPlanesCoeff;
	int *AllCoeffs;

	double y_totalSignalEnergyRefPlane, cb_totalSignalEnergyRefPlane, cr_totalSignalEnergyRefPlane;
	double y_totalSignalEnergyOtherPlanes, cb_totalSignalEnergyOtherPlanes, cr_totalSignalEnergyOtherPlanes;
	double y_samplesSum, cb_samplesSum, cr_samplesSum;
	double y_refPlaneSamplesSum, cb_refPlaneSamplesSum, cr_refPlaneSamplesSum;
	double y_otherPlanesSamplesSum, cb_otherPlanesSamplesSum, cr_otherPlanesSamplesSum;
	double energyRefPlaneCoeff, energyOtherPlanesCoeff, sumRefPlaneCoeff, sumOtherPlanesCoeff;
	double averageRefPlaneCoeff, averageOtherPlanesCoeff, stdRefPlaneCoeff, stdOtherPlanesCoeff;
	double coeffEnergy, sumCoeff;

	Statistics(int prediction);
	void saveSamplesMule(Block4D *residueBlock, Block4D *origBlock, int spectralComponent);
	void printOneBlock(Block4D *lfBlock);
	void calcReferencePlaneEnergy(Block4D *lfBlock, int spectralComponent);
	void calcOtherPlanesEnergy(Block4D *lfBlock, int spectralComponent);
	void saveDCCoeff(Block4D *lfBlock);
	void saveACCoeff(Block4D *lfBlock);
	void calcSumYRefPlaneSamples(Block4D *lfBlock);
	void calcSumCbRefPlaneSamples(Block4D *lfBlock);
	void calcSumCrRefPlaneSamples(Block4D *lfBlock);
	void calcSumYOtherPlanesSamples(Block4D *lfBlock);
	void calcSumCbOtherPlanesSamples(Block4D *lfBlock);
	void calcSumCrOtherPlanesSamples(Block4D *lfBlock);
	int getMaxRefPlane();
	int getMinRefPlane();
	int getMaxOtherPlanes();
	int getMinOtherPlanes();
	double getYRefPlaneEnergy();
	double getCbRefPlaneEnergy();
	double getCrRefPlaneEnergy();
	double getYOtherPlanesEnergy();
	double getCbOtherPlanesEnergy();
	double getCrOtherPlanesEnergy();
	double getSamplesAverage(Block4D lfBlock);
	double getYSamplesAverage(Block4D lfBlock);
	double getCbSamplesAverage(Block4D lfBlock);
	double getCrSamplesAverage(Block4D lfBlock);
	double getYRefPlaneSamplesAverage();
	double getCbRefPlaneSamplesAverage();
	double getCrRefPlaneSamplesAverage();
	double getYOtherPlanesSamplesAverage();
	double getCbOtherPlanesSamplesAverage();
	double getCrOtherPlanesSamplesAverage();
	double calcStdYRefPlane();
	double calcStdCbRefPlane();
	double calcStdCrRefPlane();
	double calcStdYOtherPlanes();
	double calcStdCbOtherPlanes();
	double calcStdCrOtherPlanes();

	void calcCoeffEnergy(Block4D *transformedBlock);
	void calcRefPlaneCoeffEnergy(Block4D *transformedBlock);
	void calcOtherPlanesCoeffEnergy(Block4D *transformedBlock);
	void calcSumCoeff(Block4D *transformedBlock);
	void calcSumRefPlaneCoeff(Block4D *transformedBlock);
	void calcSumOtherPlanesCoeff(Block4D *transformedBlock);
	double getAverageCoeff();
	double getEnergyCoeff();
	double getAverageRefPlaneCoeff();
	double getAverageOtherPlanesCoeff();
	double calcStdCoeff();
	double calcStdRefPlaneCoeff();
	double calcStdOtherPlanesCoeff();
	double getEnergyRefPlaneCoeff();
	double getEnergyOtherPlanesCoeff();
};