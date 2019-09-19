#include "Block4D.h"
#include <iostream>
#include <fstream>

using namespace std;

class Prediction 
{
public:
	ofstream yResiduesFile, cbResiduesFile, crResiduesFile, firstPlaneCoefficients, allCoefficients;
	int predictionType, counter, maxRefPlane, minRefPlane, maxOtherPlanes, minOtherPlanes;
	double totalSignalEnergyFirstPlane, totalSignalEnergyOtherPlanes;

	Prediction(int prediction);
	int getMaxRefPlane();
	int getMinRefPlane();
	int getMaxOtherPlanes();
	int getMinOtherPlanes();
	int simplePredictor(Block4D *origBlock);
	void reconstruct4DBlock(Block4D *residueBlock, int DCPredictor);
	void fourRefsPredictor(Block4D *residueBlock, Block4D *origBlock, Block4D *ref0, Block4D *ref1, Block4D *ref2, Block4D *ref3);
	void calculateResidue(Block4D *residueBlock, Block4D *origBlock, int DCPredictor);
	void differentialPredictionRaster(Block4D *residueBlock, Block4D *origBlock, int spectralComponent);
	void recDifferentialPredictionRaster(Block4D *residueBlock, Block4D *reconstructedBlock);
	void differentialPredictionCentral(Block4D *residueBlock, Block4D *origBlock, int spectralComponent);
	void recDifferentialPredictionCentral(Block4D *recBlock, Block4D *origBlock);
	void differentialPredictionRasterHalf(Block4D *residueBlock, Block4D *origBlock);
	void recDifferentialPredictionRasterHalf(Block4D *recBlock, Block4D *origBlock);
	void hierarchicalDifferentialPrediction(Block4D *residueBlock, Block4D *origBlock);
	void recHierarchicalDifferentialPrediction(Block4D *residueBlock, Block4D *reconstructedBlock);
	void hierarchicalDifferentialPrediction1Level(Block4D *residueBlock, Block4D *origBlock);
	void recHierarchicalDifferentialPrediction1Level(Block4D *recBlock, Block4D *origBlock);
	void saveSamplesMule(Block4D *residueBlock, Block4D *origBlock, int spectralComponent);
	void printOneBlock(Block4D *lfBlock);
	double calcReferencePlaneEnergy(Block4D *lfBlock);
	double calcOtherPlanesEnergy(Block4D *lfBlock);
};