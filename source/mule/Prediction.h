#include "Block4D.h"
#include <iostream>
#include <fstream>

using namespace std;

class Prediction 
{
public:
	ofstream predictionFile, yResiduesFile, cbResiduesFile, crResiduesFile;

	Prediction(int prediction);
	int counter;
	int simplePredictor(Block4D *origBlock);
	void reconstruct4DBlock(Block4D *residueBlock, int DCPredictor);
	void fourRefsPredictor(Block4D *residueBlock, Block4D *origBlock, Block4D *ref0, Block4D *ref1, Block4D *ref2, Block4D *ref3);
	void calculateResidue(Block4D *residueBlock, Block4D *origBlock, int DCPredictor);
	void differentialPredictionRaster(Block4D *residueBlock, Block4D *origBlock, int spectralComponent);
	void recDifferentialPredictionRaster(Block4D *residueBlock, Block4D *reconstructedBlock, int average);
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
};