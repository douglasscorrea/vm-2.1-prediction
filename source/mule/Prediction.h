#include "Block4D.h"

class Prediction 
{
public:
	Prediction(void);
	int simplePredictor(Block4D *origBlock);
	void reconstruct4DBlock(Block4D *residueBlock, int DCPredictor);
	void fourRefsPredictor(Block4D *residueBlock, Block4D *origBlock, Block4D *ref0, Block4D *ref1, Block4D *ref2, Block4D *ref3);
	void calculateResidue(Block4D *residueBlock, Block4D *origBlock, int DCPredictor);
	void differentialPredictionRaster(Block4D *residueBlock, Block4D *origBlock);
	void recDifferentialPredictionRaster(Block4D *recBlock, Block4D *origBlock);
	void differentialPredictionCentral(Block4D *residueBlock, Block4D *origBlock);
	void recDifferentialPredictionCentral(Block4D *recBlock, Block4D *origBlock);
	void hierarchicalDifferentialPrediction(Block4D *residueBlock, Block4D *origBlock);
	void recHierarchicalDifferentialPrediction(Block4D *recBlock, Block4D *origBlock);
	void hierarchicalDifferentialPrediction1Level(Block4D *residueBlock, Block4D *origBlock);
	void recHierarchicalDifferentialPrediction1Level(Block4D *recBlock, Block4D *origBlock);

};