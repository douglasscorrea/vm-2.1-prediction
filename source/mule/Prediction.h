#include "Block4D.h"

class Prediction 
{
public:
	Prediction(void);
	int simplePredictor(Block4D *origBlock);
	void calculateResidue(Block4D *residueBlock, Block4D *origBlock, int DCPredictor);
};