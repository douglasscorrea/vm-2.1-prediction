#include "Block4D.h"

class Prediction 
{
public:
	Prediction(void);
	int simplePredictor(Block4D *origBlock);
	void fourRefsPredictor(Block4D *residueBlock, Block4D *origBlock, Block4D *ref0, Block4D *ref1, Block4D *ref2, Block4D *ref3);
	void calculateResidue(Block4D *residueBlock, Block4D *origBlock, int DCPredictor);
	void reconstruct4DBlock(Block4D *residueBlock, int DCPredictor);
};