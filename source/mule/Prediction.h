#include "Block4D.h"

class Prediction 
{
public:
	Prediction(void);
	void simplePredictor(Block4D *targetBlock, int position_t, int position_s, int position_v, int position_u, int component);
};