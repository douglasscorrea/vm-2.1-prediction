#include "Prediction.h"
#include "Block4D.h"
#include <stdio.h>
#include <stdlib.h>

Prediction :: Prediction() {

}

int Prediction :: simplePredictor(Block4D *origBlock) {
	return origBlock->sumAllPixels()/(13*13*15*15);
}

void Prediction :: calculateResidue(Block4D *residueBlock, Block4D *origBlock, int DCPredictor) {
	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
		for(int horizontalView = 0; horizontalView < 13; horizontalView += 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
					residueBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] = 
						origBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] - DCPredictor;
				}
			}
		}
	}
	//printf("Average: %d\n", targetBlock->sumAllPixels());
}

void Prediction :: reconstruct4DBlock(Block4D *residueBlock, int DCPredictor) {
	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
		for(int horizontalView = 0; horizontalView < 13; horizontalView += 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
					residueBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] += DCPredictor;
				}
			}
		}
	}
}

void Prediction :: fourRefsPredictor(Block4D *residueBlock, Block4D *origBlock, Block4D *ref0, Block4D *ref1, Block4D *ref2, Block4D *ref3) {
	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
		for(int horizontalView = 0; horizontalView < 13; horizontalView += 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
					residueBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] =
						(0.05*ref0->mPixel[verticalView][horizontalView][viewLine][viewColumn] +
						0.50*ref1->mPixel[verticalView][horizontalView][viewLine][viewColumn] +
						0.15*ref2->mPixel[verticalView][horizontalView][viewLine][viewColumn] +
						0.30*ref3->mPixel[verticalView][horizontalView][viewLine][viewColumn]);
					//origBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] -
					//printf("res: %d\n", residueBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn]);
				}
			}
		}
	}
}           


// void Prediction :: differentialPredictionRaster(Block4D *residueBlock, Block4D *origBlock) {
// 	int residuesSum = 0;

// 	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
// 		for(int horizontalView = 0; horizontalView < 13; horizontalView += 1) {
// 			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
// 				for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
// 					if(horizontalView > 0) {
// 						for(int i = 1; i < horizontalView; i++) {
// 							residuesSum += residueBlock->mPixel[verticalView][i][viewLine][viewColumn];
// 						}
// 						residueBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] =
// 								origBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] - 
// 								(origBlock->mPixel[verticalView][0][viewLine][viewColumn] + residuesSum);
							
// 							//printf("r[%d, %d]: %d - (%d + %d)\n", verticalView, horizontalView, residueBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn],
// 							//	origBlock->mPixel[verticalView][0][viewLine][viewColumn], residuesSum);
// 					}
// 					residuesSum = 0;
// 					//if(horizontalView < 2)
// 					//	printf("v[%d, %d]: %d\n", verticalView, horizontalView, origBlock->mPixel[verticalView][0][viewLine][viewColumn]);
// 				}
// 			}
// 		}
// 	}
// }

void Prediction :: differentialPredictionRaster(Block4D *residueBlock, Block4D *origBlock) {
	int residuesSum = 0;

	for(int verticalView = 12; verticalView >= 0; verticalView -= 1) {
		for(int horizontalView = 12; horizontalView >= 0; horizontalView -= 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
					if(horizontalView > 0) {
						for(int i = 11; i > horizontalView; i--) {
							residuesSum += residueBlock->mPixel[verticalView][i][viewLine][viewColumn];
						}
						residueBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] =
								origBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] - 
								(origBlock->mPixel[verticalView][0][viewLine][viewColumn] + residuesSum);
							
							//printf("r[%d, %d]: %d - (%d + %d)\n", verticalView, horizontalView, residueBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn],
							//	origBlock->mPixel[verticalView][0][viewLine][viewColumn], residuesSum);
					}
					residuesSum = 0;
					//if(horizontalView < 2)
					//	printf("v[%d, %d]: %d\n", verticalView, horizontalView, origBlock->mPixel[verticalView][0][viewLine][viewColumn]);
				}
			}
		}
	}
}

void Prediction :: recDifferentialPredictionRaster(Block4D *recBlock, Block4D *origBlock) {
	int residuesSum = 0;
	printf("#########\n");
	for(int verticalView = 12; verticalView >= 0; verticalView -= 1) {
		for(int horizontalView = 12; horizontalView >= 0; horizontalView -= 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
					if(horizontalView > 0) {
						for(int i = 11; i > horizontalView; i--) {
							printf("#########2222222222\n");
							residuesSum += origBlock->mPixel[verticalView][i][viewLine][viewColumn];
						}
						recBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] =
								origBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] +
								recBlock->mPixel[verticalView][12][viewLine][viewColumn] + residuesSum;
							
							//printf("r[%d, %d]: %d - (%d + %d)\n", verticalView, horizontalView, residueBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn],
							//	origBlock->mPixel[verticalView][0][viewLine][viewColumn], residuesSum);
					}
					residuesSum = 0;
				}
			}
		}
	}
	
}