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