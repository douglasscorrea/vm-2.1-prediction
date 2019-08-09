#include "Prediction.h"
#include "Block4D.h"
#include <stdio.h>
#include <stdlib.h>

Prediction :: Prediction() {

}

void Prediction :: simplePredictor(Block4D *targetBlock, int position_t, int position_s, int position_v, int position_u, int component) {
	//printf("t, s, v, u: %d, %d, %d, %d\n", position_t, position_s, position_v, position_u);
	int pixelsAverageBlock = targetBlock->sumAllPixels(position_t, position_s, position_v, position_u);

	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
		for(int horizontalView = 0; horizontalView < 13; horizontalView += 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
					//printf("t, s, v, u: %d, %d, %d, %d\n", verticalView, horizontalView, viewLine, viewColumn);
					targetBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] = pixelsAverageBlock;
				}
			}
		}
	}
	printf("Average: %d\n", targetBlock->sumAllPixels(position_t, position_s, position_v, position_u));

}