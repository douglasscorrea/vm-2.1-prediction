#include "Prediction.h"
#include "Block4D.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

Prediction :: Prediction(int prediction) {
	predictionType = prediction;

	if(prediction == 0) {
		DC_coeff.open("DC_coeff_diffR.txt");
		AC_coeff.open("AC_coeff_diffR.txt");
		//firstPlaneCoefficients.open("firstPlaneCoeff_diffR.txt");
		//allCoefficients.open("allCoefficients_diffR.txt");
		//yResiduesFile.open("y_residues_diffR.txt");
		//cbResiduesFile.open("cb_residues_diffR.txt");
		//crResiduesFile.open("cr_residues_diffR.txt");
	}
	else if(prediction == 1) {
		DC_coeff.open("DC_coeff_diffC.txt");
		AC_coeff.open("AC_coeff_diffC.txt");
		//firstPlaneCoefficients.open("firstPlaneCoeff_diffC.txt");
		//allCoefficients.open("allCoefficients_diffC.txt");
		//yResiduesFile.open("y_residues_diffC.txt");
		//cbResiduesFile.open("cb_residues_diffC.txt");
		//crResiduesFile.open("cr_residues_diffC.txt");
	}
	else {
		DC_coeff.open("DC_coeff_MuLE.txt");
		AC_coeff.open("AC_coeff_MuLE.txt");
		//firstPlaneCoefficients.open("firstPlaneCoeff_MuLE.txt");
		//allCoefficients.open("allCoefficients_MuLE.txt");
		//yResiduesFile.open("y_samples_MuLE.txt");
		//cbResiduesFile.open("cb_samples_MuLE.txt");
		//crResiduesFile.open("cr_samples_MuLE.txt");
	}

	counter = 0;
	y_totalSignalEnergyFirstPlane = 0;
	cb_totalSignalEnergyFirstPlane = 0;
	cr_totalSignalEnergyFirstPlane = 0;
	y_totalSignalEnergyOtherPlanes = 0;
	cb_totalSignalEnergyOtherPlanes = 0;
	cr_totalSignalEnergyOtherPlanes = 0;
	maxRefPlane = -9999999;
	minRefPlane = 9999999;
	maxOtherPlanes = -9999999;
	minOtherPlanes = 9999999;
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
				}
			}
		}
	}
}           

void Prediction :: saveSamplesMule(Block4D *residueBlock, Block4D *origBlock, int spectralComponent) {
	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
		for(int horizontalView = 0; horizontalView < 13; horizontalView += 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
					residueBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] = 
							origBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn];

					if(spectralComponent == 0) {
						yResiduesFile << origBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] << '\n';
					}
					else if(spectralComponent == 1) {
						cbResiduesFile << origBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] << '\n';
					}
					else {
						crResiduesFile << origBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] << '\n';
					}
				}
			}
		}
	}
}

void Prediction :: differentialPredictionRaster(Block4D *residueBlock, Block4D *origBlock, int spectralComponent) {
	int residuesSum = 0;

	for(int horizontalView = 0; horizontalView < 13; horizontalView += 1) {
		for(int verticalView = 0; verticalView < 13; verticalView += 1) {
			for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
				for(int viewLine = 0; viewLine < 15; viewLine += 1) {
					// horizontal, vertical, coluna, linha
					// if(horizontalView == 0 && verticalView == 2)
					// 	printf("%d\n", origBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine]);
					if(horizontalView > 0) {
						residueBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] =
								origBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] - 
								(origBlock->mPixel[horizontalView-1][verticalView][viewColumn][viewLine]);
					}
					else {
						residueBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] = 
							origBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
					}

					if(spectralComponent == 0) {
						yResiduesFile << residueBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] << '\n';
					}
					else if(spectralComponent == 1) {
						cbResiduesFile << residueBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] << '\n';
					}
					else {
						crResiduesFile << residueBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] << '\n';
					}

					residuesSum = 0;
				}
			}
		}
	}
}

void Prediction :: differentialPredictionCentral(Block4D *residueBlock, Block4D *origBlock, int spectralComponent) {
	int residuesSum = 0;

	/* RIGHT PREDICTION */
	for(int horizontalView = 6; horizontalView < 13; horizontalView += 1) {
		for(int verticalView = 0; verticalView < 13; verticalView += 1) {
			for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
				for(int viewLine = 0; viewLine < 15; viewLine += 1) {
					if(horizontalView != 6) {
						// for(int i = 7; i < horizontalView; i++) {
						// 	residuesSum += residueBlock->mPixel[verticalView][i][viewLine][viewColumn];
						// }
						residueBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] =
								origBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] - 
								origBlock->mPixel[horizontalView-1][verticalView][viewColumn][viewLine];	
					}
					else
						residueBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] = origBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
				
					if(spectralComponent == 0) {
						yResiduesFile << residueBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] << '\n';
					}
					else if(spectralComponent == 1) {
						cbResiduesFile << residueBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] << '\n';
					}
					else {
						crResiduesFile << residueBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] << '\n';
					}
				}
			}
		}
	}

	/* LEFT PREDICTION */
	for(int horizontalView = 5; horizontalView >=0; horizontalView -= 1) {
		for(int verticalView = 0; verticalView < 13; verticalView += 1) {
			for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
				for(int viewLine = 0; viewLine < 15; viewLine += 1) {
					// for(int i = 5; i > horizontalView; i--) {
					// 	residuesSum += residueBlock->mPixel[verticalView][i][viewLine][viewColumn];
					// }
					residueBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] =
							origBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] - 
							origBlock->mPixel[horizontalView+1][verticalView][viewColumn][viewLine];	
					
					residuesSum = 0;

					if(spectralComponent == 0) {
						yResiduesFile << residueBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] << '\n';
					}
					else if(spectralComponent == 1) {
						cbResiduesFile << residueBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] << '\n';
					}
					else {
						crResiduesFile << residueBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] << '\n';
					}
				}
			}
		}
	}
}

void Prediction :: differentialPredictionRasterHalf(Block4D *residueBlock, Block4D *origBlock) {
	int residuesSum = 0;

	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
		for(int horizontalView = 0; horizontalView < 13; horizontalView += 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
					if((horizontalView % 2) != 0) {
						residueBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] =
								origBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] - 
								origBlock->mPixel[verticalView][horizontalView-1][viewLine][viewColumn];
					}
					else {
						residueBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] = origBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn];
					}
					residuesSum = 0;
				}
			}
		}
	}
}

void Prediction :: recDifferentialPredictionRaster(Block4D *residueBlock, Block4D *reconstructedBlock) {
	int residuesSum = 0;

	for(int horizontalView = 0; horizontalView < 13; horizontalView += 1) {
		for(int verticalView = 0; verticalView < 13; verticalView += 1) {
			for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
				for(int viewLine = 0; viewLine < 15; viewLine += 1) {
					//printf("Ve: %d\n", residueBlock->mPixel[verticalView][0][viewLine][viewColumn]);
					if(horizontalView > 0) {
						for(int i = 1; i < horizontalView; i++) {
							residuesSum += residueBlock->mPixel[i][verticalView][viewColumn][viewLine];
						}
						reconstructedBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] =
							residueBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] +
							residueBlock->mPixel[0][verticalView][viewColumn][viewLine] + residuesSum;
					}
					else {
						reconstructedBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] = residueBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
					}
					//printf("[rec,res,sRes]: %d, %d, %d\n", reconstructedBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn], residueBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn], residuesSum);
					residuesSum = 0;
				}
			}
		}
	}
}

void Prediction :: recDifferentialPredictionCentral(Block4D *recBlock, Block4D *origBlock) {
	int residuesSum = 0;

	/* RIGHT PREDICTION */
	for(int horizontalView = 6; horizontalView < 13; horizontalView += 1) {
		for(int verticalView = 0; verticalView < 13; verticalView += 1) {
			for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
				for(int viewLine = 0; viewLine < 15; viewLine += 1) {
					if(horizontalView != 6) {
						for(int i = 7; i < horizontalView; i++) {
							residuesSum += recBlock->mPixel[i][verticalView][viewColumn][viewLine];
						}
						origBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] =
							recBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] +
							recBlock->mPixel[6][verticalView][viewColumn][viewLine] + residuesSum;
						// if(horizontalView == 6)
						// 	printf("origR: %d - %d - %d\n", origBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn], recBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn], residuesSum);
						residuesSum = 0;
					}
					else {
						origBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] = recBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
					}	
				}
			}
		}
	}

	/* LEFT PREDICTION */
	for(int horizontalView = 5; horizontalView >=0; horizontalView -= 1) {
		for(int verticalView = 0; verticalView < 13; verticalView += 1) {
			for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
				for(int viewLine = 0; viewLine < 15; viewLine += 1) {
					for(int i = 5; i > horizontalView; i--) {
						residuesSum += recBlock->mPixel[i][verticalView][viewColumn][viewLine];
					}
					origBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] =
						recBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine] +
						recBlock->mPixel[6][verticalView][viewColumn][viewLine] + residuesSum;
					// if(horizontalView == 6)
					// 	printf("origL: %d - %d - %d\n", origBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn], recBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn], residuesSum);
					residuesSum = 0;
				}
			}
		}
	}
}

void Prediction :: recDifferentialPredictionRasterHalf(Block4D *recBlock, Block4D *origBlock) {
	int residuesSum = 0;

	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
		for(int horizontalView = 0; horizontalView < 13; horizontalView += 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
					//printf("Vc: %d\n", recBlock->mPixel[6][6][viewLine][viewColumn]);
					if(horizontalView > 0) {
						for(int i = 1; i < horizontalView; i++) {
							residuesSum += recBlock->mPixel[verticalView][i][viewLine][viewColumn];
						}
						origBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] =
							recBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] +
							recBlock->mPixel[verticalView][0][viewLine][viewColumn] + residuesSum;
					}
					else {
						origBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] = recBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn];
					}
					
					residuesSum = 0;
				}
			}
		}
	}
}

void Prediction :: hierarchicalDifferentialPrediction(Block4D *residueBlock, Block4D *origBlock) {
	for(int viewLine = 0; viewLine < 15; viewLine += 1) {
		for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
			residueBlock->mPixel[6][6][viewLine][viewColumn] = origBlock->mPixel[6][6][viewLine][viewColumn];
			residueBlock->mPixel[2][2][viewLine][viewColumn] = origBlock->mPixel[2][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn]);
			residueBlock->mPixel[2][6][viewLine][viewColumn] = origBlock->mPixel[2][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn]);
			residueBlock->mPixel[2][10][viewLine][viewColumn] = origBlock->mPixel[2][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn]);
			residueBlock->mPixel[6][2][viewLine][viewColumn] = origBlock->mPixel[6][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn]);
			residueBlock->mPixel[6][10][viewLine][viewColumn] = origBlock->mPixel[6][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn]);
			residueBlock->mPixel[10][2][viewLine][viewColumn] = origBlock->mPixel[10][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn]);
			residueBlock->mPixel[10][6][viewLine][viewColumn] = origBlock->mPixel[10][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn]);
			residueBlock->mPixel[10][10][viewLine][viewColumn] = origBlock->mPixel[10][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn]);
			residueBlock->mPixel[0][4][viewLine][viewColumn] = origBlock->mPixel[0][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[0][8][viewLine][viewColumn] = origBlock->mPixel[0][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[4][0][viewLine][viewColumn] = origBlock->mPixel[4][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[4][4][viewLine][viewColumn] = origBlock->mPixel[4][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.33 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.33 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.33);
			residueBlock->mPixel[4][8][viewLine][viewColumn] = origBlock->mPixel[4][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.33 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.33 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.33);
			residueBlock->mPixel[4][12][viewLine][viewColumn] = origBlock->mPixel[4][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.5 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[8][0][viewLine][viewColumn] = origBlock->mPixel[8][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[8][4][viewLine][viewColumn] = origBlock->mPixel[8][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.33 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.33 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.33);
			residueBlock->mPixel[8][8][viewLine][viewColumn] = origBlock->mPixel[8][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.33 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.33 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.33);
			residueBlock->mPixel[8][12][viewLine][viewColumn] = origBlock->mPixel[8][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[12][4][viewLine][viewColumn] = origBlock->mPixel[12][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[12][8][viewLine][viewColumn] = origBlock->mPixel[12][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[0][2][viewLine][viewColumn] = origBlock->mPixel[0][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[0][6][viewLine][viewColumn] = origBlock->mPixel[0][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.33 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.33 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.33);
			residueBlock->mPixel[0][10][viewLine][viewColumn] = origBlock->mPixel[0][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[2][0][viewLine][viewColumn] = origBlock->mPixel[2][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[2][4][viewLine][viewColumn] = origBlock->mPixel[2][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[2][8][viewLine][viewColumn] = origBlock->mPixel[2][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[2][12][viewLine][viewColumn] = origBlock->mPixel[2][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[4][2][viewLine][viewColumn] = origBlock->mPixel[4][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[4][6][viewLine][viewColumn] = origBlock->mPixel[4][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.33 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.33 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.33);
			residueBlock->mPixel[4][10][viewLine][viewColumn] = origBlock->mPixel[4][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[6][0][viewLine][viewColumn] = origBlock->mPixel[6][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.33 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.33 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.33);
			residueBlock->mPixel[6][4][viewLine][viewColumn] = origBlock->mPixel[6][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.33 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.33 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.33);
			residueBlock->mPixel[6][8][viewLine][viewColumn] = origBlock->mPixel[6][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.33 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.33 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.33);
			residueBlock->mPixel[6][12][viewLine][viewColumn] = origBlock->mPixel[6][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.33 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.33 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.33);
			residueBlock->mPixel[8][2][viewLine][viewColumn] = origBlock->mPixel[8][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[8][6][viewLine][viewColumn] = origBlock->mPixel[8][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.33 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.33 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.33);
			residueBlock->mPixel[8][10][viewLine][viewColumn] = origBlock->mPixel[8][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[10][0][viewLine][viewColumn] = origBlock->mPixel[10][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[10][4][viewLine][viewColumn] = origBlock->mPixel[10][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[10][8][viewLine][viewColumn] = origBlock->mPixel[10][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[10][12][viewLine][viewColumn] = origBlock->mPixel[10][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[12][2][viewLine][viewColumn] = origBlock->mPixel[12][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[12][6][viewLine][viewColumn] = origBlock->mPixel[12][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.33 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.33 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.33);
			residueBlock->mPixel[12][10][viewLine][viewColumn] = origBlock->mPixel[12][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[1][1][viewLine][viewColumn] = origBlock->mPixel[1][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][2][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][0][viewLine][viewColumn]*0.2 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[1][3][viewLine][viewColumn] = origBlock->mPixel[1][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[1][5][viewLine][viewColumn] = origBlock->mPixel[1][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[1][7][viewLine][viewColumn] = origBlock->mPixel[1][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[1][9][viewLine][viewColumn] = origBlock->mPixel[1][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[1][11][viewLine][viewColumn] = origBlock->mPixel[1][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][10][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][12][viewLine][viewColumn]*0.2 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[3][1][viewLine][viewColumn] = origBlock->mPixel[3][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[3][3][viewLine][viewColumn] = origBlock->mPixel[3][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[3][5][viewLine][viewColumn] = origBlock->mPixel[3][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[3][7][viewLine][viewColumn] = origBlock->mPixel[3][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[3][9][viewLine][viewColumn] = origBlock->mPixel[3][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[3][11][viewLine][viewColumn] = origBlock->mPixel[3][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[5][1][viewLine][viewColumn] = origBlock->mPixel[5][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[5][3][viewLine][viewColumn] = origBlock->mPixel[5][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[5][5][viewLine][viewColumn] = origBlock->mPixel[5][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[5][7][viewLine][viewColumn] = origBlock->mPixel[5][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[5][9][viewLine][viewColumn] = origBlock->mPixel[5][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[5][11][viewLine][viewColumn] = origBlock->mPixel[5][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[7][1][viewLine][viewColumn] = origBlock->mPixel[7][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[7][3][viewLine][viewColumn] = origBlock->mPixel[7][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[7][5][viewLine][viewColumn] = origBlock->mPixel[7][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[7][7][viewLine][viewColumn] = origBlock->mPixel[7][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[7][9][viewLine][viewColumn] = origBlock->mPixel[7][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[7][11][viewLine][viewColumn] = origBlock->mPixel[7][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[9][1][viewLine][viewColumn] = origBlock->mPixel[9][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[9][3][viewLine][viewColumn] = origBlock->mPixel[9][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[9][5][viewLine][viewColumn] = origBlock->mPixel[9][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[9][7][viewLine][viewColumn] = origBlock->mPixel[9][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[9][9][viewLine][viewColumn] = origBlock->mPixel[9][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[9][11][viewLine][viewColumn] = origBlock->mPixel[9][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[11][1][viewLine][viewColumn] = origBlock->mPixel[11][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][0][viewLine][viewColumn]*0.2 + residueBlock->mPixel[12][2][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.2 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[11][3][viewLine][viewColumn] = origBlock->mPixel[11][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[11][5][viewLine][viewColumn] = origBlock->mPixel[11][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[11][7][viewLine][viewColumn] = origBlock->mPixel[11][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[11][9][viewLine][viewColumn] = origBlock->mPixel[11][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[11][11][viewLine][viewColumn] = origBlock->mPixel[11][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][12][viewLine][viewColumn]*0.2 + residueBlock->mPixel[12][10][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.2 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[0][0][viewLine][viewColumn] = origBlock->mPixel[0][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[0][1][viewLine][viewColumn] = origBlock->mPixel[0][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][1][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[0][3][viewLine][viewColumn] = origBlock->mPixel[0][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][3][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[0][5][viewLine][viewColumn] = origBlock->mPixel[0][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][5][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[0][7][viewLine][viewColumn] = origBlock->mPixel[0][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][7][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[0][9][viewLine][viewColumn] = origBlock->mPixel[0][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][9][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[0][11][viewLine][viewColumn] = origBlock->mPixel[0][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][11][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[0][12][viewLine][viewColumn] = origBlock->mPixel[0][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[1][0][viewLine][viewColumn] = origBlock->mPixel[1][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][1][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[1][2][viewLine][viewColumn] = origBlock->mPixel[1][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][1][viewLine][viewColumn]*0.2 + residueBlock->mPixel[1][3][viewLine][viewColumn]*0.2 + residueBlock->mPixel[0][2][viewLine][viewColumn]*0.2 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[1][4][viewLine][viewColumn] = origBlock->mPixel[1][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[1][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[1][6][viewLine][viewColumn] = origBlock->mPixel[1][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[1][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[1][8][viewLine][viewColumn] = origBlock->mPixel[1][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[1][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[1][10][viewLine][viewColumn] = origBlock->mPixel[1][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][9][viewLine][viewColumn]*0.2 + residueBlock->mPixel[1][11][viewLine][viewColumn]*0.2 + residueBlock->mPixel[0][10][viewLine][viewColumn]*0.2 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[1][12][viewLine][viewColumn] = origBlock->mPixel[1][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][11][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[2][1][viewLine][viewColumn] = origBlock->mPixel[2][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][1][viewLine][viewColumn]*0.2 + residueBlock->mPixel[3][1][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][0][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[2][3][viewLine][viewColumn] = origBlock->mPixel[2][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[2][5][viewLine][viewColumn] = origBlock->mPixel[2][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[2][7][viewLine][viewColumn] = origBlock->mPixel[2][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[2][9][viewLine][viewColumn] = origBlock->mPixel[2][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[2][11][viewLine][viewColumn] = origBlock->mPixel[2][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][11][viewLine][viewColumn]*0.2 + residueBlock->mPixel[3][11][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][12][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[3][0][viewLine][viewColumn] = origBlock->mPixel[3][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][1][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[3][2][viewLine][viewColumn] = origBlock->mPixel[3][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[3][4][viewLine][viewColumn] = origBlock->mPixel[3][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[3][6][viewLine][viewColumn] = origBlock->mPixel[3][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[3][8][viewLine][viewColumn] = origBlock->mPixel[3][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[3][10][viewLine][viewColumn] = origBlock->mPixel[3][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[3][12][viewLine][viewColumn] = origBlock->mPixel[3][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][11][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[4][1][viewLine][viewColumn] = origBlock->mPixel[4][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[5][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[4][3][viewLine][viewColumn] = origBlock->mPixel[4][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[5][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[4][5][viewLine][viewColumn] = origBlock->mPixel[4][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][5][viewLine][viewColumn]*0.2 + residueBlock->mPixel[5][5][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[4][7][viewLine][viewColumn] = origBlock->mPixel[4][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][7][viewLine][viewColumn]*0.2 + residueBlock->mPixel[5][7][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[4][9][viewLine][viewColumn] = origBlock->mPixel[4][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[5][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[4][11][viewLine][viewColumn] = origBlock->mPixel[4][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[5][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[5][0][viewLine][viewColumn] = origBlock->mPixel[5][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][1][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[5][2][viewLine][viewColumn] = origBlock->mPixel[5][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[5][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[5][4][viewLine][viewColumn] = origBlock->mPixel[5][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][3][viewLine][viewColumn]*0.2 + residueBlock->mPixel[5][5][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[5][6][viewLine][viewColumn] = origBlock->mPixel[5][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[5][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[5][8][viewLine][viewColumn] = origBlock->mPixel[5][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][7][viewLine][viewColumn]*0.2 + residueBlock->mPixel[5][9][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[5][10][viewLine][viewColumn] = origBlock->mPixel[5][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[5][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[5][12][viewLine][viewColumn] = origBlock->mPixel[5][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][11][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[6][1][viewLine][viewColumn] = origBlock->mPixel[6][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[6][3][viewLine][viewColumn] = origBlock->mPixel[6][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[6][5][viewLine][viewColumn] = origBlock->mPixel[6][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[6][7][viewLine][viewColumn] = origBlock->mPixel[6][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[6][9][viewLine][viewColumn] = origBlock->mPixel[6][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[6][11][viewLine][viewColumn] = origBlock->mPixel[6][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[7][0][viewLine][viewColumn] = origBlock->mPixel[7][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][1][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[7][2][viewLine][viewColumn] = origBlock->mPixel[7][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[7][4][viewLine][viewColumn] = origBlock->mPixel[7][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][3][viewLine][viewColumn]*0.2 + residueBlock->mPixel[7][5][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[7][6][viewLine][viewColumn] = origBlock->mPixel[7][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[7][8][viewLine][viewColumn] = origBlock->mPixel[7][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][7][viewLine][viewColumn]*0.2 + residueBlock->mPixel[7][9][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[7][10][viewLine][viewColumn] = origBlock->mPixel[7][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[7][12][viewLine][viewColumn] = origBlock->mPixel[7][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][11][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[8][1][viewLine][viewColumn] = origBlock->mPixel[8][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[8][3][viewLine][viewColumn] = origBlock->mPixel[8][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[8][5][viewLine][viewColumn] = origBlock->mPixel[8][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][5][viewLine][viewColumn]*0.2 + residueBlock->mPixel[9][5][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[8][7][viewLine][viewColumn] = origBlock->mPixel[8][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][7][viewLine][viewColumn]*0.2 + residueBlock->mPixel[9][7][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[8][9][viewLine][viewColumn] = origBlock->mPixel[8][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[8][11][viewLine][viewColumn] = origBlock->mPixel[8][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[9][0][viewLine][viewColumn] = origBlock->mPixel[9][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][1][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[9][2][viewLine][viewColumn] = origBlock->mPixel[9][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[9][4][viewLine][viewColumn] = origBlock->mPixel[9][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[9][6][viewLine][viewColumn] = origBlock->mPixel[9][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[9][8][viewLine][viewColumn] = origBlock->mPixel[9][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[9][10][viewLine][viewColumn] = origBlock->mPixel[9][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[9][12][viewLine][viewColumn] = origBlock->mPixel[9][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][11][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[10][1][viewLine][viewColumn] = origBlock->mPixel[10][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][1][viewLine][viewColumn]*0.2 + residueBlock->mPixel[11][1][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][0][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[10][3][viewLine][viewColumn] = origBlock->mPixel[10][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[11][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[10][5][viewLine][viewColumn] = origBlock->mPixel[10][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[11][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[10][7][viewLine][viewColumn] = origBlock->mPixel[10][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[11][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[10][9][viewLine][viewColumn] = origBlock->mPixel[10][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[11][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[10][11][viewLine][viewColumn] = origBlock->mPixel[10][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][11][viewLine][viewColumn]*0.2 + residueBlock->mPixel[11][11][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][12][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[11][0][viewLine][viewColumn] = origBlock->mPixel[11][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][1][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[11][2][viewLine][viewColumn] = origBlock->mPixel[11][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][1][viewLine][viewColumn]*0.2 + residueBlock->mPixel[11][3][viewLine][viewColumn]*0.2 + residueBlock->mPixel[12][2][viewLine][viewColumn]*0.2 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[11][4][viewLine][viewColumn] = origBlock->mPixel[11][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[11][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[11][6][viewLine][viewColumn] = origBlock->mPixel[11][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[11][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[11][8][viewLine][viewColumn] = origBlock->mPixel[11][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[11][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[11][10][viewLine][viewColumn] = origBlock->mPixel[11][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][9][viewLine][viewColumn]*0.2 + residueBlock->mPixel[11][11][viewLine][viewColumn]*0.2 + residueBlock->mPixel[12][10][viewLine][viewColumn]*0.2 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.2);
			residueBlock->mPixel[11][12][viewLine][viewColumn] = origBlock->mPixel[11][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][11][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[12][0][viewLine][viewColumn] = origBlock->mPixel[12][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.17);
			residueBlock->mPixel[12][1][viewLine][viewColumn] = origBlock->mPixel[12][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][1][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[12][3][viewLine][viewColumn] = origBlock->mPixel[12][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][3][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[12][5][viewLine][viewColumn] = origBlock->mPixel[12][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][5][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[12][7][viewLine][viewColumn] = origBlock->mPixel[12][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][7][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[12][9][viewLine][viewColumn] = origBlock->mPixel[12][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][9][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[12][11][viewLine][viewColumn] = origBlock->mPixel[12][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][11][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25);
			residueBlock->mPixel[12][12][viewLine][viewColumn] = origBlock->mPixel[12][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.17);
		}
	}
}

void Prediction :: recHierarchicalDifferentialPrediction(Block4D *residueBlock, Block4D *reconstructedBlock) {
	for(int viewLine = 0; viewLine < 15; viewLine += 1) {
		for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
			reconstructedBlock->mPixel[6][6][viewLine][viewColumn] = residueBlock->mPixel[6][6][viewLine][viewColumn];
			reconstructedBlock->mPixel[2][2][viewLine][viewColumn] = residueBlock->mPixel[2][2][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn];
			reconstructedBlock->mPixel[2][6][viewLine][viewColumn] = residueBlock->mPixel[2][6][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn];
			reconstructedBlock->mPixel[2][10][viewLine][viewColumn] = residueBlock->mPixel[2][10][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn];
			reconstructedBlock->mPixel[6][2][viewLine][viewColumn] = residueBlock->mPixel[6][2][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn];
			reconstructedBlock->mPixel[6][10][viewLine][viewColumn] = residueBlock->mPixel[6][10][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn];
			reconstructedBlock->mPixel[10][2][viewLine][viewColumn] = residueBlock->mPixel[10][2][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn];
			reconstructedBlock->mPixel[10][6][viewLine][viewColumn] = residueBlock->mPixel[10][6][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn];
			reconstructedBlock->mPixel[10][10][viewLine][viewColumn] = residueBlock->mPixel[10][10][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn];
			reconstructedBlock->mPixel[0][4][viewLine][viewColumn] = residueBlock->mPixel[0][4][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.5;
			reconstructedBlock->mPixel[0][8][viewLine][viewColumn] = residueBlock->mPixel[0][8][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.5;
			reconstructedBlock->mPixel[4][0][viewLine][viewColumn] = residueBlock->mPixel[4][0][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.5;
			reconstructedBlock->mPixel[4][4][viewLine][viewColumn] = residueBlock->mPixel[4][4][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.33 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.33 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.33;
			reconstructedBlock->mPixel[4][8][viewLine][viewColumn] = residueBlock->mPixel[4][8][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.33 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.33 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.33;
			reconstructedBlock->mPixel[4][12][viewLine][viewColumn] = residueBlock->mPixel[4][12][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.5 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.5;
			reconstructedBlock->mPixel[8][0][viewLine][viewColumn] = residueBlock->mPixel[8][0][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.5;
			reconstructedBlock->mPixel[8][4][viewLine][viewColumn] = residueBlock->mPixel[8][4][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.33 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.33 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.33;
			reconstructedBlock->mPixel[8][8][viewLine][viewColumn] = residueBlock->mPixel[8][8][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.33 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.33 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.33;
			reconstructedBlock->mPixel[8][12][viewLine][viewColumn] = residueBlock->mPixel[8][12][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.5;
			reconstructedBlock->mPixel[12][4][viewLine][viewColumn] = residueBlock->mPixel[12][4][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.5;
			reconstructedBlock->mPixel[12][8][viewLine][viewColumn] = residueBlock->mPixel[12][8][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.5;
			reconstructedBlock->mPixel[0][2][viewLine][viewColumn] = residueBlock->mPixel[0][2][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.5;
			reconstructedBlock->mPixel[0][6][viewLine][viewColumn] = residueBlock->mPixel[0][6][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.33 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.33 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.33;
			reconstructedBlock->mPixel[0][10][viewLine][viewColumn] = residueBlock->mPixel[0][10][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.5;
			reconstructedBlock->mPixel[2][0][viewLine][viewColumn] = residueBlock->mPixel[2][0][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.5;
			reconstructedBlock->mPixel[2][4][viewLine][viewColumn] = residueBlock->mPixel[2][4][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[2][8][viewLine][viewColumn] = residueBlock->mPixel[2][8][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[2][12][viewLine][viewColumn] = residueBlock->mPixel[2][12][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.5;
			reconstructedBlock->mPixel[4][2][viewLine][viewColumn] = residueBlock->mPixel[4][2][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[4][6][viewLine][viewColumn] = residueBlock->mPixel[4][6][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.33 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.33 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.33;
			reconstructedBlock->mPixel[4][10][viewLine][viewColumn] = residueBlock->mPixel[4][10][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[6][0][viewLine][viewColumn] = residueBlock->mPixel[6][0][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.33 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.33 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.33;
			reconstructedBlock->mPixel[6][4][viewLine][viewColumn] = residueBlock->mPixel[6][4][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.33 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.33 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.33;
			reconstructedBlock->mPixel[6][8][viewLine][viewColumn] = residueBlock->mPixel[6][8][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.33 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.33 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.33;
			reconstructedBlock->mPixel[6][12][viewLine][viewColumn] = residueBlock->mPixel[6][12][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.33 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.33 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.33;
			reconstructedBlock->mPixel[8][2][viewLine][viewColumn] = residueBlock->mPixel[8][2][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[8][6][viewLine][viewColumn] = residueBlock->mPixel[8][6][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.33 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.33 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.33;
			reconstructedBlock->mPixel[8][10][viewLine][viewColumn] = residueBlock->mPixel[8][10][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[10][0][viewLine][viewColumn] = residueBlock->mPixel[10][0][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.5;
			reconstructedBlock->mPixel[10][4][viewLine][viewColumn] = residueBlock->mPixel[10][4][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[10][8][viewLine][viewColumn] = residueBlock->mPixel[10][8][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[10][12][viewLine][viewColumn] = residueBlock->mPixel[10][12][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.5;
			reconstructedBlock->mPixel[12][2][viewLine][viewColumn] = residueBlock->mPixel[12][2][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.5;
			reconstructedBlock->mPixel[12][6][viewLine][viewColumn] = residueBlock->mPixel[12][6][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.33 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.33 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.33;
			reconstructedBlock->mPixel[12][10][viewLine][viewColumn] = residueBlock->mPixel[12][10][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.5;
			reconstructedBlock->mPixel[1][1][viewLine][viewColumn] = residueBlock->mPixel[1][1][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][2][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][0][viewLine][viewColumn]*0.2 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[1][3][viewLine][viewColumn] = residueBlock->mPixel[1][3][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[1][5][viewLine][viewColumn] = residueBlock->mPixel[1][5][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[1][7][viewLine][viewColumn] = residueBlock->mPixel[1][7][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[1][9][viewLine][viewColumn] = residueBlock->mPixel[1][9][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[1][11][viewLine][viewColumn] = residueBlock->mPixel[1][11][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][10][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][12][viewLine][viewColumn]*0.2 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[3][1][viewLine][viewColumn] = residueBlock->mPixel[3][1][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[3][3][viewLine][viewColumn] = residueBlock->mPixel[3][3][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[3][5][viewLine][viewColumn] = residueBlock->mPixel[3][5][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[3][7][viewLine][viewColumn] = residueBlock->mPixel[3][7][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[3][9][viewLine][viewColumn] = residueBlock->mPixel[3][9][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[3][11][viewLine][viewColumn] = residueBlock->mPixel[3][11][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[5][1][viewLine][viewColumn] = residueBlock->mPixel[5][1][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[5][3][viewLine][viewColumn] = residueBlock->mPixel[5][3][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[5][5][viewLine][viewColumn] = residueBlock->mPixel[5][5][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[5][7][viewLine][viewColumn] = residueBlock->mPixel[5][7][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[5][9][viewLine][viewColumn] = residueBlock->mPixel[5][9][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[5][11][viewLine][viewColumn] = residueBlock->mPixel[5][11][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[7][1][viewLine][viewColumn] = residueBlock->mPixel[7][1][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[7][3][viewLine][viewColumn] = residueBlock->mPixel[7][3][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[7][5][viewLine][viewColumn] = residueBlock->mPixel[7][5][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[7][7][viewLine][viewColumn] = residueBlock->mPixel[7][7][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[7][9][viewLine][viewColumn] = residueBlock->mPixel[7][9][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[7][11][viewLine][viewColumn] = residueBlock->mPixel[7][11][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[9][1][viewLine][viewColumn] = residueBlock->mPixel[9][1][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[9][3][viewLine][viewColumn] = residueBlock->mPixel[9][3][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[9][5][viewLine][viewColumn] = residueBlock->mPixel[9][5][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[9][7][viewLine][viewColumn] = residueBlock->mPixel[9][7][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[9][9][viewLine][viewColumn] = residueBlock->mPixel[9][9][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[9][11][viewLine][viewColumn] = residueBlock->mPixel[9][11][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[11][1][viewLine][viewColumn] = residueBlock->mPixel[11][1][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][0][viewLine][viewColumn]*0.2 + residueBlock->mPixel[12][2][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.2 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[11][3][viewLine][viewColumn] = residueBlock->mPixel[11][3][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[11][5][viewLine][viewColumn] = residueBlock->mPixel[11][5][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[11][7][viewLine][viewColumn] = residueBlock->mPixel[11][7][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[11][9][viewLine][viewColumn] = residueBlock->mPixel[11][9][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[11][11][viewLine][viewColumn] = residueBlock->mPixel[11][11][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][12][viewLine][viewColumn]*0.2 + residueBlock->mPixel[12][10][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.2 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[0][0][viewLine][viewColumn] = residueBlock->mPixel[0][0][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[0][1][viewLine][viewColumn] = residueBlock->mPixel[0][1][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][1][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[0][3][viewLine][viewColumn] = residueBlock->mPixel[0][3][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][3][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[0][5][viewLine][viewColumn] = residueBlock->mPixel[0][5][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][5][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[0][7][viewLine][viewColumn] = residueBlock->mPixel[0][7][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][7][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[0][9][viewLine][viewColumn] = residueBlock->mPixel[0][9][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][9][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[0][11][viewLine][viewColumn] = residueBlock->mPixel[0][11][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][11][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[0][12][viewLine][viewColumn] = residueBlock->mPixel[0][12][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[1][0][viewLine][viewColumn] = residueBlock->mPixel[1][0][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][1][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[1][2][viewLine][viewColumn] = residueBlock->mPixel[1][2][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][1][viewLine][viewColumn]*0.2 + residueBlock->mPixel[1][3][viewLine][viewColumn]*0.2 + residueBlock->mPixel[0][2][viewLine][viewColumn]*0.2 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[1][4][viewLine][viewColumn] = residueBlock->mPixel[1][4][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[1][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[1][6][viewLine][viewColumn] = residueBlock->mPixel[1][6][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[1][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[1][8][viewLine][viewColumn] = residueBlock->mPixel[1][8][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[1][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[1][10][viewLine][viewColumn] = residueBlock->mPixel[1][10][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][9][viewLine][viewColumn]*0.2 + residueBlock->mPixel[1][11][viewLine][viewColumn]*0.2 + residueBlock->mPixel[0][10][viewLine][viewColumn]*0.2 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[1][12][viewLine][viewColumn] = residueBlock->mPixel[1][12][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][11][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[2][1][viewLine][viewColumn] = residueBlock->mPixel[2][1][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][1][viewLine][viewColumn]*0.2 + residueBlock->mPixel[3][1][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][0][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[2][3][viewLine][viewColumn] = residueBlock->mPixel[2][3][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[2][5][viewLine][viewColumn] = residueBlock->mPixel[2][5][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[2][7][viewLine][viewColumn] = residueBlock->mPixel[2][7][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[2][9][viewLine][viewColumn] = residueBlock->mPixel[2][9][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[2][11][viewLine][viewColumn] = residueBlock->mPixel[2][11][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][11][viewLine][viewColumn]*0.2 + residueBlock->mPixel[3][11][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][12][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[3][0][viewLine][viewColumn] = residueBlock->mPixel[3][0][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][1][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[3][2][viewLine][viewColumn] = residueBlock->mPixel[3][2][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[3][4][viewLine][viewColumn] = residueBlock->mPixel[3][4][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[3][6][viewLine][viewColumn] = residueBlock->mPixel[3][6][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[3][8][viewLine][viewColumn] = residueBlock->mPixel[3][8][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[3][10][viewLine][viewColumn] = residueBlock->mPixel[3][10][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[3][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[3][12][viewLine][viewColumn] = residueBlock->mPixel[3][12][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][11][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[4][1][viewLine][viewColumn] = residueBlock->mPixel[4][1][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[5][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[4][3][viewLine][viewColumn] = residueBlock->mPixel[4][3][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[5][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[4][5][viewLine][viewColumn] = residueBlock->mPixel[4][5][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][5][viewLine][viewColumn]*0.2 + residueBlock->mPixel[5][5][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[4][7][viewLine][viewColumn] = residueBlock->mPixel[4][7][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][7][viewLine][viewColumn]*0.2 + residueBlock->mPixel[5][7][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[4][9][viewLine][viewColumn] = residueBlock->mPixel[4][9][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[5][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[4][11][viewLine][viewColumn] = residueBlock->mPixel[4][11][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[5][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[5][0][viewLine][viewColumn] = residueBlock->mPixel[5][0][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][1][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[5][2][viewLine][viewColumn] = residueBlock->mPixel[5][2][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[5][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[5][4][viewLine][viewColumn] = residueBlock->mPixel[5][4][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][3][viewLine][viewColumn]*0.2 + residueBlock->mPixel[5][5][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[5][6][viewLine][viewColumn] = residueBlock->mPixel[5][6][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[5][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[5][8][viewLine][viewColumn] = residueBlock->mPixel[5][8][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][7][viewLine][viewColumn]*0.2 + residueBlock->mPixel[5][9][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[5][10][viewLine][viewColumn] = residueBlock->mPixel[5][10][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[5][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[5][12][viewLine][viewColumn] = residueBlock->mPixel[5][12][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][11][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[6][1][viewLine][viewColumn] = residueBlock->mPixel[6][1][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[6][3][viewLine][viewColumn] = residueBlock->mPixel[6][3][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[6][5][viewLine][viewColumn] = residueBlock->mPixel[6][5][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[6][7][viewLine][viewColumn] = residueBlock->mPixel[6][7][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[6][9][viewLine][viewColumn] = residueBlock->mPixel[6][9][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[6][11][viewLine][viewColumn] = residueBlock->mPixel[6][11][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[7][0][viewLine][viewColumn] = residueBlock->mPixel[7][0][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][1][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[7][2][viewLine][viewColumn] = residueBlock->mPixel[7][2][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[7][4][viewLine][viewColumn] = residueBlock->mPixel[7][4][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][3][viewLine][viewColumn]*0.2 + residueBlock->mPixel[7][5][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[7][6][viewLine][viewColumn] = residueBlock->mPixel[7][6][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[7][8][viewLine][viewColumn] = residueBlock->mPixel[7][8][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][7][viewLine][viewColumn]*0.2 + residueBlock->mPixel[7][9][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[7][10][viewLine][viewColumn] = residueBlock->mPixel[7][10][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[7][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[7][12][viewLine][viewColumn] = residueBlock->mPixel[7][12][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][11][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[8][1][viewLine][viewColumn] = residueBlock->mPixel[8][1][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[8][3][viewLine][viewColumn] = residueBlock->mPixel[8][3][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[8][5][viewLine][viewColumn] = residueBlock->mPixel[8][5][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][5][viewLine][viewColumn]*0.2 + residueBlock->mPixel[9][5][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[8][7][viewLine][viewColumn] = residueBlock->mPixel[8][7][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][7][viewLine][viewColumn]*0.2 + residueBlock->mPixel[9][7][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[8][9][viewLine][viewColumn] = residueBlock->mPixel[8][9][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[8][11][viewLine][viewColumn] = residueBlock->mPixel[8][11][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[9][0][viewLine][viewColumn] = residueBlock->mPixel[9][0][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][1][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[9][2][viewLine][viewColumn] = residueBlock->mPixel[9][2][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[9][4][viewLine][viewColumn] = residueBlock->mPixel[9][4][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[9][6][viewLine][viewColumn] = residueBlock->mPixel[9][6][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[9][8][viewLine][viewColumn] = residueBlock->mPixel[9][8][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[9][10][viewLine][viewColumn] = residueBlock->mPixel[9][10][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[9][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[9][12][viewLine][viewColumn] = residueBlock->mPixel[9][12][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][11][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[10][1][viewLine][viewColumn] = residueBlock->mPixel[10][1][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][1][viewLine][viewColumn]*0.2 + residueBlock->mPixel[11][1][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][0][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[10][3][viewLine][viewColumn] = residueBlock->mPixel[10][3][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[11][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[10][5][viewLine][viewColumn] = residueBlock->mPixel[10][5][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[11][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[10][7][viewLine][viewColumn] = residueBlock->mPixel[10][7][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[11][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[10][9][viewLine][viewColumn] = residueBlock->mPixel[10][9][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[11][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[10][11][viewLine][viewColumn] = residueBlock->mPixel[10][11][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][11][viewLine][viewColumn]*0.2 + residueBlock->mPixel[11][11][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][12][viewLine][viewColumn]*0.2 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[11][0][viewLine][viewColumn] = residueBlock->mPixel[11][0][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][1][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[11][2][viewLine][viewColumn] = residueBlock->mPixel[11][2][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][1][viewLine][viewColumn]*0.2 + residueBlock->mPixel[11][3][viewLine][viewColumn]*0.2 + residueBlock->mPixel[12][2][viewLine][viewColumn]*0.2 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[11][4][viewLine][viewColumn] = residueBlock->mPixel[11][4][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][3][viewLine][viewColumn]*0.17 + residueBlock->mPixel[11][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[11][6][viewLine][viewColumn] = residueBlock->mPixel[11][6][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][5][viewLine][viewColumn]*0.17 + residueBlock->mPixel[11][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[11][8][viewLine][viewColumn] = residueBlock->mPixel[11][8][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][7][viewLine][viewColumn]*0.17 + residueBlock->mPixel[11][9][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[11][10][viewLine][viewColumn] = residueBlock->mPixel[11][10][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][9][viewLine][viewColumn]*0.2 + residueBlock->mPixel[11][11][viewLine][viewColumn]*0.2 + residueBlock->mPixel[12][10][viewLine][viewColumn]*0.2 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.2 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.2;
			reconstructedBlock->mPixel[11][12][viewLine][viewColumn] = residueBlock->mPixel[11][12][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][11][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[12][0][viewLine][viewColumn] = residueBlock->mPixel[12][0][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][1][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][2][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.17;
			reconstructedBlock->mPixel[12][1][viewLine][viewColumn] = residueBlock->mPixel[12][1][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][1][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[12][3][viewLine][viewColumn] = residueBlock->mPixel[12][3][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][3][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][2][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[12][5][viewLine][viewColumn] = residueBlock->mPixel[12][5][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][5][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[12][7][viewLine][viewColumn] = residueBlock->mPixel[12][7][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][7][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][6][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[12][9][viewLine][viewColumn] = residueBlock->mPixel[12][9][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][9][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[12][11][viewLine][viewColumn] = residueBlock->mPixel[12][11][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][11][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][10][viewLine][viewColumn]*0.25 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.25 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.25;
			reconstructedBlock->mPixel[12][12][viewLine][viewColumn] = residueBlock->mPixel[12][12][viewLine][viewColumn] + residueBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][11][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][10][viewLine][viewColumn]*0.17 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.17 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.17 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.17;
		}
	}
}

void Prediction :: hierarchicalDifferentialPrediction1Level(Block4D *residueBlock, Block4D *origBlock) {
	for(int viewLine = 0; viewLine < 15; viewLine += 1) {
		for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
			residueBlock->mPixel[6][6][viewLine][viewColumn] = origBlock->mPixel[6][6][viewLine][viewColumn];
			residueBlock->mPixel[2][2][viewLine][viewColumn] = origBlock->mPixel[2][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn]);
			residueBlock->mPixel[2][6][viewLine][viewColumn] = origBlock->mPixel[2][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn]);
			residueBlock->mPixel[2][10][viewLine][viewColumn] = origBlock->mPixel[2][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn]);
			residueBlock->mPixel[6][2][viewLine][viewColumn] = origBlock->mPixel[6][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn]);
			residueBlock->mPixel[6][10][viewLine][viewColumn] = origBlock->mPixel[6][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn]);
			residueBlock->mPixel[10][2][viewLine][viewColumn] = origBlock->mPixel[10][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn]);
			residueBlock->mPixel[10][6][viewLine][viewColumn] = origBlock->mPixel[10][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn]);
			residueBlock->mPixel[10][10][viewLine][viewColumn] = origBlock->mPixel[10][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn]);
			residueBlock->mPixel[0][4][viewLine][viewColumn] = origBlock->mPixel[0][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[0][8][viewLine][viewColumn] = origBlock->mPixel[0][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[4][0][viewLine][viewColumn] = origBlock->mPixel[4][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[4][4][viewLine][viewColumn] = origBlock->mPixel[4][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][2][viewLine][viewColumn]*0.33 + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.33 + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.33);
			residueBlock->mPixel[4][8][viewLine][viewColumn] = origBlock->mPixel[4][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][6][viewLine][viewColumn]*0.33 + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.33 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.33);
			residueBlock->mPixel[4][12][viewLine][viewColumn] = origBlock->mPixel[4][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][10][viewLine][viewColumn]*0.5 + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[8][0][viewLine][viewColumn] = origBlock->mPixel[8][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[8][4][viewLine][viewColumn] = origBlock->mPixel[8][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][2][viewLine][viewColumn]*0.33 + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.33 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.33);
			residueBlock->mPixel[8][8][viewLine][viewColumn] = origBlock->mPixel[8][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.33 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.33 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.33);
			residueBlock->mPixel[8][12][viewLine][viewColumn] = origBlock->mPixel[8][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][10][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[12][4][viewLine][viewColumn] = origBlock->mPixel[12][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[12][8][viewLine][viewColumn] = origBlock->mPixel[12][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][6][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[0][2][viewLine][viewColumn] = origBlock->mPixel[0][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][4][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[0][6][viewLine][viewColumn] = origBlock->mPixel[0][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[0][10][viewLine][viewColumn] = origBlock->mPixel[0][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][8][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[2][0][viewLine][viewColumn] = origBlock->mPixel[2][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][0][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[2][4][viewLine][viewColumn] = origBlock->mPixel[2][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[2][8][viewLine][viewColumn] = origBlock->mPixel[2][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][8][viewLine][viewColumn]*0.5 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[2][12][viewLine][viewColumn] = origBlock->mPixel[2][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][12][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[4][2][viewLine][viewColumn] = origBlock->mPixel[4][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.5 + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[4][6][viewLine][viewColumn] = origBlock->mPixel[4][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[4][10][viewLine][viewColumn] = origBlock->mPixel[4][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.5 + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[6][0][viewLine][viewColumn] = origBlock->mPixel[6][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][0][viewLine][viewColumn]*0.5 + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[6][4][viewLine][viewColumn] = origBlock->mPixel[6][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[6][8][viewLine][viewColumn] = origBlock->mPixel[6][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][8][viewLine][viewColumn]*0.5 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[6][12][viewLine][viewColumn] = origBlock->mPixel[6][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][12][viewLine][viewColumn]*0.5 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[8][2][viewLine][viewColumn] = origBlock->mPixel[8][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][0][viewLine][viewColumn]*0.5 + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[8][6][viewLine][viewColumn] = origBlock->mPixel[8][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[8][10][viewLine][viewColumn] = origBlock->mPixel[8][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.5 + residueBlock->mPixel[8][12][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[10][0][viewLine][viewColumn] = origBlock->mPixel[10][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][0][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[10][4][viewLine][viewColumn] = origBlock->mPixel[10][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[10][8][viewLine][viewColumn] = origBlock->mPixel[10][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][8][viewLine][viewColumn]*0.5 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[10][12][viewLine][viewColumn] = origBlock->mPixel[10][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][12][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[12][2][viewLine][viewColumn] = origBlock->mPixel[12][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[12][4][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[12][6][viewLine][viewColumn] = origBlock->mPixel[12][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[12][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[12][8][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[12][10][viewLine][viewColumn] = origBlock->mPixel[12][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[12][8][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[1][1][viewLine][viewColumn] = origBlock->mPixel[1][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][0][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[1][3][viewLine][viewColumn] = origBlock->mPixel[1][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[1][5][viewLine][viewColumn] = origBlock->mPixel[1][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][6][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[1][7][viewLine][viewColumn] = origBlock->mPixel[1][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][6][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[1][9][viewLine][viewColumn] = origBlock->mPixel[1][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][10][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[1][11][viewLine][viewColumn] = origBlock->mPixel[1][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[0][10][viewLine][viewColumn]*0.5 + residueBlock->mPixel[2][12][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[3][1][viewLine][viewColumn] = origBlock->mPixel[3][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][0][viewLine][viewColumn]*0.5 + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[3][3][viewLine][viewColumn] = origBlock->mPixel[3][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[3][5][viewLine][viewColumn] = origBlock->mPixel[3][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[3][7][viewLine][viewColumn] = origBlock->mPixel[3][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.5 + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[3][9][viewLine][viewColumn] = origBlock->mPixel[3][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][8][viewLine][viewColumn]*0.5 + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[3][11][viewLine][viewColumn] = origBlock->mPixel[3][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[2][12][viewLine][viewColumn]*0.5 + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[5][1][viewLine][viewColumn] = origBlock->mPixel[5][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[6][0][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[5][3][viewLine][viewColumn] = origBlock->mPixel[5][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[5][5][viewLine][viewColumn] = origBlock->mPixel[5][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.5 + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[5][7][viewLine][viewColumn] = origBlock->mPixel[5][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][6][viewLine][viewColumn]*0.5 + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[5][9][viewLine][viewColumn] = origBlock->mPixel[5][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.5 + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[5][11][viewLine][viewColumn] = origBlock->mPixel[5][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[4][10][viewLine][viewColumn]*0.5 + residueBlock->mPixel[6][12][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[7][1][viewLine][viewColumn] = origBlock->mPixel[7][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][0][viewLine][viewColumn]*0.5 + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[7][3][viewLine][viewColumn] = origBlock->mPixel[7][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[7][5][viewLine][viewColumn] = origBlock->mPixel[7][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[7][7][viewLine][viewColumn] = origBlock->mPixel[7][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.5 + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[7][9][viewLine][viewColumn] = origBlock->mPixel[7][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][8][viewLine][viewColumn]*0.5 + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[7][11][viewLine][viewColumn] = origBlock->mPixel[7][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[6][12][viewLine][viewColumn]*0.5 + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[9][1][viewLine][viewColumn] = origBlock->mPixel[9][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][0][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[9][3][viewLine][viewColumn] = origBlock->mPixel[9][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][2][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[9][5][viewLine][viewColumn] = origBlock->mPixel[9][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[9][7][viewLine][viewColumn] = origBlock->mPixel[9][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][6][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[9][9][viewLine][viewColumn] = origBlock->mPixel[9][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[9][11][viewLine][viewColumn] = origBlock->mPixel[9][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[8][10][viewLine][viewColumn]*0.5 + residueBlock->mPixel[10][12][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[11][1][viewLine][viewColumn] = origBlock->mPixel[11][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][0][viewLine][viewColumn]*0.5 + residueBlock->mPixel[12][2][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[11][3][viewLine][viewColumn] = origBlock->mPixel[11][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[12][2][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[11][5][viewLine][viewColumn] = origBlock->mPixel[11][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][4][viewLine][viewColumn]*0.5 + residueBlock->mPixel[12][6][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[11][7][viewLine][viewColumn] = origBlock->mPixel[11][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.5 + residueBlock->mPixel[12][6][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[11][9][viewLine][viewColumn] = origBlock->mPixel[11][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][8][viewLine][viewColumn]*0.5 + residueBlock->mPixel[12][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[11][11][viewLine][viewColumn] = origBlock->mPixel[11][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[10][12][viewLine][viewColumn]*0.5 + residueBlock->mPixel[12][10][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[0][0][viewLine][viewColumn] = origBlock->mPixel[0][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][1][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[0][1][viewLine][viewColumn] = origBlock->mPixel[0][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][1][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[0][3][viewLine][viewColumn] = origBlock->mPixel[0][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][3][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[0][5][viewLine][viewColumn] = origBlock->mPixel[0][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][5][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[0][7][viewLine][viewColumn] = origBlock->mPixel[0][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][7][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[0][9][viewLine][viewColumn] = origBlock->mPixel[0][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][9][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[0][11][viewLine][viewColumn] = origBlock->mPixel[0][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][11][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[0][12][viewLine][viewColumn] = origBlock->mPixel[0][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][11][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[1][0][viewLine][viewColumn] = origBlock->mPixel[1][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][1][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[1][2][viewLine][viewColumn] = origBlock->mPixel[1][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][1][viewLine][viewColumn]*0.5 + residueBlock->mPixel[1][3][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[1][4][viewLine][viewColumn] = origBlock->mPixel[1][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][3][viewLine][viewColumn]*0.5 + residueBlock->mPixel[1][5][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[1][6][viewLine][viewColumn] = origBlock->mPixel[1][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][5][viewLine][viewColumn]*0.5 + residueBlock->mPixel[1][7][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[1][8][viewLine][viewColumn] = origBlock->mPixel[1][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][7][viewLine][viewColumn]*0.5 + residueBlock->mPixel[1][9][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[1][10][viewLine][viewColumn] = origBlock->mPixel[1][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][9][viewLine][viewColumn]*0.5 + residueBlock->mPixel[1][11][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[1][12][viewLine][viewColumn] = origBlock->mPixel[1][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][11][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[2][1][viewLine][viewColumn] = origBlock->mPixel[2][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][1][viewLine][viewColumn]*0.5 + residueBlock->mPixel[3][1][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[2][3][viewLine][viewColumn] = origBlock->mPixel[2][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][3][viewLine][viewColumn]*0.5 + residueBlock->mPixel[3][3][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[2][5][viewLine][viewColumn] = origBlock->mPixel[2][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][5][viewLine][viewColumn]*0.5 + residueBlock->mPixel[3][5][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[2][7][viewLine][viewColumn] = origBlock->mPixel[2][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][7][viewLine][viewColumn]*0.5 + residueBlock->mPixel[3][7][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[2][9][viewLine][viewColumn] = origBlock->mPixel[2][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][9][viewLine][viewColumn]*0.5 + residueBlock->mPixel[3][9][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[2][11][viewLine][viewColumn] = origBlock->mPixel[2][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[1][11][viewLine][viewColumn]*0.5 + residueBlock->mPixel[3][11][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[3][0][viewLine][viewColumn] = origBlock->mPixel[3][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][1][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[3][2][viewLine][viewColumn] = origBlock->mPixel[3][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][1][viewLine][viewColumn]*0.5 + residueBlock->mPixel[3][3][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[3][4][viewLine][viewColumn] = origBlock->mPixel[3][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][3][viewLine][viewColumn]*0.5 + residueBlock->mPixel[3][5][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[3][6][viewLine][viewColumn] = origBlock->mPixel[3][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][5][viewLine][viewColumn]*0.5 + residueBlock->mPixel[3][7][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[3][8][viewLine][viewColumn] = origBlock->mPixel[3][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][7][viewLine][viewColumn]*0.5 + residueBlock->mPixel[3][9][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[3][10][viewLine][viewColumn] = origBlock->mPixel[3][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][9][viewLine][viewColumn]*0.5 + residueBlock->mPixel[3][11][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[3][12][viewLine][viewColumn] = origBlock->mPixel[3][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][11][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[4][1][viewLine][viewColumn] = origBlock->mPixel[4][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][1][viewLine][viewColumn]*0.5 + residueBlock->mPixel[5][1][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[4][3][viewLine][viewColumn] = origBlock->mPixel[4][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][3][viewLine][viewColumn]*0.5 + residueBlock->mPixel[5][3][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[4][5][viewLine][viewColumn] = origBlock->mPixel[4][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][5][viewLine][viewColumn]*0.5 + residueBlock->mPixel[5][5][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[4][7][viewLine][viewColumn] = origBlock->mPixel[4][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][7][viewLine][viewColumn]*0.5 + residueBlock->mPixel[5][7][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[4][9][viewLine][viewColumn] = origBlock->mPixel[4][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][9][viewLine][viewColumn]*0.5 + residueBlock->mPixel[5][9][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[4][11][viewLine][viewColumn] = origBlock->mPixel[4][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[3][11][viewLine][viewColumn]*0.5 + residueBlock->mPixel[5][11][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[5][0][viewLine][viewColumn] = origBlock->mPixel[5][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][1][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[5][2][viewLine][viewColumn] = origBlock->mPixel[5][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][1][viewLine][viewColumn]*0.5 + residueBlock->mPixel[5][3][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[5][4][viewLine][viewColumn] = origBlock->mPixel[5][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][3][viewLine][viewColumn]*0.5 + residueBlock->mPixel[5][5][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[5][6][viewLine][viewColumn] = origBlock->mPixel[5][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][5][viewLine][viewColumn]*0.5 + residueBlock->mPixel[5][7][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[5][8][viewLine][viewColumn] = origBlock->mPixel[5][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][7][viewLine][viewColumn]*0.5 + residueBlock->mPixel[5][9][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[5][10][viewLine][viewColumn] = origBlock->mPixel[5][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][9][viewLine][viewColumn]*0.5 + residueBlock->mPixel[5][11][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[5][12][viewLine][viewColumn] = origBlock->mPixel[5][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][11][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[6][1][viewLine][viewColumn] = origBlock->mPixel[6][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][1][viewLine][viewColumn]*0.5 + residueBlock->mPixel[7][1][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[6][3][viewLine][viewColumn] = origBlock->mPixel[6][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][3][viewLine][viewColumn]*0.5 + residueBlock->mPixel[7][3][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[6][5][viewLine][viewColumn] = origBlock->mPixel[6][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][5][viewLine][viewColumn]*0.5 + residueBlock->mPixel[7][5][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[6][7][viewLine][viewColumn] = origBlock->mPixel[6][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][7][viewLine][viewColumn]*0.5 + residueBlock->mPixel[7][7][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[6][9][viewLine][viewColumn] = origBlock->mPixel[6][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][9][viewLine][viewColumn]*0.5 + residueBlock->mPixel[7][9][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[6][11][viewLine][viewColumn] = origBlock->mPixel[6][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[5][11][viewLine][viewColumn]*0.5 + residueBlock->mPixel[7][11][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[7][0][viewLine][viewColumn] = origBlock->mPixel[7][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][1][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[7][2][viewLine][viewColumn] = origBlock->mPixel[7][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][1][viewLine][viewColumn]*0.5 + residueBlock->mPixel[7][3][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[7][4][viewLine][viewColumn] = origBlock->mPixel[7][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][3][viewLine][viewColumn]*0.5 + residueBlock->mPixel[7][5][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[7][6][viewLine][viewColumn] = origBlock->mPixel[7][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][5][viewLine][viewColumn]*0.5 + residueBlock->mPixel[7][7][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[7][8][viewLine][viewColumn] = origBlock->mPixel[7][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][7][viewLine][viewColumn]*0.5 + residueBlock->mPixel[7][9][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[7][10][viewLine][viewColumn] = origBlock->mPixel[7][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][9][viewLine][viewColumn]*0.5 + residueBlock->mPixel[7][11][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[7][12][viewLine][viewColumn] = origBlock->mPixel[7][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][11][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[8][1][viewLine][viewColumn] = origBlock->mPixel[8][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][1][viewLine][viewColumn]*0.5 + residueBlock->mPixel[9][1][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[8][3][viewLine][viewColumn] = origBlock->mPixel[8][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][3][viewLine][viewColumn]*0.5 + residueBlock->mPixel[9][3][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[8][5][viewLine][viewColumn] = origBlock->mPixel[8][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][5][viewLine][viewColumn]*0.5 + residueBlock->mPixel[9][5][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[8][7][viewLine][viewColumn] = origBlock->mPixel[8][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][7][viewLine][viewColumn]*0.5 + residueBlock->mPixel[9][7][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[8][9][viewLine][viewColumn] = origBlock->mPixel[8][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][9][viewLine][viewColumn]*0.5 + residueBlock->mPixel[9][9][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[8][11][viewLine][viewColumn] = origBlock->mPixel[8][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[7][11][viewLine][viewColumn]*0.5 + residueBlock->mPixel[9][11][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[9][0][viewLine][viewColumn] = origBlock->mPixel[9][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][1][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[9][2][viewLine][viewColumn] = origBlock->mPixel[9][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][1][viewLine][viewColumn]*0.5 + residueBlock->mPixel[9][3][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[9][4][viewLine][viewColumn] = origBlock->mPixel[9][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][3][viewLine][viewColumn]*0.5 + residueBlock->mPixel[9][5][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[9][6][viewLine][viewColumn] = origBlock->mPixel[9][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][5][viewLine][viewColumn]*0.5 + residueBlock->mPixel[9][7][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[9][8][viewLine][viewColumn] = origBlock->mPixel[9][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][7][viewLine][viewColumn]*0.5 + residueBlock->mPixel[9][9][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[9][10][viewLine][viewColumn] = origBlock->mPixel[9][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][9][viewLine][viewColumn]*0.5 + residueBlock->mPixel[9][11][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[9][12][viewLine][viewColumn] = origBlock->mPixel[9][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][11][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[10][1][viewLine][viewColumn] = origBlock->mPixel[10][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][1][viewLine][viewColumn]*0.5 + residueBlock->mPixel[11][1][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[10][3][viewLine][viewColumn] = origBlock->mPixel[10][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][3][viewLine][viewColumn]*0.5 + residueBlock->mPixel[11][3][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[10][5][viewLine][viewColumn] = origBlock->mPixel[10][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][5][viewLine][viewColumn]*0.5 + residueBlock->mPixel[11][5][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[10][7][viewLine][viewColumn] = origBlock->mPixel[10][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][7][viewLine][viewColumn]*0.5 + residueBlock->mPixel[11][7][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[10][9][viewLine][viewColumn] = origBlock->mPixel[10][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][9][viewLine][viewColumn]*0.5 + residueBlock->mPixel[11][9][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[10][11][viewLine][viewColumn] = origBlock->mPixel[10][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[9][11][viewLine][viewColumn]*0.5 + residueBlock->mPixel[11][11][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[11][0][viewLine][viewColumn] = origBlock->mPixel[11][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][1][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[11][2][viewLine][viewColumn] = origBlock->mPixel[11][2][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][1][viewLine][viewColumn]*0.5 + residueBlock->mPixel[11][3][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[11][4][viewLine][viewColumn] = origBlock->mPixel[11][4][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][3][viewLine][viewColumn]*0.5 + residueBlock->mPixel[11][5][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[11][6][viewLine][viewColumn] = origBlock->mPixel[11][6][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][5][viewLine][viewColumn]*0.5 + residueBlock->mPixel[11][7][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[11][8][viewLine][viewColumn] = origBlock->mPixel[11][8][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][7][viewLine][viewColumn]*0.5 + residueBlock->mPixel[11][9][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[11][10][viewLine][viewColumn] = origBlock->mPixel[11][10][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][9][viewLine][viewColumn]*0.5 + residueBlock->mPixel[11][11][viewLine][viewColumn]*0.5);
			residueBlock->mPixel[11][12][viewLine][viewColumn] = origBlock->mPixel[11][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][11][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[12][0][viewLine][viewColumn] = origBlock->mPixel[12][0][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][1][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[12][1][viewLine][viewColumn] = origBlock->mPixel[12][1][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][1][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[12][3][viewLine][viewColumn] = origBlock->mPixel[12][3][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][3][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[12][5][viewLine][viewColumn] = origBlock->mPixel[12][5][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][5][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[12][7][viewLine][viewColumn] = origBlock->mPixel[12][7][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][7][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[12][9][viewLine][viewColumn] = origBlock->mPixel[12][9][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][9][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[12][11][viewLine][viewColumn] = origBlock->mPixel[12][11][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][11][viewLine][viewColumn]*1.0);
			residueBlock->mPixel[12][12][viewLine][viewColumn] = origBlock->mPixel[12][12][viewLine][viewColumn] - (origBlock->mPixel[6][6][viewLine][viewColumn] + residueBlock->mPixel[11][11][viewLine][viewColumn]*1.0);

		}
	}
}

void Prediction :: recHierarchicalDifferentialPrediction1Level(Block4D *recBlock, Block4D *origBlock) {
	for(int viewLine = 0; viewLine < 15; viewLine += 1) {
		for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
			origBlock->mPixel[6][6][viewLine][viewColumn] = recBlock->mPixel[6][6][viewLine][viewColumn];
			origBlock->mPixel[2][2][viewLine][viewColumn] = recBlock->mPixel[2][2][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn];
			origBlock->mPixel[2][6][viewLine][viewColumn] = recBlock->mPixel[2][6][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn];
			origBlock->mPixel[2][10][viewLine][viewColumn] = recBlock->mPixel[2][10][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn];
			origBlock->mPixel[6][2][viewLine][viewColumn] = recBlock->mPixel[6][2][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn];
			origBlock->mPixel[6][10][viewLine][viewColumn] = recBlock->mPixel[6][10][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn];
			origBlock->mPixel[10][2][viewLine][viewColumn] = recBlock->mPixel[10][2][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn];
			origBlock->mPixel[10][6][viewLine][viewColumn] = recBlock->mPixel[10][6][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn];
			origBlock->mPixel[10][10][viewLine][viewColumn] = recBlock->mPixel[10][10][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn];
			origBlock->mPixel[0][4][viewLine][viewColumn] = recBlock->mPixel[0][4][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[2][2][viewLine][viewColumn]*0.5 + recBlock->mPixel[2][6][viewLine][viewColumn]*0.5;
			origBlock->mPixel[0][8][viewLine][viewColumn] = recBlock->mPixel[0][8][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[2][6][viewLine][viewColumn]*0.5 + recBlock->mPixel[2][10][viewLine][viewColumn]*0.5;
			origBlock->mPixel[4][0][viewLine][viewColumn] = recBlock->mPixel[4][0][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[2][2][viewLine][viewColumn]*0.5 + recBlock->mPixel[6][2][viewLine][viewColumn]*0.5;
			origBlock->mPixel[4][4][viewLine][viewColumn] = recBlock->mPixel[4][4][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[2][2][viewLine][viewColumn]*0.33 + recBlock->mPixel[2][6][viewLine][viewColumn]*0.33 + recBlock->mPixel[6][2][viewLine][viewColumn]*0.33;
			origBlock->mPixel[4][8][viewLine][viewColumn] = recBlock->mPixel[4][8][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[2][6][viewLine][viewColumn]*0.33 + recBlock->mPixel[2][10][viewLine][viewColumn]*0.33 + recBlock->mPixel[6][10][viewLine][viewColumn]*0.33;
			origBlock->mPixel[4][12][viewLine][viewColumn] = recBlock->mPixel[4][12][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[2][10][viewLine][viewColumn]*0.5 + recBlock->mPixel[6][10][viewLine][viewColumn]*0.5;
			origBlock->mPixel[8][0][viewLine][viewColumn] = recBlock->mPixel[8][0][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[6][2][viewLine][viewColumn]*0.5 + recBlock->mPixel[10][2][viewLine][viewColumn]*0.5;
			origBlock->mPixel[8][4][viewLine][viewColumn] = recBlock->mPixel[8][4][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[6][2][viewLine][viewColumn]*0.33 + recBlock->mPixel[10][2][viewLine][viewColumn]*0.33 + recBlock->mPixel[10][6][viewLine][viewColumn]*0.33;
			origBlock->mPixel[8][8][viewLine][viewColumn] = recBlock->mPixel[8][8][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[6][10][viewLine][viewColumn]*0.33 + recBlock->mPixel[10][6][viewLine][viewColumn]*0.33 + recBlock->mPixel[10][10][viewLine][viewColumn]*0.33;
			origBlock->mPixel[8][12][viewLine][viewColumn] = recBlock->mPixel[8][12][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[6][10][viewLine][viewColumn]*0.5 + recBlock->mPixel[10][10][viewLine][viewColumn]*0.5;
			origBlock->mPixel[12][4][viewLine][viewColumn] = recBlock->mPixel[12][4][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[10][2][viewLine][viewColumn]*0.5 + recBlock->mPixel[10][6][viewLine][viewColumn]*0.5;
			origBlock->mPixel[12][8][viewLine][viewColumn] = recBlock->mPixel[12][8][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[10][6][viewLine][viewColumn]*0.5 + recBlock->mPixel[10][10][viewLine][viewColumn]*0.5;
			origBlock->mPixel[0][2][viewLine][viewColumn] = recBlock->mPixel[0][2][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[0][4][viewLine][viewColumn]*1.0;
			origBlock->mPixel[0][6][viewLine][viewColumn] = recBlock->mPixel[0][6][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[0][4][viewLine][viewColumn]*0.5 + recBlock->mPixel[0][8][viewLine][viewColumn]*0.5;
			origBlock->mPixel[0][10][viewLine][viewColumn] = recBlock->mPixel[0][10][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[0][8][viewLine][viewColumn]*1.0;
			origBlock->mPixel[2][0][viewLine][viewColumn] = recBlock->mPixel[2][0][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[4][0][viewLine][viewColumn]*1.0;
			origBlock->mPixel[2][4][viewLine][viewColumn] = recBlock->mPixel[2][4][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[0][4][viewLine][viewColumn]*0.5 + recBlock->mPixel[4][4][viewLine][viewColumn]*0.5;
			origBlock->mPixel[2][8][viewLine][viewColumn] = recBlock->mPixel[2][8][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[0][8][viewLine][viewColumn]*0.5 + recBlock->mPixel[4][8][viewLine][viewColumn]*0.5;
			origBlock->mPixel[2][12][viewLine][viewColumn] = recBlock->mPixel[2][12][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[4][12][viewLine][viewColumn]*1.0;
			origBlock->mPixel[4][2][viewLine][viewColumn] = recBlock->mPixel[4][2][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[4][0][viewLine][viewColumn]*0.5 + recBlock->mPixel[4][4][viewLine][viewColumn]*0.5;
			origBlock->mPixel[4][6][viewLine][viewColumn] = recBlock->mPixel[4][6][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[4][4][viewLine][viewColumn]*0.5 + recBlock->mPixel[4][8][viewLine][viewColumn]*0.5;
			origBlock->mPixel[4][10][viewLine][viewColumn] = recBlock->mPixel[4][10][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[4][8][viewLine][viewColumn]*0.5 + recBlock->mPixel[4][12][viewLine][viewColumn]*0.5;
			origBlock->mPixel[6][0][viewLine][viewColumn] = recBlock->mPixel[6][0][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[4][0][viewLine][viewColumn]*0.5 + recBlock->mPixel[8][0][viewLine][viewColumn]*0.5;
			origBlock->mPixel[6][4][viewLine][viewColumn] = recBlock->mPixel[6][4][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[4][4][viewLine][viewColumn]*0.5 + recBlock->mPixel[8][4][viewLine][viewColumn]*0.5;
			origBlock->mPixel[6][8][viewLine][viewColumn] = recBlock->mPixel[6][8][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[4][8][viewLine][viewColumn]*0.5 + recBlock->mPixel[8][8][viewLine][viewColumn]*0.5;
			origBlock->mPixel[6][12][viewLine][viewColumn] = recBlock->mPixel[6][12][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[4][12][viewLine][viewColumn]*0.5 + recBlock->mPixel[8][12][viewLine][viewColumn]*0.5;
			origBlock->mPixel[8][2][viewLine][viewColumn] = recBlock->mPixel[8][2][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[8][0][viewLine][viewColumn]*0.5 + recBlock->mPixel[8][4][viewLine][viewColumn]*0.5;
			origBlock->mPixel[8][6][viewLine][viewColumn] = recBlock->mPixel[8][6][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[8][4][viewLine][viewColumn]*0.5 + recBlock->mPixel[8][8][viewLine][viewColumn]*0.5;
			origBlock->mPixel[8][10][viewLine][viewColumn] = recBlock->mPixel[8][10][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[8][8][viewLine][viewColumn]*0.5 + recBlock->mPixel[8][12][viewLine][viewColumn]*0.5;
			origBlock->mPixel[10][0][viewLine][viewColumn] = recBlock->mPixel[10][0][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[8][0][viewLine][viewColumn]*1.0;
			origBlock->mPixel[10][4][viewLine][viewColumn] = recBlock->mPixel[10][4][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[8][4][viewLine][viewColumn]*0.5 + recBlock->mPixel[12][4][viewLine][viewColumn]*0.5;
			origBlock->mPixel[10][8][viewLine][viewColumn] = recBlock->mPixel[10][8][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[8][8][viewLine][viewColumn]*0.5 + recBlock->mPixel[12][8][viewLine][viewColumn]*0.5;
			origBlock->mPixel[10][12][viewLine][viewColumn] = recBlock->mPixel[10][12][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[8][12][viewLine][viewColumn]*1.0;
			origBlock->mPixel[12][2][viewLine][viewColumn] = recBlock->mPixel[12][2][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[12][4][viewLine][viewColumn]*1.0;
			origBlock->mPixel[12][6][viewLine][viewColumn] = recBlock->mPixel[12][6][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[12][4][viewLine][viewColumn]*0.5 + recBlock->mPixel[12][8][viewLine][viewColumn]*0.5;
			origBlock->mPixel[12][10][viewLine][viewColumn] = recBlock->mPixel[12][10][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[12][8][viewLine][viewColumn]*1.0;
			origBlock->mPixel[1][1][viewLine][viewColumn] = recBlock->mPixel[1][1][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[0][2][viewLine][viewColumn]*0.5 + recBlock->mPixel[2][0][viewLine][viewColumn]*0.5;
			origBlock->mPixel[1][3][viewLine][viewColumn] = recBlock->mPixel[1][3][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[0][2][viewLine][viewColumn]*0.5 + recBlock->mPixel[2][4][viewLine][viewColumn]*0.5;
			origBlock->mPixel[1][5][viewLine][viewColumn] = recBlock->mPixel[1][5][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[0][6][viewLine][viewColumn]*0.5 + recBlock->mPixel[2][4][viewLine][viewColumn]*0.5;
			origBlock->mPixel[1][7][viewLine][viewColumn] = recBlock->mPixel[1][7][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[0][6][viewLine][viewColumn]*0.5 + recBlock->mPixel[2][8][viewLine][viewColumn]*0.5;
			origBlock->mPixel[1][9][viewLine][viewColumn] = recBlock->mPixel[1][9][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[0][10][viewLine][viewColumn]*0.5 + recBlock->mPixel[2][8][viewLine][viewColumn]*0.5;
			origBlock->mPixel[1][11][viewLine][viewColumn] = recBlock->mPixel[1][11][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[0][10][viewLine][viewColumn]*0.5 + recBlock->mPixel[2][12][viewLine][viewColumn]*0.5;
			origBlock->mPixel[3][1][viewLine][viewColumn] = recBlock->mPixel[3][1][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[2][0][viewLine][viewColumn]*0.5 + recBlock->mPixel[4][2][viewLine][viewColumn]*0.5;
			origBlock->mPixel[3][3][viewLine][viewColumn] = recBlock->mPixel[3][3][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[2][4][viewLine][viewColumn]*0.5 + recBlock->mPixel[4][2][viewLine][viewColumn]*0.5;
			origBlock->mPixel[3][5][viewLine][viewColumn] = recBlock->mPixel[3][5][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[2][4][viewLine][viewColumn]*0.5 + recBlock->mPixel[4][6][viewLine][viewColumn]*0.5;
			origBlock->mPixel[3][7][viewLine][viewColumn] = recBlock->mPixel[3][7][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[2][8][viewLine][viewColumn]*0.5 + recBlock->mPixel[4][6][viewLine][viewColumn]*0.5;
			origBlock->mPixel[3][9][viewLine][viewColumn] = recBlock->mPixel[3][9][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[2][8][viewLine][viewColumn]*0.5 + recBlock->mPixel[4][10][viewLine][viewColumn]*0.5;
			origBlock->mPixel[3][11][viewLine][viewColumn] = recBlock->mPixel[3][11][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[2][12][viewLine][viewColumn]*0.5 + recBlock->mPixel[4][10][viewLine][viewColumn]*0.5;
			origBlock->mPixel[5][1][viewLine][viewColumn] = recBlock->mPixel[5][1][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[4][2][viewLine][viewColumn]*0.5 + recBlock->mPixel[6][0][viewLine][viewColumn]*0.5;
			origBlock->mPixel[5][3][viewLine][viewColumn] = recBlock->mPixel[5][3][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[4][2][viewLine][viewColumn]*0.5 + recBlock->mPixel[6][4][viewLine][viewColumn]*0.5;
			origBlock->mPixel[5][5][viewLine][viewColumn] = recBlock->mPixel[5][5][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[4][6][viewLine][viewColumn]*0.5 + recBlock->mPixel[6][4][viewLine][viewColumn]*0.5;
			origBlock->mPixel[5][7][viewLine][viewColumn] = recBlock->mPixel[5][7][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[4][6][viewLine][viewColumn]*0.5 + recBlock->mPixel[6][8][viewLine][viewColumn]*0.5;
			origBlock->mPixel[5][9][viewLine][viewColumn] = recBlock->mPixel[5][9][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[4][10][viewLine][viewColumn]*0.5 + recBlock->mPixel[6][8][viewLine][viewColumn]*0.5;
			origBlock->mPixel[5][11][viewLine][viewColumn] = recBlock->mPixel[5][11][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[4][10][viewLine][viewColumn]*0.5 + recBlock->mPixel[6][12][viewLine][viewColumn]*0.5;
			origBlock->mPixel[7][1][viewLine][viewColumn] = recBlock->mPixel[7][1][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[6][0][viewLine][viewColumn]*0.5 + recBlock->mPixel[8][2][viewLine][viewColumn]*0.5;
			origBlock->mPixel[7][3][viewLine][viewColumn] = recBlock->mPixel[7][3][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[6][4][viewLine][viewColumn]*0.5 + recBlock->mPixel[8][2][viewLine][viewColumn]*0.5;
			origBlock->mPixel[7][5][viewLine][viewColumn] = recBlock->mPixel[7][5][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[6][4][viewLine][viewColumn]*0.5 + recBlock->mPixel[8][6][viewLine][viewColumn]*0.5;
			origBlock->mPixel[7][7][viewLine][viewColumn] = recBlock->mPixel[7][7][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[6][8][viewLine][viewColumn]*0.5 + recBlock->mPixel[8][6][viewLine][viewColumn]*0.5;
			origBlock->mPixel[7][9][viewLine][viewColumn] = recBlock->mPixel[7][9][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[6][8][viewLine][viewColumn]*0.5 + recBlock->mPixel[8][10][viewLine][viewColumn]*0.5;
			origBlock->mPixel[7][11][viewLine][viewColumn] = recBlock->mPixel[7][11][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[6][12][viewLine][viewColumn]*0.5 + recBlock->mPixel[8][10][viewLine][viewColumn]*0.5;
			origBlock->mPixel[9][1][viewLine][viewColumn] = recBlock->mPixel[9][1][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[8][2][viewLine][viewColumn]*0.5 + recBlock->mPixel[10][0][viewLine][viewColumn]*0.5;
			origBlock->mPixel[9][3][viewLine][viewColumn] = recBlock->mPixel[9][3][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[8][2][viewLine][viewColumn]*0.5 + recBlock->mPixel[10][4][viewLine][viewColumn]*0.5;
			origBlock->mPixel[9][5][viewLine][viewColumn] = recBlock->mPixel[9][5][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[8][6][viewLine][viewColumn]*0.5 + recBlock->mPixel[10][4][viewLine][viewColumn]*0.5;
			origBlock->mPixel[9][7][viewLine][viewColumn] = recBlock->mPixel[9][7][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[8][6][viewLine][viewColumn]*0.5 + recBlock->mPixel[10][8][viewLine][viewColumn]*0.5;
			origBlock->mPixel[9][9][viewLine][viewColumn] = recBlock->mPixel[9][9][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[8][10][viewLine][viewColumn]*0.5 + recBlock->mPixel[10][8][viewLine][viewColumn]*0.5;
			origBlock->mPixel[9][11][viewLine][viewColumn] = recBlock->mPixel[9][11][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[8][10][viewLine][viewColumn]*0.5 + recBlock->mPixel[10][12][viewLine][viewColumn]*0.5;
			origBlock->mPixel[11][1][viewLine][viewColumn] = recBlock->mPixel[11][1][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[10][0][viewLine][viewColumn]*0.5 + recBlock->mPixel[12][2][viewLine][viewColumn]*0.5;
			origBlock->mPixel[11][3][viewLine][viewColumn] = recBlock->mPixel[11][3][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[10][4][viewLine][viewColumn]*0.5 + recBlock->mPixel[12][2][viewLine][viewColumn]*0.5;
			origBlock->mPixel[11][5][viewLine][viewColumn] = recBlock->mPixel[11][5][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[10][4][viewLine][viewColumn]*0.5 + recBlock->mPixel[12][6][viewLine][viewColumn]*0.5;
			origBlock->mPixel[11][7][viewLine][viewColumn] = recBlock->mPixel[11][7][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[10][8][viewLine][viewColumn]*0.5 + recBlock->mPixel[12][6][viewLine][viewColumn]*0.5;
			origBlock->mPixel[11][9][viewLine][viewColumn] = recBlock->mPixel[11][9][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[10][8][viewLine][viewColumn]*0.5 + recBlock->mPixel[12][10][viewLine][viewColumn]*0.5;
			origBlock->mPixel[11][11][viewLine][viewColumn] = recBlock->mPixel[11][11][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[10][12][viewLine][viewColumn]*0.5 + recBlock->mPixel[12][10][viewLine][viewColumn]*0.5;
			origBlock->mPixel[0][0][viewLine][viewColumn] = recBlock->mPixel[0][0][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][1][viewLine][viewColumn]*1.0;
			origBlock->mPixel[0][1][viewLine][viewColumn] = recBlock->mPixel[0][1][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][1][viewLine][viewColumn]*1.0;
			origBlock->mPixel[0][3][viewLine][viewColumn] = recBlock->mPixel[0][3][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][3][viewLine][viewColumn]*1.0;
			origBlock->mPixel[0][5][viewLine][viewColumn] = recBlock->mPixel[0][5][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][5][viewLine][viewColumn]*1.0;
			origBlock->mPixel[0][7][viewLine][viewColumn] = recBlock->mPixel[0][7][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][7][viewLine][viewColumn]*1.0;
			origBlock->mPixel[0][9][viewLine][viewColumn] = recBlock->mPixel[0][9][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][9][viewLine][viewColumn]*1.0;
			origBlock->mPixel[0][11][viewLine][viewColumn] = recBlock->mPixel[0][11][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][11][viewLine][viewColumn]*1.0;
			origBlock->mPixel[0][12][viewLine][viewColumn] = recBlock->mPixel[0][12][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][11][viewLine][viewColumn]*1.0;
			origBlock->mPixel[1][0][viewLine][viewColumn] = recBlock->mPixel[1][0][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][1][viewLine][viewColumn]*1.0;
			origBlock->mPixel[1][2][viewLine][viewColumn] = recBlock->mPixel[1][2][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][1][viewLine][viewColumn]*0.5 + recBlock->mPixel[1][3][viewLine][viewColumn]*0.5;
			origBlock->mPixel[1][4][viewLine][viewColumn] = recBlock->mPixel[1][4][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][3][viewLine][viewColumn]*0.5 + recBlock->mPixel[1][5][viewLine][viewColumn]*0.5;
			origBlock->mPixel[1][6][viewLine][viewColumn] = recBlock->mPixel[1][6][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][5][viewLine][viewColumn]*0.5 + recBlock->mPixel[1][7][viewLine][viewColumn]*0.5;
			origBlock->mPixel[1][8][viewLine][viewColumn] = recBlock->mPixel[1][8][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][7][viewLine][viewColumn]*0.5 + recBlock->mPixel[1][9][viewLine][viewColumn]*0.5;
			origBlock->mPixel[1][10][viewLine][viewColumn] = recBlock->mPixel[1][10][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][9][viewLine][viewColumn]*0.5 + recBlock->mPixel[1][11][viewLine][viewColumn]*0.5;
			origBlock->mPixel[1][12][viewLine][viewColumn] = recBlock->mPixel[1][12][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][11][viewLine][viewColumn]*1.0;
			origBlock->mPixel[2][1][viewLine][viewColumn] = recBlock->mPixel[2][1][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][1][viewLine][viewColumn]*0.5 + recBlock->mPixel[3][1][viewLine][viewColumn]*0.5;
			origBlock->mPixel[2][3][viewLine][viewColumn] = recBlock->mPixel[2][3][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][3][viewLine][viewColumn]*0.5 + recBlock->mPixel[3][3][viewLine][viewColumn]*0.5;
			origBlock->mPixel[2][5][viewLine][viewColumn] = recBlock->mPixel[2][5][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][5][viewLine][viewColumn]*0.5 + recBlock->mPixel[3][5][viewLine][viewColumn]*0.5;
			origBlock->mPixel[2][7][viewLine][viewColumn] = recBlock->mPixel[2][7][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][7][viewLine][viewColumn]*0.5 + recBlock->mPixel[3][7][viewLine][viewColumn]*0.5;
			origBlock->mPixel[2][9][viewLine][viewColumn] = recBlock->mPixel[2][9][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][9][viewLine][viewColumn]*0.5 + recBlock->mPixel[3][9][viewLine][viewColumn]*0.5;
			origBlock->mPixel[2][11][viewLine][viewColumn] = recBlock->mPixel[2][11][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[1][11][viewLine][viewColumn]*0.5 + recBlock->mPixel[3][11][viewLine][viewColumn]*0.5;
			origBlock->mPixel[3][0][viewLine][viewColumn] = recBlock->mPixel[3][0][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[3][1][viewLine][viewColumn]*1.0;
			origBlock->mPixel[3][2][viewLine][viewColumn] = recBlock->mPixel[3][2][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[3][1][viewLine][viewColumn]*0.5 + recBlock->mPixel[3][3][viewLine][viewColumn]*0.5;
			origBlock->mPixel[3][4][viewLine][viewColumn] = recBlock->mPixel[3][4][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[3][3][viewLine][viewColumn]*0.5 + recBlock->mPixel[3][5][viewLine][viewColumn]*0.5;
			origBlock->mPixel[3][6][viewLine][viewColumn] = recBlock->mPixel[3][6][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[3][5][viewLine][viewColumn]*0.5 + recBlock->mPixel[3][7][viewLine][viewColumn]*0.5;
			origBlock->mPixel[3][8][viewLine][viewColumn] = recBlock->mPixel[3][8][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[3][7][viewLine][viewColumn]*0.5 + recBlock->mPixel[3][9][viewLine][viewColumn]*0.5;
			origBlock->mPixel[3][10][viewLine][viewColumn] = recBlock->mPixel[3][10][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[3][9][viewLine][viewColumn]*0.5 + recBlock->mPixel[3][11][viewLine][viewColumn]*0.5;
			origBlock->mPixel[3][12][viewLine][viewColumn] = recBlock->mPixel[3][12][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[3][11][viewLine][viewColumn]*1.0;
			origBlock->mPixel[4][1][viewLine][viewColumn] = recBlock->mPixel[4][1][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[3][1][viewLine][viewColumn]*0.5 + recBlock->mPixel[5][1][viewLine][viewColumn]*0.5;
			origBlock->mPixel[4][3][viewLine][viewColumn] = recBlock->mPixel[4][3][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[3][3][viewLine][viewColumn]*0.5 + recBlock->mPixel[5][3][viewLine][viewColumn]*0.5;
			origBlock->mPixel[4][5][viewLine][viewColumn] = recBlock->mPixel[4][5][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[3][5][viewLine][viewColumn]*0.5 + recBlock->mPixel[5][5][viewLine][viewColumn]*0.5;
			origBlock->mPixel[4][7][viewLine][viewColumn] = recBlock->mPixel[4][7][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[3][7][viewLine][viewColumn]*0.5 + recBlock->mPixel[5][7][viewLine][viewColumn]*0.5;
			origBlock->mPixel[4][9][viewLine][viewColumn] = recBlock->mPixel[4][9][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[3][9][viewLine][viewColumn]*0.5 + recBlock->mPixel[5][9][viewLine][viewColumn]*0.5;
			origBlock->mPixel[4][11][viewLine][viewColumn] = recBlock->mPixel[4][11][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[3][11][viewLine][viewColumn]*0.5 + recBlock->mPixel[5][11][viewLine][viewColumn]*0.5;
			origBlock->mPixel[5][0][viewLine][viewColumn] = recBlock->mPixel[5][0][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[5][1][viewLine][viewColumn]*1.0;
			origBlock->mPixel[5][2][viewLine][viewColumn] = recBlock->mPixel[5][2][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[5][1][viewLine][viewColumn]*0.5 + recBlock->mPixel[5][3][viewLine][viewColumn]*0.5;
			origBlock->mPixel[5][4][viewLine][viewColumn] = recBlock->mPixel[5][4][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[5][3][viewLine][viewColumn]*0.5 + recBlock->mPixel[5][5][viewLine][viewColumn]*0.5;
			origBlock->mPixel[5][6][viewLine][viewColumn] = recBlock->mPixel[5][6][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[5][5][viewLine][viewColumn]*0.5 + recBlock->mPixel[5][7][viewLine][viewColumn]*0.5;
			origBlock->mPixel[5][8][viewLine][viewColumn] = recBlock->mPixel[5][8][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[5][7][viewLine][viewColumn]*0.5 + recBlock->mPixel[5][9][viewLine][viewColumn]*0.5;
			origBlock->mPixel[5][10][viewLine][viewColumn] = recBlock->mPixel[5][10][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[5][9][viewLine][viewColumn]*0.5 + recBlock->mPixel[5][11][viewLine][viewColumn]*0.5;
			origBlock->mPixel[5][12][viewLine][viewColumn] = recBlock->mPixel[5][12][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[5][11][viewLine][viewColumn]*1.0;
			origBlock->mPixel[6][1][viewLine][viewColumn] = recBlock->mPixel[6][1][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[5][1][viewLine][viewColumn]*0.5 + recBlock->mPixel[7][1][viewLine][viewColumn]*0.5;
			origBlock->mPixel[6][3][viewLine][viewColumn] = recBlock->mPixel[6][3][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[5][3][viewLine][viewColumn]*0.5 + recBlock->mPixel[7][3][viewLine][viewColumn]*0.5;
			origBlock->mPixel[6][5][viewLine][viewColumn] = recBlock->mPixel[6][5][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[5][5][viewLine][viewColumn]*0.5 + recBlock->mPixel[7][5][viewLine][viewColumn]*0.5;
			origBlock->mPixel[6][7][viewLine][viewColumn] = recBlock->mPixel[6][7][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[5][7][viewLine][viewColumn]*0.5 + recBlock->mPixel[7][7][viewLine][viewColumn]*0.5;
			origBlock->mPixel[6][9][viewLine][viewColumn] = recBlock->mPixel[6][9][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[5][9][viewLine][viewColumn]*0.5 + recBlock->mPixel[7][9][viewLine][viewColumn]*0.5;
			origBlock->mPixel[6][11][viewLine][viewColumn] = recBlock->mPixel[6][11][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[5][11][viewLine][viewColumn]*0.5 + recBlock->mPixel[7][11][viewLine][viewColumn]*0.5;
			origBlock->mPixel[7][0][viewLine][viewColumn] = recBlock->mPixel[7][0][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[7][1][viewLine][viewColumn]*1.0;
			origBlock->mPixel[7][2][viewLine][viewColumn] = recBlock->mPixel[7][2][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[7][1][viewLine][viewColumn]*0.5 + recBlock->mPixel[7][3][viewLine][viewColumn]*0.5;
			origBlock->mPixel[7][4][viewLine][viewColumn] = recBlock->mPixel[7][4][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[7][3][viewLine][viewColumn]*0.5 + recBlock->mPixel[7][5][viewLine][viewColumn]*0.5;
			origBlock->mPixel[7][6][viewLine][viewColumn] = recBlock->mPixel[7][6][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[7][5][viewLine][viewColumn]*0.5 + recBlock->mPixel[7][7][viewLine][viewColumn]*0.5;
			origBlock->mPixel[7][8][viewLine][viewColumn] = recBlock->mPixel[7][8][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[7][7][viewLine][viewColumn]*0.5 + recBlock->mPixel[7][9][viewLine][viewColumn]*0.5;
			origBlock->mPixel[7][10][viewLine][viewColumn] = recBlock->mPixel[7][10][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[7][9][viewLine][viewColumn]*0.5 + recBlock->mPixel[7][11][viewLine][viewColumn]*0.5;
			origBlock->mPixel[7][12][viewLine][viewColumn] = recBlock->mPixel[7][12][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[7][11][viewLine][viewColumn]*1.0;
			origBlock->mPixel[8][1][viewLine][viewColumn] = recBlock->mPixel[8][1][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[7][1][viewLine][viewColumn]*0.5 + recBlock->mPixel[9][1][viewLine][viewColumn]*0.5;
			origBlock->mPixel[8][3][viewLine][viewColumn] = recBlock->mPixel[8][3][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[7][3][viewLine][viewColumn]*0.5 + recBlock->mPixel[9][3][viewLine][viewColumn]*0.5;
			origBlock->mPixel[8][5][viewLine][viewColumn] = recBlock->mPixel[8][5][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[7][5][viewLine][viewColumn]*0.5 + recBlock->mPixel[9][5][viewLine][viewColumn]*0.5;
			origBlock->mPixel[8][7][viewLine][viewColumn] = recBlock->mPixel[8][7][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[7][7][viewLine][viewColumn]*0.5 + recBlock->mPixel[9][7][viewLine][viewColumn]*0.5;
			origBlock->mPixel[8][9][viewLine][viewColumn] = recBlock->mPixel[8][9][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[7][9][viewLine][viewColumn]*0.5 + recBlock->mPixel[9][9][viewLine][viewColumn]*0.5;
			origBlock->mPixel[8][11][viewLine][viewColumn] = recBlock->mPixel[8][11][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[7][11][viewLine][viewColumn]*0.5 + recBlock->mPixel[9][11][viewLine][viewColumn]*0.5;
			origBlock->mPixel[9][0][viewLine][viewColumn] = recBlock->mPixel[9][0][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[9][1][viewLine][viewColumn]*1.0;
			origBlock->mPixel[9][2][viewLine][viewColumn] = recBlock->mPixel[9][2][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[9][1][viewLine][viewColumn]*0.5 + recBlock->mPixel[9][3][viewLine][viewColumn]*0.5;
			origBlock->mPixel[9][4][viewLine][viewColumn] = recBlock->mPixel[9][4][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[9][3][viewLine][viewColumn]*0.5 + recBlock->mPixel[9][5][viewLine][viewColumn]*0.5;
			origBlock->mPixel[9][6][viewLine][viewColumn] = recBlock->mPixel[9][6][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[9][5][viewLine][viewColumn]*0.5 + recBlock->mPixel[9][7][viewLine][viewColumn]*0.5;
			origBlock->mPixel[9][8][viewLine][viewColumn] = recBlock->mPixel[9][8][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[9][7][viewLine][viewColumn]*0.5 + recBlock->mPixel[9][9][viewLine][viewColumn]*0.5;
			origBlock->mPixel[9][10][viewLine][viewColumn] = recBlock->mPixel[9][10][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[9][9][viewLine][viewColumn]*0.5 + recBlock->mPixel[9][11][viewLine][viewColumn]*0.5;
			origBlock->mPixel[9][12][viewLine][viewColumn] = recBlock->mPixel[9][12][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[9][11][viewLine][viewColumn]*1.0;
			origBlock->mPixel[10][1][viewLine][viewColumn] = recBlock->mPixel[10][1][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[9][1][viewLine][viewColumn]*0.5 + recBlock->mPixel[11][1][viewLine][viewColumn]*0.5;
			origBlock->mPixel[10][3][viewLine][viewColumn] = recBlock->mPixel[10][3][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[9][3][viewLine][viewColumn]*0.5 + recBlock->mPixel[11][3][viewLine][viewColumn]*0.5;
			origBlock->mPixel[10][5][viewLine][viewColumn] = recBlock->mPixel[10][5][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[9][5][viewLine][viewColumn]*0.5 + recBlock->mPixel[11][5][viewLine][viewColumn]*0.5;
			origBlock->mPixel[10][7][viewLine][viewColumn] = recBlock->mPixel[10][7][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[9][7][viewLine][viewColumn]*0.5 + recBlock->mPixel[11][7][viewLine][viewColumn]*0.5;
			origBlock->mPixel[10][9][viewLine][viewColumn] = recBlock->mPixel[10][9][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[9][9][viewLine][viewColumn]*0.5 + recBlock->mPixel[11][9][viewLine][viewColumn]*0.5;
			origBlock->mPixel[10][11][viewLine][viewColumn] = recBlock->mPixel[10][11][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[9][11][viewLine][viewColumn]*0.5 + recBlock->mPixel[11][11][viewLine][viewColumn]*0.5;
			origBlock->mPixel[11][0][viewLine][viewColumn] = recBlock->mPixel[11][0][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[11][1][viewLine][viewColumn]*1.0;
			origBlock->mPixel[11][2][viewLine][viewColumn] = recBlock->mPixel[11][2][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[11][1][viewLine][viewColumn]*0.5 + recBlock->mPixel[11][3][viewLine][viewColumn]*0.5;
			origBlock->mPixel[11][4][viewLine][viewColumn] = recBlock->mPixel[11][4][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[11][3][viewLine][viewColumn]*0.5 + recBlock->mPixel[11][5][viewLine][viewColumn]*0.5;
			origBlock->mPixel[11][6][viewLine][viewColumn] = recBlock->mPixel[11][6][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[11][5][viewLine][viewColumn]*0.5 + recBlock->mPixel[11][7][viewLine][viewColumn]*0.5;
			origBlock->mPixel[11][8][viewLine][viewColumn] = recBlock->mPixel[11][8][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[11][7][viewLine][viewColumn]*0.5 + recBlock->mPixel[11][9][viewLine][viewColumn]*0.5;
			origBlock->mPixel[11][10][viewLine][viewColumn] = recBlock->mPixel[11][10][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[11][9][viewLine][viewColumn]*0.5 + recBlock->mPixel[11][11][viewLine][viewColumn]*0.5;
			origBlock->mPixel[11][12][viewLine][viewColumn] = recBlock->mPixel[11][12][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[11][11][viewLine][viewColumn]*1.0;
			origBlock->mPixel[12][0][viewLine][viewColumn] = recBlock->mPixel[12][0][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[11][1][viewLine][viewColumn]*1.0;
			origBlock->mPixel[12][1][viewLine][viewColumn] = recBlock->mPixel[12][1][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[11][1][viewLine][viewColumn]*1.0;
			origBlock->mPixel[12][3][viewLine][viewColumn] = recBlock->mPixel[12][3][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[11][3][viewLine][viewColumn]*1.0;
			origBlock->mPixel[12][5][viewLine][viewColumn] = recBlock->mPixel[12][5][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[11][5][viewLine][viewColumn]*1.0;
			origBlock->mPixel[12][7][viewLine][viewColumn] = recBlock->mPixel[12][7][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[11][7][viewLine][viewColumn]*1.0;
			origBlock->mPixel[12][9][viewLine][viewColumn] = recBlock->mPixel[12][9][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[11][9][viewLine][viewColumn]*1.0;
			origBlock->mPixel[12][11][viewLine][viewColumn] = recBlock->mPixel[12][11][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[11][11][viewLine][viewColumn]*1.0;
			origBlock->mPixel[12][12][viewLine][viewColumn] = recBlock->mPixel[12][12][viewLine][viewColumn] + recBlock->mPixel[6][6][viewLine][viewColumn] + recBlock->mPixel[11][11][viewLine][viewColumn]*1.0;

		}
	}
}

void Prediction :: printOneBlock(Block4D *lfBlock) {
	for(int viewLine = 0; viewLine < 15; viewLine += 1) {
		for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
			printf("[coeff]: %d\n", lfBlock->mPixel[3][3][viewLine][viewColumn]);
		}
	}
}

void Prediction :: calcReferencePlaneEnergy(Block4D *lfBlock, int spectralComponent) {
	int coeff;

	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
		for(int viewLine = 0; viewLine < 15; viewLine += 1) {
			for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
				if(predictionType == 1) {
					//printf("\tdiffC\n");
					coeff = lfBlock->mPixel[verticalView][6][viewLine][viewColumn];
				}
				else {
					coeff = lfBlock->mPixel[verticalView][0][viewLine][viewColumn];
				}
				if(spectralComponent == 0)
					y_totalSignalEnergyFirstPlane += coeff*coeff;
				if(spectralComponent == 1)
					cb_totalSignalEnergyFirstPlane += coeff*coeff;
				if(spectralComponent == 2)
					cr_totalSignalEnergyFirstPlane += coeff*coeff;

				if(coeff > maxRefPlane) {
					maxRefPlane = coeff;
				}
				if(coeff < minRefPlane) {
					minRefPlane = coeff;
				}
			}
		}
	}
}

void Prediction :: calcOtherPlanesEnergy(Block4D *lfBlock, int spectralComponent) {
	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
		for(int horizontalView = 0; horizontalView < 13; horizontalView += 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
					if(predictionType == 1) {
						//printf("\tdiffC\n");
						if(horizontalView != 6) {
							int coeff = lfBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn];
							
							if(spectralComponent == 0)
								y_totalSignalEnergyOtherPlanes += coeff*coeff;
							if(spectralComponent == 1)
								cb_totalSignalEnergyOtherPlanes += coeff*coeff;
							if(spectralComponent == 2)
								cr_totalSignalEnergyOtherPlanes += coeff*coeff;

							if(coeff > maxOtherPlanes) {
								maxOtherPlanes = coeff;
							}
							if(coeff < minOtherPlanes) {
								minOtherPlanes = coeff;
							}
						}
					}
					else {
						if(horizontalView != 0) {
							int coeff = lfBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn];
							
							if(spectralComponent == 0)
								y_totalSignalEnergyOtherPlanes += coeff*coeff;
							if(spectralComponent == 1)
								cb_totalSignalEnergyOtherPlanes += coeff*coeff;
							if(spectralComponent == 2)
								cr_totalSignalEnergyOtherPlanes += coeff*coeff;

							if(coeff > maxOtherPlanes) {
								maxOtherPlanes = coeff;
							}
							if(coeff < minOtherPlanes) {
								minOtherPlanes = coeff;
							}
						}
					}
				}
			}
		}
	}
}

void Prediction :: saveDCCoeff(Block4D *lfBlock) {
	DC_coeff << lfBlock->mPixel[0][0][0][0] << '\n';
}

void Prediction :: saveACCoeff(Block4D *lfBlock) {
	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
		for(int horizontalView = 0; horizontalView < 13; horizontalView += 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
					if(horizontalView > 0 && verticalView > 0) {
						AC_coeff << lfBlock->mPixel[verticalView][horizontalView][viewLine][viewColumn] << '\n';
					}
				}
			}
		}
	}
}

int Prediction :: getMaxRefPlane() {
	return maxRefPlane;
}

int Prediction :: getMinRefPlane() {
	return minRefPlane;
}

int Prediction :: getMaxOtherPlanes() {
	return maxOtherPlanes;
}

int Prediction :: getMinOtherPlanes() {
	return minOtherPlanes;
}

double Prediction :: getYFirstPlaneEnergy() {
	return y_totalSignalEnergyFirstPlane;
}

double Prediction :: getCbFirstPlaneEnergy() {
	return cb_totalSignalEnergyFirstPlane;
}

double Prediction :: getCrFirstPlaneEnergy() {
	return cr_totalSignalEnergyFirstPlane;
}

double Prediction :: getYOtherPlanesEnergy() {
	return y_totalSignalEnergyOtherPlanes;
}

double Prediction :: getCbOtherPlanesEnergy() {
	return cb_totalSignalEnergyOtherPlanes;
}

double Prediction :: getCrOtherPlanesEnergy() {
	return cr_totalSignalEnergyOtherPlanes;
}