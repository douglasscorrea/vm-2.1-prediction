#include "Statistics.h"
#include "Block4D.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#define refPlaneSamples 3562650;
#define otherPlanesSamples 42751800;

Statistics :: Statistics(int prediction) {
	predictionType = prediction;

	if(prediction == 0) {
		//DC_coeff.open("DC_coeff_diffR.txt");
		//AC_coeff.open("AC_coeff_diffR.txt");
		//refPlaneCoefficients.open("refPlaneCoeff_diffR.txt");
		//allCoefficients.open("allCoefficients_diffR.txt");
		//yResiduesFile.open("y_residues_diffR.txt");
		//cbResiduesFile.open("cb_residues_diffR.txt");
		//crResiduesFile.open("cr_residues_diffR.txt");
	}
	else if(prediction == 1) {
		//DC_coeff.open("DC_coeff_diffC.txt");
		//AC_coeff.open("AC_coeff_diffC.txt");
		//refPlaneCoefficients.open("refPlaneCoeff_diffC.txt");
		//allCoefficients.open("allCoefficients_diffC.txt");
		//yResiduesFile.open("y_residues_diffC.txt");
		//cbResiduesFile.open("cb_residues_diffC.txt");
		//crResiduesFile.open("cr_residues_diffC.txt");
	}
	else {
		//DC_coeff.open("DC_coeff_mule.txt");
		//AC_coeff.open("AC_coeff_mule.txt");
		//refPlaneCoefficients.open("refPlaneCoeff_MuLE.txt");
		//allCoefficients.open("allCoefficients_MuLE.txt");
		//yResiduesFile.open("y_samples_MuLE.txt");
		//cbResiduesFile.open("cb_samples_MuLE.txt");
		//crResiduesFile.open("cr_samples_MuLE.txt");
	}

	yRefPlaneSamples = new int[3562650];
	cbRefPlaneSamples = new int[3562650];
	crRefPlaneSamples = new int[3562650];
	yOtherPlanesSamples = new int[42751800];
	cbOtherPlanesSamples = new int[42751800];
	crOtherPlanesSamples = new int[42751800];
	refPlaneCoeff = new int[3562650];;
	otherPlanesCoeff =  new int[42751800];
	AllCoeffs = new int[3562650+42751800];

	counterAllCoeffs = 0;
	yCounterRP = 0;
	cbCounterRP = 0;
	crCounterRP = 0;
	yCounterOP = 0;
	cbCounterOP = 0;
	crCounterOP = 0;
	counterOPCoeffs = 0;
	counterRPCoeffs = 0;
	energyRefPlaneCoeff = 0;
	energyOtherPlanesCoeff = 0;
	sumRefPlaneCoeff = 0;
	sumOtherPlanesCoeff = 0;
	coeffEnergy = 0;
	sumCoeff = 0;
	y_totalSignalEnergyRefPlane = 0;
	cb_totalSignalEnergyRefPlane = 0;
	cr_totalSignalEnergyRefPlane = 0;
	y_totalSignalEnergyOtherPlanes = 0;
	cb_totalSignalEnergyOtherPlanes = 0;
	cr_totalSignalEnergyOtherPlanes = 0;
	y_samplesSum = 0;
	cb_samplesSum = 0;
	cr_samplesSum = 0;
	y_refPlaneSamplesSum = 0;
	cb_refPlaneSamplesSum = 0;
	cr_refPlaneSamplesSum = 0;
	y_otherPlanesSamplesSum = 0;
	cb_otherPlanesSamplesSum = 0;
	cr_otherPlanesSamplesSum = 0;

	maxRefPlane = -9999999;
	minRefPlane = 9999999;
	maxOtherPlanes = -9999999;
	minOtherPlanes = 9999999;
}

void Statistics :: printOneBlock(Block4D *lfBlock) {
	for(int viewLine = 0; viewLine < 15; viewLine += 1) {
		for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
			printf("[coeff]: %d\n", lfBlock->mPixel[3][3][viewLine][viewColumn]);
		}
	}
}

void Statistics :: calcReferencePlaneEnergy(Block4D *lfBlock, int spectralComponent) {
	int samples;

	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
		for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				if(predictionType == 1) {
					//printf("\tdiffC\n");
					samples = lfBlock->mPixel[6][verticalView][viewColumn][viewLine];
				}
				else {
					samples = lfBlock->mPixel[0][verticalView][viewColumn][viewLine];
				}
				if(spectralComponent == 0)
					y_totalSignalEnergyRefPlane += samples*samples;
				if(spectralComponent == 1)
					cb_totalSignalEnergyRefPlane += samples*samples;
				if(spectralComponent == 2)
					cr_totalSignalEnergyRefPlane += samples*samples;

				if(samples > maxRefPlane) {
					maxRefPlane = samples;
				}
				if(samples < minRefPlane) {
					minRefPlane = samples;
				}
			}
		}
	}
}

void Statistics :: calcOtherPlanesEnergy(Block4D *lfBlock, int spectralComponent) {
	for(int horizontalView = 1; horizontalView < 13; horizontalView += 1) {
		for(int verticalView = 0; verticalView < 13; verticalView += 1) {
			for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
				for(int viewLine = 0; viewLine < 15; viewLine += 1) {
					if(predictionType == 1) {
						//printf("\tdiffC\n");
						if(horizontalView != 6) {
							int samples = lfBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
							
							if(spectralComponent == 0)
								y_totalSignalEnergyOtherPlanes += samples*samples;
							if(spectralComponent == 1)
								cb_totalSignalEnergyOtherPlanes += samples*samples;
							if(spectralComponent == 2)
								cr_totalSignalEnergyOtherPlanes += samples*samples;

							if(samples > maxOtherPlanes) {
								maxOtherPlanes = samples;
							}
							if(samples < minOtherPlanes) {
								minOtherPlanes = samples;
							}
						}
					}
					else {
						if(horizontalView != 0) {
							int samples = lfBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
							
							if(spectralComponent == 0)
								y_totalSignalEnergyOtherPlanes += samples*samples;
							if(spectralComponent == 1)
								cb_totalSignalEnergyOtherPlanes += samples*samples;
							if(spectralComponent == 2)
								cr_totalSignalEnergyOtherPlanes += samples*samples;

							if(samples > maxOtherPlanes) {
								maxOtherPlanes = samples;
							}
							if(samples < minOtherPlanes) {
								minOtherPlanes = samples;
							}
						}
					}
				}
			}
		}
	}
}

void Statistics :: saveDCCoeff(Block4D *lfBlock) {
	DC_coeff << lfBlock->mPixel[0][0][0][0] << '\n';
}

void Statistics :: saveACCoeff(Block4D *lfBlock) {
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

int Statistics :: getMaxRefPlane() {
	return maxRefPlane;
}

int Statistics :: getMinRefPlane() {
	return minRefPlane;
}

int Statistics :: getMaxOtherPlanes() {
	return maxOtherPlanes;
}

int Statistics :: getMinOtherPlanes() {
	return minOtherPlanes;
}

double Statistics :: getYRefPlaneEnergy() {
	return y_totalSignalEnergyRefPlane;
}

double Statistics :: getCbRefPlaneEnergy() {
	return cb_totalSignalEnergyRefPlane;
}

double Statistics :: getCrRefPlaneEnergy() {
	return cr_totalSignalEnergyRefPlane;
}

double Statistics :: getYOtherPlanesEnergy() {
	return y_totalSignalEnergyOtherPlanes;
}

double Statistics :: getCbOtherPlanesEnergy() {
	return cb_totalSignalEnergyOtherPlanes;
}

double Statistics :: getCrOtherPlanesEnergy() {
	return cr_totalSignalEnergyOtherPlanes;
}

void Statistics :: calcSumYRefPlaneSamples(Block4D *lfBlock) {
	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
		for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				yRefPlaneSamples[yCounterRP] = lfBlock->mPixel[0][verticalView][viewColumn][viewLine];
				//printf("%d\n", lfBlock->mPixel[0][verticalView][viewColumn][viewLine]);
				y_refPlaneSamplesSum += lfBlock->mPixel[0][verticalView][viewColumn][viewLine];
				yCounterRP++;
			}
		}
	}
}

void Statistics :: calcSumCbRefPlaneSamples(Block4D *lfBlock) {
	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
		for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				cbRefPlaneSamples[cbCounterRP] = lfBlock->mPixel[0][verticalView][viewColumn][viewLine];
				cb_refPlaneSamplesSum += lfBlock->mPixel[0][verticalView][viewColumn][viewLine];
				cbCounterRP++;
			}
		}
	}
}

void Statistics :: calcSumCrRefPlaneSamples(Block4D *lfBlock) {
	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
		for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				crRefPlaneSamples[crCounterRP] = lfBlock->mPixel[0][verticalView][viewColumn][viewLine];
				cr_refPlaneSamplesSum += lfBlock->mPixel[0][verticalView][viewColumn][viewLine];
				crCounterRP++;
			}
		}
	}
}

void Statistics :: calcSumYOtherPlanesSamples(Block4D *lfBlock) {
	for(int horizontalView = 1; horizontalView < 13; horizontalView += 1) {
		for(int verticalView = 0; verticalView < 13; verticalView += 1) {
			for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
				for(int viewLine = 0; viewLine < 15; viewLine += 1) {
					yOtherPlanesSamples[yCounterOP] = lfBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
					y_otherPlanesSamplesSum += lfBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
					yCounterOP++;
				}
			}
		}
	}
}

void Statistics :: calcSumCbOtherPlanesSamples(Block4D *lfBlock) {
	for(int horizontalView = 1; horizontalView < 13; horizontalView += 1) {
		for(int verticalView = 0; verticalView < 13; verticalView += 1) {
			for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
				for(int viewLine = 0; viewLine < 15; viewLine += 1) {
					cbOtherPlanesSamples[cbCounterOP] = lfBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
					cb_otherPlanesSamplesSum += lfBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
					cbCounterOP++;
				}
			}
		}
	}
}

void Statistics :: calcSumCrOtherPlanesSamples(Block4D *lfBlock) {
	for(int horizontalView = 1; horizontalView < 13; horizontalView += 1) {
		for(int verticalView = 0; verticalView < 13; verticalView += 1) {
			for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
				for(int viewLine = 0; viewLine < 15; viewLine += 1) {
					crOtherPlanesSamples[crCounterOP] = lfBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
					cr_otherPlanesSamplesSum += lfBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
					crCounterOP++;
				}
			}
		}
	}
}

double Statistics :: calcStdYRefPlane() {
	double stdDev = 0;
	double average = (double)y_refPlaneSamplesSum/3562650;

	for(int i = 0; i < 3562650; i++) {
		stdDev += pow(yRefPlaneSamples[i] - average, 2);
	}

	return sqrt(stdDev/3562650);
}

double Statistics :: calcStdCbRefPlane() {
	double stdDev = 0;
	double average = (double)cb_refPlaneSamplesSum/3562650;

	for(int i = 0; i < 3562650; i++) {
		stdDev += pow(cbRefPlaneSamples[i] - average, 2);
	}
	
	return sqrt(stdDev/3562650);
}

double Statistics :: calcStdCrRefPlane() {
	double stdDev = 0;
	double average = (double)cr_refPlaneSamplesSum/3562650;

	for(int i = 0; i < 3562650; i++) {
		stdDev += pow(crRefPlaneSamples[i] - average, 2);
	}
	
	return sqrt(stdDev/3562650);
}

double Statistics :: calcStdYOtherPlanes() {
	double stdDev = 0;
	double average = (double)y_otherPlanesSamplesSum/42751800;

	for(int i = 0; i < 42681600; i++) {
		stdDev += pow(yOtherPlanesSamples[i] - average, 2);
	}
	
	return sqrt(stdDev/42751800);
}

double Statistics :: calcStdCbOtherPlanes() {
	double stdDev = 0;
	double average = (double)cb_otherPlanesSamplesSum/42751800;

	for(int i = 0; i < 42681600; i++) {
		stdDev += pow(cbOtherPlanesSamples[i] - average, 2);
	}
	
	return sqrt(stdDev/42751800);
}

double Statistics :: calcStdCrOtherPlanes() {
	double stdDev = 0;
	double average = (double)cr_otherPlanesSamplesSum/42751800;

	for(int i = 0; i < 42681600; i++) {
		stdDev += pow(crOtherPlanesSamples[i] - average, 2);
	}
	
	return sqrt(stdDev/42751800);
}



double Statistics :: getYRefPlaneSamplesAverage() {
	return y_refPlaneSamplesSum/refPlaneSamples;
}

double Statistics :: getCbRefPlaneSamplesAverage() {
	return cb_refPlaneSamplesSum/refPlaneSamples;
}

double Statistics :: getCrRefPlaneSamplesAverage() {
	return cr_refPlaneSamplesSum/refPlaneSamples;
}

double Statistics :: getYOtherPlanesSamplesAverage() {
	return y_otherPlanesSamplesSum/otherPlanesSamples;
}

double Statistics :: getCbOtherPlanesSamplesAverage() {
	return cb_otherPlanesSamplesSum/otherPlanesSamples;
}

double Statistics :: getCrOtherPlanesSamplesAverage() {
	return cr_otherPlanesSamplesSum/otherPlanesSamples;
}

void Statistics :: calcRefPlaneCoeffEnergy(Block4D *transformedBlock) {
	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
		for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				int coeff = transformedBlock->mPixel[0][verticalView][viewColumn][viewLine];

				energyRefPlaneCoeff += coeff*coeff;
				//printf("%.2lf\n", energyRefPlaneCoeff);
				// if(coeff > maxRefPlane) {
				// 	maxRefPlane = sacoeffmples;
				// }
				// if(coeff < minRefPlane) {
				// 	minRefPlane = coeff;
				// }
			}
		}
	}
}

void Statistics :: calcOtherPlanesCoeffEnergy(Block4D *transformedBlock) {
	for(int horizontalView = 1; horizontalView < 13; horizontalView += 1) {
		for(int verticalView = 0; verticalView < 13; verticalView += 1) {
			for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
				for(int viewLine = 0; viewLine < 15; viewLine += 1) {
					int coeff = transformedBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
					
					energyOtherPlanesCoeff += coeff*coeff;
					//printf("%.2lf\n", energyOtherPlanesCoeff);
					// if(coeff > maxOtherPlanes) {
					// 	maxOtherPlanes = coeff;
					// }
					// if(coeff < minOtherPlanes) {
					// 	minOtherPlanes = coeff;
					// }
				}
			}
		}
	}
}

void Statistics :: calcSumRefPlaneCoeff(Block4D *transformedBlock) {
	for(int verticalView = 0; verticalView < 13; verticalView += 1) {
		for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
			for(int viewLine = 0; viewLine < 15; viewLine += 1) {
				int coeff = transformedBlock->mPixel[0][verticalView][viewColumn][viewLine];
				refPlaneCoeff[counterRPCoeffs] = transformedBlock->mPixel[0][verticalView][viewColumn][viewLine];
				counterRPCoeffs++;
				sumRefPlaneCoeff += coeff;

				// if(coeff > maxRefPlane) {
				// 	maxRefPlane = sacoeffmples;
				// }
				// if(coeff < minRefPlane) {
				// 	minRefPlane = coeff;
				// }
			}
		}
	}
}

void Statistics :: calcSumOtherPlanesCoeff(Block4D *transformedBlock) {
	for(int horizontalView = 1; horizontalView < 13; horizontalView += 1) {
		for(int verticalView = 0; verticalView < 13; verticalView += 1) {
			for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
				for(int viewLine = 0; viewLine < 15; viewLine += 1) {
					int coeff = transformedBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
					otherPlanesCoeff[counterOPCoeffs] = transformedBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
					counterOPCoeffs++;
					sumOtherPlanesCoeff += coeff;

					// if(coeff > maxOtherPlanes) {
					// 	maxOtherPlanes = coeff;
					// }
					// if(coeff < minOtherPlanes) {
					// 	minOtherPlanes = coeff;
					// }
				}
			}
		}
	}
}

double Statistics :: getAverageRefPlaneCoeff() {
	return sumRefPlaneCoeff/refPlaneSamples;
}

double Statistics :: getAverageOtherPlanesCoeff() {
	return sumOtherPlanesCoeff/otherPlanesSamples;
}

double Statistics :: calcStdRefPlaneCoeff() {
	double stdDev = 0;
	double average = (double)sumRefPlaneCoeff/3562650;

	for(int i = 0; i < 3562650; i++) {
		stdDev += pow(refPlaneCoeff[i] - average, 2);
	}
	
	return sqrt(stdDev/3562650);
}

double Statistics :: calcStdOtherPlanesCoeff() {
	double stdDev = 0;
	double average = (double)sumOtherPlanesCoeff/3562650;

	for(int i = 0; i < 3562650; i++) {
		stdDev += pow(otherPlanesCoeff[i] - average, 2);
	}
	
	return sqrt(stdDev/3562650);
}

double Statistics :: getEnergyRefPlaneCoeff() {
	return energyRefPlaneCoeff;
}

double Statistics :: getEnergyOtherPlanesCoeff() {
	return energyOtherPlanesCoeff;
}

void Statistics :: calcSumCoeff(Block4D *transformedBlock) {
	for(int horizontalView = 0; horizontalView < 13; horizontalView += 1) {
		for(int verticalView = 0; verticalView < 13; verticalView += 1) {
			for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
				for(int viewLine = 0; viewLine < 15; viewLine += 1) {
					int coeff = transformedBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
					AllCoeffs[counterOPCoeffs] = transformedBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
					counterAllCoeffs++;
					sumCoeff += coeff;

					// if(coeff > maxOtherPlanes) {
					// 	maxOtherPlanes = coeff;
					// }
					// if(coeff < minOtherPlanes) {
					// 	minOtherPlanes = coeff;
					// }
				}
			}
		}
	}
}

void Statistics :: calcCoeffEnergy(Block4D *transformedBlock) {
	for(int horizontalView = 0; horizontalView < 13; horizontalView += 1) {
		for(int verticalView = 0; verticalView < 13; verticalView += 1) {
			for(int viewColumn = 0; viewColumn < 15; viewColumn += 1) {
				for(int viewLine = 0; viewLine < 15; viewLine += 1) {
					int coeff = transformedBlock->mPixel[horizontalView][verticalView][viewColumn][viewLine];
					
					coeffEnergy += coeff*coeff;

					// if(coeff > maxOtherPlanes) {
					// 	maxOtherPlanes = coeff;
					// }
					// if(coeff < minOtherPlanes) {
					// 	minOtherPlanes = coeff;
					// }
				}
			}
		}
	}
}

double Statistics :: getAverageCoeff() {
	return sumCoeff/(3562650+42751800);
}

double Statistics :: getEnergyCoeff() {
	return coeffEnergy;
}

double Statistics :: calcStdCoeff() {
	double stdDev = 0;
	double average = (double)sumCoeff/3562650;

	for(int i = 0; i < 3562650; i++) {
		stdDev += pow(AllCoeffs[i] - average, 2);
	}
	
	return sqrt(stdDev/3562650);
}
