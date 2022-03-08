#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <iostream>
#include <string>
#include <Windows.h>
#define SDL_MAIN_HANDLED
#include <SDL.h>
#include <iostream>
#include <fstream>
using namespace std;

//Number of charges on disc
int N = 38;

//Charge Location Array
double Charges[2][200];

//An array to store charge locations temporarily, whilst forces are being calculated, before transfering to the Charges array
double ChargesTemp[2][200];

//Array to store forces charges are subject to
double Forces[2][200];

//Velocities of charges
double Velocities[2][200];

double PI = atan(1) * 4;

//Temperature of simulated annealing process
double T = 0;

//Frame tracker
int Frame = 0;

//Individual potential energies that charges are subject to
double ChargeEnergies[200];

//Total energy of system of charges
double TotalEnergy = 0;

//The charge locations of the configuration of charges with the lowest energies for N charges
double LowestConfig[2][200];

//Variable to track the total potential energy of the current most optimal configuration
double LowestEnergy = 100000;

//A constant used to artificially increase and decrease distances between charges in order to quickly lock them into their current local minima
int Const = 200;

//Heuristic values for number of charges on the outer edge of the disc. This value is imported using the results from https://github.com/xplosiveninja/RingConfinedChargesOnThinDiscOptimisation
int OuterRing = 0;

//Total energy of system per each frame
float Vals[20000];

//Number of Rings to be drawn as guide
int Max = 6;

//Force Calculation Function
void NetForce(int Charge) {
	for (int i = Charge; i < N; i++) {
		if (i != Charge) {
			double DistanceX = Const * (Charges[0][i] - Charges[0][Charge]);
			double DistanceY = Const * (Charges[1][i] - Charges[1][Charge]);
			double Vector[2] = { -DistanceX, -DistanceY };
			double Magnitude = sqrt((DistanceX * DistanceX) + (DistanceY * DistanceY));
			double MagCubed = pow(Magnitude, 3);
			Forces[0][Charge] += Vector[0] / MagCubed;
			Forces[1][Charge] += Vector[1] / MagCubed;
			Forces[0][i] -= Vector[0] / MagCubed;
			Forces[1][i] -= Vector[1] / MagCubed;
		}
		double RandVel = static_cast <double> ((rand()) * T) / (static_cast <double> (RAND_MAX));
		double RandAngle = 2 * PI * static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		Forces[0][Charge] += RandVel * cos(RandAngle);
		Forces[1][Charge] += RandVel * sin(RandAngle);
	}
}

//Random Initialisation Function
void Initialise() {
	for (int i = 0; i < N; i++) {
		float Dist = static_cast <float> (rand())/ static_cast <float> (RAND_MAX);
		float Angle = 2 * PI * static_cast <float> (rand())/ static_cast <float> (RAND_MAX);
		Charges[0][i] = Dist * cos(Angle);
		Charges[1][i] = Dist * sin(Angle);
	}
}

//Function to import configurations from quick search
void Initialise2() {
	string myText;
	ifstream MyReadFile("Results/" + std::to_string(N) + ".txt");
	size_t pos = 0;
	std::string token;
	std::string delimiter = " ";
	int i = 0;
	int j = 0;
	double configs[10][3];
	int count = 0;
	while (getline(MyReadFile, myText)) {
		std::cout << myText + " \n";
		while ((pos = myText.find(delimiter)) != std::string::npos) {
			token = myText.substr(0, pos);
			configs[i][j] = stoi(token);
			myText.erase(0, pos + delimiter.length());
			j += 1;
		}
		i += 1;
		j = 0;
	}
	MyReadFile.close();
	for (int i = 0; i < 10; i++) {
		for (int Nc = 0; Nc < configs[i][0]; Nc++) {
			Charges[0][count] = (configs[i][1] / 200) * cos((2 * PI * Nc / configs[i][0]) + (2 * PI * configs[i][2] / 2000));
			Charges[1][count] = (configs[i][1] / 200) * sin((2 * PI * Nc / configs[i][0]) + (2 * PI * configs[i][2] / 2000));
			ChargesTemp[0][count] = (configs[i][1] / 200) * cos((2 * PI * Nc / configs[i][0]) + (2 * PI * configs[i][2] / 2000));
			ChargesTemp[1][count] = (configs[i][1] / 200) * sin((2 * PI * Nc / configs[i][0]) + (2 * PI * configs[i][2] / 2000));
			count += 1;
		}
	}
	OuterRing = configs[0][0];
}

//Draws a circle
void Circle(SDL_Renderer* renderer) {
	for (int r = 0; r < Max; r++) {
		for (float i = 0; i < 2 * PI; i += 0.01) {
			SDL_RenderDrawLine(renderer, (int)(200 - (r*(200/ Max))) * cos(i) + 300, (int)(200 - (r * (200 / Max)))* sin(i) + 300, (int)(200 - (r * (200 / Max)))* cos(i + 0.01) + 300, (int)(200 - (r * (200 / Max)))* sin(i + 0.01) + 300);
			SDL_RenderDrawLine(renderer, (int)(200 - (r * (200 / Max))) * cos(i) + 300, (int)(200 - (r * (200 / Max))) * sin(i) + 900, (int)(200 - (r * (200 / Max))) * cos(i + 0.01) + 300, (int)(200 - (r * (200 / Max))) * sin(i + 0.01) + 900);
		}
	}
}

//Function to calculate total energy of the system
void EnergyCalc(int Charge) {
	double ChargeEnergy = 0;
	for (int i = 0; i < N; i++) {
		if (i != Charge) {
			double DistanceX = (Charges[0][i] - Charges[0][Charge]);
			double DistanceY = (Charges[1][i] - Charges[1][Charge]);
			double Magnitude = sqrt((DistanceX * DistanceX) + (DistanceY * DistanceY));
			ChargeEnergies[Charge] += 1 / Magnitude;
			ChargeEnergies[i] += 1 / Magnitude;
		}
	}
}

void main() {
	SDL_Init(SDL_INIT_EVERYTHING);
	SDL_Window *window = SDL_CreateWindow("title", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 600, 1200, SDL_WINDOW_SHOWN);
	SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, 0);
	SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
	SDL_RenderClear(renderer);
	srand(time(0));
	//Will import and analyse for charges 17 to 100, 17 is chosen as initial value as it is the first configuration with an optimal configuration that has more than 1 charge in the central area
	for (int N_Values = 36; N_Values < 101; N_Values++) {
		N = N_Values;
		LowestEnergy = 100000;
		memset(Charges, 0, sizeof(Charges));
		memset(ChargesTemp, 0, sizeof(ChargesTemp));
		memset(Forces, 0, sizeof(Forces));
		memset(Velocities, 0, sizeof(Velocities));
		memset(ChargeEnergies, 0, sizeof(ChargeEnergies));
		memset(LowestConfig, 0, sizeof(LowestConfig));
		Initialise2();
		//Each N charges runs for 20000 frames
		while (Frame < 20000) {
			int Index = 0;
			Frame++;
			//Set T to 0 if wanting to remove simulated annealing
			T = 0.005 * pow(sin((2 * PI * Frame) / 5000), 3);
			//T = 0;
			TotalEnergy = 0;
			if (abs(T) < 0.0001) {
				Const = 40;
			}
			else {
				Const = 200;
			}
			if ((Frame % 30000) == 0) {
				for (int i = 0; i < N; i++) {
					Charges[0][i] = LowestConfig[0][i];
					Charges[1][i] = LowestConfig[1][i];
				}
			}
			//Start C from Outer ring if wanting to perform heuristic analysis, else 0 to simulate all charges
			//for (int C = 0; C < N; C++) {
			for (int C = OuterRing; C < N; C++) {
				Charges[0][C] = ChargesTemp[0][C];
				Charges[1][C] = ChargesTemp[1][C];
				Forces[0][C] = 0;
				Forces[1][C] = 0;
				Velocities[0][C] /= 2;
				Velocities[1][C] /= 2;
			}
			//Run Force Calculation
			for (int C = 0; C < N; C++) {
				NetForce(C);
			}
			//Apply new forces
			for (int C = 0; C < N; C++) {
				ChargeEnergies[C] = 0;
				Velocities[0][C] += Forces[0][C];
				Velocities[1][C] += Forces[1][C];
				ChargesTemp[0][C] += Velocities[0][C];
				ChargesTemp[1][C] += Velocities[1][C];
				//simulate charge bouncing off of wall
				double Mag = sqrt(((ChargesTemp[0][C]) * (ChargesTemp[0][C])) + ((ChargesTemp[1][C]) * (ChargesTemp[1][C])));
				if (Mag > 1) {
					ChargesTemp[0][C] /= Mag;
					ChargesTemp[1][C] /= Mag;
					double NewVectorMagChange = 2 * ((Velocities[0][C] * -ChargesTemp[0][C]) + (Velocities[1][C] * -ChargesTemp[1][C])) / ((ChargesTemp[0][C] * ChargesTemp[0][C]) + (ChargesTemp[1][C] * ChargesTemp[1][C]));
					Velocities[0][C] -= NewVectorMagChange * (-ChargesTemp[0][C]);
					Velocities[1][C] -= NewVectorMagChange * (-ChargesTemp[1][C]);
					Velocities[0][C] /= 10;
					Velocities[1][C] /= 10;
				}
			}
			//Perform energy calculations and update lowest config
			if (Frame % 1 == 0) {
				for (int C = 0; C < N; C++) {
					EnergyCalc(C);
					TotalEnergy += 0.5 * ChargeEnergies[C];
				}
				if (TotalEnergy < LowestEnergy) {
					LowestEnergy = TotalEnergy;
					for (int i = 0; i < N; i++) {
						LowestConfig[0][i] = Charges[0][i];
						LowestConfig[1][i] = Charges[1][i];
					}
				}
			}
			//Update renderer
			if (Frame % 1 == 0) {
				SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
				SDL_RenderClear(renderer);
				SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
				Circle(renderer);
				for (int i = 0; i < N; i++) {
					SDL_RenderDrawLine(renderer, (int)200 * Charges[0][i] + 295, (int)200 * Charges[1][i] + 300, (int)200 * Charges[0][i] + 305, (int)200 * Charges[1][i] + 300);
					SDL_RenderDrawLine(renderer, (int)200 * Charges[0][i] + 300, (int)200 * Charges[1][i] + 295, (int)200 * Charges[0][i] + 300, (int)200 * Charges[1][i] + 305);
					SDL_RenderDrawLine(renderer, (int)200 * LowestConfig[0][i] + 295, (int)200 * LowestConfig[1][i] + 900, (int)200 * LowestConfig[0][i] + 305, (int)200 * LowestConfig[1][i] + 900);
					SDL_RenderDrawLine(renderer, (int)200 * LowestConfig[0][i] + 300, (int)200 * LowestConfig[1][i] + 895, (int)200 * LowestConfig[0][i] + 300, (int)200 * LowestConfig[1][i] + 905);
				}
				SDL_RenderPresent(renderer);
			}
			Vals[Frame] = TotalEnergy;
		}
		//Save output files
		//ofstream MyFile("SimulatedAnalysis" + std::to_string(N) + "SA.txt");
		ofstream MyFile("HeuristicAnalysis/" + std::to_string(N) + "HA.txt");

		for (int i = 0; i < N; i++) {
			MyFile << std::to_string(LowestConfig[0][i]) + " " + std::to_string(LowestConfig[1][i]) + "\n";
		}
		MyFile.close();
		Frame = 0;
	}
}