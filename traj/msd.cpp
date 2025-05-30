/*
	Code to compute the Mean-Squared-Displacement of a group of atoms
	by averaging over each particle's trajectory and over all particles
	in the group.

	Author			: Ashwin Kumar (CH23S006)
	Date created 	: 15.03.25
	Last modified	: 28.05.25
*/

#include "analysis.h"

using namespace analysis;

float** computeMeanSquaredDisplacement(Trajectory *TRAJ, System *BOX, atom_style **ATOMS, const int delFrames[], int N);

int main(int argc, char *argv[])
{
	// ----------- System Params -----------
	float Lx = 150.0, Ly = 30.0;

	// ----------- Trajectory Params -----------
	long eq_steps = long(1e7);
	long startStep = long(0e9), endStep = long(2e9);
	float dt = 5e-4;
	int frameW = int(1e5);
	int sep = 1;

	int nSample = 19;
	int delFrames[nSample] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};

	Trajectory *TRAJ = new Trajectory(dt, frameW);
	sprintf(TRAJ -> fpathI, "//media/ashwin/One Touch/ashwin_md/lane/Apr2025/lmp/Data26/traj3.xyz");
	sprintf(TRAJ -> fpathO, "//media/ashwin/One Touch/ashwin_md/lane/Apr2025/lmp/Data26/msd.dat");
	TRAJ -> openTrajectory();

	for(int i = 0; i < nSample; i++) delFrames[i] *= sep;

	int frameStart = int(startStep/frameW), frameEnd = int(endStep/frameW);
	TRAJ->totalFrames = frameEnd - frameStart + 1;

	atom_style **ATOMS = new atom_style* [TRAJ->totalFrames];
	System *BOX = new System(Lx, Ly, TRAJ->nAtoms);
	for(int i = 0; i < TRAJ->totalFrames; i++) 
		ATOMS[i] = new atom_style [TRAJ->nAtoms];

	TRAJ -> loadTrajectory(ATOMS, BOX, frameStart, frameEnd); 
	float **meanSquaredDisplacement = computeMeanSquaredDisplacement(TRAJ, BOX, ATOMS, delFrames, nSample);

	TRAJ -> write2file(meanSquaredDisplacement, delFrames, nSample);
	TRAJ -> closeTrajectory();

	delete[] ATOMS; 
	return(0);
}

float** computeMeanSquaredDisplacement(Trajectory *TRAJ, System *BOX, atom_style **ATOMS, const int delFrames[], int N)
{
	printf("\nComputing MSD...\n");

	float ***MSD = new float**[2];
	int **numSamples = new int*[2];

	for(int i = 0; i < 2; i++)
	{
		MSD[i] = new float*[5];
		for(int j = 0; j < 5; j++)
			MSD[i][j] = new float[N];

		numSamples[i] = new int[N];
	}

	for(int i = 0; i < TRAJ -> totalFrames; i++)
	{
		int ctr = 0;

		while(ctr < N)
		{
			int nextFrame = i + delFrames[ctr];
			
			if(nextFrame < TRAJ->totalFrames)
			{
				for(int j = 0; j < TRAJ->nAtoms; j++)
				{
					float dx = ATOMS[nextFrame][j].rxt1 - ATOMS[i][j].rxt1 + BOX->Lx*(ATOMS[nextFrame][j].jumpx - ATOMS[i][j].jumpx);
					float dy = ATOMS[nextFrame][j].ryt1 - ATOMS[i][j].ryt1 + BOX->Ly*(ATOMS[nextFrame][j].jumpy - ATOMS[i][j].jumpy);
					float dx2 = dx*dx;
					float dxdy = dx*dy;
					float dy2 = dy*dy;

					int pid;
					if(ATOMS[nextFrame][j].id == 'N') pid = 0;
					else if(ATOMS[nextFrame][j].id == 'O') pid = 1;

					MSD[pid][0][ctr] += dx2 + dy2;
					MSD[pid][1][ctr] += dx2;
					MSD[pid][2][ctr] += dxdy;
					MSD[pid][3][ctr] += dxdy;
					MSD[pid][4][ctr] += dy2;
					numSamples[pid][ctr] += 1;

					// meanSquaredDisplacement[0][ctr] += dx2 + dy2;
					// meanSquaredDisplacement[1][ctr] += dx2;
					// meanSquaredDisplacement[2][ctr] += dxdy;
					// meanSquaredDisplacement[3][ctr] += dxdy;
					// meanSquaredDisplacement[4][ctr] += dy2;
					// numSamples[ctr] += 1;
				}
			}

			ctr++;
		}
	}

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			MSD[0][j][i] /= numSamples[0][i];
			MSD[1][j][i] /= numSamples[1][i];
		}

		printf("%d %d %d %f %f %f %f\n", i, delFrames[i], numSamples[0][i], MSD[0][0][i], MSD[0][1][i], MSD[0][4][i], MSD[0][1][i] + MSD[0][4][i]);
		printf("%d %d %d %f %f %f %f\n\n", i, delFrames[i], numSamples[1][i], MSD[1][0][i], MSD[1][1][i], MSD[1][4][i], MSD[1][1][i] + MSD[1][4][i]);
	}



	// for(int i = 0; i < N; i++)
	// {
	// 	printf("%d %d %f %f %f %f %f %f\n", i, numSamples[i], meanSquaredDisplacement[0][i]/numSamples[i], meanSquaredDisplacement[1][i]/numSamples[i], meanSquaredDisplacement[2][i]/numSamples[i], meanSquaredDisplacement[3][i]/numSamples[i], meanSquaredDisplacement[4][i]/numSamples[i], (meanSquaredDisplacement[0][i] - meanSquaredDisplacement[1][i] - meanSquaredDisplacement[4][i])/numSamples[i]);
	// 	for(int j = 0; j < 5; j++) meanSquaredDisplacement[j][i] /= numSamples[i];
	// }

	return(nullptr);
}

void analysis::Trajectory::write2file(float **MSD, const int delFrames[], int nSample)
{
	remove(fpathO);

	fileO = fopen(fpathO, "w");
	if(fileO == NULL)
	{
		printf("Error. Cannot create new file %s. Exiting...\n", fpathO);
		exit(-1);
	}
	else
		fprintf(fileO, "time msd_all msd_xx msd_xy msd_yx msd_yy\n");

	for(int i = 0; i < nSample; i++)
		fprintf(fileO, "%f %f %f %f %f %f\n", frameWidth*delFrames[i]*timeStep, MSD[0][i], MSD[1][i], MSD[2][i], MSD[3][i], MSD[4][i]);	

	fclose(fileO);
}