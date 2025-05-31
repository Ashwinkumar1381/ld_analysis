/*
	Code to compute the Mean-Squared-Displacement of a group of atoms
	by averaging over each particle's trajectory and over all particles
	in the group.

	Author			: Ashwin Kumar
	Date created 	: 15.03.25
	Last modified	: 31.05.25
*/

#include "analysis.h"

using namespace analysis;

float** computeMeanSquaredDisplacement(Trajectory *TRAJ, System *BOX, atom_style **ATOMS, int *delFrames, int N);

int main(int argc, char *argv[])
{
	// ----------- System Params -----------
	float Lx = 150.0, Ly = 30.0;
	int nAtomTypes = 2;

	// ----------- Trajectory Params -----------
	long eq_steps = long(1e7);
	long startStep = long(0e9), endStep = long(2e9);
	float dt = 5e-4;
	int frameW = int(1e5);
	int sep = 1;

	int nSample = 100;
	int *delFrames = new int[nSample];
	for(int i = 0; i < nSample; i++)
		delFrames[i] = int(i + 1);

	Trajectory *TRAJ = new Trajectory(dt, frameW);
	sprintf(TRAJ -> fpathI, "//media/ashwin/One Touch/ashwin_md/lane/Apr2025/lmp/Data26/traj3.xyz");
	sprintf(TRAJ -> fpathO, "//media/ashwin/One Touch/ashwin_md/lane/Apr2025/lmp/Data26/msd.dat");
	TRAJ -> openTrajectory();
	TRAJ -> nAtomTypes = nAtomTypes;

	for(int i = 0; i < nSample; i++) 
		delFrames[i] = int(delFrames[i] * sep);

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

float** computeMeanSquaredDisplacement(Trajectory *TRAJ, System *BOX, atom_style **ATOMS, int *delFrames, int N)
{
	printf("\nComputing MSD...\n");

	float ***MSD = new float**[TRAJ->nAtomTypes];
	int **numSamples = new int*[TRAJ->nAtomTypes];

	for(int i = 0; i < TRAJ->nAtomTypes; i++)
	{
		MSD[i] = new float*[5];
		for(int j = 0; j < 5; j++)
		{
			MSD[i][j] = new float[N];
			for(int k = 0; k < N; k++)
				MSD[i][j][k] = 0.0;
		}

		numSamples[i] = new int[N];
		for(int k = 0; k < N; k++)
			numSamples[i][k] = 0;
	}

	float **meanSquaredDisplacement = new float*[5];
	for(int i = 0; i < 5; i++)
	{
		meanSquaredDisplacement[i] = new float[N];
		for(int j = 0; j < N; j++)
			meanSquaredDisplacement[i][j] = 0.0;
	}

	for(int currFrame = 0; currFrame < TRAJ -> totalFrames; currFrame++)
	{
		int ctr = 0;

		while(ctr < N)
		{
			int nextFrame = currFrame + delFrames[ctr];
			
			if(nextFrame < TRAJ->totalFrames)
			{
				for(int j = 0; j < TRAJ->nAtoms; j++)
				{
					float dx = ATOMS[nextFrame][j].rxt1 - ATOMS[currFrame][j].rxt1 + BOX->Lx*(ATOMS[nextFrame][j].jumpx - ATOMS[currFrame][j].jumpx);
					float dy = ATOMS[nextFrame][j].ryt1 - ATOMS[currFrame][j].ryt1 + BOX->Ly*(ATOMS[nextFrame][j].jumpy - ATOMS[currFrame][j].jumpy);
					float dx2 = dx*dx;
					float dxdy = abs(dx*dy);
					float dy2 = dy*dy;

					int pid = 0;
					if(ATOMS[nextFrame][j].id == 'N') pid = 0;
					else if(ATOMS[nextFrame][j].id == 'O') pid = 1;

					MSD[pid][0][ctr] += (dx2 + dy2);
					MSD[pid][1][ctr] += dx2;
					MSD[pid][2][ctr] += dxdy;
					MSD[pid][3][ctr] += dxdy;
					MSD[pid][4][ctr] += dy2;
					numSamples[pid][ctr] += 1;
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
			meanSquaredDisplacement[j][i] = 0.5*(MSD[0][j][i] + MSD[1][j][i]);
		}

		// printf("A: %d %d %f %f %f %f %f %d\n", i, delFrames[i], MSD[0][0][i], MSD[0][1][i], MSD[0][2][i], MSD[0][3][i], MSD[0][4][i], numSamples[0][i]);
		// printf("B: %d %d %f %f %f %f %f %d\n\n", i, delFrames[i], MSD[1][0][i], MSD[1][1][i], MSD[1][2][i], MSD[1][3][i], MSD[1][4][i], numSamples[1][i]);
	}

	return(meanSquaredDisplacement);
}

void analysis::Trajectory::write2file(float **MSD, int *delFrames, int nSample)
{
	remove(fpathO);

	fileO = fopen(fpathO, "w");
	if(fileO == NULL)
	{
		printf("Error. Cannot create new file %s. Exiting...\n", fpathO);
		exit(-1);
	}
	else
		fprintf(fileO, "time MSD_all MSD_xx MSD_xy MSD_yx MSD_yy\n");

	for(int i = 0; i < nSample; i++)
		fprintf(fileO, "%f %f %f %f %f %f\n", frameWidth*delFrames[i]*timeStep, MSD[0][i], MSD[1][i], MSD[2][i], MSD[3][i], MSD[4][i]);	

	fclose(fileO);
}