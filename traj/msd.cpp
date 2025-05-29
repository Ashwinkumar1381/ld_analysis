/*
	Code to compute the Mean-Squared-Displacement of a group of atoms
	by averaging over each particle's trajectory and over all particles
	in the group.

	Author			: Ashwin Kumar (CH23S006)
	Date created 	: 15.03.25
	Last modified	: 18.03.25
*/

#include "analysis.h"

using namespace analysis;

float* computeMeanSquaredDisplacement(Trajectory *TRAJ, System *BOX, atom_style **ATOMS, const int delFrames[], int N);

int main(int argc, char *argv[])
{
	float startTime = 0.0, endTime = 10000.0, dt = 5e-5;
	int frameW = 200000;						// (width of each frame in terms of no. of steps)
	float Lx = 100.0, Ly = 100.0;

	int nSample = 10;
	int delFrames[nSample] = {1, 2, 5, 10, 20, 50, 100, 200, 500, 1000};

	Trajectory *TRAJ = new Trajectory(dt, frameW);
	sprintf(TRAJ -> fpathI, "//media/ashwin/One Touch/ashwin_md/msd_benchmarks/Mar2025/Data278/traj2.xyz");
	sprintf(TRAJ -> fpathO, "//media/ashwin/One Touch/ashwin_md/msd_benchmarks/Mar2025/Data278/msd.dat");
	
	TRAJ -> openTrajectory(true);

	System *BOX = new System(Lx, Ly, TRAJ->nAtoms);

	atom_style **ATOMS = new atom_style* [TRAJ->totalFrames];
	for(int i = 0; i < TRAJ -> totalFrames; i++) 
		ATOMS[i] = new atom_style [TRAJ->nAtoms];

	TRAJ -> loadTrajectory(ATOMS, BOX);
	TRAJ -> closeTrajectory(); 

	float *meanSquaredDisplacement = computeMeanSquaredDisplacement(TRAJ, BOX, ATOMS, delFrames, nSample);

	TRAJ -> write2file(meanSquaredDisplacement, delFrames, nSample);

	delete[] ATOMS; 
	return(0);
}

float* computeMeanSquaredDisplacement(Trajectory *TRAJ, System *BOX, atom_style **ATOMS, const int delFrames[], int N)
{
	printf("\nComputing MSD...\n");

	float *meanSquaredDisplacement = new float[N];
	int *numSamples = new int[N];

	for(int i = 0; i < N; i++)
	{
		meanSquaredDisplacement[i] = 0.0;
		numSamples[i] = 0;
	}

	for(int i = 0; i < TRAJ -> totalFrames; i++)
	{
		int ctr = 0;

		while(ctr < N)
		{
			int nextFrame = i + delFrames[ctr];
			
			if(nextFrame < TRAJ -> totalFrames)
			{
				for(int j = 0; j < TRAJ -> nAtoms; j++)
				{
					float dx = ATOMS[nextFrame][j].rxt1 - ATOMS[i][j].rxt1 + BOX->Lx*(ATOMS[nextFrame][j].jumpx - ATOMS[i][j].jumpx);

					float dy = ATOMS[nextFrame][j].ryt1 - ATOMS[i][j].ryt1 + BOX->Ly*(ATOMS[nextFrame][j].jumpy - ATOMS[i][j].jumpy);

					meanSquaredDisplacement[ctr] += dx*dx + dy*dy;
					numSamples[ctr] += 1;
				}
			}

			ctr++;
		}
	}

	for(int i = 0; i < N; i++)
		meanSquaredDisplacement[i] /= numSamples[i];

	return(meanSquaredDisplacement);
}

void analysis::Trajectory::write2file(float *meanSquaredDisplacement, const int delFrames[], int nSamples)
{
	remove(fpathO);

	fileO = fopen(fpathO, "w");
	if(fileO == NULL)
	{
		printf("Error. Cannot create new file %s. Exiting...\n", fpathO);
		exit(-1);
	}
	else
		fprintf(fileO, "step msd\n");

	for(int i = 0; i < nSamples; i++)
		fprintf(fileO, "%d %f\n", frameWidth*delFrames[i], meanSquaredDisplacement[i]);	

	fclose(fileO);
}