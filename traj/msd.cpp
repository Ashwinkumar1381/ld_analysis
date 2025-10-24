/*
	Code to compute the Mean-Squared-Displacement of a group of atoms
	by averaging over each particle's trajectory and over all particles
	in the group.

	Author			: Ashwin Kumar
	Date created 	: 15.03.25
	Last modified	: 30.09.25
*/

#include "analysis.h"

using namespace analysis;

void computeMeanSquaredDisplacement(Trajectory *TRAJ, System *BOX, atom_style **ATOMS, char group[], int delFrames[], int nSample);

int main(int argc, char *argv[])
{
	// ----------- System Params -----------
	float Lx = 50.0, Ly = 50.0;
	int nAtomTypes = 2;

	// ----------- Trajectory Params -----------
	int frameStart = int(0e4), frameEnd = int(2e4) - 1;
	float dt = 5e-4;
	int frameW = int(9000);

	char group[5] = "1";
	int nSample = 12;
	int delFrames[nSample] = {1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000};

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ -> fpathI, "//media/ashwin/Expansion/ashwin_md/psps/tau_1e-3/test2/traj2.cfg");
	sprintf(TRAJ -> fpathO, "//media/ashwin/Expansion/ashwin_md/psps/tau_1e-3/test2/msd.dat");
	TRAJ -> openTrajectory();

	if(frameEnd < frameStart)
	{
		printf("Error: Ending frame number cannot be smaller than the starting frame. Exiting...\n");
		exit(-1);
	}
	else
	{
		TRAJ->totalFrames = frameEnd - frameStart + 1;

		atom_style **ATOMS = new atom_style* [TRAJ->totalFrames];
		System *BOX = new System(Lx, Ly, TRAJ->nAtoms, nAtomTypes);
		for(int i = 0; i < TRAJ->totalFrames; i++) 
			ATOMS[i] = new atom_style [TRAJ->nAtoms];

		TRAJ -> loadTrajectory(ATOMS, BOX, frameStart, frameEnd, true); 
		computeMeanSquaredDisplacement(TRAJ, BOX, ATOMS, group, delFrames, nSample);

		delete[] ATOMS;
	}

	TRAJ -> closeTrajectory(true, true);

	return(0);
}

void computeMeanSquaredDisplacement(Trajectory *TRAJ, System *BOX, atom_style **ATOMS, char group[], int delFrames[], int nSample)
{
	printf("\nComputing MSD...\n");

	float **MSD = new float*[4];
	int *count = new int[nSample];

	for(int i = 0; i < 4; i++)
	{
		MSD[i] = new float[nSample];
		for(int j = 0; j < nSample; j++)
			MSD[i][j] = 0.0;
	}

	for(int i = 0; i < nSample; i++)
		count[i] = 0;

	long **group_lists = TRAJ->sortAtomsByType(ATOMS[0], BOX);

	int loopMax, group_id = 0;

	if(strcmp(group, "all") == 0)
		loopMax = TRAJ->nAtoms;
	else
	{
		sscanf(group, "%d\n", &group_id);
		loopMax = group_lists[0][group_id - 1];
	}

	// printf("%d %d %d\n", group_lists[0][0], group_lists[0][1], loopMax);

	for(int currFrame = 0; currFrame < TRAJ->totalFrames; currFrame++)
	{
		int ctr = 0;

		while(ctr < nSample)
		{
			int nextFrame = currFrame + delFrames[ctr];
			
			if(nextFrame < TRAJ->totalFrames)
			{
				for(int j = 0; j < loopMax; j++)
				{
					int k = j;
					if(group_id != 0)
						k = group_lists[group_id][j] - 1;

					float dx = ATOMS[nextFrame][k].rxt1 - ATOMS[currFrame][k].rxt1 + BOX->Lx*(ATOMS[nextFrame][k].jumpx - ATOMS[currFrame][k].jumpx);
					float dy = ATOMS[nextFrame][k].ryt1 - ATOMS[currFrame][k].ryt1 + BOX->Ly*(ATOMS[nextFrame][k].jumpy - ATOMS[currFrame][k].jumpy);
					float dx2 = dx*dx;
					float dy2 = dy*dy;
					float dxdy = abs(dx)*abs(dy);

					MSD[0][ctr] += (dx2 + dy2);
					MSD[1][ctr] += dx2;
					MSD[2][ctr] += dy2;
					MSD[3][ctr] += dxdy;
					count[ctr] += 1;
				}
			}

			ctr++;
		}
	}

	for(int i = 0; i < nSample; i++)
	{
		for(int j = 0; j < 4; j++)
			MSD[j][i] /= float(count[i]);
	}

	TRAJ -> createOutputFile("time count MSD_all MSD_xx MSD_yy MSD_xy");
	TRAJ -> write2file(MSD, count, delFrames, nSample);
}

void analysis::Trajectory::write2file(float **MSD, int *count, int *delFrames, int nSample)
{
	for(int i = 0; i < nSample; i++)
		fprintf(fileO, "\n%g %d %g %g %g %g", delFrames[i]*frameWidth*timeStep, count[i], MSD[0][i], MSD[1][i], MSD[2][i], MSD[3][i]);

	fprintf(fileO, "\n");
}