/*
	Code to compute the Mean-Squared-Displacement of a group of atoms
	by averaging over each particle's trajectory and over all particles
	in the group.

	Author			: Ashwin Kumar
	Date created 	: 15.03.25
	Last modified	: 07.08.25
*/

#include "analysis.h"

using namespace analysis;

void computeMeanSquaredDisplacement(Trajectory *TRAJ, System *BOX, atom_style **ATOMS, int *delFrames, int nSample);

int main(int argc, char *argv[])
{
	// ----------- System Params -----------
	float Lx = 150.0, Ly = 30.0;
	int nAtomTypes = 2;

	// ----------- Trajectory Params -----------
	int frameStart = int(4e4), frameEnd = int(5e4);
	float dt = 5e-5;
	int frameW = int(1e5);

	int nSample = 36;

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ -> fpathI, "//media/ashwin/Expansion/ashwin_md/lane/June_July2025/Pe90/Data37/traj2.cfg");
	sprintf(TRAJ -> fpathO, "//media/ashwin/Expansion/ashwin_md/lane/June_July2025/Pe90/Data37/msd.dat");
	TRAJ -> openTrajectory();

	int *delFrames = new int[nSample];

	int ctr = 0, j = 0;
	for(int i = 0; i < nSample; i++)
	{
		j++;
		delFrames[i] = int(j * pow(10, ctr));
 		
 		if(j % 9 == 0)
		{
			j = 0;
			ctr++;
		}
	}

	TRAJ->totalFrames = frameEnd - frameStart + 1;

	atom_style **ATOMS = new atom_style* [TRAJ->totalFrames];
	System *BOX = new System(Lx, Ly, TRAJ->nAtoms, nAtomTypes);
	for(int i = 0; i < TRAJ->totalFrames; i++) 
		ATOMS[i] = new atom_style [TRAJ->nAtoms];

	TRAJ -> loadTrajectory(ATOMS, BOX, frameStart, frameEnd); 
	computeMeanSquaredDisplacement(TRAJ, BOX, ATOMS, delFrames, nSample);

	TRAJ -> closeTrajectory(true, true);

	delete[] ATOMS; 
	return(0);
}

void computeMeanSquaredDisplacement(Trajectory *TRAJ, System *BOX, atom_style **ATOMS, int *delFrames, int nSample)
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

	for(int j = 0; j < nSample; j++)
		count[j] = 0;

	for(int currFrame = 0; currFrame < TRAJ -> totalFrames; currFrame++)
	{
		int ctr = 0;

		while(ctr < nSample)
		{
			int nextFrame = currFrame + delFrames[ctr];
			
			if(nextFrame < TRAJ->totalFrames)
			{
				for(int j = 0; j < TRAJ->nAtoms; j++)
				{
					float dx = ATOMS[nextFrame][j].rxt1 - ATOMS[currFrame][j].rxt1 + BOX->Lx*(ATOMS[nextFrame][j].jumpx - ATOMS[currFrame][j].jumpx);
					float dy = ATOMS[nextFrame][j].ryt1 - ATOMS[currFrame][j].ryt1 + BOX->Ly*(ATOMS[nextFrame][j].jumpy - ATOMS[currFrame][j].jumpy);
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
			MSD[j][i] /= count[i];
	}

	TRAJ -> createOutputFile("time count MSD_all MSD_xx MSD_yy MSD_xy");
	TRAJ -> write2file(MSD, count, delFrames, nSample);
}

void analysis::Trajectory::write2file(float **MSD, int *count, int *delFrames, int nSample)
{
	for(int i = 0; i < nSample; i++)
		fprintf(fileO, "\n%d %d %g %g %g %g", int(delFrames[i]*frameWidth*timeStep), count[i], MSD[0][i], MSD[1][i], MSD[2][i], MSD[3][i]);
}