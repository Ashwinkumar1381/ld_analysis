/*
	Code to compute displacement distributions along x- and y- axis
	by averaging over all particles and over time  

	Author  	  : Ashwin Kumar
	Date created  : 21.05.25
	Last modified : 25.05.25
*/

#include "analysis.h"

using namespace analysis;

float **buildHistogram(int nBins[]);
void computeDisplacementDistribution(atom_style **ATOMS, Trajectory *TRAJ, System *BOX, float **Hist, float Rcut[], float binW[], int sep);

int main(int argc, char *argv[])
{
	// ----------- System params -----------
	float Lx = 150.0, Ly = 30.0;
	float Rcut[2] = {1e3, 1e3};
	float binW[2] = {0.1, 0.1}; 
	int nBins[2] = {int(Rcut[0]/binW[0]), int(Rcut[1]/binW[1])};

	// ----------- Trajectory params -----------
	long eq_steps = long(1e7);
	long startStep = long(8e8), endStep = long(10e8);
	float dt = 5e-4;
	int frameW = int(1e4);
	int sep = 10;

	Trajectory *TRAJ = new Trajectory(dt, frameW);
	sprintf(TRAJ->fpathI, "//media/ashwin/One Touch/ashwin_md/lane/Apr2025/lmp/Data44/traj2.xyz");
	sprintf(TRAJ->fpathO, "//media/ashwin/One Touch/ashwin_md/lane/Apr2025/lmp/Data44/dispDist.dat");
	TRAJ -> openTrajectory();

	int frameStart = int(startStep/frameW), frameEnd = int(endStep/frameW);
	TRAJ->totalFrames = frameEnd - frameStart + 1;

	atom_style **ATOMS = new atom_style*[TRAJ->totalFrames];
	for(int i = 0; i < TRAJ->totalFrames; i++)
		ATOMS[i] = new atom_style[TRAJ->nAtoms];

	System *BOX = new System(Lx, Ly, TRAJ->nAtoms);
	auto Hist = buildHistogram(nBins);

	TRAJ -> loadTrajectory(ATOMS, BOX, frameStart, frameEnd);

	computeDisplacementDistribution(ATOMS, TRAJ, BOX, Hist, Rcut, binW, sep);

	TRAJ -> write2file(Hist, nBins);
	TRAJ -> closeTrajectory();

	return(0);
}

float **buildHistogram(int nBins[])
{
	float **Hist = new float*[2];
	for(int i = 0; i < 2; i++)
	{
		Hist[i] = new float[nBins[i]];
		for(int j = 0; j < nBins[i]; j++)
			Hist[i][j] = 0.0;
	}

	return(Hist);
}

void computeDisplacementDistribution(atom_style **ATOMS, Trajectory *TRAJ, System *BOX, float **Hist, float Rcut[], float binW[], int sep)
{
	printf("\nComputing Displacement distributions...\n");

	int nBins[2] = {int(Rcut[0]/binW[0]), int(Rcut[1]/binW[1])};

	for(int i = sep; i < TRAJ->totalFrames; i++)
	{
		for(int j = 0; j < TRAJ->nAtoms; j++)
		{
			float dx = abs(ATOMS[i][j].rxt1 - ATOMS[i - sep][j].rxt1 + BOX->Lx * (ATOMS[i][j].jumpx - ATOMS[i - sep][j].jumpx));
			float dy = abs(ATOMS[i][j].ryt1 - ATOMS[i - sep][j].ryt1 + BOX->Ly * (ATOMS[i][j].jumpy - ATOMS[i - sep][j].jumpy));

			if(dx <= Rcut[0])
			{
				int bin_x = int(dx / binW[0]);
				Hist[0][bin_x]++; 
			}
			if(dy <= Rcut[1])
			{
				int bin_y = int(dy / binW[1]);
				Hist[1][bin_y]++;
			}
		}
	}

	for(int i = 0; i < 2; i++)
	{
		int total_count = 0;
		for(int j = 0; j < nBins[i]; j++) total_count += Hist[i][j];
		for(int j = 0; j < nBins[i]; j++) Hist[i][j] /= total_count;
	}
}

void analysis::Trajectory::write2file(float **Hist, int nBins[])
{
	remove(fpathO);

	fileO = fopen(fpathO, "w");
	if(fileO == NULL)
	{
		printf("Cannot create file %s. Exiting...\n", fpathO);
		exit(-1);
	}
	else
	{
		fprintf(fileO, "bin P(x)\n");
		for(int i = 0; i < nBins[0]; i++) fprintf(fileO, "%d %f\n", i + 1, Hist[0][i]);

		fprintf(fileO, "bin P(y)\n");
		for(int i = 0; i < nBins[1]; i++) fprintf(fileO, "%d %f\n", i + 1, Hist[1][i]);
	
		fclose(fileO);
	}
}