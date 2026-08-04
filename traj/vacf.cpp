#include "analysis.h"

using namespace analysis;

void computeVACF(Trajectory *TRAJ, atom_style **ATOMS, int *frame_sep, int nSamples);

int main(int argc, char* argv[])
{
	float dt = 5e-3;
	int frameW = 1;
	int frameStart = int(1e4);
	int frameEnd = int(2e4);
	int nAtomTypes = 1;

	int nSamples = 100;
	int *frame_sep = new int[nSamples];
	for(int i = 0; i < nSamples; i++)
		frame_sep[i] = int(i*1);

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/ASH_DRIVE_3/ashwin_md/lane/vacf/tau_1e-2/traj1.cfg");
	sprintf(TRAJ->fpathO, "//media/ashwin/ASH_DRIVE_3/ashwin_md/lane/vacf/tau_1e-2/vacf.dat");
	TRAJ -> openTrajectory();

	if(frameEnd < frameStart)
	{
		printf("Error: Ending frame number cannot be smaller than the starting frame. Exiting...\n");
		exit(-1);
	}
	else
	{
		TRAJ->totalFrames = frameEnd - frameStart + 1;
		
		atom_style **ATOMS = new atom_style*[TRAJ->totalFrames];
		System *BOX = new System(TRAJ->Lx, TRAJ->Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);
		for(int i = 0; i < TRAJ->totalFrames; i++)
			ATOMS[i] = new atom_style[TRAJ->nAtoms];

		TRAJ -> loadTrajectory(ATOMS, BOX, frameStart, frameEnd, false);
		computeVACF(TRAJ, ATOMS, frame_sep, nSamples);

		delete[] ATOMS;
	}

	TRAJ -> closeTrajectory(true, true);

}

void computeVACF(Trajectory *TRAJ, atom_style **ATOMS, int *frame_sep, int nSamples)
{
	float *correlations = new float[nSamples];
	long *counts = new long[nSamples];
	for(int i = 0; i < nSamples; i++)
	{
		correlations[i] = 0.0;
		counts[i] = 0;
	}

	for(int currFrame = 0; currFrame < 1; currFrame++)
	{
		for(int ctr = 0; ctr < nSamples; ctr++)
		{
			int nextFrame = currFrame + frame_sep[ctr];
			if(nextFrame < TRAJ->totalFrames)
			{
				for(int i = 0; i < TRAJ->nAtoms; i++)
				{
					float product = ATOMS[nextFrame][i].vx*ATOMS[currFrame][i].vx + ATOMS[nextFrame][i].vy*ATOMS[currFrame][i].vy + ATOMS[nextFrame][i].vz*ATOMS[currFrame][i].vz;
					
					correlations[ctr] += product;
					counts[ctr]++;
				}
			}
		}
	}

	for(int ctr = 0; ctr < nSamples; ctr++)
		correlations[ctr] /= counts[ctr];

	for(int ctr = 0; ctr < nSamples; ctr++)
	{
		correlations[ctr] /= correlations[0];
		printf("Frame_sep = %d, counts = %ld, vacf = %f\n", frame_sep[ctr], counts[ctr], correlations[ctr]);
	}

}