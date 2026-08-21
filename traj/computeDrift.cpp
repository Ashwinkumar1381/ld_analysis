

#include "analysis.h"

using namespace std;
using namespace analysis;

int main(int argc, char *argv[])
{
	// System Params
	int nAtomTypes = 2;

	// Trajectory Params
	int frameStart = int(3e4), frameEnd = int(5e4), frameW = int(1e5);
	float dt = 5e-4;

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//run/media/ethaya_lab_c1/ASH_DRIVE_3/ashwin_md/lane/Aug_Nov2025/Fd100/tau_1e-3/traj2.cfg");

	TRAJ -> openTrajectory();
	TRAJ -> totalFrames = frameEnd - frameStart + 1;

	System *BOX = new System(TRAJ->Lx, TRAJ->Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);
	atom_style **ATOMS = new atom_style*[TRAJ->totalFrames];
	for(int i = 0; i < TRAJ->totalFrames; i++)
		ATOMS[i] = new atom_style[TRAJ->nAtoms];

	TRAJ -> loadTrajectory(ATOMS, BOX, frameStart, frameEnd, true);

	float drift_vel = 0.0;
	long nSamples = long(0);

	#pragma omp parallel num_threads(5)
	{
		float partial_sum = 0.0;
		long partial_count = 0;

		#pragma omp for
		for(int currFrame = 0; currFrame < TRAJ->totalFrames; currFrame++)
		{
			if(currFrame%10 == 0)
			{
				// printf("Processing frame %d\n", currFrame);

				for(int nextFrame = currFrame + 1; nextFrame < TRAJ->totalFrames; nextFrame++)
				{
					for(int i = 0; i < TRAJ->nAtoms; i++)
					{
						float dx = ATOMS[nextFrame][i].rxt1 - ATOMS[currFrame][i].rxt1 + BOX->Lx*(ATOMS[nextFrame][i].jumpx - ATOMS[currFrame][i].jumpx);
						float delt = float((nextFrame - currFrame)*frameW*dt);

						partial_sum += float((dx*dx)/(delt*delt));
						partial_count++;
					}
				}
			}
		}

		#pragma omp critical
		{
			drift_vel += partial_sum;
			nSamples += partial_count;
		}
	}

	drift_vel /= nSamples;
	printf("\nDrift velocity averaged over %ld samples, Vd = %f\n", nSamples, drift_vel);

	TRAJ -> closeTrajectory();

	delete[] ATOMS, BOX;

	return(0);
}