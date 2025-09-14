/*
	Extract specific frames from multiple trajectories and combine into one trajectory

	Last modified : 14.09.25
*/

#include "analysis.h"

using namespace analysis;

int main(int argc, char* argv[])
{
	/* ---------- System Params ---------- */
	float Lx = 150.0, Ly = 30.0;
	int nAtomTypes = 2;

	/* ---------- Trajectory Params ---------- */
	float dt = 5e-4;
	int frameW = int(1e5);
	
	int nFiles = 3;
	long startStep[nFiles] = {1e7, 1e0, 1e0};
	long endStep[nFiles]   = {2.13e9, 2e9, 8.8e8};
	string filenames[nFiles] = {"traj2_old", "traj2_res", "traj3_res"};

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathO, "//media/ashwin/Expansion/ashwin_md/lane/Aug2025/Fd100/tau_5e-2/traj2.%s", TRAJ->format);
	TRAJ -> createOutputFile();

	atom_style *ATOMS;
	System *BOX;

	for(int i = 0; i < nFiles; i++)
	{
		sprintf(TRAJ->fpathI, "//media/ashwin/Expansion/ashwin_md/lane/Aug2025/Fd100/tau_5e-2/%s.%s", (filenames[i]).c_str(), TRAJ->format);
		TRAJ -> openTrajectory();

		int addStep = 0;

		if(i == 0)
		{
			ATOMS = new atom_style[TRAJ->nAtoms];
			BOX = new System(Lx, Ly, TRAJ->nAtoms, nAtomTypes);

			TRAJ->totalFrames = 0;
		}
		else
			addStep = endStep[i - 1];

		while( !feof(TRAJ->fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->step >= startStep[i] and TRAJ->step <= endStep[i])
				TRAJ -> writeThisFrame(ATOMS, BOX, addStep);

			if(TRAJ->step == endStep[i]) 
			{
				endStep[i] += addStep;
				TRAJ->totalFrames += TRAJ->frame_nr;

				printf("\n\nSource: %s\nDestination: %s\nFrames extracted this step: %d\nTotal frames extracted: %d\n", TRAJ->fpathI, TRAJ->fpathO, TRAJ->frame_nr, TRAJ->totalFrames);
				break;
			}
		}

		if(i < nFiles - 1)
			TRAJ -> closeTrajectory(true, false);
		else
			TRAJ -> closeTrajectory(true, true);
	}
}