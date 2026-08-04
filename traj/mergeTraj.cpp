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
	
	int nFiles = 12;
	long startStep[nFiles] = {long(1e5), long(1e5), long(1e5), long(1e5), long(1e5), long(1e5), long(1e5), long(1e5), long(1e5), long(1e5), long(1e5), long(1e5)};

	long endStep[nFiles]   = {long(1e7), long(1e7), long(1e7), long(1e7), long(1e7), long(1e7), long(1e7), long(1e7), long(1e7), long(1e7), long(1e7), long(1e7)};

	string filenames[nFiles] = {"50", "70", "80", "100", "150", "200", "250", "300", "350", "400", "450", "500"};

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathO, "//media/ashwin/ASH_DRIVE_3/ashwin_md/lane/Aug_Nov2025/traj_eq.%s", TRAJ->format);
	TRAJ -> createOutputFile();

	atom_style *ATOMS;
	System *BOX;

	for(int i = 0; i < nFiles; i++)
	{
		sprintf(TRAJ->fpathI, "//media/ashwin/ASH_DRIVE_3/ashwin_md/lane/Aug_Nov2025/Fd%s/tau_1e-3/traj1.%s", (filenames[i]).c_str(), TRAJ->format);
		TRAJ -> openTrajectory();

		long addStep = 0;

		if(i == 0)
		{
			ATOMS = new atom_style[TRAJ->nAtoms];
			BOX = new System(TRAJ->Lx, TRAJ->Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);

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

				printf("\nSource: %s\nDestination: %s\nFrames extracted this step: %d\nTotal frames extracted: %d\n", TRAJ->fpathI, TRAJ->fpathO, TRAJ->frame_nr, TRAJ->totalFrames);
				break;
			}
		}

		if(i < nFiles - 1)
			TRAJ -> closeTrajectory(true, false);
		else
			TRAJ -> closeTrajectory(true, true);
	}
}