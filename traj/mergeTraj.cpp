/*
	Extract specific frames from multiple trajectories and combine into one trajectory

	Last modified : 05.08.25
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
	long startStep1 = long(1e7);
	long endStep1 = long(4.7e9 + 1e7);
	long startStep2 = long(1);
	long endStep2 = long(3e8);

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/Expansion/ashwin_md/lane/Aug2025/Fd100/tau_5e-2/traj2_old.%s", TRAJ->format);
	sprintf(TRAJ->fpathO, "//media/ashwin/Expansion/ashwin_md/lane/Aug2025/Fd100/tau_5e-2/traj2.%s", TRAJ->format);

	TRAJ -> openTrajectory();
	TRAJ -> createOutputFile();

	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	System *BOX = new System(Lx, Ly, TRAJ->nAtoms, nAtomTypes);

	while( !feof(TRAJ->fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if(TRAJ->step >= startStep1 and TRAJ->step <= endStep1)
			TRAJ -> writeThisFrame(ATOMS, BOX, 0);

		if(TRAJ->step == endStep1) 
		{
			printf("\n\nSource: %s\nDestination: %s\nTotal frames extracted: %d\n", TRAJ->fpathI, TRAJ->fpathO, TRAJ->frame_nr);
			break;
		}
	}

	TRAJ -> closeTrajectory(true, false);

	sprintf(TRAJ->fpathI, "//media/ashwin/Expansion/ashwin_md/lane/Aug2025/Fd100/tau_5e-2/traj2_res.%s", TRAJ->format);

	TRAJ -> openTrajectory();

	while( !feof(TRAJ->fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if(TRAJ->step >= startStep2 and TRAJ->step <= endStep2)
			TRAJ -> writeThisFrame(ATOMS, BOX, endStep1);

		if(TRAJ->step == endStep2)
		{
			printf("\n\nSource: %s\nDestination: %s\nTotal frames extracted: %d\n", TRAJ->fpathI, TRAJ->fpathO, TRAJ->frame_nr);
			break;
		}
	}

	TRAJ -> closeTrajectory(true, true);
}