#include "analysis.h"

using namespace analysis;

int main(int argc, char* argv[])
{
	float dt = 5e-4;
	int frameW = int(1e5);
	long startStep1 = long(1e7);
	long endStep1 = long(997e7);
	long startStep2 = long(0);
	long endStep2 = long(4e7);

	Trajectory *TRAJ = new Trajectory(dt , frameW, "xyz");
	sprintf(TRAJ->fpathI, "//media/ashwin/One Touch/ashwin_md/lane/May2025/lmp/Data26/traj2.xyz");
	sprintf(TRAJ->fpathO, "//media/ashwin/One Touch/ashwin_md/lane/May2025/lmp/Data26/traj2_merged.xyz");

	TRAJ -> openTrajectory();

	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];

	while( !feof(TRAJ->fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if(TRAJ->step == startStep1)
			TRAJ -> writeThisFrame(ATOMS, -1);
		else
			TRAJ -> writeThisFrame(ATOMS);

		if(TRAJ->step == endStep1) 
		{
			printf("\n\nSource: %s\nDestination: %s\nTotal frames extracted: %d\n", TRAJ->fpathI, TRAJ->fpathO, TRAJ->frame_nr);
			break;
		}
	}

	TRAJ -> closeTrajectory();

	sprintf(TRAJ->fpathI, "//media/ashwin/One Touch/ashwin_md/lane/May2025/lmp/Data26/traj3.xyz");

	TRAJ -> openTrajectory();

	while( !feof(TRAJ->fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if(TRAJ->step >= startStep2 and TRAJ->step <= endStep2)
			TRAJ -> writeThisFrame(ATOMS, endStep1);

		if(TRAJ->step == endStep2)
		{
			printf("\n\nSource: %s\nDestination: %s\nTotal frames extracted: %d\n", TRAJ->fpathI, TRAJ->fpathO, TRAJ->frame_nr);
			break;
		}
	}

	TRAJ -> closeTrajectory();
	fclose(TRAJ->fileO);
}