#include "analysis.h"

using namespace std;
using namespace analysis;

int main(int argc, char *argv[])
{
	int frameStart = int(0), frameEnd = int(5e3);
	int frameW = int(1);
	float dt = 5e-4;

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//run/media/ethaya_lab_c1/ASH_DRIVE_3/ashwin_md/binary_collisions/tests/traj.cfg");
	sprintf(TRAJ->fpathO, "//run/media/ethaya_lab_c1/ASH_DRIVE_3/ashwin_md/binary_collisions/tests/coords.dat");

	TRAJ -> openTrajectory();
	TRAJ -> createOutputFile("step x0 y0 jumpx0 jumpy0 vx vy x1 y1 jumpx1 jumpy1");

	TRAJ -> totalFrames = int(frameEnd - frameStart + 1);
	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];

	for(int currFrame = 0; currFrame < TRAJ->totalFrames; currFrame++)
	{
		TRAJ -> readThisFrame(ATOMS);
		if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
			fprintf(TRAJ->fileO, "\n%ld %f %f %d %d %f %f %f %f %d %d", TRAJ->step, ATOMS[0].rxt1, ATOMS[0].ryt1, ATOMS[0].jumpx, ATOMS[0].jumpy, ATOMS[0].vx, ATOMS[0].vy, ATOMS[1].rxt1, ATOMS[1].ryt1, ATOMS[1].jumpx, ATOMS[1].jumpy);
	}

	TRAJ -> closeTrajectory(true, true);

	return(0);
}