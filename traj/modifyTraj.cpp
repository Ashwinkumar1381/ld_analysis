/*
	modifyTraj.cpp : Read a current trajectory file and modify certain field values. Then output to a new trajectory file

	Date created  : 01.07.25
	Last modified : 01.07.25
*/

#include "analysis.h"

using namespace analysis;

int main(int argc, char *argv[])
{
	float Lx = 150.0, Ly = 30.0;
	int nAtomTypes = 2;

	float dt = 5e-4;
	int frameW = int(1e3);
	int frameStart = int(0), frameEnd = int(1e5);

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/One Touch/ashwin_md/lane/June2025/Pe30/Data10/traj2.cfg");
	sprintf(TRAJ->fpathO, "//media/ashwin/One Touch/ashwin_md/lane/June2025/Pe30/Data10/traj2_new.cfg");
	TRAJ -> openTrajectory();

	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	System *BOX = new System(Lx, Ly, TRAJ->nAtoms, nAtomTypes);

	TRAJ -> createOutputFile();

	while( !feof(TRAJ->fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <=frameEnd)
			TRAJ -> writeThisFrame(ATOMS, BOX);

		if(TRAJ->frame_nr == frameEnd)
			break;
	}

	TRAJ -> closeTrajectory();
	return(0);
}