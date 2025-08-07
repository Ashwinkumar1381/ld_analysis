// computeInteractions.cpp

#include "analysis.h"

using namespace program;
using namespace analysis;

int main(int argc, char *argv[])
{
	// System Params
	float Lx = 150.0, Ly = 30.0;
	int nAtomTypes = 2;

	// Trajectory Params
	int frameStart = int(1e3), frameEnd = int(1e3); 
	float dt = 5e-4;
	int frameW = int(1e5);

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/Expansion/ashwin_md/lane/June_July2025/Pe90/Data38/traj2.cfg");
	sprintf(TRAJ->fpathO, "//media/ashwin/Expansion/ashwin_md/lane/June_July2025/Pe90/Data38/energy.dat");
	TRAJ -> openTrajectory();

	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	System *BOX = new System(Lx, Ly, TRAJ->nAtoms, nAtomTypes);
	WCA_2P *INTERACTIONS = new WCA_2P(); 

	TRAJ -> createOutputFile("step PotEng");
	while( !feof(TRAJ->fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
		{
			printf("Processing frame %d\n", TRAJ->frame_nr);
			float pe = computeNonBondedInteractions(ATOMS, BOX, INTERACTIONS);
			// printf("\nPE = %f", pe);
			TRAJ -> write2file(pe);
		}

		if(TRAJ->frame_nr == frameEnd)
			break;
	}

	TRAJ -> closeTrajectory();
	return(0);
}

void analysis::Trajectory::write2file(float pe)
{
	fprintf(fileO, "%ld %f\n", step, pe);
}