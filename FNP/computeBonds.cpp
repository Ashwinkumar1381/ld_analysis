#include "../traj/analysis.h"

using namespace analysis;

int main(int argc, char *argv[])
{
	int nAtomTypes = 2;

	int frameStart = int(0), frameEnd = int(1e3);
	float dt = 5e-2;
	int frameW = 100; 

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/ASH_DRIVE_3/ashwin_md/FNP/polymer_benchmarks/spring/K_200/traj1.cfg");
	sprintf(TRAJ->fpathO, "//media/ashwin/ASH_DRIVE_3/ashwin_md/FNP/polymer_benchmarks/spring/K_200/bonds.dat");

	TRAJ -> openTrajectory();
	TRAJ -> createOutputFile("atom_i atom_j bond_lengths");

	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	System *BOX = new System(TRAJ->Lx, TRAJ->Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);

	while( !feof(TRAJ -> fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
		{
			for(int i = 1; i < TRAJ->nAtoms; i++)
			{
				float dx = ATOMS[i].rxt1 - ATOMS[i-1].rxt1;
				float dy = ATOMS[i].ryt1 - ATOMS[i-1].ryt1;
				float dz = ATOMS[i].rzt1 - ATOMS[i-1].rzt1;

				BOX -> checkMinImage(&dx, &dy, &dz);

				float dr = sqrt(dx*dx + dy*dy + dz*dz);

				fprintf(TRAJ->fileO, "\n%d %d %f", ATOMS[i-1].atom_id, ATOMS[i].atom_id, dr);
			}
		}
	}

	TRAJ -> closeTrajectory();
}