// computeInteractions.cpp

#include "analysis.h"

using namespace program;
using namespace analysis;

float computePotentialEnergy(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, WCA_2P *INTERACTIONS);

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
			float pe = computePotentialEnergy(TRAJ, ATOMS, BOX, INTERACTIONS);
			printf("\nPE = %f", pe);
			// TRAJ -> write2file(pe);
		}

		if(TRAJ->frame_nr == frameEnd)
			break;
	}

	TRAJ -> closeTrajectory();
	return(0);
}

float computePotentialEnergy(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, WCA_2P *INTERACTIONS)
{
	float pe_total = 0.0;

	#pragma omp parallel num_threads(6)
	{
		float pe = 0.0;

		#pragma omp for
		for(int i = 0; i < TRAJ->nAtoms; i++)
		{
			for(int j = i + 1; j < TRAJ->nAtoms; j++)
			{
				float dxij = ATOMS[j].rxt1 - ATOMS[i].rxt1;
				float dyij = ATOMS[j].ryt1 - ATOMS[i].ryt1;

				BOX -> checkMinImage(&dxij, &dyij);

				float drij2 = dxij*dxij + dyij*dyij;

				if(drij2 < INTERACTIONS->rcut*INTERACTIONS->rcut)
				{
					float *energies = INTERACTIONS->get_forces(drij2);
					pe += energies[0];
				}
			}
		}

		#pragma omp critical
		pe_total = pe_total + pe;		
	}

	return(pe_total);
}

void analysis::Trajectory::write2file(float pe)
{
	fprintf(fileO, "%ld %f\n", step, pe);
}