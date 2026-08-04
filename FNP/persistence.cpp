#include "../traj/analysis.h"

using namespace analysis;

void computeAngularCorrelations1(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, int *separations, float *correlations, long *counts, int nMonomers);

void computeAngularCorrelations2(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, int *separations, float *correlations, long *counts, int nMonomers);

int main(int argc, char *argv[])
{
	int nAtomTypes = 2;
	int nMonomers = 50;

	int nReplicas = 3;
	int replicas[nReplicas] = {1, 2, 3};
	int frameStart[nReplicas] = {0, 0, 0}; 
	int frameEnd[nReplicas] = {10000, 10000, 10000};

	float dt = 5e-2;
	int frameW = 50;
	char fname[25] = "persistence.dat";

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "/media/ashwin/ASH_DRIVE_3/ashwin_md/FNP/single_chain/K_theta_0/scaling/N_50/a_25");
	sprintf(TRAJ->fpathO, "%s", TRAJ->fpathI);

	int *separations = new int[nMonomers-1];
	float *correlations = new float[nMonomers-1];
	long *counts = new long[nMonomers-1];

	for(int i = 0; i < nMonomers-1; i++)
	{
		separations[i] = i;
		correlations[i] = 0.0;
		counts[i] = 0;
	}

	for(int i = 0; i < nReplicas; i++)
	{
		sprintf(TRAJ->fpathI, "%s/R%d/traj1.cfg", TRAJ->fpathO, replicas[i]);

		TRAJ -> openTrajectory();

		atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
		System *BOX = new System(TRAJ->Lx, TRAJ->Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);

		if(frameEnd[i] < frameStart[i])
		{
			printf("Error: Ending frame number cannot be smaller than the starting frame. Exiting...\n");
			exit(-1);
		}
		while( !feof(TRAJ->fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->frame_nr >= frameStart[i] and TRAJ->frame_nr <= frameEnd[i])
				computeAngularCorrelations2(TRAJ, ATOMS, BOX, separations, correlations, counts, nMonomers);
		}

		delete[] ATOMS;
		delete BOX;

		if(i == nReplicas - 1)
		{
			sprintf(TRAJ->fpathO, "%s/%s", TRAJ->fpathO, fname);
			TRAJ -> createOutputFile("separation counts correlations");

			for(int k = 0; k < nMonomers-1; k++)
			{
				correlations[k] /= counts[k];
				fprintf(TRAJ->fileO, "\n%d %d %f", separations[k], counts[k], correlations[k]);
			}

			TRAJ -> closeTrajectory(true, true);
		}
		else
			TRAJ -> closeTrajectory(true, false);
	}

	delete[] separations, counts, correlations;
	return(0);
}

void computeAngularCorrelations1(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, int *separations, float *correlations, long *counts, int nMonomers)
{
	for(int i = 0; i < TRAJ->nMols; i++)
	{
		int mol_ID = i + 1;

		for(int j = 1; j < TRAJ->nAtoms; j++)
		{
			if(ATOMS[j].mol_id == mol_ID)
			{
				int atom_j = ATOMS[j].atom_id;

				for(int k = 0; k < nMonomers-3; k++)
				{
					int atom_k = atom_j + separations[k];

					if(atom_k < int(nMonomers*mol_ID))
					{
						float delx_j = ATOMS[atom_j].rxt1 - ATOMS[atom_j-2].rxt1 + TRAJ->Lx*(ATOMS[atom_j].jumpx - ATOMS[atom_j-2].jumpx);
						float dely_j = ATOMS[atom_j].ryt1 - ATOMS[atom_j-2].ryt1 + TRAJ->Ly*(ATOMS[atom_j].jumpy - ATOMS[atom_j-2].jumpy);

						float delx_k = ATOMS[atom_k].rxt1 - ATOMS[atom_k-2].rxt1 + TRAJ->Lx*(ATOMS[atom_k].jumpx - ATOMS[atom_k-2].jumpx);
						float dely_k = ATOMS[atom_k].ryt1 - ATOMS[atom_k-2].ryt1 + TRAJ->Ly*(ATOMS[atom_k].jumpy - ATOMS[atom_k-2].jumpy);

						float mod_delr_j = sqrt(delx_j*delx_j + dely_j*dely_j);
						float mod_delr_k = sqrt(delx_k*delx_k + dely_k*dely_k);

						float product = (delx_j*delx_k + dely_j*dely_k)/(mod_delr_j*mod_delr_k);

						correlations[k] += product;
						counts[k] += 1;
					}
				}
			}
		}
	}
}

void computeAngularCorrelations2(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, int *separations, float *correlations, long *counts, int nMonomers)
{
	for(int j = 0; j < nMonomers; j++)
	{
		int atom_j = ATOMS[j].atom_id;

		for(int k = 0; k < nMonomers-1; k++)
		{
			int atom_k = atom_j + separations[k];

			if(atom_k < nMonomers)
			{
				float delx_j = ATOMS[atom_j].rxt1 - ATOMS[atom_j-1].rxt1 + TRAJ->Lx*(ATOMS[atom_j].jumpx - ATOMS[atom_j-1].jumpx);
				float dely_j = ATOMS[atom_j].ryt1 - ATOMS[atom_j-1].ryt1 + TRAJ->Ly*(ATOMS[atom_j].jumpy - ATOMS[atom_j-1].jumpy);

				float delx_k = ATOMS[atom_k].rxt1 - ATOMS[atom_k-1].rxt1 + TRAJ->Lx*(ATOMS[atom_k].jumpx - ATOMS[atom_k-1].jumpx);
				float dely_k = ATOMS[atom_k].ryt1 - ATOMS[atom_k-1].ryt1 + TRAJ->Ly*(ATOMS[atom_k].jumpy - ATOMS[atom_k-1].jumpy);
				
				float mod_delr_j = sqrt(delx_j*delx_j + dely_j*dely_j);
				float mod_delr_k = sqrt(delx_k*delx_k + dely_k*dely_k);

				float product = (delx_j*delx_k + dely_j*dely_k)/(mod_delr_j*mod_delr_k);

				correlations[k] += product;
				counts[k] += 1;
			}
		}
	}	
}