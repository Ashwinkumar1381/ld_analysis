#include "../traj/analysis.h"

using namespace analysis;

void computeMeanSquaredDistance(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, int *separations, float *msd, long *counts, int nMonomers);

void computeEndDistances(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, int nMonomers);

int main(int argc, char *argv[])
{
	int nAtomTypes = 2;
	int nMonomers = 50;

	float dt = 5e-2;
	int frameW = 50;
	char option[50] = "End_Distance";
	char fname[50] = "end_distances.dat";

	int nReplicas = 3;
	int replicas[nReplicas] = {1, 2, 3};
	int frameStart[nReplicas] = {0, 0, 0}; 
	int frameEnd[nReplicas] = {int(3e4), int(3e4), int(3e4)};

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/ASH_DRIVE_3/ashwin_md/FNP/polymer_benchmarks/angles/potential_WLC/N50/K_100.0/R1");

	if(strcmp(option, "MSD") == 0)
	{
		sprintf(TRAJ->fpathO, "%s", TRAJ->fpathI);

		int *separations = new int[nMonomers-1];
		float *msd = new float[nMonomers-1];
		long *counts = new long[nMonomers-1];

		for(int i = 0; i < nMonomers-1; i++)
		{
			separations[i] = i + 1;
			msd[i] = 0.0;
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
					computeMeanSquaredDistance(TRAJ, ATOMS, BOX, separations, msd, counts, nMonomers);
			}

			delete[] ATOMS;
			delete BOX;

			if(i == nReplicas - 1)
			{
				sprintf(TRAJ->fpathO, "%s/%s", TRAJ->fpathO, fname);
				TRAJ -> createOutputFile("separation counts msd");

				for(int k = 0; k < nMonomers-1; k++)
				{
					msd[k] /= counts[k];
					fprintf(TRAJ->fileO, "\n%d %d %f", separations[k], counts[k], msd[k]);
				}

				TRAJ -> closeTrajectory(true, true);
			}
			else
				TRAJ -> closeTrajectory(true, false);
		}

		delete[] separations, counts, msd;
	}

	else if(strcmp(option, "End_Distance") == 0)
	{
		sprintf(TRAJ->fpathO, "%s/%s", TRAJ->fpathI, fname);
		sprintf(TRAJ->fpathI, "%s/traj1.cfg", TRAJ->fpathI);

		TRAJ -> openTrajectory();

		char headers[int(5*TRAJ->nMols+100)] = "Step";
		for(int i = 0; i < TRAJ->nMols; i++)
			sprintf(headers, "%s R%d", headers, i + 1);

		TRAJ -> createOutputFile(headers);

		atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
		System *BOX = new System(TRAJ->Lx, TRAJ->Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);

		while( !feof(TRAJ->fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->frame_nr >= frameStart[0] and TRAJ->frame_nr <= frameEnd[0])
				computeEndDistances(TRAJ, ATOMS, BOX, nMonomers);
		}

		fprintf(TRAJ->fileO, "\n");

		TRAJ -> closeTrajectory(true, true);
	}

	return(0);
}

void computeMeanSquaredDistance(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, int *separations, float *msd, long *counts, int nMonomers)
{
	for(int i = 0; i < TRAJ->nMols; i++)
	{
		int mol_ID = i + 1;

		for(int j = 0; j < TRAJ->nAtoms; j++)
		{
			if(ATOMS[j].mol_id == mol_ID)
			{
				int atom_j = ATOMS[j].atom_id;

				for(int k = 0; k < nMonomers-1; k++)
				{
					int atom_k = atom_j + separations[k];

					if(atom_k <= nMonomers*mol_ID)
					{
						float dx = ATOMS[atom_k-1].rxt1 - ATOMS[atom_j-1].rxt1 + TRAJ->Lx*(ATOMS[atom_k-1].jumpx - ATOMS[atom_j-1].jumpx);
						float dy = ATOMS[atom_k-1].ryt1 - ATOMS[atom_j-1].ryt1 + TRAJ->Ly*(ATOMS[atom_k-1].jumpy - ATOMS[atom_j-1].jumpy);
						float dz = ATOMS[atom_k-1].rzt1 - ATOMS[atom_j-1].rzt1 + TRAJ->Lz*(ATOMS[atom_k-1].jumpz - ATOMS[atom_j-1].jumpz);
						float dr2 = dx*dx + dy*dy + dz*dz;

						msd[k] += dr2;
						counts[k] += 1;
					}
				}
			}
		}
	}
}

void computeEndDistances(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, int nMonomers)
{
	float *end_distances = new float[TRAJ->nMols];

	for(int i = 0; i < TRAJ->nMols; i++)
	{
		int mol_ID = i + 1;

		for(int j = 0; j < TRAJ->nAtoms; j++)
		{
			if(ATOMS[j].mol_id == mol_ID and ATOMS[j].atom_id == int(i*nMonomers + 1))
			{
				float dx = ATOMS[j].rxt1 - ATOMS[j + nMonomers - 1].rxt1 + BOX->Lx*(ATOMS[j].jumpx - ATOMS[j + nMonomers - 1].jumpx);
				float dy = ATOMS[j].ryt1 - ATOMS[j + nMonomers - 1].ryt1 + BOX->Ly*(ATOMS[j].jumpy - ATOMS[j + nMonomers - 1].jumpy);
				float dz = ATOMS[j].rzt1 - ATOMS[j + nMonomers - 1].rzt1 + BOX->Lz*(ATOMS[j].jumpz - ATOMS[j + nMonomers - 1].jumpz);
				float dr2 = dx*dx + dy*dy + dz*dz;

				end_distances[i] = sqrt(dr2);
			}
		}
	}

	fprintf(TRAJ->fileO, "\n%ld", TRAJ->step);

	for(int i = 0; i < TRAJ->nMols; i++)
		fprintf(TRAJ->fileO, " %f", end_distances[i]);
}