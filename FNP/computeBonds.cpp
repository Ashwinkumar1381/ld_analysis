#include "../traj/analysis.h"

using namespace analysis;

int main(int argc, char *argv[])
{
	int nAtomTypes = 2;

	int frameStart = int(0), frameEnd = int(1000);
	float dt = 2e-2;
	int frameW = 50;
	char option[50] = "bonds_angles"; 

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/ASH_DRIVE_3/ashwin_md/FNP/polymer_benchmarks/dpd_large/aSS_200/K_theta_10/a_400/R1");

	if(strcmp(option, "bonds") == 0 or strcmp(option, "bonds_angles") == 0)
	{
		sprintf(TRAJ->fpathO, "%s", TRAJ->fpathI);
		sprintf(TRAJ->fpathI, "%s/traj1.cfg", TRAJ->fpathO);

		TRAJ -> openTrajectory();

		sprintf(TRAJ->fpathI, "%s", TRAJ->fpathO);
		sprintf(TRAJ->fpathO, "%s/bonds.dat", TRAJ->fpathI);

		TRAJ -> createOutputFile("Step Atom_i Atom_j Bond_length");

		atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
		System *BOX = new System(TRAJ->Lx, TRAJ->Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);

		printf("\nComputing Bond separations...\n");

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

					fprintf(TRAJ->fileO, "\n%ld %d %d %f", TRAJ->step, ATOMS[i-1].atom_id, ATOMS[i].atom_id, dr);
				}
			}
		}

		TRAJ -> closeTrajectory();
	}

	if(strcmp(option, "angles") == 0 or strcmp(option, "bonds_angles") == 0)
	{
		sprintf(TRAJ->fpathO, "%s", TRAJ->fpathI);
		sprintf(TRAJ->fpathI, "%s/traj1.cfg", TRAJ->fpathO);

		TRAJ -> openTrajectory();

		sprintf(TRAJ->fpathI, "%s", TRAJ->fpathO);
		sprintf(TRAJ->fpathO, "%s/angles.dat", TRAJ->fpathI);

		TRAJ -> createOutputFile("Step Atom_i Atom_j Atom_k Cos_theta Sin_theta");

		atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
		System *BOX = new System(TRAJ->Lx, TRAJ->Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);

		printf("\nComputing Bond angles...\n");

		while ( !feof(TRAJ -> fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
			{
				for(int i = 1; i < TRAJ->nAtoms - 1; i++)
				{
					float bond_A[4], bond_B[4], bond_C[4];
					for(int j = 0; j < 4; j++)
					{
						bond_A[j] = 0.0;
						bond_B[j] = 0.0;
						bond_C[j] = 0.0;
					}
					
					float dx = ATOMS[i-1].rxt1 - ATOMS[i].rxt1;
					float dy = ATOMS[i-1].ryt1 - ATOMS[i].ryt1;
					float dz = ATOMS[i-1].rzt1 - ATOMS[i].rzt1;
					BOX -> checkMinImage(&dx, &dy, &dz);

					bond_A[0] = dx; bond_A[1] = dy; bond_A[2] = dz;
					bond_A[3] = sqrt(dx*dx + dy*dy + dz*dz);

					dx = ATOMS[i+1].rxt1 - ATOMS[i].rxt1;
					dy = ATOMS[i+1].ryt1 - ATOMS[i].ryt1;
					dz = ATOMS[i+1].rzt1 - ATOMS[i].rzt1;
					BOX -> checkMinImage(&dx, &dy, &dz);

					bond_B[0] = dx; bond_B[1] = dy; bond_B[2] = dz;
					bond_B[3] = sqrt(dx*dx + dy*dy + dz*dz);

					float proj_B_A = (bond_B[0]*bond_A[0] + bond_B[1]*bond_A[1] + bond_B[2]*bond_A[2])/bond_A[3];

					bond_C[0] = bond_B[0] - proj_B_A*(bond_A[0]/bond_A[3]);
					bond_C[1] = bond_B[1] - proj_B_A*(bond_A[1]/bond_A[3]);
					bond_C[2] = bond_B[2] - proj_B_A*(bond_A[2]/bond_A[3]);
					bond_C[3] = sqrt(bond_C[0]*bond_C[0] + bond_C[1]*bond_C[1] + bond_C[2]*bond_C[2]);

					float Cos_theta = proj_B_A/bond_B[3];
					float Sin_theta = bond_C[3]/bond_B[3];

					fprintf(TRAJ->fileO, "\n%ld %d %d %d %f %f", TRAJ->step, i-1, i, i+1, Cos_theta, Sin_theta);
				}
			}
		}

		TRAJ -> closeTrajectory();
	}
}