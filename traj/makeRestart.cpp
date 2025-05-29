/*
	makeRestart.cpp

	Generates a data file for LAMMPS restart simulations by importing coordinates and velocoties 
	from a .cfg file or just coordinates from a .xyz file  

	Date created  : 05.05.25
	Last modified : 29.05.25
*/

#include "analysis.h"

using namespace analysis;

int main(int argc, char* argv[])
{
	// ---------- System Params ----------
	float Lx = 150.0, Ly = 30.0;
	int nAtomTypes = 2;
	char typeLabels[2] = {'N', 'O'};

	// ---------- Trajectory Params ----------
	float dt = 5e-4;
	int frameW = int(1e2);
	long eqStep = long(1e7);
	long chooseStep = long(8e9);

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/One Touch/ashwin_md/lane/May2025/lmp/Data30/frames/frame.10010000000.cfg");
	sprintf(TRAJ->fpathO, "//media/ashwin/One Touch/ashwin_md/lane/May2025/lmp/Data30/lane.10010000000.res");
	TRAJ -> openTrajectory();

	atom_style *ATOMS = new atom_style[TRAJ -> nAtoms]; 

	if(strcmp(TRAJ->format, "xyz") == 0)
	{
		while( !feof(TRAJ -> fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ -> step == (eqStep + chooseStep))
			{
				printf("Step %ld extracted for restart simulations\n", TRAJ->step);
				break;
			}
		}	
	}

	else if(strcmp(TRAJ->format, "cfg") == 0)
	{
		TRAJ -> readThisFrame(ATOMS);
		printf("Step %ld extracted for restart simulations\n", TRAJ->step);
	}

	TRAJ -> closeTrajectory();

	remove(TRAJ->fpathO);
	TRAJ->fileO = fopen(TRAJ->fpathO, "w");

	fprintf(TRAJ->fileO, "makeRestart.cpp : LAMMPS data file for restart simulation generated using %ld step from trajectory file %s\n\n", TRAJ->step, TRAJ->fpathI);

	fprintf(TRAJ->fileO, "%d atoms\n", TRAJ->nAtoms);
	fprintf(TRAJ->fileO, "%f %f xlo xhi\n", 0.0, Lx);
	fprintf(TRAJ->fileO, "%f %f ylo yhi\n", 0.0, Ly);
	fprintf(TRAJ->fileO, "%d atom types\n\n", nAtomTypes);
	
	fprintf(TRAJ->fileO, "Atom Type Labels\n\n");

	for(int i = 0; i < nAtomTypes; i++)
		fprintf(TRAJ->fileO, "%d %c\n", i + 1, typeLabels[i]);

	fprintf(TRAJ->fileO, "\nMasses\n\n");

	for(int i = 0; i < nAtomTypes; i++)
		fprintf(TRAJ->fileO, "%d %f\n", i + 1, 1.0);

	fprintf(TRAJ->fileO, "\nAtoms # atomic\n\n");

	for(int i = 0; i < TRAJ->nAtoms; i++)
		fprintf(TRAJ->fileO, "%d %c %f %f %f\n", i + 1, ATOMS[i].id, ATOMS[i].rxt1, ATOMS[i].ryt1, 0.0);

	if(strcmp(TRAJ->format, "cfg") == 0)
	{
		fprintf(TRAJ->fileO, "\nVelocities # atomic\n\n");

		for(int i = 0; i < TRAJ->nAtoms; i++)
			fprintf(TRAJ->fileO, "%d %f %f %f\n", i + 1, ATOMS[i].vx, ATOMS[i].vy, ATOMS[i].vz);
	}
	
	fclose(TRAJ->fileO);
}