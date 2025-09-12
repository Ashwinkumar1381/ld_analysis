/*
	makeRestart.cpp

	Generates a data file for LAMMPS restart simulations by importing coordinates and velocoties 
	from a .cfg file or just coordinates from a .xyz file  

	Date created  : 05.05.25
	Last modified : 31.07.25
*/

#include "analysis.h"

using namespace analysis;

int main(int argc, char* argv[])
{
	// ---------- System Params ----------
	float Lx = 150.0, Ly = 30.0;
	int nAtomTypes = 2;
	char typeLabels[2] = {'O', 'N'};

	// ---------- Trajectory Params ----------
	float dt = 5e-4;
	int frameW = int(1e5);
	long eqStep = long(1e7);
	long chooseStep = long(248e7);

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/Expansion/ashwin_md/lane/Aug_Sept2025/Fd80/tau_5e-1/traj2.cfg");
	sprintf(TRAJ->fpathO, "//media/ashwin/Expansion/ashwin_md/lane/Aug_Sept2025/Fd80/tau_5e-1/lane.%ld.res", eqStep + chooseStep);
	TRAJ -> openTrajectory();

	atom_style *ATOMS = new atom_style[TRAJ -> nAtoms]; 

	while( !feof(TRAJ -> fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if(TRAJ -> step == (eqStep + chooseStep))
		{
			printf("Step %ld extracted for restart simulations\n", TRAJ->step);
			break;
		}
	}	

	TRAJ -> createOutputFile("");

	fprintf(TRAJ->fileO, "makeRestart.cpp : LAMMPS data file for restart simulation generated using %ld step from trajectory file %s\n", TRAJ->step, TRAJ->fpathI);

	fprintf(TRAJ->fileO, "\n%d atoms", TRAJ->nAtoms);
	fprintf(TRAJ->fileO, "\n%f %f xlo xhi", 0.0, Lx);
	fprintf(TRAJ->fileO, "\n%f %f ylo yhi", 0.0, Ly);
	fprintf(TRAJ->fileO, "\n%d atom types\n", nAtomTypes);
	
	fprintf(TRAJ->fileO, "\nAtom Type Labels\n");

	for(int i = 0; i < nAtomTypes; i++)
		fprintf(TRAJ->fileO, "\n%d %c", i + 1, typeLabels[i]);

	fprintf(TRAJ->fileO, "\n\nMasses\n");

	for(int i = 0; i < nAtomTypes; i++)
		fprintf(TRAJ->fileO, "\n%d %f", i + 1, 1.0);

	fprintf(TRAJ->fileO, "\n\nAtoms # atomic\n");

	for(int i = 0; i < TRAJ->nAtoms; i++)
		fprintf(TRAJ->fileO, "\n%d %c %f %f %f", i + 1, ATOMS[i].id, ATOMS[i].rxt1, ATOMS[i].ryt1, 0.0);

	if(strcmp(TRAJ->format, "cfg") == 0)
	{
		fprintf(TRAJ->fileO, "\n\nVelocities # atomic\n");

		for(int i = 0; i < TRAJ->nAtoms; i++)
			fprintf(TRAJ->fileO, "\n%d %f %f %f", i + 1, ATOMS[i].vx, ATOMS[i].vy, ATOMS[i].vz);
	}
	
	TRAJ -> closeTrajectory();
}