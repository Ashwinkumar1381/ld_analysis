/* power.cpp
 
 Performs the power balance on a set of particle configurations

 Date created  : 21.07.25
 Last modified : 21.07.25 
 
*/

#include "analysis.h"

using namespace analysis;
using namespace program;

float *computePowerDistribution(atom_style *ATOMS, System *BOX, WCA_2P *INTERACTIONS, float Fd, float tau);

int main(int argc, char* argv[])
{
	// ---------- Trajectory Params ----------
	float dt = 5e-4;
	int frameW = int(1e5);
	int frameStart = int(4e4), frameEnd = int(5e4);

	// ---------- System Params ----------
	float Lx = 150.0, Ly = 30.0;
	int nAtomTypes = 2;
	float Fd = 90.0;
	float tau = 1e-3;

	char *option = new char[20];
	// sprintf(option, "time_evolve_traj");
	sprintf(option, "time_avg_traj");

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/Expansion/ashwin_md/lane/June_July2025/Pe90/Data38/traj2.cfg");
	sprintf(TRAJ->fpathO, "//media/ashwin/Expansion/ashwin_md/lane/June_July2025/Pe90/Data38/power.dat");
	TRAJ -> openTrajectory();

	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	System *BOX = new System(Lx, Ly, TRAJ->nAtoms, nAtomTypes);
	WCA_2P *INTERACTIONS = new WCA_2P();

	BOX -> buildCellMaps();

	if(strcmp(option, "time_evolve_traj") == 0)
	{
		TRAJ -> createOutputFile("step Wext Qdiss Qint Epot");

		int ctr = 0;
		while( !feof(TRAJ->fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
			{
				printf("\nProcessing frame %d", TRAJ->frame_nr);

				float *power = computePowerDistribution(ATOMS, BOX, INTERACTIONS, Fd, tau);
				TRAJ -> write2file(power);
				ctr++;

				delete[] power;
			}

			if(TRAJ->frame_nr == frameEnd)
			{
				printf("\n\nPower balance computed for %d frames.", ctr);
				break;
			}
		}
	}

	else if(strcmp(option, "time_avg_traj") == 0)
	{
		float power_avg[3];

		int ctr = 0;
		while( !feof(TRAJ->fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->frame_nr == frameStart)
			{
				for(int i = 0; i < 3; i++)
					power_avg[i] = 0.0;
			}

			if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
			{
				float *power = computePowerDistribution(ATOMS, BOX, INTERACTIONS, Fd, tau);

				for(int i = 0; i < 3; i++)
					power_avg[i] += power[i];

				ctr++;

				delete[] power;
			}

			if(TRAJ->frame_nr == frameEnd)
			{
				for(int i = 0; i < 3; i++)
					power_avg[i] /= ctr;

				break;
			}

		}

		printf("\n\nAveraged Power Statistics: %g %g %g", power_avg[0], power_avg[1], power_avg[2]);
	}

	TRAJ -> closeTrajectory();
}

float *computePowerDistribution(atom_style *ATOMS, System *BOX, WCA_2P *INTERACTIONS, float Fd, float tau)
{
	float Wext = 0.0, Qdiss = 0.0, Qint = 0.0;

	float pe = computeNonBondedInteractions(ATOMS, BOX, INTERACTIONS);

	for(int i = 0; i < BOX->nAtoms; i++)
	{
		if(ATOMS[i].id == 'O') Wext += ATOMS[i].vx * Fd;
		else if(ATOMS[i].id == 'N') Wext += ATOMS[i].vx * (-1.0 * Fd);

		Qdiss += (ATOMS[i].vx * ATOMS[i].vx + ATOMS[i].vy * ATOMS[i].vy) / tau;
		Qint += (ATOMS[i].vx * ATOMS[i].fx + ATOMS[i].vy * ATOMS[i].fy);
	}

	float *power = new float[4];

	power[0] = Wext / BOX->nAtoms;
	power[1] = Qdiss / BOX->nAtoms;
	power[2] = -1.0 * Qint / BOX->nAtoms;
	power[3] = pe;

	return(power);
}

void analysis::Trajectory::write2file(float *power)
{
	fprintf(fileO, "\n%ld %g %g %g %g", step, power[0], power[1], power[2], power[3]);
}