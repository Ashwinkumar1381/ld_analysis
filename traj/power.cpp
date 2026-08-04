/* power.cpp
 
 Performs the power balance on a set of particle configurations

 Date created  : 21.07.25
 Last modified : 21.08.25 
 
*/

#include "analysis.h"

using namespace analysis;
using namespace program;

float *computePowerDistribution(atom_style *ATOMS, System *BOX, WCA_2P *INTERACTIONS, float Fd, float tau, char compute_option[5] = "full");

int main(int argc, char* argv[])
{
	// ---------- Trajectory Params ----------
	float dt = 5e-4;
	int frameW = int(1e5);
	int frameStart = int(0), frameEnd = int(5e4);

	// ---------- System Params ----------
	float Lx = 150.0, Ly = 30.0;
	int nAtomTypes = 2;
	float Fd = 100.0;
	float tau = 1e-3;

	char *option = new char[20];
	sprintf(option, "time_evolve_traj");
	// sprintf(option, "time_avg_traj");

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/ASH_DRIVE_3/ashwin_md/lane/Aug_Nov2025/Fd100/tau_1e-3/traj2.cfg");
	sprintf(TRAJ->fpathO, "//media/ashwin/ASH_DRIVE_3/ashwin_md/lane/Aug_Nov2025/Fd100/tau_1e-3/power.dat");
	TRAJ -> openTrajectory();

	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	System *BOX = new System(Lx, Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);
	WCA_2P *INTERACTIONS = new WCA_2P();

	BOX -> buildCellMaps();

	if(strcmp(option, "time_evolve_traj") == 0)
	{
		TRAJ -> createOutputFile("step Wext Qdiss Qint Wfluc");

		int ctr = 0;
		while( !feof(TRAJ->fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
			{
				printf("\nProcessing frame %d", TRAJ->frame_nr);

				float *power = computePowerDistribution(ATOMS, BOX, INTERACTIONS, Fd, tau, "full");
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
		float power_avg[4];

		int ctr = 0;
		while( !feof(TRAJ->fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->frame_nr == frameStart)
			{
				for(int i = 0; i < 4; i++)
					power_avg[i] = 0.0;
			}

			if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
			{
				float *power = computePowerDistribution(ATOMS, BOX, INTERACTIONS, Fd, tau);

				for(int i = 0; i < 4; i++)
					power_avg[i] += power[i];

				ctr++;

				delete[] power;
			}

			if(TRAJ->frame_nr == frameEnd)
			{
				for(int i = 0; i < 4; i++)
					power_avg[i] /= ctr;

				break;
			}

		}

		printf("\n\nAveraged Power Statistics: %g %g %g %g", power_avg[0], power_avg[1], power_avg[2], power_avg[3]);
	}

	TRAJ -> closeTrajectory();
}

float *computePowerDistribution(atom_style *ATOMS, System *BOX, WCA_2P *INTERACTIONS, float Fd, float tau, char compute_option[5])
{
	float pe = computeNonBondedInteractions(ATOMS, BOX, INTERACTIONS);

	float Vxcom = 0.0, Vycom = 0.0;
	for(int i = 0; i < BOX->nAtoms; i++)
	{
		Vxcom += ATOMS[i].vx;
		Vycom += ATOMS[i].vy; 
	}

	Vxcom /= (BOX->nAtoms);
	Vycom /= (BOX->nAtoms);

	for(int i = 0; i < BOX->nAtoms; i++)
	{
		ATOMS[i].vxth = ATOMS[i].vx - Vxcom;
		ATOMS[i].vyth = ATOMS[i].vy - Vycom;
	}

	float Wext = 0.0, Qdiss = 0.0, Qint = 0.0, Wfluc = 0.0;

	if(strcmp(compute_option, "full") == 0)
	{
		for(int i = 0; i < BOX->nAtoms; i++)
		{
			float FBx = float(ATOMS[i].fx + ATOMS[i].vxth/tau - ATOMS[i].fx_int - Fd*ATOMS[i].si);
			float FBy = float(ATOMS[i].fy + ATOMS[i].vyth/tau - ATOMS[i].fy_int);

			Wext  += ATOMS[i].vxth*Fd*ATOMS[i].si;
			Qdiss += (ATOMS[i].vxth*ATOMS[i].vxth + ATOMS[i].vyth*ATOMS[i].vyth)/tau;
			Qint  += (ATOMS[i].vxth*ATOMS[i].fx_int + ATOMS[i].vyth*ATOMS[i].fy_int);
			Wfluc += (ATOMS[i].vxth*FBx + ATOMS[i].vyth*FBy);
		}	
	}

	else if(strcmp(compute_option, "x") == 0)
	{
		for(int i = 0; i < BOX->nAtoms; i++)
		{
			float FBx = float(ATOMS[i].fx + ATOMS[i].vxth/tau - ATOMS[i].fx_int - Fd*ATOMS[i].si);

			Wext  += ATOMS[i].vxth*Fd*ATOMS[i].si;
			Qdiss += (ATOMS[i].vxth*ATOMS[i].vxth)/tau;
			Qint  += (ATOMS[i].vxth*ATOMS[i].fx_int);
			Wfluc += (ATOMS[i].vxth*FBx);
		}
	}

	float *power = new float[4];

	power[0] = Wext / BOX->nAtoms;
	power[1] = Qdiss / BOX->nAtoms;
	power[2] = -1.0 * Qint / BOX->nAtoms;
	power[3] = Wfluc / BOX->nAtoms;

	return(power);
}

void analysis::Trajectory::write2file(float *power)
{
	fprintf(fileO, "\n%ld %g %g %g %g", step, power[0], power[1], power[2], power[3]);
}