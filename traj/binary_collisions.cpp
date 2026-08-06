#include "analysis.h"

using namespace std;
using namespace analysis;

int main(int argc, char *argv[])
{
	float Lx = 10.0, Ly = 10.0;
	int nAtomTypes = 2;

	int frameStart = int(0), frameEnd = int(2e3);
	float dt = 5e-4;
	int frameW = int(100);

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//run/media/ethaya_lab_c1/ASH_DRIVE_3/ashwin_md/binary_collisions/tests/traj.cfg");
	sprintf(TRAJ->fpathO, "//run/media/ethaya_lab_c1/ASH_DRIVE_3/ashwin_md/binary_collisions/tests/collision_stats.dat");

	int run_nr = 0;
	sscanf(argv[1], "%d", &run_nr);

	TRAJ -> openTrajectory();

	if(run_nr == 1)
		TRAJ -> createOutputFile("run_nr tau_col dely_init delx_dist dely_dist delx_disp dely_disp vel_corr", true);
	else
		TRAJ -> createOutputFile("", false);

	TRAJ->totalFrames = frameEnd - frameStart + 1;
	atom_style **ATOMS = new atom_style*[TRAJ->totalFrames];
	System *BOX = new System(Lx, Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);
	for(int i = 0; i < TRAJ->totalFrames; i++)
		ATOMS[i] = new atom_style[TRAJ->nAtoms];

	TRAJ -> loadTrajectory(ATOMS, BOX, frameStart, frameEnd);

	int nFrames = 0;
	float dely_init = 0.0, dist_x = 0.0, dist_y = 0.0, disp_x = 0.0, disp_y = 0.0;
	bool track = false;

	float *u_init = new float[2];
	float *u_final = new float[2];
	for(int i = 0; i < 2; i++)
	{
		u_init[i] = 0.0;
		u_final[i] = 0.0;
	}
	
	for(int currFrame = 0; currFrame < TRAJ->totalFrames; currFrame++)
	{
		float dxij = ATOMS[currFrame][1].rxt1 - ATOMS[currFrame][0].rxt1 + BOX->Lx*(ATOMS[currFrame][1].jumpx - ATOMS[currFrame][0].jumpx);
		float dyij = ATOMS[currFrame][1].ryt1 - ATOMS[currFrame][0].ryt1 + BOX->Ly*(ATOMS[currFrame][1].jumpy - ATOMS[currFrame][0].jumpy);

		if(track == true)
		{
			float dxTracer = ATOMS[currFrame][0].rxt1 - ATOMS[currFrame - 1][0].rxt1 + BOX->Lx*(ATOMS[currFrame][0].jumpx - ATOMS[currFrame - 1][0].jumpx);
			float dyTracer = ATOMS[currFrame][0].ryt1 - ATOMS[currFrame - 1][0].ryt1 + BOX->Ly*(ATOMS[currFrame][0].jumpy - ATOMS[currFrame - 1][0].jumpy);

			dist_x += abs(dxTracer);
			dist_y += abs(dyTracer);
			nFrames += 1;

			if(dxij < 0.0 and abs(dyij) >= 1.0)
			{
				track = false;

				disp_x = ATOMS[currFrame][0].rxt1 + BOX->Lx*ATOMS[currFrame][0].jumpx - disp_x;
				disp_y = ATOMS[currFrame][0].ryt1 + BOX->Ly*ATOMS[currFrame][0].jumpy - disp_y;
				float tau_col = float(nFrames*frameW*dt);

				float vx_Tracer = ATOMS[currFrame][0].vx;
				float vy_Tracer = ATOMS[currFrame][0].vy;
				float vel_Tracer = sqrt(vx_Tracer*vx_Tracer+vy_Tracer*vy_Tracer);
				u_final[0] = float(vx_Tracer/vel_Tracer);
				u_final[1] = float(vy_Tracer/vel_Tracer);

				float vel_corr = u_init[0]*u_final[0] + u_init[1]*u_final[1];

				fprintf(TRAJ->fileO, "\n%d %g %g %g %g %g %g %g", run_nr, tau_col, dely_init, dist_x, dist_y, disp_x, disp_y, vel_corr);

				delete[] ATOMS;

				TRAJ->closeTrajectory(true, true);
				break;
			}
		}

		if(track == false)
		{
			if(dxij < 1.0 and abs(dyij) < 1.0)
			{
				track = true;
				disp_x = ATOMS[currFrame][0].rxt1 + BOX->Lx*ATOMS[currFrame][0].jumpx;
				disp_y = ATOMS[currFrame][0].ryt1 + BOX->Ly*ATOMS[currFrame][0].jumpy;
				dely_init = dyij;

				float vx_Tracer = ATOMS[currFrame][0].vx;
				float vy_Tracer = ATOMS[currFrame][0].vy;
				float vel_Tracer = sqrt(vx_Tracer*vx_Tracer+vy_Tracer*vy_Tracer);
				u_init[0] = float(vx_Tracer/vel_Tracer);
				u_init[1] = float(vy_Tracer/vel_Tracer);
			}
		}
	}

	return(0);
}