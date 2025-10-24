/*
	collisions.cpp

	Tracks particle collisions in a binary fluid mixture of oppositely 
	driven colliods and calculates the x- and y-displacements of the 
	tracer particle as collisions take place. 

	Date created  : 15.10.25
	Last modified : 18.10.25
*/

#include "analysis.h"

using namespace std;
using namespace analysis;

int main(int argc, char *argv[])
{
	// -------- System params --------
	float Lx = 50.0, Ly = 50.0;
	int nAtomTypes = 2;

	// -------- Trajectory params --------
	int frameStart = int(1e4), frameEnd = int(5e4)-1;
	float dt = 5e-4;
	int frameW = int(90);

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ -> fpathI, "//media/ashwin/Expansion/ashwin_md/psps/tau_5e-2/Fd500/traj2.cfg");
	sprintf(TRAJ -> fpathO, "//media/ashwin/Expansion/ashwin_md/psps/tau_5e-2/Fd500/collisions.dat");

	if(frameEnd < frameStart)
	{
		printf("Error: Ending frame cannot be smaller than the starting frame. Exiting...\n");
		exit(-1);
	}
	else
	{
		TRAJ -> openTrajectory();
		TRAJ -> createOutputFile("frame_col tau_col delx dely");

		TRAJ->totalFrames = frameEnd - frameStart + 1;
		atom_style **ATOMS = new atom_style*[TRAJ->totalFrames];
		System *BOX = new System(Lx, Ly, TRAJ->nAtoms, nAtomTypes, true, 1.5, 1.5);
		for(int i = 0; i < TRAJ->totalFrames; i++)
			ATOMS[i] = new atom_style[TRAJ->nAtoms];

		TRAJ -> loadTrajectory(ATOMS, BOX, frameStart, frameEnd, true);

		long **group_lists = TRAJ -> sortAtomsByType(ATOMS[0], BOX);		
		int tracerID = int(group_lists[1][0] - 1);

		int pairID = 0;
		float delx = 0.0, dely = 0.0;
		int steps = 0, frame_col = 0;
		bool track = false;

		for(int currFrame = 0; currFrame < TRAJ->totalFrames; currFrame++)
		{
			printf("\nProcessing frame %d", currFrame);

			if(track == true)
			{
				float dxij = ATOMS[currFrame][pairID].rxt1 - ATOMS[currFrame][tracerID].rxt1;
				float dxTracer = ATOMS[currFrame][tracerID].rxt1 - ATOMS[currFrame - 1][tracerID].rxt1;
				float dyTracer = ATOMS[currFrame][tracerID].ryt1 - ATOMS[currFrame - 1][tracerID].ryt1; 

				BOX -> checkMinImage(&dxij, NULL);
				BOX -> checkMinImage(&dxTracer, &dyTracer);

				delx += abs(dxTracer);
				dely += abs(dyTracer);
				steps += 1;

				if(dxij <= 0)
				{
					fprintf(TRAJ->fileO, "\n%d %f %f %f", frame_col, float(steps*TRAJ->timeStep*TRAJ->frameWidth), delx, dely);
					// printf("Particles %d and %d underwent a collision event that spanned %d frames during which the tracer particle traversed dx=%f and dy=%f\n", tracerID, pairID, steps, delx, dely);
					delx = 0.0;
					dely = 0.0;
					steps = 0;	
					track = false;
				}

			}

			if(track == false)
			{
				BOX -> buildCellList(ATOMS[currFrame]);

				float rxTracer = ATOMS[currFrame][tracerID].rxt1;
				float ryTracer = ATOMS[currFrame][tracerID].ryt1;
				int icellx = int(rxTracer/BOX->rcellx);
				int icelly = int(ryTracer/BOX->rcelly);

				int searchCells[6];
				searchCells[0] = BOX->cellindex(icellx, icelly);
				searchCells[1] = BOX->cellindex(icellx + 1, icelly);
				searchCells[2] = BOX->cellindex(icellx + 1, icelly + 1);
				searchCells[3] = BOX->cellindex(icellx + 1, icelly - 1);
				searchCells[4] = BOX->cellindex(icellx, icelly + 1);
				searchCells[5] = BOX->cellindex(icellx, icelly - 1);

				for(int nbor = 0; nbor < 6; nbor++)
				{
					int cell_id = searchCells[nbor];

					int i = BOX->HEAD[cell_id];
					while(i != 0)
					{
						int ii = i - 1;
						if(ii != tracerID)
						{
							float dxij = ATOMS[currFrame][ii].rxt1 - rxTracer;
							float dyij = ATOMS[currFrame][ii].ryt1 - ryTracer;

							BOX -> checkMinImage(&dxij, &dyij);

							if(dxij <= 1.0 and abs(dyij) < 0.6)
							{
								pairID = ii;
								track = true;
								frame_col = frameStart + currFrame;
								// printf("Warning from frame %d! Particles %d and %d are about to collide at distance of dx=%f, dy=%f\n", frameStart + currFrame, tracerID, pairID, dxij, dyij);
								break;
							}
						}

						i = BOX->LIST[i];
					}
				}
			}
		}

		printf("\n\nCollision tracking complete!\n");

		delete[] ATOMS;

		TRAJ -> closeTrajectory(true, true);
	}

	return(0);
}