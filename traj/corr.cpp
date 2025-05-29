// corr.cpp

#include "analysis.h"

using namespace analysis;

Hist2D ***buildHistogram(System *BOX, int nAtomTypes, float rcut_x, float rcut_y);
void computeCorrelations(Hist2D ***RDF, atom_style *ATOMS, System *BOX, float rcut_x, float rcut_y);
void computeCorrelations_NSQ(Hist2D ***RDF, atom_style *ATOMS, System *BOX, float rcut_x, float rcut_y);

int main(int argc, char *argv[])
{
	// Trajectory params
	float timeStart = 0.9e3, timeEnd = 1.0e3, dt = 5e-4;
	int frameW = 1000;

	// System params
	int nAtomTypes = 2;
	float Lx = 120.0, Ly = 30.0, rho = 0.45;

	// Histogram params
	float binWidth_x = 0.01, binWidth_y = 0.01;
	float rcut_x = 15.0, rcut_y = 15.0;

	Trajectory *TRAJ = new Trajectory(dt, frameW);
	sprintf(TRAJ->fpathI, "//media/ashwin/One Touch/ashwin_md/Apr2025/lane/cpp/Data11/traj1.xyz");
	sprintf(TRAJ->fpathO, "//media/ashwin/One Touch/ashwin_md/Apr2025/lane/cpp/Data11/rdf.dat");
	
	TRAJ -> openTrajectory();
	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	System *BOX = new System(Lx, Ly, TRAJ->nAtoms, binWidth_x, binWidth_y);
	Hist2D ***RDF = buildHistogram(BOX, nAtomTypes + 1, rcut_x, rcut_y);

	int nFrames = 0;

	while( !feof(TRAJ->fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		TRAJ->time = TRAJ->step*TRAJ->timeStep;

		if(TRAJ->time >= timeStart and TRAJ->time <= timeEnd)
		{
			nFrames++;
			printf("Processing timeframe %f, frame %d\n", TRAJ->time, nFrames);
			computeCorrelations_NSQ(RDF, ATOMS, BOX, rcut_x, rcut_y);
		}
	}

	for(int i = 0; i <= nAtomTypes; i++)
	{
		for(int j = 0; j <= nAtomTypes; j++)
		{
			if((i == j) or (i > 0 and j > 0))
				RDF[i][j] -> normalize(rho * TRAJ->nAtoms * nFrames);
		}
	}

	TRAJ -> write2file(RDF, nAtomTypes + 1);
	TRAJ -> closeTrajectory();

	return(0);
}

Hist2D ***buildHistogram(System *BOX, int nAtomTypes, float rcut_x, float rcut_y)
{
	Hist2D ***HIST = new Hist2D**[nAtomTypes];

	for(int i = 0; i < nAtomTypes; i++)
	{
		HIST[i] = new Hist2D*[nAtomTypes];
		for(int j = 0; j < nAtomTypes; j++)
			HIST[i][j] = new Hist2D(rcut_x, rcut_y, BOX->rcellx, BOX->rcelly);
	}

	return(HIST);
}

void computeCorrelations_NSQ(Hist2D ***RDF, atom_style *ATOMS, System *BOX, float rcut_x, float rcut_y)
{
	for(int i = 0; i < BOX->nAtoms; i++)
	{
		float rxi = ATOMS[i].rxt1;
		float ryi = ATOMS[i].ryt1;

		for(int j = i + 1; j < BOX->nAtoms; j++)
		{
			float rxj = ATOMS[j].rxt1;
			float ryj = ATOMS[j].ryt1;

			float dxij = rxj - rxi;
			float dyij = ryj - ryi;

			BOX -> checkMinImage(&dxij, &dyij);

			if(abs(dxij) < rcut_x and abs(dyij) < rcut_y)
			{
				// if(abs(dxij) < BOX->rcellx and abs(dyij) < BOX->rcelly)
				// 	printf("Atoms %d and %d are separated by dx = %f, dy = %f\n", i, j, dxij, dyij);

				RDF[0][0]->addToBin('X', 'X', dxij, dyij);

				if(ATOMS[i].id == 'N')
				{
					if(ATOMS[j].id == 'N')
						RDF[1][1]->addToBin('N', 'N', dxij, dyij); 

					else if(ATOMS[j].id == 'O')
					{
						RDF[1][2]->addToBin('N', 'O', dxij, dyij);
						RDF[2][1]->addToBin('O', 'N', -dxij, -dyij);
					}
				}

				else if(ATOMS[i].id == 'O')
				{
					if(ATOMS[j].id == 'N')
					{
						RDF[2][1]->addToBin('O', 'N', dxij, dyij); 
						RDF[1][2]->addToBin('N', 'O', -dxij, -dyij);
					}

					else if(ATOMS[j].id == 'O')
						RDF[2][2]->addToBin('O', 'O', dxij, dyij);
				}
			}
		}
	}
}

void computeCorrelations(Hist2D ***RDF, atom_style *ATOMS, System *BOX, float rcut_x, float rcut_y)
{
	BOX -> buildCellList(ATOMS);

	int ncellx = BOX -> Ncellx;
	int ncelly = BOX -> Ncelly;
	int maxBin_x = int(rcut_x/BOX->rcellx);
	int maxBin_y = int(rcut_y/BOX->rcelly);

	for(int iy = 0; iy < ncelly; iy++)
	{
		for(int ix = 0; ix < ncellx; ix++)		
		{
			int icell = BOX -> cellindex(ix, iy);

			// intracellular pairs
			int i = BOX->HEAD[icell];
			while (i != 0)
			{
				int ii = i - 1;
				float rxi = ATOMS[ii].rxt1;
				float ryi = ATOMS[ii].ryt1;

				int j = BOX->LIST[i];
				while(j != 0)
				{
					int jj = j - 1;
					float rxj = ATOMS[jj].rxt1;
					float ryj = ATOMS[jj].ryt1;

					float dxij = rxj - rxi;
					float dyij = ryj - ryi;

					BOX -> checkMinImage(&dxij, &dyij);

					if(abs(dxij) < rcut_x and abs(dyij) < rcut_y)
					{

						if(ATOMS[ii].id == 'N')
						{
							if(ATOMS[jj].id == 'N')
								RDF[0][0]->addToBin('N', 'N', dxij, dyij); 

							else if(ATOMS[jj].id == 'O')
							{
								RDF[0][1]->addToBin('N', 'O', dxij, dyij);
								RDF[1][0]->addToBin('O', 'N', -dxij, -dyij);
							}
						}

						else if(ATOMS[ii].id == 'O')
						{
							if(ATOMS[jj].id == 'N')
							{
								RDF[1][0]->addToBin('O', 'N', dxij, dyij); 
								RDF[0][1]->addToBin('N', 'O', -dxij, -dyij);
							}

							else if(ATOMS[jj].id == 'O')
								RDF[1][1]->addToBin('O', 'O', dxij, dyij);
						}
					}
					
					j = BOX->LIST[j];
				}

				i = BOX->LIST[i];
			}

			// intercellular pairs
			for(int jy = 0; jy <= maxBin_y + 1; jy++)
			{
				int jx_start = -(maxBin_x + 1);
				if(jy == 0) jx_start = 1;

				for(int jx = jx_start; jx <= maxBin_x + 1; jx++)
				{
					int jcell = BOX -> cellindex(ix + jx, iy + jy);

					int j = BOX->HEAD[jcell];
					while(j != 0)
					{
						int jj = j - 1;
						float rxj = ATOMS[jj].rxt1;
						float ryj = ATOMS[jj].ryt1;

						int i = BOX->HEAD[icell];
						while(i != 0)
						{
							int ii = i - 1;
							float rxi = ATOMS[ii].rxt1;
							float ryi = ATOMS[ii].ryt1;

							float dxij = rxj - rxi;
							float dyij = ryj - ryi;

							BOX -> checkMinImage(&dxij, &dyij);

							if(abs(dxij) < rcut_x and abs(dyij) < rcut_y)
							{
								
								if(ATOMS[ii].id == 'N')
								{
									if(ATOMS[jj].id == 'N')
										RDF[0][0]->addToBin('N', 'N', dxij, dyij); 

									else if(ATOMS[jj].id == 'O')
									{
										RDF[0][1]->addToBin('N', 'O', dxij, dyij);
										RDF[1][0]->addToBin('O', 'N', -dxij, -dyij);
									}
								}

								else if(ATOMS[ii].id == 'O')
								{
									if(ATOMS[jj].id == 'N')
									{
										RDF[1][0]->addToBin('O', 'N', dxij, dyij); 
										RDF[0][1]->addToBin('N', 'O', -dxij, -dyij);
									}

									else if(ATOMS[jj].id == 'O')
										RDF[1][1]->addToBin('O', 'O', dxij, dyij);
								}
							}

							i = BOX->LIST[i];
						}

						j = BOX->LIST[j];
					}
				}
			}
		}
	}
}

void analysis::Trajectory::write2file(Hist2D ***RDF, int nAtomTypes)
{
	remove(fpathO);

	fileO = fopen(fpathO, "w");
	if(fileO == NULL)
	{
		printf("Couldn't create new data file %s. Exiting...\n", fpathO);
		exit(-1);
	}

	char ***name = new char**[nAtomTypes];
	for(int i = 0; i < nAtomTypes; i++)
	{
		name[i] = new char*[nAtomTypes];
		for(int j = 0; j < nAtomTypes; j++)
			name[i][j] = new char[10];
	}

	sprintf(name[0][0], "all");
	sprintf(name[1][1], "AA");
	sprintf(name[1][2], "AB");
	sprintf(name[2][1], "BA");
	sprintf(name[2][2], "BB");

	for(int i = 0; i < nAtomTypes; i++)
	{
		for(int j = 0; j < nAtomTypes; j++)
		{
			if((i == j) or (i > 0 and j > 0))
			{	
				fprintf(fileO, "g%s\n", name[i][j]);

				for(int ix = 0; ix < RDF[i][j]->nBin_x; ix++)
				{
					for(int iy = 0; iy < RDF[i][j]->nBin_y; iy++)
						fprintf(fileO, "%f ", RDF[i][j]->bin[ix][iy]);

					if(int iy = RDF[i][j]->nBin_y) 
						fprintf(fileO, "\n");
				}
			}
		}
	}	

	fclose(fileO);
}