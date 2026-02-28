/*
	computeRDF_2D.cpp

	Calculates 2D full and partial pair distribution functions - g(x, y), on the 2D plane 

	Author		  : Ashwin Kumar
	Date created  : 25.04.25
	Last modified : 20.02.26
*/

#include "analysis.h"

using namespace analysis;

float ***buildHistogram(int nRDFtypes, float Rcut[], float binW[]);
void computeRDF_2D(float ***RDF, atom_style *ATOMS, System *BOX, int nRDFtypes, float Rcut[], float binW[]);

int main(int argc, char *argv[])
{
	// ---------- Trajectory params ----------
	long eq_steps = long(1e7);
	long startSteps = long(6e9), endSteps = long(8e9);
	float dt = 5e-4;
	int frameW = int(1e5);

	// ---------- System params ----------
	int nAtomTypes = 2;
	float xA = 0.5, xB = 0.5;
	float Lx = 150.0, Ly = 30.0, rho = 0.45;

	// ---------- Histogram params ----------
	float Rcut[2] = {10.0, 10.0};
	float binW[2] = {0.05, 0.05};

	Trajectory *TRAJ = new Trajectory(dt, frameW);
	sprintf(TRAJ->fpathI, "//media/ashwin/One Touch/ashwin_md/Apr2025/lane/lmp/Data25/traj2.xyz");
	sprintf(TRAJ->fpathO, "//media/ashwin/One Touch/ashwin_md/Apr2025/lane/lmp/Data25/laneRDF_2D.dat");
	
	TRAJ -> openTrajectory();
	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	System *BOX = new System(Lx, Ly, TRAJ->nAtoms, nAtomTypes);

	// ---------- RDF initialization ----------
	int nRDFtypes = 1;
	if(nAtomTypes > 1) nRDFtypes = nAtomTypes*nAtomTypes + 1;

	int Nbins[2] = {2*int(Rcut[0]/binW[0]), 2*int(Rcut[1]/binW[1])}; 

	float ***RDF = buildHistogram(nRDFtypes, Rcut, binW);

	int nFrames = 0;
	while( !feof(TRAJ->fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);
		if(TRAJ->step >= (startSteps + eq_steps) and TRAJ->step <= (endSteps + eq_steps))
		{
			nFrames++;
			printf("Processing timestep %ld, frame %d\n", TRAJ->step - eq_steps, nFrames);
			computeRDF_2D(RDF, ATOMS, BOX, nRDFtypes, Rcut, binW);
		}
	}

	// ---------- RDF normalization ----------
	float FAC[nRDFtypes] = {1.0, xA*xA, xA*xB, xB*xA, xB*xB};
	float norm = rho * BOX->nAtoms * binW[0] * binW[1] * nFrames;

	for(int i = 0; i < nRDFtypes; i++)
	{
		for(int j = 0; j < Nbins[0]; j++)
		{
			for(int k = 0; k < Nbins[1]; k++)
			{
				RDF[i][j][k] /= (norm * FAC[i]);
			}
		}
	}

	TRAJ -> write2file(RDF, nRDFtypes, Nbins);
	TRAJ -> closeTrajectory();

	return(0);
}

float ***buildHistogram(int nRDFtypes, float Rcut[], float binW[])
{
	int nBin_x = 2*int(Rcut[0] / binW[0]);
	int nBin_y = 2*int(Rcut[1] / binW[1]);

	float ***HIST = new float**[nRDFtypes];

	for(int i = 0; i < nRDFtypes; i++)
	{
		HIST[i] = new float*[nBin_x];

		for(int j = 0; j < nBin_x; j++)
		{
			HIST[i][j] = new float[nBin_y];

			for(int k = 0; k < nBin_y; k++)
			{
				HIST[i][j][k] = 0.0;  	
			}
		}
	}

	return(HIST);
}

void computeRDF_2D(float ***RDF, atom_style *ATOMS, System *BOX, int nRDFtypes, float Rcut[], float binW[])
{
	int Nbins[2] = {2*int(Rcut[0]/binW[0]), 2*int(Rcut[1]/binW[1])};

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

			if(abs(dxij) < Rcut[0] and abs(dyij) < Rcut[1])
			{
				int bin_xi = int((Rcut[0] + dxij) / binW[0]);
				int bin_xj = int((Rcut[0] - dxij) / binW[0]);
				int bin_yi = int((Rcut[1] + dyij) / binW[1]);
				int bin_yj = int((Rcut[1] - dyij) / binW[1]);

				if(bin_xi == Nbins[0]) bin_xi--;
				if(bin_xj == Nbins[0]) bin_xj--;
				if(bin_yi == Nbins[1]) bin_yi--;
				if(bin_yj == Nbins[1]) bin_yj--;

				RDF[0][bin_xi][bin_yi] += 1.0;
				RDF[0][bin_xj][bin_yj] += 1.0;

				if(nRDFtypes > 1)
				{
					if(ATOMS[i].element == ATOMS[j].element)
					{
						if(ATOMS[i].element == 'O')
						{
							RDF[1][bin_xi][bin_yi] += 1.0;
							RDF[1][bin_xj][bin_yj] += 1.0;
						}

						else if(ATOMS[i].element == 'N')
						{
							RDF[4][bin_xi][bin_yi] += 1.0;
							RDF[4][bin_xj][bin_yj] += 1.0;
						}
					}

					else
					{
						RDF[2][bin_xi][bin_yi] += 1.0;
						RDF[3][bin_xj][bin_yj] += 1.0;
					}
				} 
			}
		}
	}
}

void analysis::Trajectory::write2file(float ***RDF, int nRDFtypes, int Nbins[])
{
	remove(fpathO);

	fileO = fopen(fpathO, "w");
	if(fileO == NULL)
	{
		printf("Couldn't create new data file %s. Exiting...\n", fpathO);
		exit(-1);
	}

	char **name = new char*[nRDFtypes];
	for(int i = 0; i < nRDFtypes; i++)
		name[i] = new char[10];

	sprintf(name[0], "gall");
	sprintf(name[1], "gAA");
	sprintf(name[2], "gAB");
	sprintf(name[3], "gBA");
	sprintf(name[4], "gBB");

	for(int i = 0; i < nRDFtypes; i++)
	{
		fprintf(fileO, "%s\n", name[i]);

		int ix, iy;
		for(ix = 0; ix < Nbins[0]; ix++)
		{
			for(iy = 0; iy < Nbins[1]; iy++)
				fprintf(fileO, "%f ", RDF[i][ix][iy]);

			if(iy == Nbins[1])
				fprintf(fileO, "\n");
		}
	}

	fclose(fileO);
}