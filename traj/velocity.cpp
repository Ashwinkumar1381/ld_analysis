// velocity.cpp

#include "analysis.h"

using namespace analysis;

int main(int argc, char *argv[])
{
	// ------------ Trajectory Params ------------
	float dt = 5e-4;
	int frameW = int(1e5);
	int frameStart = int(4e4), frameEnd = int(5e4);

	// ------------ System Params ------------
	int nAtomTypes = 2;

	// ------------ Bins Params ------------
	int nBins = 100;

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/Expansion/ashwin_md/lane/June_July2025/Pe90/Data44/traj2.cfg");
	sprintf(TRAJ->fpathO, "//media/ashwin/Expansion/ashwin_md/lane/June_July2025/Pe90/Data44/velocityDist.dat");
	TRAJ -> openTrajectory();

	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	velocityDist *Dist = new velocityDist(2, nBins); 

	printf("\nScanning particle velocities to set distribution bounds...\n");

	while( !feof(TRAJ->fileI) )
	{
		if(TRAJ->frame_nr == frameStart - 1)
			fgetpos(TRAJ->fileI, &TRAJ->curr_pos);

		TRAJ -> readThisFrame(ATOMS);

		if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
			Dist -> scanVelocities(ATOMS, TRAJ->nAtoms);

		if(TRAJ->frame_nr == frameEnd)
		{
			TRAJ->frame_nr = frameStart - 1;
			fsetpos(TRAJ->fileI, &TRAJ->curr_pos);
			break;
		}
	}

	printf("\nIdentified bounds. Creating velocity distribution...\n");

	Dist -> createBins();

	printf("\nBinning velocities...\n");

	while ( !feof(TRAJ->fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if(TRAJ->frame_nr <= frameEnd)
			Dist -> binVelocities(ATOMS, TRAJ->nAtoms);

		if(TRAJ->frame_nr == frameEnd)
			break;
	}

	printf("\nVelocities binned for vx:%d frames, vy:%d frames.\n", int(Dist->ctr[0]/TRAJ->nAtoms), int(Dist->ctr[1]/TRAJ->nAtoms));

	Dist -> normalize();

	TRAJ -> write2file(Dist);
	TRAJ -> closeTrajectory();

	delete TRAJ;
	delete Dist;
	delete[] ATOMS;
}

analysis::velocityDist::velocityDist(int nDims, int nBins)
{
	this->nDims = nDims;
	this->nBins = nBins;

	Bins = new float*[nDims];
	min_val = new float[nDims];
	max_val = new float[nDims];
	binW = new float[nDims];
	ctr = new long[nDims];

	for(int i = 0; i < nDims; i++)
	{
		min_val[i] = 1e5;
		max_val[i] = 0.0;
		binW[i] = 0.0;
		ctr[i] = 0;
	}
}

void analysis::velocityDist::scanVelocities(atom_style *ATOMS, int nAtoms)
{
	for(int i = 0; i < nAtoms; i++)
	{
		if(ATOMS[i].vx > max_val[0]) max_val[0] = ATOMS[i].vx;
		if(ATOMS[i].vy > max_val[1]) max_val[1] = ATOMS[i].vy;

		if(ATOMS[i].vx < min_val[0]) min_val[0] = ATOMS[i].vx;
		if(ATOMS[i].vy < min_val[1]) min_val[1] = ATOMS[i].vy;
	}
}

void analysis::velocityDist::createBins()
{	
	for(int i = 0; i < nDims; i++)
	{
		max_val[i] = bounds(max_val[i]);
		min_val[i] = bounds(min_val[i]);

		if(abs(max_val[i]) < abs(min_val[i])) max_val[i] = abs(min_val[i]); 

		binW[i] = 2 * max_val[i] / nBins;

		Bins[i] = new float[nBins];
		for(int j = 0; j < nBins; j++)
			Bins[i][j] = 0.0;
	}
}

void analysis::velocityDist::binVelocities(atom_style *ATOMS, int nAtoms)
{
	for(int i = 0; i < nAtoms; i++)
	{
		int bin_vx = int((ATOMS[i].vx + max_val[0]) / binW[0]);
		int bin_vy = int((ATOMS[i].vy + max_val[1]) / binW[1]);

		Bins[0][bin_vx] += 1.0;
		Bins[1][bin_vy] += 1.0;

		ctr[0] += 1;
		ctr[1] += 1;
	}
}

void analysis::velocityDist::normalize()
{
	for(int i = 0; i < nBins; i++)
	{
		Bins[0][i] /= (ctr[0] * binW[0]);
		Bins[1][i] /= (ctr[1] * binW[1]);
	}
}

void analysis::Trajectory::write2file(velocityDist *Dist)
{
	char string[100];
	sprintf(string, "bin vx_dist vy_dist\nVxmax %g Vymax %g VxbinW %g VybinW %g", Dist->max_val[0], Dist->max_val[1], Dist->binW[0], Dist->binW[1]);

	createOutputFile(string);

	for(int i = 0; i < Dist->nBins; i++)
		fprintf(fileO, "\n%d %g %g", i + 1, Dist->Bins[0][i], Dist->Bins[1][i]);
}