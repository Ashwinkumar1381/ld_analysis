// structure.cpp

#define PI 3.14159265359 
#include "analysis.h"

using namespace analysis;

float **computeStructureFactor(atom_style *ATOMS, Trajectory *TRAJ, float kx[], float ky[], float delk[]);

int main(int argc, char *argv[])
{
	// -------------- Trajectory Params --------------
	long eq_steps = long(1e7);
	long startStep = long(0e9), endStep = long(2e9);
	float dt = 5e-4;
	int frameW = int(1e5);

	// -------------- System Params --------------
	float Lx = 150, Ly = 30; 
	int nAtomTypes = 2;
	float delk[2] = {float(2*PI/Lx), float(2*PI/Ly)};

	int fac[2] = {75, 60};
	float kx[2] = {-fac[0] * delk[0], fac[0] * delk[0]};
	float ky[2] = {-fac[1] * delk[1], fac[1] * delk[1]};
	int Nk[2] = {int((kx[1] - kx[0])/delk[0]), int((ky[1] - ky[0])/delk[1])};

	Trajectory *TRAJ = new Trajectory(dt, frameW);
	sprintf(TRAJ->fpathI, "//media/ashwin/One Touch/ashwin_md/lane/Apr2025/lmp/Data27/traj3.xyz");
	sprintf(TRAJ->fpathO, "//media/ashwin/One Touch/ashwin_md/lane/Apr2025/lmp/Data27");
	TRAJ -> openTrajectory();

	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	System *BOX = new System(Lx, Ly, TRAJ->nAtoms, nAtomTypes);

	TRAJ -> write2file(NULL, kx, ky, Nk, 1);

	int ctr = 0;
	while( !feof(TRAJ -> fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if((TRAJ->step >= eq_steps + startStep) and (TRAJ->step <= eq_steps + endStep))
		{
			ctr++;
			printf("Processing step %ld, frame %d\n", TRAJ->step - eq_steps, ctr);
			computeStructureFactor(ATOMS, TRAJ, kx, ky, delk);
		}
	}

	TRAJ -> closeTrajectory();

	return(0);
}

float **computeStructureFactor(atom_style *ATOMS, Trajectory *TRAJ, float kx[], float ky[], float delk[])
{
	int Nk[2] = {int((kx[1] - kx[0])/delk[0]), int((ky[1] - ky[0])/delk[1])};

	float **Sk_x_y = new float*[2];
	for(int i = 0; i < 2; i++)
	{
		Sk_x_y[i] = new float[Nk[i] + 1];
		for(int j = 0; j <= Nk[i]; j++)
			Sk_x_y[i][j] = 0.0;
	}

	for(int i = 0; i <= Nk[0]; i++)
	{
		float kn = kx[0] + i*delk[0];
		float ax = 0.0, bx = 0.0;

		for(int j = 0; j < TRAJ->nAtoms; j++)
		{
			if(ATOMS[j].element == 'O')
			{
				float cx = kn * ATOMS[j].rxt1;
				ax += cos(cx);
				bx += sin(cx); 
			}
		}

		Sk_x_y[0][i] = 2 * (ax*ax + bx*bx) / TRAJ->nAtoms;
	}

	for(int i = 0; i <= Nk[1]; i++)
	{
		float kn = ky[0] + i*delk[1];
		float ay = 0.0, by = 0.0;

		for(int j = 0; j < TRAJ->nAtoms; j++)
		{
			if(ATOMS[j].element == 'O')
			{
				float cy = kn * ATOMS[j].ryt1;
			
				ay += cos(cy);
				by += sin(cy);
			}	
		}

		Sk_x_y[1][i] = 2 * (ay*ay + by*by) / TRAJ->nAtoms;		
	}

	TRAJ -> write2file(Sk_x_y, kx, ky, Nk);
	delete[] Sk_x_y; 
}

void analysis::Trajectory::write2file(float **Sk_x_y, float kx[], float ky[], int Nk[], int tag)
{
	char *fname_x = new char[500];
	char *fname_y = new char[500];

	sprintf(fname_x, "%s/Sk_x.dat", fpathO);
	sprintf(fname_y, "%s/Sk_y.dat", fpathO);

	if(tag == 1)
	{
		remove(fname_x);
		remove(fname_y);

		fileO = fopen(fname_x, "w");
		if(fileO == NULL)
		{
			printf("Could not create file %s. Exiting...\n", fname_x);
			exit(-1);
		}
		else
			fprintf(fileO, "time kmin %f kmax %f Nk %d\n", kx[0], kx[1], Nk[0] + 1);

		fclose(fileO);

		fileO = fopen(fname_y, "w");
		if(fileO == NULL)
		{
			printf("Could not create file %s. Exiting...\n", fname_y);
			exit(-1);
		}
		else
			fprintf(fileO, "time kmin %f kmax %f Nk %d\n", ky[0], ky[1], Nk[1] + 1);

		fclose(fileO);
	}
	else
	{
		int i = 0;

		fileO = fopen(fname_x, "a+");
		if(fileO == NULL)
		{
			printf("Could not open file %s. Exiting...\n", fname_x);
			exit(-1);
		}
		else
		{
			fprintf(fileO, "%ld ", step);

			for(i = 0; i <= Nk[0]; i++) fprintf(fileO, "%f ", Sk_x_y[0][i]);

			fprintf(fileO, "\n");
		}

		fclose(fileO);

		fileO = fopen(fname_y, "a+");
		if(fileO == NULL)
		{
			printf("Could not open file %s. Exiting...\n", fname_y);
			exit(-1);
		}
		else
		{
			fprintf(fileO, "%ld ", step);

			for(i = 0; i <= Nk[1]; i++) fprintf(fileO, "%f ", Sk_x_y[1][i]);

			fprintf(fileO, "\n");
		}

		fclose(fileO);	
	}
}