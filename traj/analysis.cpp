// analysis.cpp

#include "analysis.h"

using namespace analysis;

/* ----------------- atomsXYZ members -----------------*/

analysis::atomsXYZ::atomsXYZ()
{
	rxt1 = ryt1 = rzt1 = 0.0;
	jumpx = jumpy = 0;
}

analysis::atomsXYZ::~atomsXYZ(){}

/* ----------------- System members -----------------*/

analysis::System::System(float Lx, float Ly, float nAtoms, float rcellx, float rcelly)
{
	this -> Lx = Lx;
	this -> Ly = Ly;
	this -> nAtoms = nAtoms;
	this -> rcellx = rcellx;
	this -> rcelly = rcelly;

	Ncellx = int(this->Lx/this->rcellx);
	Ncelly = int(this->Ly/this->rcelly);
	ncells = int(Ncellx*Ncelly);

	this->rcellx = this->Lx/Ncellx;
	this->rcelly = this->Ly/Ncelly; 

	for(int i = 0; i < MAXCELL; i++)
	{
		MAPS[i] = 0;
		LIST[i] = 0;
		HEAD[i] = 0;
	}

	//buildCellMaps();
}

analysis::System::~System(){}

int analysis::System::cellindex(int ix, int iy)
{   
    if (ix >= Ncellx) ix = ix - Ncellx; 
	else if (ix <= -1) ix = ix + Ncellx; 
	if (iy >= Ncelly) iy = iy - Ncelly;
	else if (iy <= -1) iy = iy + Ncelly;

	if(Ncellx >= Ncelly)
		return (1 + ix + iy*Ncellx);
	else
		return (1 + ix + iy*Ncelly);
}

void analysis::System::buildCellMaps()
{
	for (int ix = 0; ix < Ncellx; ix++)
	{
		for (int iy = 0; iy < Ncelly; iy++)
		{
			int imap = 4*(cellindex(ix,iy)-1);
			MAPS[imap+1] = cellindex(ix-1,iy);
			MAPS[imap+2] = cellindex(ix-1,iy+1);
			MAPS[imap+3] = cellindex(ix,iy+1);
			MAPS[imap+4] = cellindex(ix+1,iy+1);
		}
	}
    
    printf("Successfully constructed MAPS array with %d cells.\n", int(Ncellx * Ncelly));
}

void analysis::System::buildCellList(atom_style *ATOMS)
{
	for(int i = 0; i < MAXCELL; i++)
	{
		LIST[i] = 0;
		HEAD[i] = 0;
	}

	for (int i = 1; i <= nAtoms; i++)
	{
		int ii = i - 1;
		int ix = int(ATOMS[ii].rxt1/rcellx);
		int iy = int(ATOMS[ii].ryt1/rcelly);
		
		int icell = cellindex(ix, iy);	        
		LIST[i] = HEAD[icell];
		HEAD[icell] = i;		    
	}	
}

void analysis::System::checkMinImage(float *dx, float *dy)
{
	if(dx != NULL)
	{
		if(*dx >= 0.5*Lx) *dx -= Lx;
		else if(*dx <= -0.5*Lx) *dx += Lx;
	}
	if(dy != NULL)
	{
		if(*dy >= 0.5*Ly) *dy -= Ly;
		else if(*dy <= -0.5*Ly) *dy += Ly; 
	}
}

/* ----------------- Bin1D members -----------------*/

analysis::Bin1D::Bin1D(float Lx, float Ly, float binWidth)
{
	this -> Lx = Lx;
	this -> Ly = Ly;
	this -> binWidth = binWidth;

	nBins = int(Lx/binWidth);

	bin = new float[nBins];
	binavg = new float[nBins];

	zero(bin);
	zero(binavg);
}

analysis::Bin1D::~Bin1D(){}

void analysis::Bin1D::addToBin(float rx, float val)
{
	if(val == 1.0)
		bin[int(rx/binWidth)] += 1;
	else
		bin[int(rx/binWidth)] += val;	
}

void analysis::Bin1D::normalize(float *binTonorm, float *normv, float normc)
{
	if(normv != NULL)
	{
		for(int i = 0; i < nBins; i++)
			binTonorm[i] /= normv[i];
	}
	else
	{
		if(normc == 1.0) normc = Ly*binWidth;

		for(int i = 0; i < nBins; i++)
			binTonorm[i] /= normc;	
	}
}

void analysis::Bin1D::zero(float *binTozero)
{
	for(int i = 0; i < nBins; i++)
		binTozero[i] = 0.0;
}

void analysis::Bin1D::addBins(float *bin1, float *bin2)
{
	for(int i = 0; i < nBins; i++)
		bin1[i] += bin2[i];
}

/* ----------------- Hist2D members -----------------*/

analysis::Hist2D::Hist2D(float rcut_x, float rcut_y, float binW_x, float binW_y)
{
	this -> rcut_x = rcut_x;
	this -> rcut_y = rcut_y;
	this -> binW_x = binW_x;
	this -> binW_y = binW_y;

	nBin_x = 2*int(rcut_x / binW_x);
	nBin_y = 2*int(rcut_y / binW_y);

	bin = new float*[nBin_x];
	for(int i = 0; i < nBin_x; i++)
	{
		bin[i] = new float[nBin_y];

		for(int j = 0; j < nBin_y; j++)
			bin[i][j] = 0.0;
	}
} 

analysis::Hist2D::~Hist2D(){}

void analysis::Hist2D::addToBin(char iid, char jid, float dxij, float dyij)
{
	int bin_xi = int((rcut_x + dxij)/binW_x);
	int bin_xj = int((rcut_x - dxij)/binW_x);
	int bin_yi = int((rcut_y + dyij)/binW_y);
	int bin_yj = int((rcut_y - dyij)/binW_y);

	if(bin_xi == nBin_x) bin_xi--;
	if(bin_xj == nBin_x) bin_xj--;
	if(bin_yi == nBin_y) bin_yi--;
	if(bin_yj == nBin_y) bin_yj--;

	if(iid == jid)
	{
		bin[bin_xi][bin_yi] += 1.0;
		bin[bin_xj][bin_yj] += 1.0;
	}
	else
		bin[bin_xi][bin_yi] += 1.0;	
}

void analysis::Hist2D::normalize(float fac)
{
	for(int i = 0; i < nBin_x; i++)
	{
		for(int j = 0; j < nBin_y; j++)
			bin[i][j] /= (fac * binW_x * binW_y);
	}
}

/* ----------------- Trajectory members -----------------*/

analysis::Trajectory::Trajectory(float timeStep, int frameWidth)
{
	fpathI = new char [500];
	fpathO = new char [500];
	pipeString = new char [500];
	pipeChar = new char [500];
	fileI = nullptr;
	fileO = NULL;

	frame_nr = -1;
	this->timeStep = timeStep;
	this->frameWidth = frameWidth;
	xCom = yCom = zCom = 0.0;
}

analysis::Trajectory::~Trajectory(){}

void analysis::Trajectory::openTrajectory(bool count)
{
	fileI = fopen(fpathI, "r");
	if(fileI == NULL)
	{
		printf("Could not open %s. Exiting ...\n", fpathI);
		exit(-1);
	}
	else
	{
		fgets(pipeString, 500, fileI);
		sscanf(pipeString, "%d", &nAtoms);
		rewind(fileI);
		printf("\nTrajectory %s file opened and ready to be read...\n", fpathI);
	}

	if(count == true)
		countFrames();
}

void analysis::Trajectory::closeTrajectory()
{
	printf("\nClosing trajectory file.\n");
	fclose(fileI);
}

void analysis::Trajectory::importTrajectory(atom_style **ATOMS, System *BOX, int frameStart, int frameEnd)
{
	atom_style *tempATOMS = new atom_style();

	frame_nr = -1;
	while( !feof(fileI) )
	{
		readThisFrame(tempATOMS);

		if(frame_nr >= frameStart and frame_nr <= frameEnd)
		{
			int new_frame = frame_nr - frameStart;

			ATOMS[new_frame] = tempATOMS;

			if(new_frame == 0) 
				printf("First frame: Step %ld\n", step);

			if(new_frame > 0)
			{
				for(int i = 0; i < nAtoms; i++)
				{
					ATOMS[new_frame][i].jumpx = ATOMS[new_frame - 1][i].jumpx;
					ATOMS[new_frame][i].jumpy = ATOMS[new_frame - 1][i].jumpy;

					float dx = ATOMS[new_frame][i].rxt1 - ATOMS[new_frame - 1][i].rxt1;
					float dy = ATOMS[new_frame][i].ryt1 - ATOMS[new_frame - 1][i].ryt1;

					if(dx <= -0.5*BOX->Lx) ATOMS[new_frame][i].jumpx++;
					else if(dx >= 0.5*BOX->Lx) ATOMS[new_frame][i].jumpx--;

					if(dy <= -0.5*BOX->Ly) ATOMS[new_frame][i].jumpy++;
					else if(dy >= 0.5*BOX->Ly) ATOMS[new_frame][i].jumpy--;
				}
			}

			if(new_frame == frameEnd - frameStart)
			{
				printf("Last frame: Step %ld\n", step);
				break;
			} 
		}
	}

	printf("\nCoordinates imported successfully!\n");
	rewind(fileI);

	delete tempATOMS;
}

void analysis::Trajectory::countFrames()
{
	atom_style *ATOMS = new atom_style [nAtoms];

	frame_nr = 0;
	while( !(feof(fileI)) )
		readThisFrame(ATOMS);

	totalFrames = frame_nr - 1;

	printf("Counting %d frames, %d atoms\n", totalFrames, nAtoms);

	frame_nr = -1;

	rewind(fileI);
	delete[] ATOMS;
}

void analysis::Trajectory::readThisFrame(atom_style *ATOMS)
{
	fgets(pipeString, 500, fileI);
	
	fgets(pipeString, 500, fileI);
	// sscanf(pipeString, "%*s %d", &step);
	sscanf(pipeString, "%*s %*s %ld", &step);

	for(int i = 0; i < nAtoms; i++)
	{
		fgets(pipeString, 500, fileI);
		sscanf(pipeString, "%c %f %f %*f %*f %*f", &ATOMS[i].id, &ATOMS[i].rxt1, &ATOMS[i].ryt1);
	}

	frame_nr++;
}

void analysis::Trajectory::computeCom(atom_style *ATOMS)
{
	xCom = yCom = zCom = 0.0;

	for(int i = 0; i < nAtoms; i++)
	{
		xCom += ATOMS[i].rxt1;
		yCom += ATOMS[i].ryt1;
		zCom += ATOMS[i].rzt1;
	}	

	xCom /= nAtoms;
	yCom /= nAtoms;
	zCom /= nAtoms;
}