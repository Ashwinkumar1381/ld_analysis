// analysis.h

#ifndef ANALYSIS_H
#define ANALYSIS_H

#define MAXCELL 9000000
#define rcut 1.122462048

#include "library.h"

namespace analysis {

	typedef class atomsXYZ {

	public:

	char id;
	int type;
	float rxt1, ryt1, rzt1; 		// *) Used in MSD calculations 
	float rxt2, ryt2, rzt2; 
	float vx, vy, vz;
	float dx, dy, dz;
	int jumpx, jumpy;

	atomsXYZ();
	~atomsXYZ();
	
	} atom_style;

	class System {

	public:

	float Lx, Ly, rcellx, rcelly;
	int Ncellx, Ncelly, ncells;
	int nAtoms;

	int MAPS[MAXCELL], HEAD[MAXCELL], LIST[MAXCELL];

	System(float Lx, float Ly, float nAtoms, float rcellx = rcut, float rcelly = rcut);
	~System();

	int cellindex(int ix, int iy);
	void buildCellMaps();
	void buildCellList(atom_style *ATOMS);
	void checkMinImage(float *dx = NULL, float *dy = NULL);

	};

	class Bin1D {

	public:

	float Lx, Ly, binWidth;
	int nBins;
	float *bin;
	float *binavg;

	Bin1D(float Lx, float Ly, float binWidth);
	Bin1D();
	~Bin1D();

	void addToBin(float rx, float val = 1.0);
	void normalize(float *binTonorm, float *normv = NULL, float normc = 1.0);
	void zero(float *binTozero);
	void addBins(float *bin1, float *bin2);

	};

	class Hist2D {

		public:

		float rcut_x, rcut_y;
		float binW_x, binW_y;
		int nBin_x, nBin_y;
		float **bin;

		Hist2D(float rcut_x, float rcut_y, float binW_x, float binW_y);
		~Hist2D();

		void addToBin(char iid, char jid, float dxij, float dyij);
		void normalize(float fac = 1.0);
	};

	class Trajectory {

	public:

	long step; 
	int nAtoms, nAtomTypes;
	int frame_nr, totalFrames, frameWidth;
	float timeStep, time;
	float xCom, yCom, zCom; 
	char *format;

	char *fpathI, *fpathO, *pipeString, *pipeChar;
	FILE *fileI, *fileO;

	Trajectory(float timeStep = 1.0, int frameWidth = 1, char fileFormat[10] = "xyz");
	~Trajectory();

	void openTrajectory(bool count = false);
	void closeTrajectory();
	void loadTrajectory(atom_style **ATOMS, System *BOX, int frameStart, int frameEnd);
	void countFrames();
	void readThisFrame(atom_style *ATOMS);
	void readNextFrame(atom_style *ATOMS);
	void writeThisFrame(atom_style *ATOMS, long add_step = 0);
	void copyThisFrame(atom_style *From, atom_style *To);
	void computeCom(atom_style *ATOMS);

	void write2file();
	void write2file(Bin1D *bin1, int ctr = 0);
	void write2file(Bin1D *binA, Bin1D *binB, int ctr = 0, float timeToAvg = 0.0);
	void write2file(Bin1D **Pin, Bin1D **Pkin, Bin1D *Pswim, int ctr = 0);
	void write2file(float **MSD, int *delFrames, int nSample);
	void write2file(float time, float order, int step = -1);
	void write2file(float *rn, float **RDF, int nRDF, int nBins);
	void write2file(float ***RDF_x_y, int nRDF, int nBins[], float binW[], char option[]);
	void write2file(float ***RDFxy, int nRDFtypes, int nBins[]);
	void write2file(float **Sk_x_y, float kx[], float ky[], int Nk[], int tag = 0);
	void write2file(float **Hist, int nBins[]);
	};
}

#endif /*ANALYSIS_H*/