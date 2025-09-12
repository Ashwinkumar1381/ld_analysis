// computeInteractions.cpp

#include "analysis.h"

using namespace program;
using namespace analysis;

namespace analysis {

	float computeNonBondedInteractions(atom_style *ATOMS, System *BOX, WCA_2P *INTERACTIONS);
}

float analysis::computeNonBondedInteractions(atom_style *ATOMS, System *BOX, WCA_2P *INTERACTIONS)
{
	float rcut2 = INTERACTIONS->rcut * INTERACTIONS->rcut;

	for(int i = 0; i < BOX->nAtoms; i++)
	{
		ATOMS[i].fx_int = 0.0;
		ATOMS[i].fy_int = 0.0;
	}

	float pe = 0.0;
	
	BOX -> buildCellList(ATOMS);

	for(int icell = 1; icell <= BOX->ncells; icell++)
	{
		int i = BOX->HEAD[icell];

		while(i != 0)
		{
			int ii = i - 1;
			float rxi = ATOMS[ii].rxt1;
			float ryi = ATOMS[ii].ryt1;

			int j = BOX->LIST[i];
			while(j != 0)
			{
				int jj = j - 1;
				float dxij = ATOMS[jj].rxt1 - rxi;
				float dyij = ATOMS[jj].ryt1 - ryi;

				float r2ij = dxij*dxij + dyij*dyij;

				if(r2ij <= rcut2)
				{
					float *pairs = INTERACTIONS -> get_forces(r2ij);

					pe += pairs[0];

					ATOMS[ii].fx_int += dxij * pairs[1];
					ATOMS[jj].fx_int += -dxij * pairs[1];
					ATOMS[ii].fy_int += dyij * pairs[1];
					ATOMS[jj].fy_int += -dyij * pairs[1];

					delete[] pairs;	
				}

				j = BOX->LIST[j];
			}

			i = BOX->LIST[i];
		} 
	}

	int nNbors = 4;
	for(int icell = 1; icell <= BOX->ncells; icell++)
	{
		int icell_index = nNbors*(icell - 1);

		int i = BOX->HEAD[icell];
		while(i != 0)
		{
			int ii = i - 1;
			float rxi = ATOMS[ii].rxt1;
			float ryi = ATOMS[ii].ryt1;

			for(int nbor = 1; nbor <= nNbors; nbor++)
			{
				int jcell = BOX->MAPS[icell_index + nbor];

				int j = BOX->HEAD[jcell];
				while(j != 0)
				{
					int jj = j - 1;

					float dxij = rxi - ATOMS[jj].rxt1;
					float dyij = ryi - ATOMS[jj].ryt1;

					BOX -> checkMinImage(&dxij, &dyij);

					float r2ij = dxij*dxij + dyij*dyij;

					if(r2ij <= rcut2)
					{
						float *pairs = INTERACTIONS -> get_forces(r2ij);

						pe += pairs[0];

						ATOMS[ii].fx_int += dxij * pairs[1];
						ATOMS[jj].fx_int += -dxij * pairs[1];
						ATOMS[ii].fy_int += dyij * pairs[1];
						ATOMS[jj].fy_int += -dyij * pairs[1];

						delete[] pairs;		
					}

					j = BOX->LIST[j];
				}
			}

			i = BOX->LIST[i];
		}
	}

	return(pe);
}