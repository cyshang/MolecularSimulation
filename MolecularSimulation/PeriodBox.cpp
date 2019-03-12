#include <cmath>
#include "PeriodBox.h"

void PeriodBox::Init()
{
	double half_period[3];

	for (int i = 0; i < 3; ++i) {
		half_period[i] = 0.5 * period[i];
	}

	int iBox = 0;	
	for (int k = -1; k < 3; ++k) {
		for (int j = -1; j < 3; ++j) {
			for (int i = -1; i < 3; ++i) {
				if ((i == 0 || i == 1) && (j == 0 || j == 1) && (k == 0 || k == 1))
					continue;

				int x = (i == -1) ? 1 : ((i == 2) ? 0 : i);
				int y = (j == -1) ? 1 : ((j == 2) ? 0 : j);
				int z = (k == -1) ? 1 : ((k == 2) ? 0 : k);

				imgBox[iBox].octant = &(octant[x + (y << 1) + (z << 2)]);
				imgBox[iBox].offset[0] = (i - x) * half_period[0];
				imgBox[iBox].offset[1] = (j - y) * half_period[1];
				imgBox[iBox].offset[2] = (k - z) * half_period[2];

				for (int ii = -1; ii < 2; ++ii) {
					for (int jj = -1; jj < 2; ++jj) {
						for (int kk = -1; kk < 2; ++kk) {
							x = i + ii;
							y = j + jj;
							z = k + kk;
							if ((x == 0 || x == 1) && (y == 0 || y == 1) && (z == 0 || z == 1)) {
								octant[x + (y << 1) + (z << 2)].imgBox.push_back(&(imgBox[iBox]));
							}
						}
					}
				}

				++iBox;
			}
		}
	}
}

void PeriodBox::MolcPartition()
{
	double	trans[3];
	double	pos[3];
	int		oct[3];
	bool	IfTrans;

	for (int i = 0; i < 8; ++i) {
		octant[i].molc.clear();
	}

	for (size_t iMolc = 0; iMolc < molc.size(); ++iMolc) {

		IfTrans = false;
		
		for (int i = 0; i < 3; ++i) {
			pos[i] = molc[iMolc]->centroid[0];

			if (pos[i] < lattice[i][0]) {
				trans[i] = period[i];
				IfTrans = true;
				oct[i] = 1;
			}
			else if (pos[i] < lattice[i][1]) {
				oct[i] = 0;
			}
			else if (pos[i] < lattice[i][2]) {
				oct[i] = 1;
			}
			else {
				trans[i] = -period[i];
				IfTrans = true;
				oct[i] = 0;
			}
		}

		if (IfTrans) {
			molc[iMolc]->Translation(trans);
		}

		octant[oct[0] + (oct[1] << 1) + (oct[2] << 2)].molc.push_back(molc[iMolc]);
	}
}

void PeriodBox::BuildNeighborList()
{
	for (size_t i = 0; i < molc.size(); ++i) {
		molc[i]->neighbors.clear();
	}

	double Rij = 0;
	double tmp = 0;
	for (size_t i = 0; i < molc.size() - 1; ++i) {
		for (size_t j = i + 1; j < molc.size(); ++j) {
			Rij = 0;

			tmp = (molc[i]->centroid[0] - molc[j]->centroid[0]);
			Rij += tmp * tmp;
			tmp = (molc[i]->centroid[1] - molc[j]->centroid[1]);
			Rij += tmp * tmp;
			tmp = (molc[i]->centroid[2] - molc[j]->centroid[2]);
			Rij += tmp * tmp;

			Rij = sqrt(Rij);

			if (Rij < period_cutoff) {
				molc[i]->neighbors.emplace_back(*(molc[j]), Rij);
				molc[j]->neighbors.emplace_back(*(molc[i]), Rij);
			}
		}
	}

	for (int iBox = 0; iBox < 8; ++iBox) {
		
		size_t nMolc = octant[iBox].molc.size();
		for (size_t iMolc = 0; iMolc < nMolc; ++iMolc) {
			for (int jBox = 0; jBox < 19; ++jBox) {
				
				size_t mMolc = imgBox[jBox].imgMolc.size();
				for (size_t jMolc = 0; jMolc < mMolc; ++jMolc) {
					
					Molecule & molcA = *(octant[iBox].molc[iMolc]);
					ImgMolecule & molcB = imgBox[jBox].imgMolc[jMolc];

					Rij = 0;
					tmp = molcA.centroid[0] - molcB.centroid[0];					
					Rij += tmp * tmp;
					tmp = molcA.centroid[1] - molcB.centroid[1];
					Rij += tmp * tmp;
					tmp = molcA.centroid[2] - molcB.centroid[2];
					Rij += tmp * tmp;
					Rij = sqrt(Rij);

					if (Rij < period_cutoff) {
						molcA.neighbors.emplace_back(molcB, Rij);
					}
				}
			}
		}
	}
}

void PeriodBox::ImgBox::BuildImgMolc() 
{
	imgMolc.clear();

	for (size_t i = 0; i < octant->molc.size(); ++i) {
		imgMolc.emplace_back(octant->molc[i], offset);
	}
}