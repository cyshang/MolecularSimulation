#ifndef PERIODBOX_H_
#define PERIODBOX_H_

#include <vector>
#include <array>
#include "Container.h"
#include "Atom.h"
#include "Molecule.h"
#include "ImgMolecule.h"

struct PeriodBox :public Container
{
	struct Octant;
	struct ImgBox;

	double lattice[3][3];				// lattice = [[x_left, x_mid, x_right],
										//            [y_left, y_mid, y_right],
										//            [z_left, z_mid, z_right]]

	double period[3];					// period = {X, Y, Z}
	double period_cutoff;


	std::array<Octant, 8>	octant;
	std::array<ImgBox, 56>	imgBox;
	std::vector<Molecule*>	molc;

	void Init();
	void MolcPartition();
	void BuildNeighborList();

};

struct PeriodBox::Octant
{
	std::vector<Molecule*>	molc;
	std::vector<ImgBox*>	imgBox;
};

struct PeriodBox::ImgBox
{
	const Octant*	octant;
	double			offset[3];

	std::vector<ImgMolecule> imgMolc;

	void BuildImgMolc();
};

#endif // !PERIODBOX_H_

