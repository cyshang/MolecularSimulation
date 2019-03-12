#ifndef IMGMOLECULE_H_
#define IMGMOLECULE_H_

#include <vector>

struct Molecule;

struct ImgMolecule
{
	int type;
	std::vector<double> atomXYZ;
	double centroid[3];

	ImgMolecule(const Molecule * molc, double offset[3]);
};

#endif // !IMGMOLECULE_H_
