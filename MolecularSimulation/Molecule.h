#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <vector>

struct Atom;
struct ImgMolecule;

struct Molecule
{
	struct NeighborMolecule;
	static std::vector<size_t>	nAtom;
	static std::vector<int>		cAtom;

	int type;	
	std::vector<Atom*> atoms;
	std::vector<NeighborMolecule> neighbors;

	double centroid[3];

	Molecule(const int & type, const size_t & offset, std::vector<Atom> & atom_list);
	void CalCentroid();
	void Translation(const double trans[3]);
};

struct Molecule::NeighborMolecule
{	
	int type;
	double Rij;
	std::vector<double> atomXYZ;

	NeighborMolecule(const Molecule & molc, const double & distance);
	NeighborMolecule(const ImgMolecule & molc, const double & distance);
};

#endif // !MOLECULE_H_

