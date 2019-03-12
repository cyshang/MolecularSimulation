#include "Molecule.h"
#include "Atom.h"
#include "ImgMolecule.h"
#include <cstring>

Molecule::Molecule(const int & molc_type, const size_t & offset, std::vector<Atom> & atom_list)
	:type(molc_type)
{
	for (size_t iAtom = 0; iAtom < nAtom[type]; ++iAtom) {
		atoms.push_back(&(atom_list[offset + iAtom]));
	}
}

void Molecule::CalCentroid()
{
	if (cAtom[type] < 0 || cAtom[type] > nAtom[type] - 1) {

		centroid[0] = 0;
		centroid[1] = 0;
		centroid[2] = 0;

		double totMass = 0;
		double mass = 0;
		for (int iAtom = 0; iAtom < nAtom[type]; ++iAtom) {
			mass = atoms[iAtom]->mass;
			centroid[0] += mass * atoms[iAtom]->position[0];
			centroid[1] += mass * atoms[iAtom]->position[1];
			centroid[2] += mass * atoms[iAtom]->position[2];
			totMass += mass;
		}

		centroid[0] /= totMass;
		centroid[1] /= totMass;
		centroid[2] /= totMass;
	}
	else {
		int center = cAtom[type];
		centroid[0] = atoms[center]->position[0];
		centroid[1] = atoms[center]->position[1];
		centroid[2] = atoms[center]->position[2];
	}
}

void Molecule::Translation(const double trans[3])
{
	for (int iAtom = 0; iAtom < nAtom[type]; ++iAtom) {
		atoms[iAtom]->position[0] += trans[0];
		atoms[iAtom]->position[1] += trans[1];
		atoms[iAtom]->position[2] += trans[2];
	}

}

Molecule::NeighborMolecule::NeighborMolecule(const Molecule & molc, const double & distance)
	:type(molc.type), Rij(distance)
{
	size_t nAtom = Molecule::nAtom[type];
	atomXYZ.resize(3 * nAtom);

	size_t pos = 0;
	for (size_t iAtom = 0; iAtom < nAtom; ++iAtom) {
		atomXYZ[pos] = molc.atoms[iAtom]->position[0];
		atomXYZ[pos + 1] = molc.atoms[iAtom]->position[1];
		atomXYZ[pos + 2] = molc.atoms[iAtom]->position[2];

		pos += 3;
	}
}

Molecule::NeighborMolecule::NeighborMolecule(const ImgMolecule & molc, const double & distance)
	:type(molc.type), Rij(distance)
{
	atomXYZ = molc.atomXYZ;
}