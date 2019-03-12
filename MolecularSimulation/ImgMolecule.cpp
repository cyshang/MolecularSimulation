#include "ImgMolecule.h"
#include "Molecule.h"
#include "Atom.h"

ImgMolecule::ImgMolecule(const Molecule * molc, double offset[3])
	:type(molc->type)
{
	size_t nAtom = Molecule::nAtom[type];
	atomXYZ.resize(3 * nAtom);

	centroid[0] = molc->centroid[0] + offset[0];
	centroid[1] = molc->centroid[1] + offset[1];
	centroid[2] = molc->centroid[2] + offset[2];

	size_t pos = 0;
	for (size_t iAtom = 0; iAtom < nAtom; ++iAtom) {
		atomXYZ[pos] = molc->atoms[iAtom]->position[0] + offset[0];
		atomXYZ[pos + 1] = molc->atoms[iAtom]->position[1] + offset[1];
		atomXYZ[pos + 2] = molc->atoms[iAtom]->position[2] + offset[2];

		pos += 3;
	}
}