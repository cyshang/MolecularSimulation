#ifndef SIMULATOR_H_
#define SIMULATOR_H_

#include <vector>
#include "Container.h"
#include "Integrator.h"
#include "Atom.h"
#include "Molecule.h"

struct Integrator;
struct Container;

struct Simulator
{
	static size_t totAtom;
	static size_t totMolecule;

	std::vector<Atom> atom_list;
	std::vector<Molecule> molecule_list;

	Container *container;
	Integrator *integrator;

	Simulator();
	~Simulator();
};

#endif // !SIMULATOR_H_

