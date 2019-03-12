#include "VelocityVerlet.h"
#include "Simulator.h"
#include "Atom.h"

VelocityVerlet::VelocityVerlet() {}

void VelocityVerlet::Init()
{
	last_force.resize(3 * Simulator::totAtom);
}

void VelocityVerlet::FirstStep(std::vector<Atom> & atoms)
{
	double tmp = 0;
	for (size_t iAtom = 0; iAtom < Simulator::totAtom; ++iAtom) {
		Atom & atom = atoms[iAtom];
		size_t pos = iAtom * 3;

		tmp = 0.5 * time_step * time_step / atom.mass;
		atom.position[0] += atom.velocity[0] * time_step + atom.force[0] * tmp;
		atom.position[1] += atom.velocity[1] * time_step + atom.force[1] * tmp;
		atom.position[2] += atom.velocity[2] * time_step + atom.force[2] * tmp;

		last_force[pos] = atom.force[0];
		last_force[pos + 1] = atom.force[1];
		last_force[pos + 2] = atom.force[2];
	}
}

void VelocityVerlet::Integrate(std::vector<Atom> & atoms)
{	
	double tmp = 0;
	for (size_t iAtom = 0; iAtom < Simulator::totAtom; ++iAtom) {
		Atom & atom = atoms[iAtom];
		size_t pos = iAtom * 3;

		tmp = time_step * 0.5 / atom.mass;
		atom.velocity[0] += (atom.force[0] + last_force[pos]) * tmp;
		atom.velocity[1] += (atom.force[1] + last_force[pos + 1]) * tmp;
		atom.velocity[2] += (atom.force[2] + last_force[pos + 2]) * tmp;

		tmp *= time_step;  // tmp = time_step**2 / (2 * mass)
		atom.position[0] += atom.velocity[0] * time_step + atom.force[0] * tmp;
		atom.position[1] += atom.velocity[1] * time_step + atom.force[1] * tmp;
		atom.position[2] += atom.velocity[2] * time_step + atom.force[2] * tmp;

		last_force[pos] = atom.force[0];
		last_force[pos + 1] = atom.force[1];
		last_force[pos + 2] = atom.force[2];
	}
}