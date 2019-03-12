#ifndef VELOCITYVERLET_H_
#define VELOCITYVERLET_H_

#include <vector>
#include "Integrator.h"

struct Atom;

struct VelocityVerlet :public Integrator
{
	static double time_step;
	std::vector<double> last_force;

	VelocityVerlet();
	void Init();
	void FirstStep(std::vector<Atom> &);
	void Integrate(std::vector<Atom> &);
};

#endif // !INTEGRATOR_H_

