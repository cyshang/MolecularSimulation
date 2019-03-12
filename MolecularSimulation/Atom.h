#ifndef ATOM_H_
#define ATOM_H_

struct Atom
{
	
	double mass;
	double position[3];
	double velocity[3];
	double force[3];
	
	Atom();
};

#endif // !ATOM_H_

