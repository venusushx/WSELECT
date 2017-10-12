#ifndef _WSDECLARE_
#define _WSDECLARE_

#include <iostream>
#include <exception>
#include <stdexcept>
#include <string>
#include <map>
#include <vector>

class PDBSpecs
{
public:
	std::string atm_name;
	std::string res_name;
	int res_number;
	std::vector<double> coord;
	std::string element;
	double hyd_c;
public:
	PDBSpecs();
	PDBSpecs(const PDBSpecs &x);
	PDBSpecs(std::vector<double> coord);
	PDBSpecs& operator=(const PDBSpecs& rhs);
	PDBSpecs& operator=(std::vector<double> coord);
	~PDBSpecs();
	double dis_sqrd(const PDBSpecs &x);
	void pdb_read(std::string buffer);
	bool close_to_backbone(int x, int y, std::vector<PDBSpecs> &p, double cutoff);
	PDBSpecs find_mid(const PDBSpecs &x);
	double mid_rad_sqrd(const PDBSpecs &x, const PDBSpecs &y);
};

double hb_cosine(double dis_sqrd_1, double dis_sqrd_2, double dis_sqrd_3);
double gaus_mod(double height, double dis, double u, double sigma);
double rmsd(double *dis_sqrd, int dis_sqrd_size);

void ws_help();
bool ascendSorting(double p, double q);
bool descendSorting(double p, double q);

extern const std::map<std::string, double> hydp_map;

#endif
