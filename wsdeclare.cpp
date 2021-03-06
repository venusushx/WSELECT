#include "wsdeclare.h"
#include "general.h"
#include <math.h>

/* PDBSpecs class constructors, operators, and destructor */
PDBSpecs::PDBSpecs():coord(3),hyd_c(0.0) {};
PDBSpecs::PDBSpecs(const PDBSpecs &x)  
{
	this->atm_name=x.atm_name;
	this->res_name=x.res_name;
	this->res_number=x.res_number;
	this->coord=x.coord;
	this->element=x.element;
	this->hyd_c=x.hyd_c;
}
PDBSpecs::PDBSpecs(std::vector<double> coord)
{
	this->coord=coord;
}

PDBSpecs& PDBSpecs::operator=(const PDBSpecs& rhs)
{
	this->atm_name=rhs.atm_name;
	this->res_name=rhs.res_name;
	this->res_number=rhs.res_number;
	this->coord=rhs.coord;
	this->element=rhs.element;
	this->hyd_c=rhs.hyd_c;
}
PDBSpecs& PDBSpecs::operator=(std::vector<double> coord)
{
	this->coord=coord;
}

PDBSpecs::~PDBSpecs() {};

/* PDBSpecs class member funtions */
void PDBSpecs::pdb_read(std::string buffer)
{
	this->atm_name=buffer.substr(12,4);
	this->atm_name.erase(remove_if(this->atm_name.begin(),this->atm_name.end(),::isspace),this->atm_name.end()); // remove spaces
	this->res_name=buffer.substr(17,3);
	this->res_number=atoi(buffer.substr(22,4).c_str());
	this->coord[0]=atof(buffer.substr(30,8).c_str());
	this->coord[1]=atof(buffer.substr(38,8).c_str());
	this->coord[2]=atof(buffer.substr(46,8).c_str());
	this->element=buffer.substr(76,2);
	this->element.erase(remove_if(this->element.begin(),this->element.end(),::isspace),this->element.end()); // remove spaces
}

double PDBSpecs::dis_sqrd(const PDBSpecs &x)
{
	return ((this->coord[0] - x.coord[0]) * (this->coord[0] - x.coord[0]) +
			(this->coord[1] - x.coord[1]) * (this->coord[1] - x.coord[1]) +
			(this->coord[2] - x.coord[2]) * (this->coord[2] - x.coord[2]));
}

bool PDBSpecs::close_to_backbone(int x, int y, std::vector<PDBSpecs> &p, double cutoff)
{
	double opdis_sqrd = this->dis_sqrd(p[x]);
	if (opdis_sqrd <= cutoff && p[x].atm_name == "O")
	{
		double COdis_sqrd = p[x].dis_sqrd(p[x-1]);
		double CWdis_sqrd = this->dis_sqrd(p[x-1]);
		double cosO = hb_cosine(opdis_sqrd, COdis_sqrd, CWdis_sqrd);
		if (cosO <= -0.4226)
			return true;
	}
	else if (opdis_sqrd <= cutoff && p[x].atm_name == "H")
	{
		double NHdis_sqrd = p[x].dis_sqrd(p[y]);
		double NWdis_sqrd = this->dis_sqrd(p[y]);
		double cosH = hb_cosine(opdis_sqrd, NHdis_sqrd, NWdis_sqrd);
		if (cosH <= -0.5)
			return true;
	}
	else
		return false;
}

PDBSpecs PDBSpecs::find_mid(const PDBSpecs &x)
{
	PDBSpecs mid_point;
	mid_point.coord[0]=(this->coord[0]+x.coord[0])/2;
	mid_point.coord[1]=(this->coord[1]+x.coord[1])/2;
	mid_point.coord[2]=(this->coord[2]+x.coord[2])/2;
	return mid_point;
}

double PDBSpecs::mid_rad_sqrd(const PDBSpecs &x, const PDBSpecs &y)
{
	PDBSpecs mid_prt;
	mid_prt=this->find_mid(x);
	if ( mid_prt.dis_sqrd(y) <= (4.0 / 9.0 * mid_prt.dis_sqrd(x)))
		return mid_prt.dis_sqrd(y);
	else
		return -1.0;
}

/* global funtions */
bool ascendSorting(double p, double q)
{
    return p < q;
}

bool descendSorting(double p, double q)
{
	return p > q;
}


double hb_cosine(double dis_sqrd_1, double dis_sqrd_2, double dis_sqrd_3)
{
	return ((dis_sqrd_1+dis_sqrd_2-dis_sqrd_3)/(2*sqrt(dis_sqrd_1)*sqrt(dis_sqrd_2)));
}

double gaus_mod(double height, double dis, double u, double sigma)
{
	return ((height*exp(-((dis-u)*(dis-u)/(2.0*sigma*sigma)))));
}

double rmsd(double *dis_sqrd, int dis_sqrd_size)
{
	double sum=0.0;
	for (int i=0; i<dis_sqrd_size; ++i)
		sum += dis_sqrd[i];
	return (sqrt(sum/dis_sqrd_size));
}

void ws_help()
{
	std::cout << "Usage: \n";
	std::cout << "    WSELECT -pdb complex.pdb -opdb water.pdb [ -cpdb crystal.pdb ] [ -t ( l LIG | p 1 10) | -center x y z [ -rad r ]] [ -cutoff cutoff ] [-hydcut hydcut]\n\n";
	std::cout << "Command line info:\n\n";
	std::cout << "    -pdb    : The modified PDB file for the complex. Ignore if analysis is about orginal crystal water.\n\n";
	std::cout << "    -opdb   : The PDB file for oxygen positions generated by 3D-RISM.\n\n";
	std::cout << "    -cpdb   : The reference PDB of crystal water. Ignore if not available.\n\n";
	std::cout << "    -o      : The prefix of output files. Default: \"WS\". \n\n";
	std::cout << "    -cutoff : The distance cut-off for water selection. Default 4.0 Angstrom.\n\n";
	std::cout << "    -hydcut : The hydropathic character cut-off value. Default -3.9 \n\n";
	std::cout << "  For holo complex, ligand specification should be provided:\n";
	std::cout << "    -t [l/p] [ligname | begres endres]:\n";
	std::cout << "      [l/p]: 'l' for normal organic ligand. 'p' for peptide-type ligand.\n";
	std::cout << "      [ligname | begres endres]:\n";
	std::cout << "        If 'l' option was set, enter the ligand name in PDB file.\n";
	std::cout << "        If 'p' option was set, enter the begin and end residue number separeted by a space.\n";
	std::cout << "        e.g. -t l LIG : normal organic ligand with name \"LIG\".\n";
	std::cout << "             -t p 1 10 : peptide as ligand, sequence number from 1 to 10.\n\n";
	std::cout << "  For apo complex, binding site center should be provided:\n";
	std::cout << "    -center [x y z]: the coordinate for the center of the binding site.\n\n";
	std::cout << "    -rad    : The radius for water selection from the center. Default: 7.0.\n\n";
}

/* hydropathic map to each residue name */
const std::map<std::string, double> hydp_map = {
	{"ILE",4.5},
	{"VAL",4.2},
	{"LEU",3.8},
	{"PHE",2.8},
	{"CYS",2.5},
	{"CYX",2.5},
	{"MET",1.9},
	{"ALA",1.8},
	{"GLY",-0.4},
	{"THR",-0.7},
	{"TRP",-0.9},
	{"SER",-0.8},
	{"TYR",-1.3},
	{"PRO",-1.6},
	{"HIE",-3.2},
	{"HID",-3.2},
	{"HIP",-3.2},
	{"HIS",-3.2},
	{"GLU",-3.5},
	{"GLN",-3.5},
	{"ASP",-3.5},
	{"ASN",-3.5},
	{"LYS",-3.9},
	{"ARG",-4.5},
};
