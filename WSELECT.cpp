#include "general.h"
#include "wsdeclare.h"

using namespace std;

int main(int argc, char *argv[])
{
	char *cmplxFile, *opdbFile, *crysFile;
	string o_prefix = "WS";
	char lgndType;
	double cutoff = 4.0;
	PDBSpecs center;
	double rad = 7.0;
	int apo = 0;
	double hydP_cutoff=-3.9;
	string lgndName;
	int begres, endres;
	int i, j;

	cout << "\n-------WSELECT Program-------\n";
	cout << "Author: X.hu (xiao.hu@unimi.it)\n";

	// error message when there is no arguments
	if ( argc <= 1 )
	{
		cout << "\nERROR: Not enough or invalid arguments, \'-help\' for usage info.\n";
		cout << "\n-------Program terminated with ERROR-------\n\n";
		return 0;
	}

	// argument parsing and reading files
	string argv_buffer(argv[1]);
	if ( argv_buffer == "-help" )
	{
		ws_help();
		return 0;
	}
	else
	{
		// initial output
		cout << "\n\'-help\' for detailed info.\n";
		time_t now = time(0);
		char* time = ctime(&now);
		cout << "Started on " << time;

		// read augments into variables
		for (i = 1; i < argc; i++) 
		{
			argv_buffer = argv[i];
			if (argv_buffer == "-pdb")
				cmplxFile = argv[i+1];
			else if (argv_buffer == "-opdb")
				opdbFile = argv[i+1];
			else if (argv_buffer == "-cpdb")
				crysFile = argv[i+1];
			else if (argv_buffer == "-t")
			{
				lgndType = *argv[i+1];
				switch (lgndType)
				{
					case 'l':
						lgndName = argv[i+2];
						break;
					case 'p':
						begres = atoi (argv[i+2]);
						endres = atoi (argv[i+3]);
						break;
				}
			}
			else if (argv_buffer == "-center")
			{
				center.coord[0] = atof (argv[i+1]);
				center.coord[1] = atof (argv[i+2]);
				center.coord[2] = atof (argv[i+3]);
			}
			else if (argv_buffer == "-rad")
				rad = atof (argv[i+1]);
			else if (argv_buffer == "-cutoff")
				cutoff = atof (argv[i+1]);
			else if (argv_buffer == "-o")
				o_prefix = argv[i+1];
			else if (argv_buffer == "-hydcut")
				hydP_cutoff = atof (argv[i+1]);
		}

	}
	
	// temporary pdb for protein and ligand
	vector<PDBSpecs> temp;
	PDBSpecs temp_var;
	ifstream File;
	string line;
	
	File.open(cmplxFile);
	
	if (File.good())
	{
		while (getline(File, line))
		{
			if (line.compare(0,6, "ATOM  ")==0)
			{
				temp_var.pdb_read(line); 
				// read hydropathic index to temp 
				temp_var.hyd_c = (hydp_map.find(temp_var.res_name))->second;
				
				temp.push_back(temp_var);
			}
		}
	}
	else
	{
		cout << "\nERROR: PDB file for complex missing or erred.\n";
		File.close();
		cout << "\n-------Program terminated with ERROR-------\n\n";
		exit (0);
	}

	File.close();

	// read data into protein and ligand array from temporary array
	vector<PDBSpecs> prtn; 
	vector<PDBSpecs> lgnd;
	int temp_size = (int)temp.size();
	int prtn_size;
	int lgnd_size;

	switch (lgndType)
	{
		case 'l':
			for ( i = 0; i < temp_size; ++i)
			{
				string temp_res(temp[i].res_name);
				string lgnd_name(lgndName);
				if (temp_res != lgnd_name)
					prtn.push_back(temp[i]);
				else
					lgnd.push_back(temp[i]);
			}
			prtn_size = (int)prtn.size();
			lgnd_size = (int)lgnd.size();
			break;

		case 'p':
			for ( i = 0; i < temp_size; ++i)
			{
				if (temp[i].res_number < begres || temp[i].res_number > endres)
					prtn.push_back(temp[i]);
				else
					lgnd.push_back(temp[i]);
			}
			prtn_size = (int)prtn.size();
			lgnd_size = (int)lgnd.size();
			break;

		default:
			cout << "\nNo ligand defined. Proceed with apo structure.\n\n";
			for ( i = 0; i < temp_size; ++i)
				prtn.push_back(temp[i]);
			prtn_size = (int)prtn.size();
			apo = 1;
			break;
	}

	// water oxygen postion
	vector<PDBSpecs> opdb;
	File.open(opdbFile);
	i = 0;
	
	if (File.good())
	{
		while (getline(File, line))
    	{
			if (line.compare(0,6, "ATOM  ")==0)
			{
				temp_var.pdb_read(line); 
				opdb.push_back(temp_var);
			}
		}
	}
	else
	{
		cout << "\nERROR: Water oxygen PDB file missing or erred.\n";
		File.close();
		File.clear();
		cout << "\n-------Program terminated with ERROR-------\n\n";
		exit (0);
	}

	File.close();
	File.clear();
	int opdb_size = (int)opdb.size();

	// Water postion from crystal structure
	vector<PDBSpecs> cwat;
	File.open(crysFile);

	if (File.good())
	{
		while (getline(File, line))
		{
			if (line.compare(0,6, "HETATM")==0 && line.compare(13,1, "O")==0 && line.compare(17,3,"HOH")==0)
			{
				temp_var.pdb_read(line);
				cwat.push_back(temp_var);
			}
		}
	}
	else
	{
		cout << "\nWARNING: No reference PDB. RMSD and relative distance calculations will be skipped.\n";
	}

	File.close();
	File.clear();

	// calculate minimum radius cut-off for apo structures
	double prad_sqrd = ((rad+2) * (rad+2));
	double rad_sqrd = (rad * rad);
	double minrad_sqrd;
	if (apo==1)
	{
		vector<double>mid_prt_rad;
		vector<PDBSpecs> prtn_inrad;
		for (auto &pn : prtn)
			if (pn.dis_sqrd(center) <= prad_sqrd)
				prtn_inrad.push_back(pn);
		int prtn_inrad_size = (int)prtn_inrad.size();
		if (prtn_inrad_size > 0)
		{
			for (i = 0; i < prtn_inrad_size; ++i)
			{
				for (j = i+1; j < prtn_inrad_size; ++j)
				{
					double mid_radsqrd=prtn_inrad[i].mid_rad_sqrd(prtn_inrad[j], center);
					if (mid_radsqrd >= 0)
						mid_prt_rad.push_back(mid_radsqrd);
				}
			}
			sort(mid_prt_rad.begin(), mid_prt_rad.end(), descendSorting);
			minrad_sqrd=mid_prt_rad[0];
			cout << "Minimum radius at binding site:" << sqrt(minrad_sqrd) << "\n";
			if (rad_sqrd <= minrad_sqrd+1)
			{
				cout << "The binding site radius is too small. Better increase the radius cut-off. \n";
				exit (0);
			}
		}
	}

	// select oxygen according to multiple criteria
	vector<PDBSpecs> oslct;
	double cutoff_sqrd = (cutoff * cutoff);
	int oslct_size = 0;
	int heav_atom = 0;

	string ofile, pfile, dumhead;
	dumhead = o_prefix;
	ofile = dumhead.append("_selected.pdb");
	dumhead = o_prefix;
	pfile = dumhead.append("_cmd_input");
	ofstream oFile, pFile;
	oFile.open(ofile, ios::out);
	pFile.open(pfile, ios::out);
	
	for ( i = 0; i < opdb_size; ++i )
	{

		// initiate distance array for squared distances to PROTEIN and calculate the squared distances
		vector<double> pdis_sqrd;
		for ( j = 0; j < prtn_size; ++j )
			if ( prtn[j].element != "H" && prtn[j].element != "C" )
				pdis_sqrd.push_back(opdb[i].dis_sqrd(prtn[j]));
		sort(pdis_sqrd.begin(), pdis_sqrd.end(), ascendSorting);
		
		// initiate distance array for squared distances to LIGAND and calculate the squared distances
		vector<double> ldis_sqrd;
		if (apo == 0)
		{
			for ( j = 0; j < lgnd_size; ++j )
				if ( lgnd[j].element != "H" && lgnd[j].element != "C" )
				{
					ldis_sqrd.push_back(opdb[i].dis_sqrd(lgnd[j]));
					heav_atom = 1;
				}
			if ( heav_atom == 1 )
				sort(ldis_sqrd.begin(), ldis_sqrd.end(), ascendSorting);
			else
			{
				cout << "\nThe ligand inside the pocket is likely hydrophobic."
					<<"\nPlease proceed without water or use an alternative structure.\n";
				exit (0);
			}
		}
		else
		{
			ldis_sqrd.push_back(opdb[i].dis_sqrd(center));
		}

		if ( pdis_sqrd[0] <= cutoff_sqrd && ((apo == 0 && ldis_sqrd[0] <= cutoff_sqrd) || (apo == 1 && ldis_sqrd[0] <= rad_sqrd && ldis_sqrd[0] >= minrad_sqrd )) )
		{
			double O_hydP = 0.0;
			int sc_count = 0;

			for ( j = 0; j < prtn_size; ++j )
			{
				int c = 0;
				int a = j;

				// count number of atoms (int c) for each protein residue
				while (true)
				{
					if (a+1 == prtn_size)
						break;
					int x = prtn[a].res_number;
					int y = prtn[a+1].res_number;
					if (y > x)
						break;
					++c;
					++a;
				}

				// output oxygen into "O_selected.pdb" if the atom is close to backbone
				int k;
				a = j;
				for (k = 0; k <= c; ++k, ++a)
				{
					if (opdb[i].close_to_backbone(a, j, prtn, cutoff_sqrd))
					{
						oslct.push_back(opdb[i]);
						O_hydP = 0.0;

						oFile << "ATOM  "
						<< right << setw(5) << i
						<< right << setw(5) << "O"
						<< right << setw(4) << "WAT" << " 1"
						<< right << setw(4) << opdb[i].res_number
						<< right << setw(12) << fixed << setprecision(3) << opdb[i].coord[0]
						<< right << setw(8) << opdb[i].coord[1]
						<< right << setw(8) << opdb[i].coord[2]
						<< "  1.00  0.00           O  " << "\n";

						pFile << "water_molecule\t" << opdb[i].coord[0] << "\t"
						<< opdb[i].coord[1] << "\t" << opdb[i].coord[2] << "\t3\n";

						j = prtn_size;
						k = c + 1;
						break;
					}
				}

				// calculate smallest distance for each protein residue to water molecule
				if (j < prtn_size)
				{
					vector<double> resdis_sqrd(c+1);
					for (k = 0; k <= c; ++k, ++j)
					{
						if (prtn[j].atm_name == "H" || prtn[j].atm_name == "O") // safeguard backbone
							resdis_sqrd[k] = 100;
						else if ( prtn[j].element != "H" && prtn[j].element != "C" )
							resdis_sqrd[k] = opdb[i].dis_sqrd(prtn[j]);
						else
							resdis_sqrd[k] = 100;
					}
					sort(resdis_sqrd.begin(), resdis_sqrd.end(), ascendSorting);

					/* the smallest distance shorter than strong hydrogen bond cutoff: 
					    add hydropathic index */
					if (resdis_sqrd[0] <= 9.0)
					{
						O_hydP += prtn[j-1].hyd_c;
						sc_count += 1;
					}
					/* the smallest distance between strong hydrogen bond cutoff and 
					    vdw interaction cutoff: 
					    add gaussian weighted hydropathic index */
					else if (resdis_sqrd[0] > 9.0 && resdis_sqrd[0] <= cutoff_sqrd)
					{
						O_hydP += gaus_mod(prtn[j-1].hyd_c, sqrt(resdis_sqrd[0]), 3.0, 0.6667); // 2 as half width, 3 sigma confidence
						sc_count += 1;
					}
				}
			}

			// output oxygen into "O_selected.pdb" if hydropathic character is not bigger than -3.9
			if (O_hydP <= hydP_cutoff && sc_count >= 2)
			{
				oslct.push_back(opdb[i]);

				oFile << "ATOM  "
				<< right << setw(5) << i
				<< right << setw(5) << "O"
				<< right << setw(4) << "WAT" << " 1"
				<< right << setw(4) << opdb[i].res_number
				<< right << setw(12) << fixed << setprecision(3) << opdb[i].coord[0]
				<< right << setw(8) << opdb[i].coord[1]
				<< right << setw(8) << opdb[i].coord[2]
				<< "  0.00"
				<< right << setw(6) << setprecision(2) << O_hydP
				<< "           O  " << "\n";

				pFile << "water_molecule\t" << opdb[i].coord[0] << "\t"
				<< opdb[i].coord[1] << "\t" << opdb[i].coord[2] << "\t3\n";
			}
		}
	}
	
	oFile.close();
	pFile.close();
	oslct_size = (int)oslct.size();

	// calculate distances of selected water molecules to crystal water and output values into file
	File.open(crysFile);
	if (oslct_size != 0 && File.good())
	{
		double ocdis_sqrd[oslct_size];
		int ocdis_count = 0; // count number of selected oxygen that is close to crystal water
		
		string disfile;
		dumhead = o_prefix;
		disfile = dumhead.append("_selected_info.dat");
		ofstream disFile(disfile, ios::out);

		disFile << "Distance cut-off: " << cutoff << "\n";
		disFile << "Hydropathic cut-off: " << hydP_cutoff << "\n";
		if (apo!=0)
		{
			disFile << "Radius to center: " << rad << "\n";
			disFile << "Center coordinates: " << center.coord[0] << " " << center.coord[1]
				<< " " << center.coord[2] << "\n\n";
		}
		disFile << "No." << "\t" << "Distance" << "\n";
		for (i = 0; i < oslct_size; ++i)
		{
		// distance array for distances to all crystal water oxygens
			int cwat_size = (int)cwat.size();
			double oodis_sqrd[cwat_size];
			for ( j = 0; j < cwat_size; ++j)
				oodis_sqrd[j] = oslct[i].dis_sqrd(cwat[j]);
			sort(oodis_sqrd, oodis_sqrd+cwat_size, ascendSorting);
			ocdis_sqrd[i] = oodis_sqrd[0];
			disFile << i+1 << "\t" << sqrt(ocdis_sqrd[i]) << "\n";
			if ( ocdis_sqrd[i] <= 9.0 )
				++ocdis_count;
		}

		// calculate RMSD and output into file
		disFile << "\n";
		double oc_rmsd=0.0;
		oc_rmsd=rmsd(ocdis_sqrd, oslct_size);
		disFile << "RMSD to crystal water molecule: " << oc_rmsd << "\n";
		disFile << "Number of selected water molecules that are close to crystal water: " << ocdis_count << "\n";
		disFile.close();
	}
	else if (oslct_size == 0)
		cout << "\nWARNING:\n\tNo water molecules are selected.\n";

	File.close();
	File.clear();

	cout << "\n-------Program terminated successfully-------\n\n";

	return (0);
}
