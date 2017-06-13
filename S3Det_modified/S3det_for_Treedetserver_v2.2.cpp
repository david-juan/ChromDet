// NOTICE: TO compile the program:: g++ -o S3det_for_Treedetserver_v2.exe S3det_for_Treedetserver_v2.cpp cluster.o subset.o /usr/lib/libboost_regex.a -L /home/arausell/Newmat10B/ -lnewmat -lm -I /home/arausell/Newmat10B/ -B /home/arausell/Newmat10B/
// NOTICE: To execute the program: S3det_for_Treedetserver_v2.exe -i Multiple_alignment_Fasta_file -o Output_file


#include <stdio.h> //Necessary for the sprintf function used in converting the "int type" into "char type"
#include <string.h>
// #include <iostream>  //Already included in the " #define WANT_STREAM " of the newmat10B library
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <math.h>
# include <set>
#include <ctype.h> //Necesasry for the toupper() function. Converts char to uppercase in function loading fasta file
#include <cstring>
//These ones are necessary for the "subset" functions
# include <cstdlib>
# include <iomanip>
// # include <cmath>
# include <ctime>
# include <limits>
# include "subset.H"


//This one is necessary for the int to binary conversion
#include <bitset>

using namespace std;

//The following statements are necessary for the NEWMAT10 LIBRARY - a Matrix Library with "eigen tools" 
#define WANT_STREAM
#include "newmatap.h"
#include "newmatio.h"
#include "myexcept.h"
#ifdef use_namespace
using namespace RBD_LIBRARIES;
#endif

#include "S3det_for_Treedetserver_v2.2.h"//Declaration & definition of structures and functions

//_______________________________________________________________________________________________________________________
//                                             STARTS MAIN
int main (int argc, char *argv[]) {

// system("clear");
try{
// string tmp_directory="/local/arausell/tmp_C_trials/";
// string exec_directory="/home/arausell/MCdet/Professional/";//Directory where R excutable files are hosted


int Analysis;
char Functional_study;
char * Function_coding_file;
char Number_axes_option;
int Fixed_number_of_axes;
int Maximum_number_of_axes;
char Number_groups_option;
int Fixed_number_of_groups;
int Maximum_number_of_groups;
int Number_of_axes;
int npass;
int npass_user_provided=0;
char * fasta_file;
char * out_file;

//default parameters
double percentage_of_gaps=0.10;
double percentage_of_similarity=0.4;
double selected_pertentage_of_residues=0.1;  //NOTICE: just considered in the unsupervised analysis
int averagerank_cutoff=10;
double Wilcoxon_cutoff=0.010;                 //NOTICE: For the moment it selects the last axes fulfiling the selected cut-off
double minimum_group_size=3;//Groups with less than minimum_group_size are disregarded when predinting residues but not in clustering output //NOTICE: just considered in the unsupervised analysis
string All_axes_in_column_distances_assessments="NO";// "YES" or "NO"
Analysis=2;//Analysis type: 1.MCA 2.MCA-Greenacre 3.Non-centered-PCA (Notice: non identical to FASS) 4.PCA
int conservation_threshold=98;
Functional_study='U';
string Verbose="NO";
double npass_threshold=0.05;//Percentage stablishing the threshold for "ifound" in kmeans (ifound=how often the optimal clustering solution was found out of "npass" runs)
int kmeans_repeats=5;
Number_axes_option='I';
Maximum_number_of_axes=10;
Number_groups_option='I';
Maximum_number_of_groups=1000;//Don't worry, it's just a hipersized limit that will be changed

string Index_type="CH_Index";// CH_Index_squared, CH_Index, DB_Index, DB_Index_squared, C_Index (Selection indeces for the number of clusters)
string combinations="YES";
string Interacton_option="NO";// "YES" or "NO"

//Parameters list for Treedetserver:
// if (parameter_list_for_Treedet_Server=="YES"){
// 	string fasta_file_string=fasta_file;
// 	srand((unsigned)time(0));
// 	int random_integer;
// 	int lowest=2000000, highest=20000000;
// 	int range=(highest-lowest)+1;
// 	random_integer = rand();
	
// 	char position_char[200];
// 	sprintf (position_char,"%d",random_integer);
	
// 	tmp_directory+="S3DET";
// 	tmp_directory+=position_char;
// 	tmp_directory+="/";
	
// 	string new_order_to_system="mkdir ";new_order_to_system+=tmp_directory;	
// 	const char* order_to_system=new_order_to_system.c_str();
// 	system(order_to_system);
//}

//string exec_directory="/gpfs_home/treedet/SoftTreedetV3/LIB/Methods/S3DET/";//Directory where R excutable files are hosted
// double selected_pertentage_of_residues=0.1;  //NOTICE: just considered in the unsupervised analysis
// double Wilcoxon_cutoff=0.010;                //NOTICE: For the moment it selects the last axes fulfiling the selected cut-off
// double percentage_of_gaps=0.10;
// double percentage_of_similarity=0.4;
// double minimum_group_size=3;//Groups with less than minimum_group_size are disregarded when predinting residues but not in clustering output //NOTICE: just considered in the unsupervised analysis

//Rest of parameters list:

// cout << "Choose Functional study type: S.Supervised U.Unsupervised : "; cin >> Functional_study;
// if(Functional_study=='S'){cout << "Type Functional coding file (absolute path) : "; cin >> Function_coding_file;}
// if(Functional_study=='U'){
// 	cout << "Choose option: F.Fixed number of axes or I.Indefinite number of axes (allow the program to choose it) : "; cin >> Number_axes_option;
// 	if(Number_axes_option=='F'){cout << "Select Number of axes (fixed) : "; cin >> Fixed_number_of_axes;}
// 	if(Number_axes_option=='I'){cout << "Select Maximum Number of axes to study : "; cin >> Maximum_number_of_axes;}	
// 	cout << "Choose option: F.Fixed number of groups or I.Indefinite number of groups (allow the program to choose it) : "; cin >> Number_groups_option;
// 	if(Number_groups_option=='F'){cout << "Select Number of groups (fixed) : "; cin >> Fixed_number_of_groups;}
// 	if(Number_groups_option=='I'){cout << "Select Maximum Number of groups to study : "; cin >> Maximum_number_of_groups;}	
// }
// cout << "Choose Analysis type: 1.MCA 2.MCA-Greenacre 3.Non-centered-PCA (Notice: non identical to FASS) 4.PCA : ";cin >> Analysis;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
if(argc<=2){
	cout << endl;
	cout << "   Required arguments: " << endl<< endl;
	
	cout << "     -i   infile                   Multiple Sequence Alignment in fasta format"      << endl;
	cout << "     -o   outfile                  Text file where program output will be printed"   << endl << endl;

	cout << "   Optional arguments regarding the analysis variants" << endl<< endl;
	
	cout << "     -f   Function coding file     Text file where a binary matrix codes for a super-"<< endl;
	cout << "                                     vised classification."                           << endl;
	cout << "                                     Default value is unsupervised"                   << endl << endl;
	
/*	cout << "     -t   Analysis type            Integer value that can take the following options"<< endl;
	cout << "                                     1 - MCA "                                       << endl;
	cout << "                                     2 - MCA with Greenacre corrections (Default)"   << endl;
	cout << "                                     3 - Non Centered PCA"                           << endl;
	cout << "                                     Default value is 2"                            << endl << endl;*/
	
	cout << "   Optional arguments regarding the processing of the Multiple sequence alignment: " << endl<< endl;
	
	cout << "     -g   % of gaps                Columns with a number of gaps greater than this"  << endl;
	cout << "                                     value are considered gappy columns and disregarded"<< endl;
	cout << "                                     from the analysis"                            << endl;
	cout << "                                     Default value is 10"                          << endl << endl;
	
	cout << "     -s   % of similarity          Sequences with a similarity to anyother less than"<< endl;
	cout << "          (for outliers)             this value are considered outliers and disregarded" << endl;
	cout << "                                     from the analysis"                            << endl;
	cout << "                                     Default value is 40"                          << endl<< endl;
	

	cout << "   Optional arguments regarding the prediction of Specificity Determining Positions: "  << endl<< endl;

	cout << "     -r   % of residues            Percentage of residues within those attributed "<< endl;
	cout << "                                     to each sequence cluster (or its combinations)"<< endl;
	cout << "                                     that is considered when predicting SDPs"       << endl;
	cout << "                                     Default value is 10"                          << endl<< endl;

	cout << "     -c   Average rank cut-off     Maximum average rank for a position being"       << endl;
	cout << "                                     considered as Specificity Determining (SDP)"   << endl;
	cout << "                                     Default value is 10"                            << endl<<endl;

// 	cout << "     -q (No Factorial Combination) This option deactivates the factorial combination" << endl;
// 	cout << "                                     of the initial clustering solution"              << endl<<endl;

	
	cout << "     -m   Minimum group size       Groups with a number of sequences less than this" << endl;
	cout << "                                     value are disregarded when predinting SDPs" << endl;
	cout << "                                     (but not in clustering output)"                 << endl;
	cout << "                                     Default value is 3"                             << endl<< endl;

	
	cout << "   Optional arguments regarding the unsupervised k-means clustering process: " << endl<< endl;
	
	cout << "     -k   Fixed number of groups   Integer value representing a fixed number of groups"  << endl;
	cout << "                                     upon which the k-means clustering will be performed" << endl<<endl;

	cout << "     -n   Number of iterationsps   Integer value representing the number of times the"  << endl;
	cout << "          in k-means algorithm       Expectation-Maximization algorithm performing"       << endl;
	cout << "                                     k-means clustering should be run"                    <<endl;
	cout << "                                     Default value is (number_of_sequences*10) with a minimum of 500"<< endl<<endl;
	
	cout << "   Optional arguments regarding the selection of the number of axes: " << endl<< endl;

	cout << "     -w   Wilcoxon cut-off         The last axis fulfiling the chosen cut-off will"  << endl;
	cout << "                                     be selected"                                    << endl;
	cout << "                                     Default value is 0.01"                          << endl<< endl;
	
	cout << "     -x   Fixed number of axes     Integer value representing a fixed number of axes "  << endl;
	cout << "                                     upon which the analysis will be performed."       << endl;
	cout << "                                     In such a case no selection of the number of "   << endl;
	cout << "                                     informative axes (nor Wilcoxon test) will be done"   << endl<< endl;
	
// 	cout << "     -a   All axes option          Option to consider distances from residues to"    << endl;
// 	cout << "                                     cluster centroids based just on the selected axes"<< endl<< endl;
	cout << "   Other optional arguments: " << endl<< endl;
	
	cout << "     -v   (Verbose)                Prints all additional outputs"                    << endl<<endl;
	
	//cout << "     -h   (Help)                   Prints help and program is not executed"          << endl<<endl;

	return 1;
}

static const boost::regex nonvalidsyntax("^(-)([iogsrcwmvfxkn])([^\\s])(.+)");
static const boost::regex nonvalidoption("^(-)([^iogsrcwmvfxkn]*)(.*)");
static const boost::regex intnumber("^(\\d)+");
static const boost::regex doublenumbergreaterthan0andlessthan1("^(0+)(\\.)(\\d*)([1-9]+)(\\d*)");
static const boost::regex matchsomecharacterotherthanhyphen("^[^-](.*)(\\w+)(.*)");

int myflag_iset=0;
int myflag_oset=0;
int myflag_fset=0;

for(int i=0;i<argc;i++){
	if(strcmp("-i",argv[i])==0){
		myflag_iset=1;
		if(boost::regex_match(argv[i+1],matchsomecharacterotherthanhyphen)){
			fasta_file=strdup(argv[i+1]);
		}else{
			cout << "\t" << "Error: Option "  << argv[i] << " requires an input file" << endl;
			cout << "\t" << "Please execute program without arguments to see \"help\" on valid options" << endl;
			cout << "\t" << "Program not executed" << endl;
			return 1;
		}
	}else if(strcmp("-o",argv[i])==0){
		myflag_oset=1;
		if(boost::regex_match(argv[i+1],matchsomecharacterotherthanhyphen)){
			out_file=strdup(argv[i+1]);
		}else{
			cout << "\t" << "Error: Option "  << argv[i] << " requires a file name to write the output" << endl;
			cout << "\t" << "Please execute program without arguments to see \"help\" on valid options" << endl;
			cout << "\t" << "Program not executed" << endl;
			return 1;
		}
	}else if(strcmp("-f",argv[i])==0){
		if(boost::regex_match(argv[i+1],matchsomecharacterotherthanhyphen)){
			Functional_study='S';
			Function_coding_file=strdup(argv[i+1]);
		}else{
			cout << "\t" << "Error: Option "  << argv[i] << " requires a function coding file" << endl;
			cout << "\t" << "Please execute program without arguments to see \"help\" on valid options" << endl;
			cout << "\t" << "Program not executed" << endl;
			return 1;
		}
	}else if(strcmp("-g",argv[i])==0){
		if((boost::regex_match(argv[i+1],intnumber))and(atof(argv[i+1])>=0)and(atof(argv[i+1])<=100)){
			percentage_of_gaps=atof(argv[i+1]);
			percentage_of_gaps=percentage_of_gaps/100;
		}else{
			cout << "\t" << argv[i+1] << " is not a valid argument for option " << argv[i] << endl;
			cout << "\t" << "Option "  << argv[i] << " requires an integer number between 0 and 100" << endl;
			cout << "\t" << "Program not executed" << endl;
			return 1;
		}	
	}else if(strcmp("-s",argv[i])==0){
		if((boost::regex_match(argv[i+1],intnumber))and(atof(argv[i+1])>=0)and(atof(argv[i+1])<=100)){
			percentage_of_similarity=atof(argv[i+1]);
			percentage_of_similarity=percentage_of_similarity/100;
		}else{
			cout << "\t" << argv[i+1] << " is not a valid argument for option " << argv[i] << endl;
			cout << "\t" << "Option "  << argv[i] << " requires an integer number between 0 and 100" << endl;
			cout << "\t" << "Program not executed" << endl;
			return 1;
		}
		
	}else if(strcmp("-r",argv[i])==0){
		if((boost::regex_match(argv[i+1],intnumber))and(atof(argv[i+1])>=0)and(atof(argv[i+1])<=100)){
			selected_pertentage_of_residues=atof(argv[i+1]);
			selected_pertentage_of_residues=selected_pertentage_of_residues/100;
		}else{
			cout << "\t" << argv[i+1] << " is not a valid argument for option " << argv[i] << endl;
			cout << "\t" << "Option "  << argv[i] << " requires an integer number between 0 and 100" << endl;
			cout << "\t" << "Program not executed" << endl;
			return 1;
		}
	}else if(strcmp("-c",argv[i])==0){
		if((boost::regex_match(argv[i+1],intnumber))and(atof(argv[i+1])>=0)){
                averagerank_cutoff=atoi(argv[i+1]);
		}else{
			cout << "\t" << argv[i+1] << " is not a valid argument for option " << argv[i] << endl;
			cout << "\t" << "Option "  << argv[i] << " requires an integer number greater or equal to 0" << endl;
			cout << "\t" << "Program not executed" << endl;
			return 1;
		}
	}else if(strcmp("-w",argv[i])==0){
		if(boost::regex_match(argv[i+1],doublenumbergreaterthan0andlessthan1)){
			Wilcoxon_cutoff=atof(argv[i+1]);
		}else{
			cout << "\t" << argv[i+1] << " is not a valid argument for option " << argv[i] << endl;
			cout << "\t" << "Option "  << argv[i] << " requires a value greater than 0 and lesser than 1 represented in floating point" << endl;
			cout << "\t" << "Program not executed" << endl;
			return 1;
		}
	}else if(strcmp("-m",argv[i])==0){
		if((boost::regex_match(argv[i+1],intnumber))and(atof(argv[i+1])>=0)){
			minimum_group_size=atof(argv[i+1]);
		}else{
			cout << "\t" << argv[i+1] << " is not a valid argument for option " << argv[i] << endl;
			cout << "\t" << "Option "  << argv[i] << " requires an integer number greater or equal to 0" << endl;
			cout << "\t" << "Program not executed" << endl;
			return 1;
		}
	}
	else if(strcmp("-a",argv[i])==0){
		All_axes_in_column_distances_assessments="YES";
	}
	else if(strcmp("-v",argv[i])==0){
		Verbose="YES";
	}
	//else if(strcmp("-t",argv[i])==0){
	//	Analysis=atoi(argv[i+1]);
	//}
	else if(strcmp("-x",argv[i])==0){
		if((boost::regex_match(argv[i+1],intnumber))and(atof(argv[i+1])>0)){
			Number_axes_option='F';
			Fixed_number_of_axes=atoi(argv[i+1]);
		}else{
			cout << "\t" << argv[i+1] << " is not a valid argument for option " << argv[i] << endl;
			cout << "\t" << "Option "  << argv[i] << " requires an integer number greater than 0" << endl;
			cout << "\t" << "Program not executed" << endl;
			return 1;
		}
	}else if(strcmp("-k",argv[i])==0){
		if((boost::regex_match(argv[i+1],intnumber))and(atof(argv[i+1])>=2)){
			Number_groups_option='F';
			Fixed_number_of_groups=atoi(argv[i+1]);
		}else{
			cout << "\t" << argv[i+1] << " is not a valid argument for option " << argv[i] << endl;
			cout << "\t" << "Option "  << argv[i] << " requires an integer number greater than 1" << endl;
			cout << "\t" << "Program not executed" << endl;
			return 1;
		}
	}else if(strcmp("-n",argv[i])==0){
		if((boost::regex_match(argv[i+1],intnumber))and(atof(argv[i+1])>0)){
			npass_user_provided=atoi(argv[i+1]);
		}else{
			cout << "\t" << argv[i+1] << " is not a valid argument for option " << argv[i] << endl;
			cout << "\t" << "Option "  << argv[i] << " requires an integer number greater than 0" << endl;
			cout << "\t" << "Program not executed" << endl;
			return 1;
		}
	}
	//else if(strcmp("-q",argv[i])==0){
	//	combinations="NO";
	//}
	else if(boost::regex_match(argv[i],nonvalidsyntax)){
		cout << "\t" << string (argv[i]) << " is not a valid syntax" << endl;
		cout << "\t" << "Please separate with a white space the options from their arguments (e.g. -i dir/infile -g 20)" << endl;
		cout << "\t" << "Program not executed" << endl;
		return 1;
	}else if(boost::regex_match(argv[i],nonvalidoption)){
		cout << "\t" << string (argv[i]) << " is not a valid option" << endl;
		cout << "\t" << "Please execute program without arguments to see \"help\" on valid options" << endl;
		cout << "\t" << "Program not executed" << endl;
		return 1;
	}
}
if(myflag_iset==0){
		cout << "\t" << "Error: No infile provided" << endl;
		cout << "\t" << "Please indicate an infile using \"-i dir/infile\"" << endl;
		cout << "\t" << "Execute program without arguments to see \"help\" on required arguments" << endl;
		cout << "\t" << "Program not executed" << endl;
		return 1;
}
if(myflag_oset==0){
		cout << "\t" << "Error: No outfile provided" << endl;
		cout << "\t" << "Please indicate an outfile using \"-o dir/outfile\"" << endl;
		cout << "\t" << "Execute program without arguments to see \"help\" on required arguments" << endl;
		cout << "\t" << "Program not executed" << endl;
		return 1;
}


ofstream fs(out_file); //Create output file with the chosen name
if (fs.fail()){
	string out_file_string=out_file;
	cout << "\t" << "\tError opening file " <<  out_file_string << endl << "\t\tFile not found or couldn't be opened" << endl;
	cout << "\t" << "\tExecution aborted" << endl;
	fs.close();
	return 1;
}
//***********************************************************************************************************************
 string conf_file_string2="conf.h";
 char * conf_file = new char [200];
 strcpy (conf_file,conf_file_string2.c_str());


 string tmp_directory;
 string exec_directory;
 string order;
 int error_loading_conf_file;
 int Conf_error_couldnt_open;

loading_conf_file(
		  conf_file,
		  tmp_directory,
		  exec_directory,
		  order,
		  error_loading_conf_file,
		  Conf_error_couldnt_open
);
if(Conf_error_couldnt_open==1){
	string conf_file_string=conf_file;
	cout << "\t" << "\tError loading file " <<  conf_file_string << endl << "\t\tFile not found or couldn't be opened" << endl;
	fs << "\t" << "\tError loading file " <<  conf_file_string << endl << "\t\tFile not found or couldn't be opened" << endl;
	cout << "\t" << "\tExecution aborted" << endl;
	fs << "\tExecution aborted" << endl;
	fs.close();
	return 1;
}
if(error_loading_conf_file==1){
	string conf_file_string=conf_file;
	cout << "\t" << "\tError loading file " <<  conf_file_string << endl << "\t\tRequired arguments are not in a proper format" << endl;
	fs << "\tError loading file " <<  conf_file_string << endl << "\t\tRequired arguments are not in a proper format" << endl;
	cout << "\t" << "\tExecution aborted" << endl;
	fs << "\tExecution aborted" << endl;
	fs.close();
	return 1;
}



char** letter_matrix;
vector<string> protein_name_vector; //it will contain the names of the proteins taking part in the alignment
vector<string> position_vector;//vector with aminoacid positions which have at least one non-null element
int number_of_sequences; //number of sequences in the alignment
int number_of_positions; //number of positions in the alignment
int error_loading_fasta;
int error_size_of_sequences;
int error_couldnt_open;
loading_fasta_file(
	fasta_file,
	letter_matrix,
	protein_name_vector,
	position_vector,
	number_of_sequences,
	number_of_positions,
	error_loading_fasta,
	error_size_of_sequences,
	error_couldnt_open
);
if(error_couldnt_open==1){
	string fasta_file_string=fasta_file;
	cout << "\t" << "\tError loading file " <<  fasta_file_string << endl << "\t\tFile not found or couldn't be opened" << endl;
	fs << "\t" << "\tError loading file " <<  fasta_file_string << endl << "\t\tFile not found or couldn't be opened" << endl;
	cout << "\t" << "\tExecution aborted" << endl;
	fs << "\tExecution aborted" << endl;
	fs.close();
	return 1;
}
if(error_loading_fasta==1){
	string fasta_file_string=fasta_file;
	cout << "\t" << "\tError loading file " <<  fasta_file_string << endl << "\t\tInput alignment is not in .mfa format" << endl;
	fs << "\tError loading file " <<  fasta_file_string << endl << "\t\tInput alignment is not in .mfa format" << endl;
	cout << "\t" << "\tExecution aborted" << endl;
	fs << "\tExecution aborted" << endl;
	fs.close();
	return 1;
}
if(error_size_of_sequences==1){
	string fasta_file_string=fasta_file;
	cout << "\t" << "\tError loading file " <<  fasta_file_string << endl << "\t\tDifferent sequence sizes within the alignment" << endl;
	fs << "\tError loading file " <<  fasta_file_string << endl << "\t\tDifferent sequence sizes within the alignment" << endl;
	cout << "\t" << "\tExecution aborted" << endl;
	fs << "\tExecution aborted" << endl;
	fs.close();
	return 1;
}

cout << "\t" << "Initial number of sequences: " << number_of_sequences << "\n";
cout << "\t" << "Initial number of positions: " << number_of_positions << "\n";
fs << "UI: Initial number of sequences: " << number_of_sequences << "\n";
fs << "UI: Initial number of positions: " << number_of_positions << "\n";
int initial_number_of_positions=number_of_positions;
int number_of_gappy_columns;
removing_gappy_columns(
	letter_matrix,
	percentage_of_gaps,
	position_vector,
	number_of_sequences,
	number_of_positions,
	number_of_gappy_columns
);

// fs << "Number of gappy columns: " << number_of_gappy_columns << endl;
// fs << "New number of positions: " << number_of_positions << "\n";

int number_of_outliers;
vector<string> outliers_name_vector;
removing_outliers(
	letter_matrix,
	percentage_of_similarity,
	number_of_sequences, 
	number_of_positions,
	number_of_outliers,
	protein_name_vector,
	outliers_name_vector
);
cout << "\t" << "Number of outliers: " << number_of_outliers << "\n";
fs << "UI: Number of outliers: " << number_of_outliers << "\n";
for(int i=0;i<outliers_name_vector.size();i++){
	cout << "\t" << "OUTLIER REMOVED: "<< outliers_name_vector[i] << endl;
	fs << "OUT: "<< outliers_name_vector[i] << endl;
}
fs << "UI: Final number of sequences after removing outliers: " << number_of_sequences << "\n";
fs << "Final number of sequences after removing outliers: " << number_of_sequences << "\n";

removing_gappy_columns(
	letter_matrix,
	percentage_of_gaps,
	position_vector,
	number_of_sequences,
	number_of_positions,
	number_of_gappy_columns
);
// fs << "New number of gappy columns: " << number_of_gappy_columns << endl;
fs << "UI: Final number of positions after removing gappy columns: " << number_of_positions << "\n";
cout << "\t" << "Final number of positions after removing gappy columns: " << number_of_positions << "\n";


Matrix A;
vector<string> position_aminoacid_vector;
vector<double> sum_of_elements_by_column_vector;
vector<double> sum_of_elements_by_row_vector;
double total_sum;
vector<int> Conserved_positions;

disjunctive_coding(
	A,
	letter_matrix,
	protein_name_vector,
	position_vector,
	number_of_sequences,
	number_of_positions,
	position_aminoacid_vector,
	sum_of_elements_by_column_vector,
	sum_of_elements_by_row_vector,
	total_sum,
	Conserved_positions,
	conservation_threshold
);
cout << "\t" << "Disjunctive coding finished" << "\n";
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int number_of_rows=number_of_sequences;
int number_of_columns=position_aminoacid_vector.size();
Matrix Z(number_of_rows,number_of_columns);
fs << "II: Final number of columns in the disjunctive matrix: " << number_of_columns << endl;
DiagonalMatrix D_Reduced,D_adjusted;
Matrix V_Reduced;
DiagonalMatrix Dc(number_of_columns);
Matrix Column_Coordinates;
//Matrix Row_Coordinates_standard;
Matrix Row_Coordinates_principal_indirect;

eigen_descomposition(
	Z,
	D_Reduced,
	D_adjusted,//Only used in the MCA-Greenacre analysis (Analysis 2)
	V_Reduced,
	Dc,
	Column_Coordinates,
	//Row_Coordinates_standard,
	Row_Coordinates_principal_indirect,
	A,
	number_of_rows,
	number_of_positions,
	number_of_columns,
	sum_of_elements_by_column_vector,
	sum_of_elements_by_row_vector,
	total_sum,
	Analysis
);
cout << "\t" << "Eigen decomposition finished" << "\n";
if(Verbose=="YES"){
/*	fs << "Row_Coordinates_principal_indirect: " << endl << setw(8) << setprecision(10) << Row_Coordinates_principal_indirect << "\n" << "\n";
	//fs << "Row_Coordinates_standard: " << endl << setw(8) << setprecision(10) << Row_Coordinates_standard << "\n" << "\n";
	fs << "Column_Coordinates: " << endl << setw(8) << setprecision(10) << Column_Coordinates << "\n" << "\n";*/

	//fs << "Row_Coordinates_principal_indirect: " << endl;
	fs << "Sequence_Coordinates: " << endl;
	fs << "Axes\t1st\t2nd\t3rd\t4th\t5th\t6th\t7th\t8th\t9th\t10th" << endl;
	for (int i=0; i<Row_Coordinates_principal_indirect.Nrows();i++){
		fs << "SeqCoord:" << "\t" << protein_name_vector.at(i) << "\t";
		//for (int j=Row_Coordinates_principal_indirect.Ncols()-1;j>Row_Coordinates_principal_indirect.Ncols()-1-2; j--){
		for (int j=Row_Coordinates_principal_indirect.Ncols()-1;j>Row_Coordinates_principal_indirect.Ncols()-1-9; j--){
			fs << Row_Coordinates_principal_indirect[i][j] << "\t";
		}
		fs << Row_Coordinates_principal_indirect[i][Row_Coordinates_principal_indirect.Ncols()-1-9] << endl;
		//fs << Row_Coordinates_principal_indirect[i][Row_Coordinates_principal_indirect.Ncols()-1-2] << endl;
	}
	fs << endl;
	fs << "Residue_Coordinates: " << endl;
	fs << "Axes\t1st\t2nd\t3rd\t4th\t5th\t6th\t7th\t8th\t9th\t10th" << endl;
	for (int i=0; i<Column_Coordinates.Nrows();i++){
		fs << "ResCoord:" << "\t" << position_aminoacid_vector.at(i) << "\t";
		for (int j=Column_Coordinates.Ncols()-1;j>Column_Coordinates.Ncols()-1-9; j--){
			fs << Column_Coordinates[i][j] << "\t";
		}
		fs << Column_Coordinates[i][Column_Coordinates.Ncols()-1-9] << endl;
	}
	fs << endl;
}
DiagonalMatrix D_Reduced_sqrt(D_Reduced.Ncols());//Only used in the calculation of function coordinates in the MCA-Greenacre Analysis
DiagonalMatrix D_Reduced_sqrt_inv(D_Reduced.Ncols());
DiagonalMatrix D_adjusted_sqrt(D_Reduced.Ncols());
if(Analysis==2){
	for (int i=D_Reduced.Ncols();i>0; i--){
		D_Reduced_sqrt(i,i)=sqrt(D_Reduced(i,i));
		D_Reduced_sqrt_inv(i,i)=(1/sqrt(D_Reduced(i,i)));
		D_adjusted_sqrt(i,i)=sqrt(D_adjusted(i,i));
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Assessing and printing eigenvalues
double sum_eigenvalues=0;
DiagonalMatrix D_percentage;
DiagonalMatrix D_percentage_accumulated(D_Reduced.Ncols());
if(Analysis==1||Analysis==3||Analysis==4){
	for (int i=D_Reduced.Ncols();i>0; i--){
		 sum_eigenvalues+=D_Reduced(i,i);
	}
	D_percentage=(1/sum_eigenvalues)*D_Reduced;
	D_percentage_accumulated(D_percentage.Ncols(),D_percentage.Ncols())=D_percentage(D_percentage.Ncols(),D_percentage.Ncols());
	for (int i=D_Reduced.Ncols()-1;i>0; i--){
		D_percentage_accumulated(i,i)=D_percentage_accumulated(i+1,i+1)+D_percentage(i,i);
	}
	//for (int i=D_Reduced.Ncols();i>0; i--){
	//	fs << D_Reduced(i,i)  << "\t" << D_percentage(i,i) << "\t" << D_percentage_accumulated(i,i) << "\n";
	//}
}
if(Analysis==2){
	for (int i=D_adjusted.Ncols();i>0; i--){
		sum_eigenvalues+=D_adjusted(i,i);
	}
	D_percentage=(1/sum_eigenvalues)*D_adjusted;
	D_percentage_accumulated(D_percentage.Ncols(),D_percentage.Ncols())=D_percentage(D_percentage.Ncols(),D_percentage.Ncols());
	for (int i=D_adjusted.Ncols()-1;i>0; i--){
		D_percentage_accumulated(i,i)=D_percentage_accumulated(i+1,i+1)+D_percentage(i,i);
	}
	//for (int i=D_adjusted.Ncols();i>0; i--){
	//	fs << D_adjusted(i,i)  << "\t" << D_percentage(i,i) << "\t" << D_percentage_accumulated(i,i) << "\n";
	//}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Matrix Function_matrix;
vector<string> function_names_vector;
int number_of_functions;
double flag_for_interaction;

vector<function_information_structure> function_information;
int number_of_groups_greater_than_minimum_group_size=0;
int number_of_sequences_within_groups_greater_than_minimum_group_size=0;
int number_of_groups_with_less_than_minimum_group_size=0;


int initial_number_of_functions=0;
if(Functional_study=='S'){
	vector<string> function_names_vector_DUMMY;
	char * function_file=Function_coding_file;
	string function_file_string=function_file;
	int myerror_couldnt_open;
	int myerror_file_is_not_in_proper_format;
	int myerror_sequence_not_found_in_function_file;
	int myflag_sequence_with_not_assigned_function;
	int myflag_sequence_with_sum_of_assigned_functions_different_from_1;
	string myprotein_name_causing_error;
	int flag_FUZZYClassification;
	//cout << "\t" << "function_file " << function_file << endl;
	loading_supervised_function_coding(
		function_file,
		protein_name_vector,
		Function_matrix,
		function_names_vector_DUMMY,
		number_of_functions,
		flag_for_interaction,
		myerror_couldnt_open,
		myerror_file_is_not_in_proper_format,
		myerror_sequence_not_found_in_function_file,
		myflag_sequence_with_not_assigned_function,
		myflag_sequence_with_sum_of_assigned_functions_different_from_1,
		myprotein_name_causing_error,
		flag_FUZZYClassification
	);
	if(myerror_couldnt_open==1){
		string function_file_string=function_file;
		cout << "\t" << "\tError loading file " <<  function_file_string << endl << "\t\tFile not found or couldn't be opened" << endl;
		fs << "\t" << "\tError loading file " <<  function_file_string << endl << "\t\tFile not found or couldn't be opened" << endl;
		cout << "\t" << "\tExecution aborted" << endl;
		fs << "\t"<< "\tExecution aborted" << endl;
		fs.close();
		return 1;
	}
	if(myerror_file_is_not_in_proper_format==1){
		cout << "\t" << "\tError loading file: " <<  function_file_string << " has not a proper format" << endl;
		fs << "\t"<< "\tError loading file " <<  function_file_string << " has not a proper format" << endl;
		cout << "\t" << "\tExecution aborted" << endl;
		fs << "\t"<< "\tExecution aborted" << endl;
		fs.close();
		return 1;
	}
	if(myflag_sequence_with_not_assigned_function==1){
		cout << "\t" << "\tError loading file " <<  function_file_string << endl;
		cout << "\t" << "\tSequence " << myprotein_name_causing_error << " has not assigned function in " << function_file_string << endl;
		fs << "\t"<< "\tError loading file " <<  function_file_string << endl;
		fs << "\t"<< "\tSequence " << myprotein_name_causing_error << " has not assigned function in " << function_file_string << endl;
		cout << "\t" << "\tExecution aborted" << endl;
		fs << "\t"<< "\tExecution aborted" << endl;
		fs.close();
		return 1;
	}
	if(myflag_sequence_with_sum_of_assigned_functions_different_from_1==1){
		cout << "\t" << "\tError loading file " <<  function_file_string << endl;
		cout << "\t" << "\tSequence " << myprotein_name_causing_error << " assigned functions values sum up a number different from 1 " << function_file_string << endl;
		fs << "\t"<< "\tError loading file " <<  function_file_string << endl;
		fs << "\t"<< "\tSequence " << myprotein_name_causing_error << " assigned functions values sum up a number different from 1 " << function_file_string << endl;
		cout << "\t" << "\tExecution aborted" << endl;
		fs << "\t"<< "\tExecution aborted" << endl;
		fs.close();
		return 1;
	}
	if(myerror_sequence_not_found_in_function_file==1){
		cout << "\t" << "\tError loading file " <<  function_file_string << endl;
		cout << "\t" << "\tSequence " << myprotein_name_causing_error << " not found in " << function_file_string << endl;
		fs << "\t"<< "\tError loading file " <<  function_file_string << endl;
		fs << "\t"<< "\tSequence " << myprotein_name_causing_error << " not found in " << function_file_string << endl;
		cout << "\t" << "\tExecution aborted" << endl;
		fs << "\t"<< "\tExecution aborted" << endl;
		fs.close();
		return 1;
	}
	cout << "\t" << "Number of functions: " << number_of_functions << endl;
	initial_number_of_functions=number_of_functions;

	
	for(int i=1;i<=Function_matrix.Ncols();i++){
		char number_of_functions_char[200];
		sprintf (number_of_functions_char,"%d",i);//converts the "int type" into "char type"
		string number_of_functions_string= number_of_functions_char;//converts the "char type" into "string type"
		//string function_name_tmp="Cluster_"+number_of_functions_string;
		string function_name_tmp=function_names_vector_DUMMY.at(i-1);
		vector<int> clusters_belonging_vector_tmp;
		clusters_belonging_vector_tmp.push_back(i);
		function_names_vector.push_back(function_name_tmp);
		
		function_information_structure tmp;
		tmp.name=function_name_tmp;
		tmp.size= (int) Function_matrix.Column(i).Sum();
		tmp.clusters_belonging_vector=clusters_belonging_vector_tmp;
		
		function_information.push_back(tmp);
		if (tmp.size>=minimum_group_size){
			number_of_sequences_within_groups_greater_than_minimum_group_size+=tmp.size;
		}
		if (tmp.size<minimum_group_size){
			number_of_groups_with_less_than_minimum_group_size++;
		}
		if (flag_for_interaction > tmp.size){
			flag_for_interaction = tmp.size;
		}
	}
	Fixed_number_of_groups=number_of_functions;
	number_of_groups_greater_than_minimum_group_size=Fixed_number_of_groups-number_of_groups_with_less_than_minimum_group_size;
	
	if(flag_FUZZYClassification==1){
	  fs << "UI: NOTICE: Supervised Fuzzy Classification is active" << endl;
	  cout << "\t" << "NOTICE: Supervised Fuzzy Classification is active" << endl;
	}
	fs << "Supervised classification" << endl;
	for(int i=0; i<number_of_rows; i++){
	  fs << "CL: "<< protein_name_vector[i];
	  int labelsup=0;
	  for(int j=1;j<=Function_matrix.Ncols();j++){
	    if(Function_matrix(i+1,j)==1){
	      fs << "\t" << j;
	    }
	    else{ 
	      if(Function_matrix(i+1,j)>0){;
		fs << "\t" << j << "(Prob:" << Function_matrix(i+1,j) << ")";
	      }
	    }
	  }
	  fs << endl;
	}
	fs << endl;
}
if(Functional_study=='U'){
	npass=number_of_rows*10;
	if(npass<500){npass=500;}
	if(npass_user_provided>npass){npass=npass_user_provided;}
	Number_of_axes=Column_Coordinates.Ncols();
	if(All_axes_in_column_distances_assessments!="YES")
	{
		if(Number_axes_option=='F'){
			if(Fixed_number_of_axes>V_Reduced.Ncols()){Number_of_axes=V_Reduced.Ncols();}
			else{Number_of_axes=Fixed_number_of_axes;}
		}
		else
		{
			if(Number_axes_option=='I'){	// Selection of the number of axes by the Wilcoxon test
				Number_of_axes=0;
				double* Wilcoxon_pvalues;
				int Wilcox_error_couldnt_open;
				Wilcoxon_test(
					order,
					Row_Coordinates_principal_indirect,
					tmp_directory,
					exec_directory,
					Wilcoxon_pvalues,
					Maximum_number_of_axes,
					Wilcox_error_couldnt_open
				);
				if(Wilcox_error_couldnt_open==1){
				  cout << "\t" << "\tError executing S3det_Wilcoxon_test.R ; its temporary output file couldn't be found or opened" << endl;
				  fs << "\t" << "Error executing S3det_Wilcoxon_test.R ; its temporary output file couldn't be found or opened" << endl;
				  cout << "\t" << "\tWARNING: Please, check that R-package is installed and running according to the conf.h file" << endl;
				  cout << "\t" << "\tWARNING: You might also consider changing the conf.h file and then recompiling S3det" << endl;
				  fs << "\t" << "WARNING: Please, check that R-package is installed and running according to the conf.h file" << endl;
				  fs << "\t" << "WARNING: You might also consider changing the conf.h file and then recompiling S3det" << endl;
				  cout << "\t" << "\tExecution aborted" << endl;
				  fs << "\tExecution aborted" << endl;
				  fs.close();
				  return 1;
				}
				if(Verbose=="YES"){ fs << "Wilcoxon values:" << endl;}
				for(int i=0;i<Maximum_number_of_axes-2;i++){
					if(Verbose=="YES"){ fs << Wilcoxon_pvalues[i] << endl;}
					if(Wilcoxon_pvalues[i]>Wilcoxon_cutoff && Number_of_axes==0){Number_of_axes=i+1;}//NOTICE: For the moment it selects the last axes fulfiling the selected cut-off
					//cout << Wilcoxon_pvalues[i] << endl;
				}
			}
		}
	}
	if(All_axes_in_column_distances_assessments=="YES")
	{
		cout << "\t" << "Number of axes selected: " << Number_of_axes << "(all)" <<endl;
		cout << "\t" << "Percentage of initial variance considered informative: 100 %" << endl;
	}else
	{
		cout << "\t" << "Number of axes selected: " << Number_of_axes <<endl;
			cout << "\t" << "Percentage of initial variance considered informative: " << (D_percentage_accumulated(D_percentage_accumulated.Ncols()+1-Number_of_axes,D_percentage_accumulated.Ncols()+1-Number_of_axes))*100 << " %" << endl;
		for(int i=1;i<=Number_of_axes;i++){
			cout << "\t\t" << "Percentage of variance explained by the axis "<<  i  <<  "  (selected): " << (D_percentage_accumulated(D_percentage_accumulated.Ncols()+1-i,D_percentage_accumulated.Ncols()+1-i))*100 << " %" << endl;
		}		
	}
	fs << "UI: Number of axes selected: " << Number_of_axes <<endl;
	fs << "UI: Percentage of initial variance considered informative: " << (D_percentage_accumulated(D_percentage_accumulated.Ncols()+1-Number_of_axes,D_percentage_accumulated.Ncols()+1-Number_of_axes))*100 << " %" << endl;
	for(int i=1;i<=Number_of_axes;i++){
		fs << "UI: Percentage of variance explained by the axis "<<  i   <<  "  (selected): " << (D_percentage_accumulated(D_percentage_accumulated.Ncols()+1-i,D_percentage_accumulated.Ncols()+1-i))*100 << " %" << endl;
	}		

	//if(Verbose=="YES"){
	//	fs << "D_percentage_accumulated :" << endl << setw(8) << setprecision(10) << D_percentage_accumulated << endl << endl;
	//}
	//*************************************************************************************************************************
	// Clustering process
	const int nrows = number_of_rows;
	const int ncols = Number_of_axes;

	int* clusterid = new int [nrows];
	double Index_value;
	int ifound;
	double F_Fisher;

	double **data= new double *[number_of_rows];
	for (int i=0; i<number_of_rows; i++){
		data[i]=new double [Number_of_axes];
	}
	int **mask= new int *[number_of_rows];
	for (int i=0; i<number_of_rows; i++){
		mask[i]=new int [Number_of_axes];
	}
	double** distmatrix;
	for (int i=0; i<number_of_rows; i++){
		for (int j=0; j<ncols; j++){
			data[i][j]=Row_Coordinates_principal_indirect[i][Row_Coordinates_principal_indirect.Ncols()-1-j];
			mask[i][j]=1;
		}
	}
	if(Number_groups_option=='I'){	// Selection of the number of groups
		Maximum_number_of_groups=(number_of_rows/4);
		if((number_of_rows%4)>0){Maximum_number_of_groups++;}
		if (Maximum_number_of_groups>50){Maximum_number_of_groups=50;}
		cout << "\t" << "Performing k-means clustering from 2 to a maximum number of groups equal to: " << Maximum_number_of_groups << endl;
		fs << "UI: Checking a maximum number of groups equal to: " << Maximum_number_of_groups << endl;
		int** Kvectors_clusterid   =new int * [Maximum_number_of_groups-1];
		for(int i=0; i<Maximum_number_of_groups-1;i++){Kvectors_clusterid[i]=new int [nrows];}
		double*  Kvectors_ifound      =new double [Maximum_number_of_groups-1];
		double*  Kvectors_Index_value =new double [Maximum_number_of_groups-1];
		double*  Kvectors_F_Fisher    =new double [Maximum_number_of_groups-1];
		double*  Kvectors_NumGroups   =new double [Maximum_number_of_groups-1];
		int num_successful_kmeans=-1;
		for (int n=2;n<=Maximum_number_of_groups;n++){
		  const int nclusters = n;
		  int* clusterid_inner = new int [nrows];
		  int ifound_inner;
		  double F_Fisher_inner;
		  double Index_value_inner;
		  int flag_unstable_clustering;
		  kmeans_v2_2(nrows,ncols,data,mask,nclusters,Analysis,Index_type,npass,npass_threshold,kmeans_repeats,flag_unstable_clustering,clusterid_inner,ifound_inner,Index_value_inner,F_Fisher_inner);
		  if(flag_unstable_clustering==0){
		    num_successful_kmeans++;
		    for(int i=0; i<nrows;i++){Kvectors_clusterid[num_successful_kmeans][i]=clusterid_inner[i];}
		    Kvectors_ifound[num_successful_kmeans]=ifound_inner;
		    Kvectors_Index_value[num_successful_kmeans]=Index_value_inner;
		    Kvectors_F_Fisher[num_successful_kmeans]=F_Fisher_inner;
		    Kvectors_NumGroups[num_successful_kmeans]=n;
		  }
		  delete [] clusterid_inner;
		}
		if(num_successful_kmeans==-1){
		  cout << "\t" << "No stable cluster results have been found" << endl;
		  fs << "UI: No stable cluster results have been found" << endl;
		  fs << "SC: NO" << endl;
		  return 0;
		}
		if(Verbose=="YES"){
		  fs << "UI: Stable cluster results have been found for the following number of groups:" << endl;;
		  for(int k=0;k<=num_successful_kmeans;k++){
		    fs << "CS: When looking for " << Kvectors_NumGroups[k] << " groups, solution is found " << Kvectors_ifound[k]<< " times out of " << npass << " iterations with a " <<  Index_type << "=" << Kvectors_Index_value[k] << "; Clustering solution:\t";
		    for(int i=0; i<nrows; i++){fs << Kvectors_clusterid[k][i]+1 << " ";}
		    fs << endl;
		  }
		}
		int best_solution_index=-1;
		double best_Index_value=Kvectors_Index_value[0];
		Fixed_number_of_groups=2;
		for(int i=0;i<=num_successful_kmeans;i++){
		  if(Index_type=="CH_Index"||Index_type=="CH_Index_squared"){
		    if (Kvectors_Index_value[i]>=best_Index_value){
		      best_Index_value=Kvectors_Index_value[i];
		      Fixed_number_of_groups=Kvectors_NumGroups[i];
		      best_solution_index=i;
		    }
		  }
		  if(Index_type=="DB_Index"||Index_type=="DB_Index_squared"||Index_type=="C_Index"){
		    if (Kvectors_Index_value[i]<=best_Index_value){//NOTICE: WARNING: If There were several C_Indices equal to zero, for the moment it is choosing the last one!!!
		      best_Index_value=Kvectors_Index_value[i];
		      Fixed_number_of_groups=Kvectors_NumGroups[i];
		      best_solution_index=i;
		    }
		  }
		}
		cout << "\t" << "Stable cluster results have been found" << endl;
		cout << "\t" << "Number of groups selected: " << Fixed_number_of_groups <<endl;
		cout << "\t" << "Best solution found " << Kvectors_ifound[best_solution_index] << " times out of " << npass << " iterations searching for " << Fixed_number_of_groups << " groups"  << endl;

		fs << "UI: Stable cluster results have been found" << endl;
		fs << "SC: YES" << endl;
		fs << "UI: Number of groups selected: " << Fixed_number_of_groups <<endl;
		fs << "II: Clustering index value: " << Index_type << " = " << best_Index_value << endl;
		fs << "II: F-Fisher: " << Kvectors_F_Fisher[best_solution_index] << endl;
		fs << "II: Best solution found " << Kvectors_ifound[best_solution_index] << " times out of " << npass << " iterations searching for " << Fixed_number_of_groups << " groups"  << endl;

		for(int i=0; i<nrows;i++){clusterid[i]=Kvectors_clusterid[best_solution_index][i];}
		Index_value= best_Index_value;
		ifound= Kvectors_ifound[best_solution_index];
		F_Fisher= Kvectors_F_Fisher[best_solution_index];
	}
	if(Number_groups_option=='F'){

		Maximum_number_of_groups=(number_of_rows/4);
		if((number_of_rows%4)>0){Maximum_number_of_groups++;}
		if(Maximum_number_of_groups>50){Maximum_number_of_groups=50;}
		if(Fixed_number_of_groups>Maximum_number_of_groups){Fixed_number_of_groups=Maximum_number_of_groups;}


		cout << "\t" << "Performing k-means clustering for a fixed number of groups equal to: " << Fixed_number_of_groups << endl;
		fs << "UI: Checking a fixed number of groups equal to: " << Fixed_number_of_groups << endl;
		const int nclusters = Fixed_number_of_groups;
		int* clusterid_inner = new int [nrows];
		int ifound_inner;
		double F_Fisher_inner;
		double Index_value_inner;
		int flag_unstable_clustering;
		kmeans_v2_2(nrows,ncols,data,mask,nclusters,Analysis,Index_type,npass,npass_threshold,kmeans_repeats,flag_unstable_clustering,clusterid_inner,ifound_inner,Index_value_inner,F_Fisher_inner);
		if(flag_unstable_clustering==1){
		  cout << "\t" << "No stable cluster results have been found" << endl;
		  fs << "UI: No stable cluster results have been found" << endl;
		  fs << "SC: NO" << endl;
		  return 0;
		}
		if(flag_unstable_clustering==0){
		  cout << "\t" << "Stable cluster results have been found" << endl;
		  cout << "\t" << "Number of groups selected: " << Fixed_number_of_groups <<endl;
		  cout << "\t" << "Best solution found " << ifound_inner << " times out of " << npass << " iterations searching for " << Fixed_number_of_groups << " groups"  << endl;

		  fs << "UI: Stable cluster results have been found" << endl;
		  fs << "SC: YES" << endl;
		  fs << "II: Clustering index value: " << Index_type << " = " << Index_value_inner << endl;
		  fs << "II: F-Fisher: " << F_Fisher_inner << endl;
		  fs << "II: Best solution found " << ifound_inner << " times out of " << npass << " iterations searching for " << Fixed_number_of_groups << " groups"  << endl;

		  for(int i=0; i<nrows;i++){clusterid[i]=clusterid_inner[i];}
		  Index_value= Index_value_inner;
		  ifound= ifound_inner;
		  F_Fisher= F_Fisher_inner;

		}
		delete [] clusterid_inner;
	}

	fs << "Clustering of sequences: " << endl;
	for(int i=0; i<number_of_rows; i++){
		fs << "CL: "<< protein_name_vector[i] << "\t"<< clusterid[i]+1 << endl;
	}
	fs << endl;

	for (int i=0; i<number_of_rows; i++){
		delete [] data[i];
		delete [] mask[i];
	}
	delete  [] data;
	delete  [] mask;
	
	number_of_functions=Fixed_number_of_groups;
	Matrix Function_matrix_aux(number_of_rows,number_of_functions);
	Function_matrix_aux=0;
	for(int j=0; j<number_of_functions; j++){
		for(int i=0; i<number_of_rows; i++){
			if(j==clusterid[i]){Function_matrix_aux[i][j]=1;}
		}
	}

	flag_for_interaction=Function_matrix_aux.Column(1).Sum();
	for(int i=1;i<=Function_matrix_aux.Ncols();i++){
		char number_of_functions_char[200];
		sprintf (number_of_functions_char,"%d",i);//converts the "int type" into "char type"
		string number_of_functions_string= number_of_functions_char;//converts the "char type" into "string type"
		string function_name_tmp="Cluster_"+number_of_functions_string;
		vector<int> clusters_belonging_vector_tmp;
		clusters_belonging_vector_tmp.push_back(i);
		function_names_vector.push_back(function_name_tmp);
		
		function_information_structure tmp;
		tmp.name=function_name_tmp;
		tmp.size= (int) Function_matrix_aux.Column(i).Sum();
		tmp.clusters_belonging_vector=clusters_belonging_vector_tmp;
		
		function_information.push_back(tmp);
		if (tmp.size>=minimum_group_size){
			number_of_sequences_within_groups_greater_than_minimum_group_size+=tmp.size;
		}
		if (tmp.size<minimum_group_size){
			number_of_groups_with_less_than_minimum_group_size++;
		}
		if (flag_for_interaction > tmp.size){
			flag_for_interaction = tmp.size;
		}
	}
		cout << "\t" << "Number of groups with less than "<< minimum_group_size << " sequences: " << number_of_groups_with_less_than_minimum_group_size <<endl;
	fs << "UI: Number of groups with less than "<< minimum_group_size << " sequences: " << number_of_groups_with_less_than_minimum_group_size <<endl;
	number_of_groups_greater_than_minimum_group_size=Fixed_number_of_groups-number_of_groups_with_less_than_minimum_group_size;
	fs << "Is there enough group diversity to perform the analysis?:" << endl << "GD: ";
	if (Fixed_number_of_groups-number_of_groups_with_less_than_minimum_group_size< 2){
		fs << "NO" << endl;
		fs << "UI: There is not enough group diversity to perform the analysis" << endl;
		cout << "\t" << "There is not enough group diversity to perform the analysis" << endl;
		return 0;
	}else{
		cout << "\t" << "There is enough group diversity to perform the analysis" << endl;
		fs << "YES" << endl << endl;
	}
	initial_number_of_functions=number_of_functions;
	Function_matrix=Function_matrix_aux;
	Function_matrix_aux.Release();
}

if((combinations=="YES")&&(Fixed_number_of_groups>1)){
	int number_of_combinations=0;//including groups with less sequences than the minimum_group_size; AND also including the "total"
	for(int n=2;n<=number_of_functions;n++){
		number_of_combinations += i4_choose (number_of_functions,n);//function from the subset.o object: i4_choose=binomialCoefficient = Factorial(m)/(Factorial(n) * Factorial(m-n))
	}
	Matrix Function_matrix_aux_second(number_of_rows,number_of_combinations);
	Function_matrix_aux_second=0;
	
	int k=0;
	for(int i=1;i<=(number_of_functions+number_of_combinations);i++){
		
		bitset<50> mybitset(i);
		if(mybitset.count()>1){
			k++;
			string function_name_tmp="";
			if(Functional_study=='U'){
				function_name_tmp="Cluster_";
			}
			if(Functional_study=='S'){
				function_name_tmp="";
			}

			vector<int> clusters_belonging_vector_tmp;
			for (size_t j=0; j<number_of_functions; j++){
				if(mybitset.test(j)){
					for(int r=1;r<=number_of_rows;r++){
						Function_matrix_aux_second(r,k)+=Function_matrix(r,(int)j+1);
					}
					char number_of_function_char[200];
					sprintf (number_of_function_char,"%d",(int)j+1);//converts the "int type" into "char type"
					string number_of_function_string="";
					if(Functional_study=='U'){
						number_of_function_string=number_of_function_char;//converts the "char type" into "string type"
					}
					if(Functional_study=='S'){
						number_of_function_string=function_names_vector.at((int)j);
					}
						function_name_tmp+=number_of_function_string+"_";
					clusters_belonging_vector_tmp.push_back((int)j+1);
				}
			}
			function_names_vector.push_back(function_name_tmp);
			function_information_structure tmp;
			tmp.name=function_name_tmp;
			tmp.size= (int) Function_matrix_aux_second.Column(k).Sum();
			tmp.clusters_belonging_vector=clusters_belonging_vector_tmp;
			function_information.push_back(tmp);
		}
	}
	Function_matrix=Function_matrix|Function_matrix_aux_second;
	Function_matrix_aux_second.Release();
	number_of_functions+=number_of_combinations;
}
// copy(function_names_vector.begin(), function_names_vector.end(), ostream_iterator<string>(fs,"\n"));
// fs << endl;
// fs << "Function_matrix: " << endl << setw(2) << setprecision(0) << Function_matrix;


// copy(protein_name_vector.begin(), protein_name_vector.end(), ostream_iterator<string>(fs,"\n"));
// fs << endl;
Matrix Function_Coordinates_indirect(number_of_functions,V_Reduced.Ncols());
Function_Coordinates_indirect=0;
if(Analysis==1||Analysis==2){
	DiagonalMatrix Dr_f(Function_matrix.Nrows());
	for (int i=1;i<=Function_matrix.Nrows();i++ ){
		Dr_f(i,i)=1/(sqrt(sum_of_elements_by_row_vector[i-1]/total_sum));
	}
	DiagonalMatrix Dc_f(Function_matrix.Ncols());
	for (int i=1;i<=Function_matrix.Ncols();i++ ){
		Dc_f(i,i)=1/Function_matrix.Column(i).Sum();
	}
	Matrix Function_matrix_tranf=(Dr_f*(Function_matrix*Dc_f));
	if(Analysis==1){Function_Coordinates_indirect=Function_matrix_tranf.t()*V_Reduced;}
	if(Analysis==2){Function_Coordinates_indirect=((Function_matrix_tranf.t()*V_Reduced)*D_adjusted_sqrt)*D_Reduced_sqrt_inv;}
}
if(Analysis==3){Function_Coordinates_indirect=Function_matrix.t()*V_Reduced;}

if(Verbose=="YES"){
	fs << "Cluster_Coordinates: " << endl;
	fs << "Axes\t1st\t2nd\t3rd\t4th\t5th\t6th\t7th\t8th\t9th\t10th" << endl;
	for (int i=0; i<Function_Coordinates_indirect.Nrows()-1;i++){
		fs << "ClusCoord:" << "\t" << function_names_vector.at(i) << "\t";
		for (int j=Function_Coordinates_indirect.Ncols()-1;j>Function_Coordinates_indirect.Ncols()-1-9; j--){
			fs << Function_Coordinates_indirect[i][j] << "\t";
		}
		fs << Function_Coordinates_indirect[i][Function_Coordinates_indirect.Ncols()-1-9] << endl;
	}
	fs << endl;
}
// copy(function_names_vector.begin(), function_names_vector.end(), ostream_iterator<string>(fs,"\n"));
// fs << endl;
// fs << "Function_matrix: " << endl << setw(2) << setprecision(0) << Function_matrix;

// fs << "Function Coordinates_indirect: " << endl << setw(10) << setprecision(5) << Function_Coordinates_indirect;

// Matrix Centroid_Row_Coordinates_standard(number_of_functions,Row_Coordinates_standard.Ncols());
// Centroid_Row_Coordinates_standard=0;
// for (int i=0;i<number_of_functions;i++){
//     for (int n=0;n<number_of_sequences;n++ ){
// 	if(Function_matrix[n][i]==1){
// 		for (int j=0;j<Row_Coordinates_standard.Ncols();j++){
// 			Centroid_Row_Coordinates_standard[i][j]+=Row_Coordinates_standard[n][j]/Function_matrix.Column(i+1).Sum();
// 		}
// 	}
//     }
// }
// fs << "Centroid_Row_Coordinates_standard" << endl << setw(10) << setprecision(5) << Centroid_Row_Coordinates_standard;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Matrix B;
vector<int> position_aminoacid_interaction_vector_first_element;
vector<int> position_aminoacid_interaction_vector_second_element;
vector<double> sum_of_elements_by_column_interaction_vector;
Matrix Interaction_Column_Coordinates;
if(Interacton_option=="YES"){
	interaction_coding(
		A,
		position_aminoacid_vector,
		sum_of_elements_by_column_vector,
		flag_for_interaction,
		B,
		position_aminoacid_interaction_vector_first_element,
		position_aminoacid_interaction_vector_second_element,
		sum_of_elements_by_column_interaction_vector
	);
	A.Release();//The memory used by A is returned to the system, to improve efficiency
	
	// CONSTRUCTING THE MATRIX WITH THE COORDINATES OF EXTRA INTERACTION COLUMNS (by rows) for each EIGENFACTOR (by columns)
	Interaction_Column_Coordinates.ReSize(position_aminoacid_interaction_vector_first_element.size(),V_Reduced.Ncols());
	Interaction_Column_Coordinates=0;
	if(Analysis==1||Analysis==2){	
		DiagonalMatrix Dr_i(B.Nrows());
		for (int i=1;i<=B.Nrows();i++ ){
			Dr_i(i,i)=1/(sqrt(sum_of_elements_by_row_vector[i-1]/total_sum));
		}
		DiagonalMatrix Dc_i(B.Ncols());
		for (int i=1;i<=B.Ncols();i++ ){
			Dc_i(i,i)=1/sum_of_elements_by_column_interaction_vector[i-1];
		}
		Matrix B_tranf=(Dr_i*(B*Dc_i));
		if(Analysis==1){Interaction_Column_Coordinates=B_tranf.t()*V_Reduced;}
		if(Analysis==2){Interaction_Column_Coordinates=((B_tranf.t()*V_Reduced)*D_adjusted_sqrt)*D_Reduced_sqrt_inv;}
	}
	if(Analysis==3){Interaction_Column_Coordinates=B.t()*V_Reduced;}
}
//***********************************************************************************************************************
if(Functional_study=='S')
{
	Number_of_axes=Column_Coordinates.Ncols();
	if(All_axes_in_column_distances_assessments!="YES")
	{
		if(Number_axes_option=='F'){
			if(Fixed_number_of_axes>V_Reduced.Ncols()){Number_of_axes=V_Reduced.Ncols();}
			else{Number_of_axes=Fixed_number_of_axes;}
		}
		else
		{
			if(Number_axes_option=='I'){	// Selection of the number of axes by the Wilcoxon test
				Number_of_axes=1;
				double* Wilcoxon_pvalues;
				int Wilcox_error_couldnt_open;
				Wilcoxon_test(
					order,
					Row_Coordinates_principal_indirect,
					tmp_directory,
					exec_directory,
					Wilcoxon_pvalues,
					Maximum_number_of_axes,
					Wilcox_error_couldnt_open
				);
				if(Wilcox_error_couldnt_open==1){
				  cout << "\t" << "\tError executing S3det_Wilcoxon_test.R ; its temporary output file couldn't be found or opened" << endl;
				  fs << "\t" << "Error executing S3det_Wilcoxon_test.R ; its temporary output file couldn't be found or opened" << endl;
				  cout << "\t" << "\tWARNING: Please, check that R-package is installed and running according to the conf.h file" << endl;
				  cout << "\t" << "\tWARNING: You might also consider changing the conf.h file and then recompiling S3det" << endl;
				  fs << "\t" << "WARNING: Please, check that R-package is installed and running according to the conf.h file" << endl;
				  fs << "\t" << "WARNING: You might also consider changing the conf.h file and then recompiling S3det" << endl;
				  cout << "\t" << "\tExecution aborted" << endl;
				  fs << "\tExecution aborted" << endl;
				  fs.close();
				  return 1;
				}
				if(Verbose=="YES"){ fs << "Wilcoxon values:" << endl;}
				for(int i=0;i<Maximum_number_of_axes-1;i++){
					if(Verbose=="YES"){ fs << Wilcoxon_pvalues[i] << endl;}
					if(Wilcoxon_pvalues[i]<Wilcoxon_cutoff){Number_of_axes=i+2;}//NOTICE: For the moment it selects the last axes fulfiling the selected cut-off
					//cout << Wilcoxon_pvalues[i] << endl;
				}
			}
		}
	}
	if(All_axes_in_column_distances_assessments=="YES")
	{
		cout << "\t" << "Number of axes selected: " << Number_of_axes << "(all)" <<endl;
		cout << "\t" << "Percentage of initial variance considered informative: 100 %" << endl;
	}else
	{
		cout << "\t" << "Number of axes selected: " << Number_of_axes <<endl;
		cout << "\t" << "Percentage of initial variance considered informative: " << (D_percentage_accumulated(D_percentage_accumulated.Ncols()+1-Number_of_axes,D_percentage_accumulated.Ncols()+1-Number_of_axes))*100 << " %" << endl;
		for(int i=0;i<Number_of_axes;i++){
			cout << "\t\t" << "Percentage of variance explained by the axis "<<  i  + 1 <<  "  (selected): " << (D_percentage_accumulated(D_percentage_accumulated.Ncols()-i,D_percentage_accumulated.Ncols()-i))*100 << " %" << endl;
		}		
		for(int i=Number_of_axes;i<10;i++){
			cout << "\t\t" << "Percentage of variance explained by the axis "<<  i  + 1 <<  "  (not selected): " << (D_percentage_accumulated(D_percentage_accumulated.Ncols()-i,D_percentage_accumulated.Ncols()-i))*100 << " %" << endl;
		}	
	}
} //NOTICE: Tambi\E9n se puede optar aqu\ED por restringir el numero de ejes

int Final_size=Column_Coordinates.Nrows();
if(Interacton_option=="YES"){Final_size+=Interaction_Column_Coordinates.Nrows();}
// cout << "Final_size: " << Final_size << endl;

for (int j=0; j<Column_Coordinates.Nrows();j++){
	double minimum_distance=100000000;
	int assigned_function=-1;
	double* Distance_to_Functions_Vector= new double [Function_Coordinates_indirect.Nrows()];
	for (int i=0; i<Function_Coordinates_indirect.Nrows();i++){
		double square_distance_sum=0;
		for (int k=Column_Coordinates.Ncols()-1;k>(Column_Coordinates.Ncols()-1-Number_of_axes);k--){
			double diference= Column_Coordinates[j][k]-Function_Coordinates_indirect[i][k];
			square_distance_sum+=diference*diference;
		}
		double distance=sqrt(square_distance_sum);if(distance<0.000000001){distance=0;}
		if(distance<minimum_distance){
			minimum_distance=distance;
			assigned_function=i;
		}
		Distance_to_Functions_Vector[i]=distance;
	}
	int Num_of_functions_at_minimum_distance=0;
	for (int i=0; i<Function_Coordinates_indirect.Nrows();i++){
	  if(Distance_to_Functions_Vector[i]==minimum_distance){Num_of_functions_at_minimum_distance++;}
	}
	delete [] Distance_to_Functions_Vector;
	if(Num_of_functions_at_minimum_distance==1){
	  vector<string> tmp_aminoacid_position;
	  Tokenize(position_aminoacid_vector.at(j),tmp_aminoacid_position,"/");
	  position_aa_distance tmp;
	  tmp.distance=minimum_distance;
	  tmp.position_1st= atoi(tmp_aminoacid_position[0].c_str()); ;
	  tmp.aminoacid_1st=tmp_aminoacid_position[1];
	  function_information[assigned_function].elements.push_back(tmp);
	}
}
Matrix Position_rank_per_fingerprint(initial_number_of_positions,number_of_functions-1);Position_rank_per_fingerprint=0;
Matrix Position_coverage_per_group(initial_number_of_positions,initial_number_of_functions);Position_coverage_per_group=0;

for(int i=0; i<function_information.size()-1;i++){
	
	function_information[i].elements_array = new position_aa_distance[function_information[i].elements.size()];
	for(int k=0;k<function_information[i].elements.size();k++){
		function_information[i].elements_array[k]=function_information[i].elements[k];
	}
	combsort_position_aa_distance(function_information[i].elements_array,function_information[i].elements.size());//Sorting the array
	if(Verbose=="YES"){
		fs << "Cluster_name" <<"\t"<< " Num_of_seqs" <<"\t"<< "Total_num_of_residues" << "\n";
		fs << function_information[i].name<<"\t"<<function_information[i].size<<"\t"<<function_information[i].elements.size()<< endl;
		fs << "Rank\tPosition\tAminoAcid\tChi-squared_distance_to_cluster\n" ;
	}
	for(int j=0;j<function_information[i].elements.size();j++){
		if(j==0){function_information[i].elements_array[j].rank=1;}
		else{
			if(float (function_information[i].elements_array[j].distance)==float(function_information[i].elements_array[j-1].distance)){
				function_information[i].elements_array[j].rank=function_information[i].elements_array[j-1].rank;
			}else{
				function_information[i].elements_array[j].rank=function_information[i].elements_array[j-1].rank+1;
			}
		}
		if(Verbose=="YES"){
			fs<< function_information[i].elements_array[j].rank << "\t"<<function_information[i].elements_array[j].position_1st<<"\t"<<function_information[i].elements_array[j].aminoacid_1st<<"\t"<< setprecision(8) <<function_information[i].elements_array[j].distance<<endl;
		}
	}
	
	int mysize= int(floor(function_information[i].elements.size()*selected_pertentage_of_residues));
	//int mysize=int(selected_pertentage_of_residues);
	
	char gap_char[1]; gap_char[1]='#';
	string gap=gap_char;
	
	vector<int> elements_array_selected;
	vector<int> elements_array_selected_rank;
	vector<string> elements_array_selected_amin;
	vector<double> elements_array_selected_dist;

	vector<int> no_gap_covered_ranks;
	int next=0;
	for(int j=0;(j<function_information[i].elements.size() && elements_array_selected.size()<mysize);j++){
		next++;
		if(function_information[i].elements_array[j].aminoacid_1st!=gap){
			elements_array_selected.push_back(function_information[i].elements_array[j].position_1st);
			elements_array_selected_rank.push_back(function_information[i].elements_array[j].rank);
			elements_array_selected_amin.push_back(function_information[i].elements_array[j].aminoacid_1st);
			elements_array_selected_dist.push_back(function_information[i].elements_array[j].distance);
		}
	}
	if(next<function_information[i].elements.size()) {
		for(int j=next;j<function_information[i].elements.size();j++) {
			if(function_information[i].elements_array[j-1].rank==function_information[i].elements_array[j].rank){
				if(function_information[i].elements_array[j].aminoacid_1st!=gap){
					elements_array_selected.push_back(function_information[i].elements_array[j].position_1st);
					elements_array_selected_rank.push_back(function_information[i].elements_array[j].rank);
					elements_array_selected_amin.push_back(function_information[i].elements_array[j].aminoacid_1st);
					elements_array_selected_dist.push_back(function_information[i].elements_array[j].distance);
				}
			}else{ break;}
		}
	}
	
	function_information[i].number_of_selected_residues=elements_array_selected.size();
	function_information[i].elements_array_selected = new int[elements_array_selected.size()];
	function_information[i].elements_array_selected_rank = new int[elements_array_selected.size()];
	function_information[i].elements_array_selected_amin = new string[elements_array_selected.size()];
	function_information[i].elements_array_selected_dist = new double[elements_array_selected.size()];

	if(function_information[i].number_of_selected_residues>0){	
		function_information[i].elements_array_selected[0]=elements_array_selected.at(0);
		function_information[i].elements_array_selected_amin[0]=elements_array_selected_amin.at(0);
		function_information[i].elements_array_selected_dist[0]=elements_array_selected_dist.at(0);
		function_information[i].elements_array_selected_rank[0]=1;
		for(int k=1;k<elements_array_selected.size();k++){
			function_information[i].elements_array_selected[k]=elements_array_selected.at(k);
			function_information[i].elements_array_selected_amin[k]=elements_array_selected_amin.at(k);
			function_information[i].elements_array_selected_dist[k]=elements_array_selected_dist.at(k);
			if(elements_array_selected_rank.at(k)==elements_array_selected_rank.at(k-1)) {
				function_information[i].elements_array_selected_rank[k]=function_information[i].elements_array_selected_rank[k-1];
			}else{
				function_information[i].elements_array_selected_rank[k]=function_information[i].elements_array_selected_rank[k-1]+1;
			}
		}
		for(int k=0;k<elements_array_selected.size();k++){
			if(Position_rank_per_fingerprint(function_information[i].elements_array[k].position_1st,i+1)==0){
				Position_rank_per_fingerprint(function_information[i].elements_array[k].position_1st,i+1)=function_information[i].elements_array[k].rank;
				for(int m=0;m<function_information[i].clusters_belonging_vector.size();m++){
					Position_coverage_per_group(function_information[i].elements_array[k].position_1st,function_information[i].clusters_belonging_vector.at(m))+=1;
				}
			}else{
				if(function_information[i].elements_array[k].rank<Position_rank_per_fingerprint(function_information[i].elements_array[k].position_1st,i+1)){
					Position_rank_per_fingerprint(function_information[i].elements_array[k].position_1st,i+1)=function_information[i].elements_array[k].rank;
					for(int m=0;m<function_information[i].clusters_belonging_vector.size();m++){
						Position_coverage_per_group(function_information[i].elements_array[k].position_1st,function_information[i].clusters_belonging_vector.at(m))+=1;
					}
				}
			}
		}
	}
	//combsort_int(function_information[i].elements_array_selected,function_information[i].number_of_selected_residues);//Sorting the array
	if(Verbose=="YES"){
		fs << "RCS: Cluster_name" <<"\t"<< " Num_of_seqs" <<"\t"<< "Total_num_of_residues" << "\t"<< "Num_of_residues_selected"<< "\n";
		fs << "RCS: " << function_information[i].name<<"\t"<<function_information[i].size<<"\t"<<function_information[i].elements.size()<<"\t"<< function_information[i].number_of_selected_residues<< "\n";
		//copy(function_information[i].clusters_belonging_vector.begin(), function_information[i].clusters_belonging_vector.end(), ostream_iterator<int>(fs,"\t"));fs << endl;
		fs << "RCL: Cluster_name\tRank\tPosition_selected\tAminoAcid\tChi-squared_distance_to_cluster\n" ;
		for(int j=0;j<function_information[i].number_of_selected_residues;j++){
		  fs<< "RCL: " << function_information[i].name <<"\t" << function_information[i].elements_array_selected_rank[j] << "\t" << function_information[i].elements_array_selected[j]  << "\t" << function_information[i].elements_array_selected_amin[j]  << "\t" << function_information[i].elements_array_selected_dist[j] <<endl;
		}
	}
}
if(Verbose=="YES"){
// 	fs << "Position_rank_per_fingerprint " << endl << setw(2) << setprecision(0) << Position_rank_per_fingerprint;
// 	fs << "Position_coverage_per_group "   << endl << setw(2) << setprecision(0) << Position_coverage_per_group;
	fs << "RCM: Position_rank_per_cluster " << endl;
	fs << "RCM: AlignPosition\t";
	for (int j=0; j<number_of_functions-1;j++){
		fs << "    " << function_information[j].name;
	}
	fs << endl;
	int mypointer=0;
	for (int i=0; i<initial_number_of_positions;i++){
	  fs << "RCM: "  << i+1 << "\t";
		int empty=0;
		for (int j=mypointer; j<Column_Coordinates.Nrows();j++){
			vector<string> tmp_aminoacid_position;
			Tokenize(position_aminoacid_vector.at(j),tmp_aminoacid_position,"/");
			int tmp_pos=atoi(tmp_aminoacid_position[0].c_str());
			if(i+1==tmp_pos){
				empty=1;
				mypointer++;
				break;
			}
		}
		for (int j=0; j<number_of_functions-1;j++){
			if(empty==1){
				fs << setw(4) << setprecision(0) << (int)Position_rank_per_fingerprint[i][j];
			}else{
				fs << setw(4) << setprecision(0) << "-";
			}
		}
		fs << endl;
	}
	fs << endl;
}

vector<int> InitialFunctions2check;
vector<int> WholeFunctions2check;
for(int j=0;j<initial_number_of_functions;j++){
	if(function_information[j].size>=minimum_group_size){
		InitialFunctions2check.push_back(j+1);
	}
}
WholeFunctions2check=InitialFunctions2check;
for(int j=initial_number_of_functions;j<number_of_functions-1;j++){
	if(function_information[j].size>=minimum_group_size){
		WholeFunctions2check.push_back(j+1);
	}
}
vector<global_position_structure> Global_selected_positions;
for(int i=1;i<=initial_number_of_positions;i++){
	int complete_partition=1;
	int overlapping_partition=0;
	for(int j=0;j<InitialFunctions2check.size();j++){
		if(Position_coverage_per_group(i,InitialFunctions2check.at(j))==0){complete_partition=0;break;}
		if(Position_coverage_per_group(i,InitialFunctions2check.at(j))>1) {overlapping_partition=1;}
	}
	if((complete_partition==1)&&(overlapping_partition==0)){
		double denominator=0;
		double average_rank=0;
		//fs << "Position:" << i << "\t";
		for(int j=0;j<WholeFunctions2check.size();j++){
			if(Position_rank_per_fingerprint(i,WholeFunctions2check.at(j))>0){
				average_rank+=Position_rank_per_fingerprint(i,WholeFunctions2check.at(j));
				denominator+=1;
				//fs << Position_rank_per_fingerprint(i,WholeFunctions2check.at(j)) << "\t";
			}
		}
		average_rank=average_rank/denominator;
		//fs << setprecision(5) << average_rank << "\t" << denominator << endl;
		if((denominator>1)&&(average_rank<=averagerank_cutoff)){
			global_position_structure tmp;
			tmp.position=i;
			tmp.rank=average_rank;
			tmp.number_of_groups_based=denominator;
			Global_selected_positions.push_back(tmp);
			//fs << "Position:" << i << "\t" << setprecision(5) << average_rank << "\t" << denominator << endl;
		}
	}
	// if((complete_partition==1)&&(overlapping_partition==1)){
	// 	fs << "Overlapping Position:" << i << endl;
	// }
}


int number_of_Global_predicted_positions=Global_selected_positions.size();
fs << "Are there any position predicted as being important for the whole family?: "<< endl << "PR: ";
if(number_of_Global_predicted_positions==0){
	fs << "NO" << endl;
	fs << "UI: There are no positions predicted to be important for the group segregation" << endl;
	cout << "\t" << "There are no positions predicted to be important for the group segregation (SDPs)" << endl;
	return 0;
}
else{fs << "YES" << endl;}

global_position_structure *Global_selected_positions_array;
Global_selected_positions_array=new  global_position_structure [Global_selected_positions.size()];
for(int k=0;k<Global_selected_positions.size();k++){
	Global_selected_positions_array[k]=Global_selected_positions[k];
}
combsort_global_position_structure(Global_selected_positions_array,Global_selected_positions.size());//Sorting the array

fs << setprecision(2);
fs << "Positions predicted (SDPs): " << endl;
fs << "Position\tAverage_rank\tNumber_of_groups_within_the_complete_partition\tList_of_clusters(rank_within_cluster)" << endl;
for(int i=0;i<Global_selected_positions.size();i++){
	fs << "RE: " << Global_selected_positions_array[i].position << "\t" << Global_selected_positions_array[i].rank << "\t" << Global_selected_positions_array[i].number_of_groups_based << "\t" << "Partition_based_on:";
	for (int j=0; j<number_of_functions-1;j++){
		if(Position_rank_per_fingerprint[Global_selected_positions_array[i].position-1][j]>0){
		fs << " " << function_information[j].name <<"(rank:"<<Position_rank_per_fingerprint[Global_selected_positions_array[i].position-1][j]<<")";
		}
	}
	fs << endl;

	//fs << "RE: " << Global_selected_positions_array[i].position << " " << setprecision(5) << Global_selected_positions_array[i].rank << " " << Global_selected_positions_array[i].number_of_groups_based << endl;
	//fs << "RE: " << Global_selected_positions_array[i].position  << endl;
}
fs << endl;
fs << "UI: Number of positions predicted to be important for the group segregation: " << number_of_Global_predicted_positions << endl << endl;
cout << "\t" << "Number of positions predicted to be important for the group segregation: " << number_of_Global_predicted_positions << endl;

// cout << "\t" << "Number of Globally important positions predicted: " << number_of_Global_predicted_positions << endl;



fs << "Are there any conserved position at "<< conservation_threshold <<" % of identity?: "<< endl << "CN: ";
int number_of_Conserved_positions=Conserved_positions.size();
if(number_of_Conserved_positions==0){
	fs << "NO" << endl;
	fs << "UI: There are no conserved positions at "<< conservation_threshold <<" % of identity" << endl;
	return 0;
}
else{fs << "YES" << endl;}
// sort(Conserved_positions.begin(),Conserved_positions.end());
fs << "Conserved positions: " << endl;
for(int i=0;i<Conserved_positions.size();i++){
	fs << "CP: " << Conserved_positions[i] << endl;
}
fs << endl;
fs << "UI: Number of conserved positions (at "<< conservation_threshold <<" % of identity): " << number_of_Conserved_positions << endl;
cout << "\t" << "Program finished OK. All outputs written to outfile" << endl;


//____________________________________________________________________________________________________________________
// position_aminoacid_structure *position_aminoacid_checking;
// position_aminoacid_checking=new position_aminoacid_structure[position_aminoacid_vector.size()];
// for(int k=0;k<position_aminoacid_vector.size();k++){
// 	position_aminoacid_checking[k].position_aminoacid=position_aminoacid_vector[k];
// }
// 
// for (int i=0; i<Function_Coordinates_indirect.Nrows();i++){
// 
// 	for(int k=0;k<position_aminoacid_vector.size();k++){
// 		position_aminoacid_checking[k].checked=false;
// 	}
// 
// 	cout << "FUNCTION NUMBER " << i+1 << endl;
// 	distance_aaposition_pair *array;
// 	array = new distance_aaposition_pair[Final_size];
// 	for (int j=0; j<Column_Coordinates.Nrows();j++){
// 		double square_distance_sum=0;
// 		for (int k=Column_Coordinates.Ncols()-1;k>(Column_Coordinates.Ncols()-1-Number_of_axes);k--){
// 			double diference= Column_Coordinates[j][k]-Function_Coordinates_indirect[i][k];
// 			square_distance_sum+=diference*diference;
// 		}
// 		array[j].distance=sqrt(square_distance_sum);if(array[j].distance<0.000000001){array[j].distance=0;}
// 		array[j].aminoacid_position_1st=j;
// 		array[j].aminoacid_position_2nd=0;
// 		array[j].interaction=false;
// 	}
// 	if(Interacton_option=="YES"){
// 		for (int j=0; j<Interaction_Column_Coordinates.Nrows();j++){
// 			double square_distance_sum=0;
// 			for (int k=Interaction_Column_Coordinates.Ncols()-1;k>(Interaction_Column_Coordinates.Ncols()-1-Number_of_axes);k--){
// 				double diference= Interaction_Column_Coordinates[j][k]-Function_Coordinates_indirect[i][k];
// 				square_distance_sum+=diference*diference;
// 			}
// 			array[j+Column_Coordinates.Nrows()].distance=sqrt(square_distance_sum);if(array[j+Column_Coordinates.Nrows()].distance<0.000000001){array[j+Column_Coordinates.Nrows()].distance=0;}
// 			array[j+Column_Coordinates.Nrows()].aminoacid_position_1st=position_aminoacid_interaction_vector_first_element.at(j);
// 			array[j+Column_Coordinates.Nrows()].aminoacid_position_2nd=position_aminoacid_interaction_vector_second_element.at(j);
// 			array[j+Column_Coordinates.Nrows()].interaction=true;
// 		}
// 	}
// 	combsort_distance_aaposition_pair(array,Final_size);//Sorting the array
// 	fs << "Best scores for function " << function_names_vector.at(i) << endl;
// 	
// 	int number_of_different_clusters=0;
// 	int cluster_start=0;
// 	int flag=0;
// 	while (cluster_start<Final_size-1){
// 		number_of_different_clusters++;
// 		vector<int> ISOdistance_cluster;
// 		ISOdistance_cluster.push_back(cluster_start);
// 		
// 		int current=cluster_start;
// 		if((current<Final_size-1 && (flag==0))){
// 			while( (array[current].distance==array[current+1].distance) && (flag==0)){
// 				ISOdistance_cluster.push_back(current+1);
// 				current=current+1;
// 				if(current==Final_size-1){flag=1;}
// 			}
// 			cluster_start=current+1;
// 		}
// 		for(int i=0;i<ISOdistance_cluster.size();i++){
// 			if(array[ISOdistance_cluster[i]].interaction==false){position_aminoacid_checking[array[ISOdistance_cluster[i]].aminoacid_position_1st].checked=true;}
// 		}
// 		for(int i=0;i<ISOdistance_cluster.size();i++){
// 			if(array[ISOdistance_cluster[i]].interaction==true){
// 				if(
// 				(position_aminoacid_checking[array[ISOdistance_cluster[i]].aminoacid_position_1st].checked==false) &&
// 				(position_aminoacid_checking[array[ISOdistance_cluster[i]].aminoacid_position_2nd].checked==false)
// 				){
// 				 fs <<  number_of_different_clusters << "\t" << setw(10) << setprecision(8) << array[ISOdistance_cluster[i]].distance
// 				 << "\t"<< setw(10) << setprecision(8) << position_aminoacid_vector[array[ISOdistance_cluster[i]].aminoacid_position_1st]
// 				 << "\t"<< setw(10) << setprecision(8) << position_aminoacid_vector[array[ISOdistance_cluster[i]].aminoacid_position_2nd]
// 				 << endl;
// 				 }
// 			}else{
// 				fs <<  number_of_different_clusters << "\t"  << setw(10) << setprecision(8) << array[ISOdistance_cluster[i]].distance
// 				<<  "\t"<< setw(10) << setprecision(8) << position_aminoacid_vector[array[ISOdistance_cluster[i]].aminoacid_position_1st]
// 				<<  endl;
// 			}
// 		}
// 	}
// 
// /*	for (int j=0;j<Final_size;j++){
// 	//for (int j=0;j<15;j++){
// 		if(array[j].interaction==true){
// 			 fs     << setw(10) << setprecision(8) << array[j].distance
// 			 << "\t"<< setw(10) << setprecision(8) << position_aminoacid_vector[array[j].aminoacid_position_1st]
// 			 << "\t"<< setw(10) << setprecision(8) << position_aminoacid_vector[array[j].aminoacid_position_2nd]
// 			 << endl;
// 		}else{
// 			fs      << setw(10) << setprecision(8) << array[j].distance
// 			<<  "\t"<< setw(10) << setprecision(8) << position_aminoacid_vector[array[j].aminoacid_position_1st]
// 			<<  endl;
// 		}
// 	}
// 	fs << endl << endl << endl;	*/
// }

// cout << "Fin: " << endl;
// Closing output file
fs.close();

}
catch(Exception){cout << Exception::what() << endl;}

return 0;

}
