 //NOTICE!!! HAY QUE COMPLETAR CON LA PLASTICIDAD QUE PUEDE TENER EL FORMATO FASTA (VARIOS NOMBRES PARA UN MISMO AMINOACIDO, letras minúsculas O INTRODUCCIÓN DE COMENTARIOS) ASÍ COMO PARA OTROS FORMATOS.

//Hay que quitarle el ">" al nombre de las secuencias



//NOTICE: TO compile the program: g++ MCdet_scores.cpp -o MCdet_scores.exe -lboost_regex  -L. -lnewmat -lm  /usr/lib64/libkdeinit_krandrinithack.so
//NOTICE: To compile the program statically: g++ MCdet_scores.cpp -o MCdet_scores.exe  /usr/lib/libboost_regex.a -L. -lnewmat -lm  /usr/lib64/libkdeinit_krandrinithack.so
//NOTICE: To execute the program: ./MCdet.exe    Multiple_alignment_Fasta_file    Output_file_chosen
//____________________________________________________________________________________________________________________________________________

// #define FECHA_INICIAL __DATE__ 
// #define HORA_INICIAL __TIME__


#include <stdio.h> //Necessary for the sprintf function used in converting the "int type" into "char type"
//#include <iostream>  //Already included in the " #define WANT_STREAM " of the newmat10B library
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <math.h>
#include <ctime>
#include <cstdlib>

using namespace std;


#include <boost/regex.hpp> //Necessary for the BOOST_REGEX LIBRARY - C Library dealing with Regular Expressions

//The following statements are necessary for the NEWMAT10 LIBRARY - a Matrix Library with "eigen tools" 
#define WANT_STREAM
#include "newmatap.h"
#include "newmatio.h"              
#ifdef use_namespace
using namespace RBD_LIBRARIES;
#endif

//____________________________________________________________________________________________________________________________________________

struct distance_aaposition_pair{
	double distance;
	string aminoacid_position;
	double pvalue;
	double pvalue_separated;
	
};

//________________________________________________________________________________________________________________________
/*
                                               FUNCTION DECLARATION
*/

//This function converts a given string in a vector called "tokens" splitting by a given delimiter
void Tokenize(const string& str, vector<string>& tokens, const string& delimiters =" ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}


//Boolean funtion who returns TRUE if the string s starts with a > (this will be used to indicate if a line in a FASTA input file contains the name of a protein 
bool fasta_protein_name (const std::string& s)
{
   static const boost::regex e("^>(.*)");
   return boost::regex_match(s, e);
}

bool fasta_protein_sequence (const std::string& s)
{
	//static const boost::regex e("^(\w|.|-)(.*)");
	static const boost::regex e0("(.*)([ABCDEFGHIKLMNPQRSTUVWYZXabcdefghiklmnpqrstuvwyzx\\*\\-\\.])(.*)");
	static const boost::regex e1("(.*)([^ABCDEFGHIKLMNPQRSTUVWYZXabcdefghiklmnpqrstuvwyzx\\*\\-\\.])(.*)");
	bool re0= boost::regex_match(s, e0);
	bool re1= boost::regex_match(s, e1);
	if(re0 and not re1){return 1;}else{return 0;}
}


int setsrand()
{
FILE *f;
char l[200];

char tim[3];
f=popen("date","r");
if (f==NULL) return -1;
if(fgets(l,999,f)==0) return -2;
pclose(f);
tim[0]=l[15]; tim[1]=l[18]; tim[2]='\0';
srand(atoi(tim));
return 0;
}

static int newGap(int gap) {
  gap = (gap * 10) / 13;
  if (gap == 9 || gap == 10)
    gap = 11;
  if (gap < 1)
    gap = 1;
  return gap;
}

static void combsort(double a[], int aSize) {
  int gap = aSize;
  for (;;) {
    gap = newGap(gap);
    bool swapped = false;
    for (int i = 0; i < aSize - gap; i++) {
      int j = i + gap;
      if (a[i] > a[j]) {
        std::swap(a[i], a[j]);
        swapped = true;
      }
    }
    if (gap == 1 && !swapped)
      break;
  }
}

//_______________________________________________________________________________________________________________________
//                                             STARTS MAIN
//_______________________________________________________________________________________________________________________


int main (int argc, char *argv[]) {
  
//***********************************************************************************************************************
/*
                                           VARIABLE DECLARATION
*/

vector<string> protein_name_vector; //it will contain the names of the proteins taking part in the alignment
int number_of_sequences; //number of sequences of the alignment
vector<string> protein_sequence_vector; //it will contain the sequences of the proteins taking part in the alignment
int number_of_positions; //number of positions of the alignment

vector<string> position_aminoacid_vector;//vector with elements the 21xL position-aminoacid flags, where L is the length of the alignment

//char letter_matrix[number_of_sequences][number_of_positions];letter_matrix: a bidimensional char array with the "amino acid letters" of the alignment. It will be declared later, otherwise compilation errors will abort the program


char aminoacid[]={'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','#'};

string line;//it will be used in the input file reading process
string protein_sequence_part="";//it will be used in constructing the protein_sequence_vector


//***********************************************************************************************************************
// cout << "Inicio: " << endl;
// cout << FECHA_INICIAL << endl;
// cout << HORA_INICIAL << endl;

//Create output file with the chosen name

int help=0;
char * fasta_file;
char * out_file;
char * Function_coding_file;
if(argc<=2){
	cout << endl;
	cout << "   Required arguments: " << endl<< endl;

	cout << "     -i   infile                   Multiple Sequence Alignment in fasta format"      << endl;
	cout << "     -o   outfile                  Text file where program output will be printed"   << endl ;
	cout << "     -f   Function coding file     Text file where a binary matrix codes for a super-"<< endl;
	cout << "                                     vised classification."                           << endl<<endl;
	
	cout << "   Optional arguments:" << endl<< endl;
	
	cout << "     -p                            activates p-value calculation"<< endl<< endl;


	return 1;
}
string score="";//converts the "char type" position into "string type"

static const boost::regex nonvalidsyntax("^(-)([iofp])([^\\s])(.+)");
static const boost::regex nonvalidoption("^(-)([^iofp]*)(.*)");
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
			cout << "\t" << "Error: Option "  << argv[i] << " requires a file to write the output" << endl;
			cout << "\t" << "Please execute program without arguments to see \"help\" on valid options" << endl;
			cout << "\t" << "Program not executed" << endl;
			return 1;
		}
	}else if(strcmp("-f",argv[i])==0){
		myflag_fset=1;
		if(boost::regex_match(argv[i+1],matchsomecharacterotherthanhyphen)){
			Function_coding_file=strdup(argv[i+1]);
		}else{
			cout << "\t" << "Error: Option "  << argv[i] << " requires a function coding file" << endl;
			cout << "\t" << "Please execute program without arguments to see \"help\" on valid options" << endl;
			cout << "\t" << "Program not executed" << endl;
			return 1;
		}
	}else if(strcmp("-p",argv[i])==0){
		score="-p";
	}else if(boost::regex_match(argv[i],nonvalidsyntax)){
		cout << "\t" << string (argv[i]) << " is not a valid syntax" << endl;
		cout << "\t" << "Please separate with a white space the options from their arguments (e.g. -i dir/infile -o dir/outfile -f dir/functionfile)" << endl;
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
if(myflag_fset==0){
		cout << "\t" << "Error: No function coding file provided" << endl;
		cout << "\t" << "Please indicate a function file  \"-f dir/functionfile\"" << endl;
		cout << "\t" << "Execute program without arguments to see \"help\" on required arguments" << endl;
		cout << "\t" << "Program not executed" << endl;
		return 1;
}


//***********************************************************************************************************************
//                                     READING FASTA FORMAT INPUT
ofstream fs(out_file); 
if (fs.fail()){
	string out_file_string=out_file;
	cout << "\t" << "\tError opening file " <<  out_file_string << endl << "\t\tFile not found or couldn't be opened" << endl;
	cout << "\t" << "\tExecution aborted" << endl;
	fs.close();
	return 1;
}

//Open fasta format input file wich is indicated as the fist argument in the prompt
ifstream infile (fasta_file);
if(infile.fail()){
	string fasta_file_string=fasta_file;
	cout << "\t" << "\tError loading file " <<  fasta_file_string << endl << "\t\tFile not found or couldn't be opened" << endl;
	fs << "\t" << "\tError loading file " <<  fasta_file_string << endl << "\t\tFile not found or couldn't be opened" << endl;
	cout << "\t" << "\tExecution aborted" << endl;
	fs << "\tExecution aborted" << endl;
	fs.close();
	return 1;
}
static const boost::regex e5("^$");
int error_loading_fasta=0;
int number_of_positions_tmp=0;
int error_size_of_sequences=0;
//Reading line by line the input file (NOTICE: Getline eliminates the 'endline' character)
while (getline(infile,line,'\n')) {
	//if the line starts with a ">" the whole line is taken as the protein name an added to the sequence_name_vector:
	if ( fasta_protein_name(line)  ){
		if ( fasta_protein_sequence(protein_sequence_part) ) {
			protein_sequence_vector.push_back(protein_sequence_part);
			
			if(protein_sequence_vector.size()==1){number_of_positions_tmp=protein_sequence_part.length();}
			if(protein_sequence_vector.size()>1){
				if(protein_sequence_part.length()!=number_of_positions_tmp){
					error_size_of_sequences=1;
					string fasta_file_string=fasta_file;
					cout << "\t" << "\tError loading file " <<  fasta_file_string << endl << "\t\tDifferent sequence sizes within the alignment" << endl;
					fs << "\tError loading file " <<  fasta_file_string << endl << "\t\tDifferent sequence sizes within the alignment" << endl;
					cout << "\t" << "\tExecution aborted" << endl;
					fs << "\tExecution aborted" << endl;
					fs.close();
					return 1;
				}
			}
		}
		protein_sequence_part="";//reinitializing the construction of the protein sequence
		line.erase(0,1);
		protein_name_vector.push_back(line);
	}	
	else if ( fasta_protein_sequence(line)  ){
		protein_sequence_part=protein_sequence_part+line;
	} else if(boost::regex_match(line, e5)){
		continue;
	}else {
		error_loading_fasta=1;
	}
}
// just to load into the protein_sequence_vector the last sequence
if ( fasta_protein_sequence(protein_sequence_part) ) {
	protein_sequence_vector.push_back(protein_sequence_part);
	if(protein_sequence_vector.size()==1){number_of_positions_tmp=protein_sequence_part.length();}
	if(protein_sequence_vector.size()>1){
		if(protein_sequence_part.length()!=number_of_positions_tmp){
			error_size_of_sequences=1;
			string fasta_file_string=fasta_file;
			cout << "\t" << "\tError loading file " <<  fasta_file_string << endl << "\t\tDifferent sequence sizes within the alignment" << endl;
			fs << "\tError loading file " <<  fasta_file_string << endl << "\t\tDifferent sequence sizes within the alignment" << endl;
			cout << "\t" << "\tExecution aborted" << endl;
			fs << "\tExecution aborted" << endl;
			fs.close();
			return 1;
		}
	}
}
//Closing input file
infile.close();

if(error_loading_fasta==1){
	string fasta_file_string=fasta_file;
	cout << "\t" << "\tError loading file " <<  fasta_file_string << endl << "\t\tInput alignment is not in .mfa format" << endl;
	fs << "\tError loading file " <<  fasta_file_string << endl << "\t\tInput alignment is not in .mfa format" << endl;
	cout << "\t" << "\tExecution aborted" << endl;
	fs << "\tExecution aborted" << endl;
	fs.close();
	return 1;
}


number_of_sequences = protein_name_vector.size();
number_of_positions = number_of_positions_tmp;

cout << "\t"<< "Number of sequences: " << number_of_sequences << "\n";
cout << "\t"<< "Number of positions: " << number_of_positions << "\n";


char letter_matrix[number_of_sequences][number_of_positions];//	letter_matrix: a bidimensional char array with the "amino acid letters" of the alignment
//Constructing a bidimensional char array called "letter_matrix" with the "amino acid letters" of the alignment 
for (int i=0; i<(number_of_sequences);i++ ){
	
	const char* puntero= protein_sequence_vector[i].c_str();//it converts a string in a char array
	for (int j=0; j<(number_of_positions);j++){
		
		if ((puntero[j]=='B')|(puntero[j]=='U')|(puntero[j]=='Z')|(puntero[j]=='X')|(puntero[j]=='*')|(puntero[j]=='.')|(puntero[j]=='-')){
			letter_matrix[i][j]='#';
		} else {
			letter_matrix[i][j]=puntero[j];
		}
		//fs << letter_matrix[i][j] << " ";
		
	}
	//fs << endl;
}

/*
'B','U','Z','X' and '*' are treated as gap positions

    *A  alanine                         *P  proline
    B  aspartate or asparagine         *Q  glutamine
    *C  cystine                         *R  arginine
    *D  aspartate                       *S  serine
    *E  glutamate                       *T  threonine
    
    *F  phenylalanine                   U  selenocysteine
    *G  glycine                         *V  valine
    *H  histidine                       *W  tryptophan
    *I  isoleucine                      *Y  tyrosine
    *K  lysine                          Z  glutamate or glutamine
    *L  leucine                         X  any
    *M  methionine                      *  translation stop
    *N  asparagine                      -  gap of indeterminate length
*/


//**************************************************************************************************************************

//Constructing the "position_aminoacid_vector" with elements the 21xL position-aminoacid flags, where L is the length of the alignment
for (int p=0; p<(number_of_positions);p++){
	
	char position_char[200];
	sprintf (position_char,"%d",p+1);//converts the "int type" position into "char type"
	string position_string= position_char;//converts the "char type" position into "string type"
	
	for (int i=0;i<21;i++ ){
		string position_aminoacid_string=position_string+aminoacid[i];//fusion "position" + "aminoacid" strings
		position_aminoacid_vector.push_back(position_aminoacid_string);
		//fs << position_aminoacid_vector.back() << " ";
	}
	//fs << endl;
}

//*********************************************************************************************************************************
//                                         READING FUNCTION MATRIX
//Open the input file indicated as the second argument in the prompt which will contain the binary classification of the proteins of the alignment (by rows) corresponding a diferent functions (arranged by columns)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	Matrix Function_matrix;
	vector<string> function_names_vector;
	int number_of_functions;
	double flag_for_interaction;
	int error_file_is_not_in_proper_format=0;
	int error_sequence_not_found_in_Function_coding_file=0;
	int flag_sequence_with_not_assigned_function=0;
	int flag_sequence_with_sum_of_assigned_functions_different_from_1=0;
	string protein_name_causing_error="";
	line="";
	ifstream infile2 (Function_coding_file);
	if(infile2.fail()){
		string Function_coding_file_string=Function_coding_file;
		cout << "\t" << "\tError loading file " <<  Function_coding_file_string << endl << "\t\tFile not found or couldn't be opened" << endl;
		fs << "\t" << "\tError loading file " <<  Function_coding_file_string << endl << "\t\tFile not found or couldn't be opened" << endl;
		cout << "\t" << "\tExecution aborted" << endl;
		fs << "\t"<< "\tExecution aborted" << endl;
		fs.close();
		return 1;
	}

	getline(infile2,line,'\n');
	Tokenize(line,function_names_vector,"\t");
	number_of_functions=function_names_vector.size();
	infile2.close();
	
	Matrix Function_matrix_aux(protein_name_vector.size(),number_of_functions);
	Function_matrix_aux=0;
	
	static const boost::regex e0("0");
	static const boost::regex e1("1");
	static const boost::regex e2("(^[a-zA-Z]|^[0-9]|^>)(.*)");
	static const boost::regex e3("(^[0-9]+)(\\.*)([0-9]*)");
	
	for (int i=0;i<protein_name_vector.size();i++){
		ifstream infile2 (Function_coding_file);
		int flag=0;
		getline(infile2,line,'\n');
		while ( (getline(infile2,line,'\n')) && (flag==0) ){
			if (boost::regex_match(line,e2)){
				vector<string> tokens;
				Tokenize(line,tokens,"\t");
				if(tokens.size()!=number_of_functions+1){
					error_file_is_not_in_proper_format=1;
				}
				if(protein_name_vector.at(i)==tokens.at(0)){
					flag=1;
					double sum_per_row=0;
					for (int j = 1; j < tokens.size(); j++) {
						if(boost::regex_match(tokens.at(j),e3)){
							const char* temp= tokens.at(j).c_str();
							double temp_double=atof(temp);
							if((temp_double>=0)and(temp_double<=1)){
								Function_matrix_aux[i][j-1]=temp_double;
								sum_per_row+=temp_double;
								/*if (boost::regex_match(tokens[j],e1)){
									Function_matrix_aux[i][j-1]=1;
								}else if (boost::regex_match(tokens[j],e0)){
									Function_matrix_aux[i][j-1]=0;
								}*/
							}else{
								error_file_is_not_in_proper_format=1;
							}
						}else{
							error_file_is_not_in_proper_format=1;
						}
					}
					if(sum_per_row==0){
						flag_sequence_with_not_assigned_function=1;
						protein_name_causing_error=protein_name_vector.at(i);
					}
					if(sum_per_row!=1){
						flag_sequence_with_sum_of_assigned_functions_different_from_1=1;
						protein_name_causing_error=protein_name_vector.at(i);
					}
				}
			}
		}
		infile2.close();
		if(flag==0){
			protein_name_causing_error=protein_name_vector.at(i);
			error_sequence_not_found_in_Function_coding_file=1;
		}
	}	
	flag_for_interaction=Function_matrix_aux.Column(1).Sum();
	for(int i=1;i<Function_matrix_aux.Ncols();i++){
		if (flag_for_interaction > Function_matrix_aux.Column(i+1).Sum()){
			flag_for_interaction = Function_matrix_aux.Column(i+1).Sum();
		}
	}
	Function_matrix=Function_matrix_aux;
	Function_matrix_aux.Release();	
	string function_file_string=Function_coding_file;
	if(error_file_is_not_in_proper_format==1){
		cout << "\t" << "\tError loading file: " <<  function_file_string << " has not a proper format" << endl;
		fs << "\t"<< "\tError loading file " <<  function_file_string << " has not a proper format" << endl;
		cout << "\t" << "\tExecution aborted" << endl;
		fs << "\t"<< "\tExecution aborted" << endl;
		fs.close();
		return 1;
	}
// 	if(flag_sequence_with_not_assigned_function==1){
// 		cout << "\t" << "\tError loading file " <<  function_file_string << endl;
// 		cout << "\t" << "\tSequence " << myprotein_name_causing_error << " has not assigned function in " << function_file_string << endl;
// 		fs << "\tError loading file " <<  function_file_string << endl;
// 		fs << "\tSequence " << myprotein_name_causing_error << " has not assigned function in " << function_file_string << endl;
// 		cout << "\t" << "\tExecution aborted" << endl;
// 		fs << "\tExecution aborted" << endl;
// 		fs.close();
// 		return 1;
// 	}
// 	if(myflag_sequence_with_sum_of_assigned_functions_different_from_1==1){
// 		cout << "\t" << "\tError loading file " <<  function_file_string << endl;
// 		cout << "\t" << "\tSequence " << myprotein_name_causing_error << " assigned functions values sum up a number different from 1 " << function_file_string << endl;
// 		fs << "\tError loading file " <<  function_file_string << endl;
// 		fs << "\tSequence " << myprotein_name_causing_error << " assigned functions values sum up a number different from 1 " << function_file_string << endl;
// 		cout << "\t" << "\tExecution aborted" << endl;
// 		fs << "\tExecution aborted" << endl;
// 		fs.close();
// 		return 1;
// 	}
	if(error_sequence_not_found_in_Function_coding_file==1){
		cout << "\t" << "\tError loading file " <<  function_file_string << endl;
		cout << "\t" << "\tSequence " << protein_name_causing_error << " not found in " << function_file_string << endl;
		fs << "\t"<< "\tError loading file " <<  function_file_string << endl;
		fs << "\t"<< "\tSequence " << protein_name_causing_error << " not found in " << function_file_string << endl;
		cout << "\t" << "\tExecution aborted" << endl;
		fs << "\t"<< "\tExecution aborted" << endl;
		fs.close();
		return 1;
	}

	
	
	
	
	
fs << "List of functions: "<< endl;
copy(function_names_vector.begin(), function_names_vector.end(), ostream_iterator<string>(fs,"\t"));
fs << endl <<  "______________________________________________________________________________" << endl << endl << endl;
    
// fs << "Number_of_functions: " << number_of_functions << endl;
cout << "\t"<< "Number of functions: " << number_of_functions << endl;



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
int flag_for_interaction=Function_matrix.Column(1).Sum();
for(int i=1;i<Function_matrix.Ncols();i++){
	if (flag_for_interaction > Function_matrix.Column(i+1).Sum()){
		flag_for_interaction = Function_matrix.Column(i+1).Sum();
	}
}
*/

// fs << "Function matrix: "<< endl << setw(1) << setprecision(5) << Function_matrix;
//fs << Function_matrix.Nrows() << "\n";
//fs << Function_matrix.Ncols() << "\n";




//************************************************************************************************************************

cout << "\t"<< "Performing Eigenvalue descomposition..." << endl;

//Constructing the disjunctive matrix A
//________________________________________________
//SOME MORE VARIABLE DECLARATION
Matrix A(number_of_sequences,1);
SymmetricMatrix X;
DiagonalMatrix D,D_Reduced;
Matrix V,V_Reduced;
vector<string> cleaned_position_aminoacid_vector;//vector with aminoacid positions which has at least one non-null element
vector<double> sum_of_elements_by_column_vector;
vector<double> sum_of_elements_by_row_vector(number_of_sequences);
double total_sum=0;

int position_aminoacid_vector_MARKER=-1;
int columns_asigned =0;
//________________________________________________

A=0;

for (int j=0; j<(number_of_positions);j++){
	for (int k=0; k<21; k++){
		
		position_aminoacid_vector_MARKER += 1;
		Matrix C(number_of_sequences,1);//new column to add to A if has at least one non-null element
		C=0;
		int number_of_ones=0;
		for (int i=0;i<(number_of_sequences);i++ ){	
			
			if (letter_matrix[i][j]== aminoacid[k]){ //codifying for each aminoacid type k
				C[i][0]=1;
				sum_of_elements_by_row_vector[i]+=1;
				number_of_ones +=1;
				
			} 
		}
		if (number_of_ones!=0){
			columns_asigned +=1;
			if (columns_asigned > 1){
				 //because the first time you need to overwrite the first default null column (honestly this statement is a shoddy piece of work)
				A=A|C;
			} else {	
				A << C;
			}
			//constructing the vector with aminoacid positions which has at least one non-null element
			
			cleaned_position_aminoacid_vector.push_back(position_aminoacid_vector.at(position_aminoacid_vector_MARKER));
			sum_of_elements_by_column_vector.push_back(number_of_ones);
			total_sum+=number_of_ones;
//			fs << cleaned_position_aminoacid_vector.back() << " ";
		}
	}
//	fs << endl;
}

//************************************************************************************************************************
//							INTERACTION OPTION

// if (argv[4]) {
// 
// 	for (int i=0; i<columns_asigned-1;i++){
// 		//if (sum_of_elements_by_column_vector[i]*4 > flag_for_interaction*3 ){	
// 			for (int j=i+1;j<columns_asigned;j++ ){	
// 				//if (sum_of_elements_by_column_vector[j]*4 > flag_for_interaction*3 ){	
// 					Matrix C(number_of_sequences,1);
// 					C=0;
// 					int number_of_ones=0;
// 					
// 					for (int k=0;k<(number_of_sequences);k++ ){	
// 						C[k][0]= A[k][i]*A[k][j];
// 						if (C[k][0]!=0){
// 							sum_of_elements_by_row_vector[k]+=1;
// 							number_of_ones +=1;
// 						}
// 					}
// 				
// 					if(number_of_ones!=0){
// 						A=A|C;
// 						string interaction_position=cleaned_position_aminoacid_vector.at(i)+cleaned_position_aminoacid_vector.at(j);
// 						cleaned_position_aminoacid_vector.push_back(interaction_position);
// 						sum_of_elements_by_column_vector.push_back(number_of_ones);
// 						total_sum+=number_of_ones;
// 			//			fs << cleaned_position_aminoacid_vector.back() << " ";
// 					}
// 				//}	
// 			}
// 			//fs << endl;
// 		//}
// 	}
// }
//************************************************************************************************************************/


Matrix Z(number_of_sequences,cleaned_position_aminoacid_vector.size()); //Constructing Matrix Z
for (int i=0;i<(number_of_sequences);i++ ){
	for (int j=0;j<cleaned_position_aminoacid_vector.size();j++){
		Z[i][j]=A[i][j]/(sqrt(sum_of_elements_by_row_vector[i])*sqrt(sum_of_elements_by_column_vector[j]));
	}
}
A.Release();//The memory used by A is returned to the system, to improve efficiency

DiagonalMatrix Dc(cleaned_position_aminoacid_vector.size()); //Constructing DiagonalMatrix Dc
for (int i=1;i<(cleaned_position_aminoacid_vector.size()+1);i++ ){
		Dc(i,i)=1/(sqrt(sum_of_elements_by_column_vector[i-1]/total_sum));
}



X << Z*Z.t();
EigenValues(X,D,V);
X.Release();

int number_of_non_null_eigen_values=0;
for (int i=number_of_sequences;((i>0) and (D(i,i)>0.0000000001));i--){
	number_of_non_null_eigen_values+=1;
}
// cout << "Number of non null eigen values: "<< number_of_non_null_eigen_values << "\n";

V_Reduced << V.SubMatrix(1,number_of_sequences,(number_of_sequences+1)-number_of_non_null_eigen_values,number_of_sequences-1);
V.Release();
D_Reduced << D.SubMatrix((number_of_sequences+1)-number_of_non_null_eigen_values,number_of_sequences-1,(number_of_sequences+1)-number_of_non_null_eigen_values,number_of_sequences-1);
D.Release();

Matrix Column_Coordenates=Dc*(Z.t()*V_Reduced);//Matrix wich each file i has the coordinates of position_aminoacid i in factor j (columns arranging the coordinates for each factorj)
Z.Release();

//fs << "V_Reduced: " << endl << setw(8) << setprecision(5) << V_Reduced << "\n" << "\n";
//fs << "D_Reduced: " << endl << setw(8) << setprecision(5) << D_Reduced << "\n" << "\n";
//fs << "Column Coordinates: " << endl << setw(8) << setprecision(5) << Column_Coordenates << "\n" << "\n";
//cout << Column_Coordenates.Nrows() << "\n";
//cout << Column_Coordenates.Ncols() << "\n";
//*********************************************************************************************************************************
/*
SymmetricMatrix Xt;
DiagonalMatrix Dt,Dt_Reduced;
Matrix Vt,Vt_Reduced;

Xt << Z.t()*Z;
EigenValues(Xt,Dt,Vt);

int number_of_non_null_eigen_values_transp=0;
for (int i=cleaned_position_aminoacid_vector.size();((i>0) and (Dt(i,i)>0.0000000001));i--){
	number_of_non_null_eigen_values_transp+=1;
}
cout << "Number of non null eigen values_transp: "<< number_of_non_null_eigen_values_transp << "\n";

Vt_Reduced << Vt.SubMatrix(1,cleaned_position_aminoacid_vector.size(),(cleaned_position_aminoacid_vector.size()+1)-number_of_non_null_eigen_values_transp,cleaned_position_aminoacid_vector.size()-1);
Vt.Release();
Dt_Reduced << Dt.SubMatrix((cleaned_position_aminoacid_vector.size()+1)-number_of_non_null_eigen_values_transp,cleaned_position_aminoacid_vector.size()-1,(cleaned_position_aminoacid_vector.size()+1)-number_of_non_null_eigen_values_transp,cleaned_position_aminoacid_vector.size()-1);
Dt.Release();

fs << "Vt_Reduced: " << endl << setw(8) << setprecision(5) << Vt_Reduced << "\n" << "\n";
fs << "Dt_Reduced: " << endl << setw(8) << setprecision(5) << Dt_Reduced << "\n" << "\n";
*/

//**************************************************************************************************************************************
//             CONSTRUCTING THE MATRIX WITH THE COORDINATES OF FUNCTIONS (by rows) for each EIGENFACTOR (by columns)


Matrix Function_Coordinates(number_of_functions,V_Reduced.Ncols());
Function_Coordinates=0;
for (int i=0;i<number_of_functions;i++){
    for (int j=0;j<V_Reduced.Ncols();j++){
	for (int n=0;n<number_of_sequences;n++ ){
		Function_Coordinates[i][j]+=(Function_matrix[n][i]/(Function_matrix.Column(i+1).Sum()))*(V_Reduced[n][j]/(sqrt(sum_of_elements_by_row_vector[n]/total_sum)));
	}	
    }
    	//
	//
	//fs << Function_matrix.Column(i+1).Sum()<< "\t";
}
//fs << endl; 

//fs << "Function Coordinates: " << endl << setw(10) << setprecision(5) << Function_Coordinates;
//fs << Function_Coordinates.Nrows() << "\n";
//fs << Function_Coordinates.Ncols() << "\n";

//****************************************************************************************************************************************



//***************************************************************************************************************************************
//              CALCULATING THE DISTANCE BETWEEN EACH AMINOACID_position_VECTOR WITH EACH PROJECTED_FUNCTION_VECTOR

cout << "\t"<< "Calculating chi-squared distances from residues to functions..." << endl;


if (score=="-p") {
	cout << "\t"<< "p-value calculation option selected (it could result in long run time...) " << endl;
	
	for (int i=0; i<Function_Coordinates.Nrows();i++) {
		
		Matrix Random_Function_matrix(number_of_sequences,1000);
		Random_Function_matrix=0;
		int Number_of_ones=(int)Function_matrix.Column(i+1).Sum();
		
		
		int random_integer,max,min;
		if(setsrand()<0) {
			printf("** Couldn't set seed for RND using 'date'.\n");
			exit(1);
		}
		min=0; max=number_of_sequences-1;
	
		int index_array[Number_of_ones];
		for (int d=0; d<1000; d++){
			for(int k=0; k<Number_of_ones; k++){
				index_array[k]=number_of_sequences;
			}
			for(int l=0; l<Number_of_ones; l++){
				int flag=0;
				while (flag<Number_of_ones){
					flag=0;
					random_integer=min+(int)((double)(max-min+1)*rand()/(double)RAND_MAX);
					//fs << "random integer: " << random_integer <<endl;
					for(int k=0;k<Number_of_ones; k++){
						if (random_integer != index_array[k]){
							flag++;
						}
					}
				}
				
				index_array[l]=random_integer;
				//cout << index_array[l] << endl;
				
			}
			//cout << endl;
			//int number_of_ones_assignation=0;
			for(int k=0; k<Number_of_ones; k++){
				int index=index_array[k];
				//fs << "random index array: " << index <<endl;
				Random_Function_matrix[index][d]=1;
				//number_of_ones_assignation++;
				
			}
			
			//fs << endl<< "\t" << number_of_ones_assignation<< endl << endl;
		}
		
	
		//fs << "Vt_Reduced: " << endl << setw(1) << setprecision(0) << Random_Function_matrix << "\n" << "\n";
		
		Matrix Random_Function_Coordinates(1000,V_Reduced.Ncols());
		Random_Function_Coordinates=0;
		for (int k=0;k<1000;k++){
			for (int j=0;j<V_Reduced.Ncols();j++){
				for (int n=0;n<number_of_sequences;n++){
				Random_Function_Coordinates[k][j]+=(Random_Function_matrix[n][k]/(Random_Function_matrix.Column(k+1).Sum()))*(V_Reduced[n][j]/(sqrt(sum_of_elements_by_row_vector[n]/total_sum)));
				}	
			}
		}
		
		//fs << "Random_Function_Coordinates: " << endl << setw(8) << setprecision(5) << Random_Function_Coordinates << "\n" << "\n";
		
	
		distance_aaposition_pair *array;
		array = new distance_aaposition_pair[Column_Coordenates.Nrows()];
		
		double *random_distances_distribution;
		random_distances_distribution = new double[1000*Column_Coordenates.Nrows()];
		int l=-1;
		
		for (int j=0; j<Column_Coordenates.Nrows();j++){
			double square_distance_sum=0;
			for (int k=0; k<Column_Coordenates.Ncols();k++){
				
				double diference= Column_Coordenates[j][k]-Function_Coordinates[i][k];
				//fs << diference << "\t"  << square_distance_sum << endl;
				square_distance_sum+=diference*diference;
				
				//cout << array[j].distance;
			}
			//fs << endl;
			array[j].distance=sqrt(square_distance_sum);
			array[j].aminoacid_position=cleaned_position_aminoacid_vector.at(j);
			
			double *random_distances_distribution_separated;
			random_distances_distribution_separated = new double[1000];
	
			
			for (int m=0;m<1000;m++){
				
				double square_distance_sum=0;
				for (int k=0; k<Column_Coordenates.Ncols();k++){
					
					double diference= Column_Coordenates[j][k]-Random_Function_Coordinates[m][k];
					//fs << diference << "\t"  << square_distance_sum << endl;
					square_distance_sum+=diference*diference;
					
					//cout << array[j].distance;
				}
				//fs << endl;
				random_distances_distribution_separated[m]=sqrt(square_distance_sum);
				l++;
				random_distances_distribution[l]=random_distances_distribution_separated[m];
				//fs << random_distances_distribution[m]<< endl;
			}
			//fs << endl;
			
			combsort(random_distances_distribution_separated,1000);
		
			int f=0;
			while ((array[j].distance>=random_distances_distribution_separated[f]) and (f<1000)) {
			//while (array[j].distance>random_distances_distribution[f])  {
				//fs << " random_distances_distribution[f] "<< "\t" << setw(0) << setprecision(20)<< array[j].distance << "\t" << f << "\t" << setw(10) << setprecision(20) << random_distances_distribution[f] << endl;
				f++;
			}
			//fs << endl;
			
			double pvalue_separated=((double)f/(double)1000);
			array[j].pvalue_separated=pvalue_separated;
		
			delete [] random_distances_distribution_separated;
		}
	
		combsort(random_distances_distribution,1000*Column_Coordenates.Nrows());
		
		for (int j=0;j<Column_Coordenates.Nrows()-1;j++){
			for (int k=j+1;k<Column_Coordenates.Nrows();k++){
				//if (array[j].distance > array[k].distance){
				if (array[j].distance > array[k].distance){	
				
					distance_aaposition_pair temp=array[j];
					array[j]=array[k];
					array[k]=temp;
				}
			}
		}
		
	
		int total=1000*Column_Coordenates.Nrows();
		for (int j=0;j<Column_Coordenates.Nrows();j++){
			
			
			
		
			/*
			int imax_inicial=total-1;
			int imin_inicial=0;
			
			int imax=total-1;
			int imin=0;
			
			while (imax!=imin){
				if (array[j].distance>random_distances_distribution[imax_inicial]) {
					imax=imax_inicial+1;
					imin=imax;
				}else if (array[j].distance<random_distances_distribution[imin_inicial]){
					imax=imin_inicial;
					imin=imax;
				}else if (imax-imin==1){
					imin++;
					imax=imin;
				} else if (array[j].distance>random_distances_distribution[(int)((imax-imin)/2)+imin]) {
					imin=(int)(((imax-imin)/2))+imin;
				} else if (array[j].distance<random_distances_distribution[(int)((imax-imin)/2)+imin]){
					imax=(int)(((imax-imin)/2))+imin;
				} else if (array[j].distance == random_distances_distribution[(int)((imax-imin)/2)+imin]){
					imax=(int)((imax-imin)/2)+imin;
					while (random_distances_distribution[imax]==random_distances_distribution[imax-1]) {
						imax=imax-1;
					}
					imin=imax;
				}
			}
			*/
			
			int f=0;
			while ((array[j].distance>=random_distances_distribution[f]) and (f<1000*Column_Coordenates.Nrows())) {
				f++;
			}
			double pvalue=((double)f/(double)total);
			
			
			//double pvalue=((double)imax/(double)total);
	
			array[j].pvalue=pvalue;
		}
		
	
		
		
		
		
		
		fs << "Best scores for function " << function_names_vector.at(i) << endl << endl;
		fs  << setw(10) << "Residue" << "\t"<< setw(10) << "X2-dist"  << "\t"<< setw(10) << "global-pvalue" << "\t"<< setw(10) << "local-pvalue" << endl;
		for (int j=0;j<cleaned_position_aminoacid_vector.size();j++){
			
			fs  << setw(10) << setprecision(8) << array[j].aminoacid_position<< "\t"<< setw(10) << setprecision(8) << array[j].distance  << "\t"<< setw(10) << setprecision(8) << array[j].pvalue << "\t"<< setw(10) << setprecision(8) << array[j].pvalue_separated << endl;
		}
		fs << "______________________________________________________________________________" << endl << endl << endl;
		
		delete [] random_distances_distribution;
		delete [] array;
	
	}
}  else {
	cout << "\t"<< "p-value calculation option NOT selected" << endl;
	for (int i=0; i<Function_Coordinates.Nrows();i++) {
		
		
			//fs << "Random_Function_Coordinates: " << endl << setw(8) << setprecision(5) << Random_Function_Coordinates << "\n" << "\n";
			
		
		distance_aaposition_pair *array;
		array = new distance_aaposition_pair[Column_Coordenates.Nrows()];
		
		int l=-1;
		
		for (int j=0; j<Column_Coordenates.Nrows();j++){
			double square_distance_sum=0;
			for (int k=0; k<Column_Coordenates.Ncols();k++){
				
				double diference= Column_Coordenates[j][k]-Function_Coordinates[i][k];
				//fs << diference << "\t"  << square_distance_sum << endl;
				square_distance_sum+=diference*diference;
				
				//cout << array[j].distance;
			}
			//fs << endl;
			array[j].distance=sqrt(square_distance_sum);
			array[j].aminoacid_position=cleaned_position_aminoacid_vector.at(j);
			
		
			
		}
	
		
		for (int j=0;j<Column_Coordenates.Nrows()-1;j++){
			for (int k=j+1;k<Column_Coordenates.Nrows();k++){
				//if (array[j].distance > array[k].distance){
				if (array[j].distance > array[k].distance){	
				
					distance_aaposition_pair temp=array[j];
					array[j]=array[k];
					array[k]=temp;
				}
			}
		}
			
		
					
			
			
			
			
		fs << "Best scores for function " << function_names_vector.at(i) << endl << endl;
		fs  << setw(10) << "Residue" << "\t"<< setw(10) << "X2-dist"  << endl;
		for (int j=0;j<cleaned_position_aminoacid_vector.size();j++){
			
			fs  << setw(10) << setprecision(8) << array[j].aminoacid_position<< "\t"<< setw(10) << setprecision(8) << array[j].distance  << endl;
		}
		fs << "______________________________________________________________________________" << endl << endl << endl;
		
		
		
		delete [] array;
	}
	
}

//************************************************************************************************************************

// #define FECHA_FINAL __DATE__ 
// #define HORA_FINAL __TIME__
// 
// cout << "Fin: " << endl;
// cout << FECHA_FINAL << endl;
// cout << HORA_FINAL << endl;

cout << "\t" << "Program finished OK. All outputs written to outfile" << endl;

// Closing output file
fs.close();

return 0;
}








