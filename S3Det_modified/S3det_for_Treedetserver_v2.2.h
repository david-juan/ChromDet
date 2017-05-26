//General
#include "./boost_1_42_0/boost/regex.hpp" //Necessary for the BOOST_REGEX LIBRARY - C Library dealing with Regular Expressions
extern "C"{
#include "./cluster-1.49/src/cluster.h"
//#include "Cluster_library/cluster.h"
//#include "./conf.h"
}

//For eixe
// #include <boost/regex.hpp> //Necessary for the BOOST_REGEX LIBRARY - C Library dealing with Regular Expressions
// extern "C"{
// #include </home/arausell/Newmat10B/cluster.h>
// #include "./conf.h"
// }


// #include <ctype.h> //Necesasry for the toupper() function. Converts char to uppercase in function loading fasta file
// #include <cstring>


struct position_aminoacid_structure{
	string position_aminoacid;
	bool checked;
};

struct distance_aaposition_pair{
	double distance;
	int aminoacid_position_1st;
	int aminoacid_position_2nd;
	bool interaction;
	double pvalue;
};

struct position_aa_distance{
	double distance;
	int rank;
	int position_1st;
	string aminoacid_1st;
};

struct function_information_structure{
	string name;
	vector<int> clusters_belonging_vector;
	int size;
	vector<position_aa_distance> elements;
	position_aa_distance *elements_array; //same as previous but it can be sorted by combsort (combsort doesn't work with vectors[...])
	int *elements_array_selected;
	int *elements_array_selected_rank;
	string *elements_array_selected_amin;
	double *elements_array_selected_dist;
	int number_of_selected_residues;
};

struct global_position_structure{
	int position;
	double rank;
	double number_of_groups_based;
};

//////////////////////////////////////////////
// Fast sorting algorithm
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

static void combsort_int(int a[], int aSize) {
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
/////////////////////////////////////////////

static void combsort_distance_aaposition_pair(distance_aaposition_pair a[], int aSize) {
  int gap = aSize;
  for (;;) {
    gap = newGap(gap);
    bool swapped = false;
    for (int i = 0; i < aSize - gap; i++) {
      int j = i + gap;
      if (a[i].distance > a[j].distance) {
        std::swap(a[i], a[j]);
        swapped = true;
      }
    }
    if (gap == 1 && !swapped)
      break;
  }
}

static void combsort_position_aa_distance(position_aa_distance a[], int aSize) {
  int gap = aSize;
  for (;;) {
    gap = newGap(gap);
    bool swapped = false;
    for (int i = 0; i < aSize - gap; i++) {
      int j = i + gap;
      if (a[i].distance > a[j].distance) {
        std::swap(a[i], a[j]);
        swapped = true;
      }
    }
    if (gap == 1 && !swapped)
      break;
  }
}

static void combsort_global_position_structure(global_position_structure a[], int aSize) {
  int gap = aSize;
  for (;;) {
    gap = newGap(gap);
    bool swapped = false;
    for (int i = 0; i < aSize - gap; i++) {
      int j = i + gap;
      if (a[i].rank > a[j].rank) {
        std::swap(a[i], a[j]);
        swapped = true;
      }
    }
    if (gap == 1 && !swapped)
      break;
  }
}

//////////////////////////////////////////////
int factorial(int n)
 {
    int fn=1;

    if(n>1)
       fn=(n*factorial(n-1));

    return fn;
 }
//////////////////////////////////////////////

bool fasta_protein_name (const std::string&);
bool fasta_protein_sequence (const std::string&);

//Boolean funtion which returns TRUE if the string s starts with a > (it will be used to indicate if a line in a FASTA input file contains the name of a protein 
bool fasta_protein_name (const std::string& s)
{
	static const boost::regex e("^>(.*)");
	return boost::regex_match(s, e);
}

bool fasta_protein_sequence (const std::string& s)
{
	//static const boost::regex e("^(\w|.|-)(.*)");
	static const boost::regex e0("(.*)([ABCDEFGHIJKLMNOPQRSTUVWYZXabcdefghijklmnopqrstuvwyzx\\*\\-\\.])(.*)");
	static const boost::regex e1("(.*)([^ABCDEFGHIJKLMNOPQRSTUVWYZXabcdefghijklmnopqrstuvwyzx\\*\\-\\.])(.*)");
	bool re0= boost::regex_match(s, e0);
	bool re1= boost::regex_match(s, e1);
	if(re0 and not re1){return 1;}else{return 0;}
}

//This function converts a given string in a vector called "tokens" splitting it by a given delimiter
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


// Funtion that reads a fasta file with amino acid letters in uppercase and gaps as "-"
// returns a letter_matrix as a bidimensional char array with the "amino acid letters" of the alignment
// Dependencies: fasta_protein_name; fasta_protein_sequence; #include <boost/regex.hpp>; #include <ctype.h> //Necesasry for the toupper() function. Converts char to uppercase in function loading fasta file

void loading_conf_file(
	char * conf_file,
        string & tmp_directory,
        string & exec_directory,
        string & order,
	int & error_loading_conf_file,
	int & error_couldnt_open
)
	
{
	string line;//it will be used in the input file reading process
	ifstream infile (conf_file);
	if(infile.fail()){error_couldnt_open=1;return;}else{error_couldnt_open=0;}
	static const boost::regex tmp_regex("(.*)(string tmp_directory=)(.*)");
	static const boost::regex comment_regex("^//(.*)");
	static const boost::regex exec_regex("(.*)(string exec_directory=)(.*)");
	static const boost::regex order_regex("(.*)(string order=)(.*)");
	//Reading line by line the input file (NOTICE: Getline eliminates the 'endline' character) 
	while (getline(infile,line,'\n')) {
	  if      ( boost::regex_match(line,tmp_regex)   and not boost::regex_match(line,comment_regex) ) {
	    vector<string> tokens;
	    Tokenize(line,tokens,"\"");
	    if(tokens.size()!=3){
	      error_loading_conf_file=1;
	      return;
	    }
	    tmp_directory=tokens.at(1);
	    //cout << "Here is tmp: " << tmp_directory << endl;
	  }
	  else if ( boost::regex_match(line,exec_regex)  and not boost::regex_match(line,comment_regex) ) {
	    vector<string> tokens;
	    Tokenize(line,tokens,"\"");
	    if(tokens.size()!=3){
	      error_loading_conf_file=1;
	      return;
	    }
	    exec_directory=tokens.at(1);
	    //cout << "Here is exec: " << exec_directory << endl;
	  }
	  else if ( boost::regex_match(line,order_regex) and not boost::regex_match(line,comment_regex) ) {
	    vector<string> tokens;
	    Tokenize(line,tokens,"\"");
	    if(tokens.size()!=3){
	      error_loading_conf_file=1;
	      return;
	    }
	    order=tokens.at(1);
	    //cout << "Here is order: " << order << endl;
	  }
	  else{
	    continue;
	  }
	}
	//Closing input file
	infile.close();
}

void loading_fasta_file(
	char * fasta_file,
	char** & letter_matrix,
	vector<string> & protein_name_vector, //it will contain the names of the proteins taking part in the alignment
	vector<string> & position_vector,//vector with aminoacid positions which have at least one non-null element
	int & number_of_sequences, //number of sequences of the alignment
	int & number_of_positions, //number of positions of the alignment
	int & error_loading_fasta,
	int & error_size_of_sequences,
	int & error_couldnt_open
)
	
{
	vector<string> protein_sequence_vector; //it will contain the sequences of the proteins taking part in the alignment
	int number_of_positions_tmp=0;
	string line;//it will be used in the input file reading process
	string protein_sequence_part="";//it will be used while constructing the protein_sequence_vector
	ifstream infile (fasta_file);
	if(infile.fail()){error_couldnt_open=1;return;}else{error_couldnt_open=0;}
	static const boost::regex e3("^$");
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
						return;
					}
				}
			}
			protein_sequence_part="";//reinitializing the construction of the protein sequence
			line.erase(0,1);
			protein_name_vector.push_back(line);
		}	
		else if ( fasta_protein_sequence(line)  ){
			protein_sequence_part=protein_sequence_part+line;
		} else if(boost::regex_match(line, e3)){
			continue;
		}else {
			error_loading_fasta=1;
			return;
		}
	}
	error_loading_fasta=0;
	// just to load into the protein_sequence_vector the last sequence
	if ( fasta_protein_sequence(protein_sequence_part) ) {
		protein_sequence_vector.push_back(protein_sequence_part);
		if(protein_sequence_vector.size()==1){number_of_positions_tmp=protein_sequence_part.length();}
		if(protein_sequence_vector.size()>1){
			if(protein_sequence_part.length()!=number_of_positions_tmp){
				error_size_of_sequences=1;
				return;
			}
		}
	}
	//Closing input file
	infile.close();
	
	number_of_sequences = protein_name_vector.size();
	number_of_positions = number_of_positions_tmp;
		
	for (int p=0; p<(number_of_positions);p++){
		char position_char[200];
		sprintf (position_char,"%d",p+1);//converts the "int type" position into "char type"
		string position_string= position_char;//converts the "char type" position into "string type"
		position_vector.push_back(position_string);
	}
	
	//Constructing a bidimensional char array called "letter_matrix" with the "amino acid letters" of the alignment 
	//letter_matrix: a bidimensional char array with the "amino acid letters" of the alignment
	letter_matrix=new char *[number_of_sequences];
	for(int i=0;i<number_of_sequences;i++){letter_matrix[i]=new char[number_of_positions];}
	
	for (int i=0; i<(number_of_sequences);i++ ){
		char * puntero;
		puntero = new char [protein_sequence_vector[i].size()+1];
		strcpy (puntero,protein_sequence_vector[i].c_str());
 		//const char* puntero= protein_sequence_vector[i].c_str();//it converts a string in a char array
		for (int j=0; j<(number_of_positions);j++){
//			puntero[j]=(char)toupper(puntero[j]);
//			if ((puntero[j]=='B')|(puntero[j]=='U')|(puntero[j]=='Z')|(puntero[j]=='X')|(puntero[j]=='*')|(puntero[j]=='.')|(puntero[j]=='-')){
//			if ((puntero[j]=='Z')|(puntero[j]=='X')|(puntero[j]=='*')|(puntero[j]=='.')|(puntero[j]=='-')){
                      if ((puntero[j]=='*')|(puntero[j]=='.')|(puntero[j]=='-')){
				letter_matrix[i][j]='#';
			} else {
				letter_matrix[i][j]=puntero[j];
			}
		}
 		delete [] puntero;
	}
	/*
	'B','U','Z','X' and '*' are treated as gap positions
	
	    *A  alanine                        *P  proline
	     B  aspartate or asparagine        *Q  glutamine
	    *C  cystine                        *R  arginine
	    *D  aspartate                      *S  serine
	    *E  glutamate                      *T  threonine
	    *F  phenylalanine                   U  selenocysteine
	    *G  glycine                        *V  valine
	    *H  histidine                      *W  tryptophan
	    *I  isoleucine                     *Y  tyrosine
	    *K  lysine                          Z  glutamate or glutamine
	    *L  leucine                         X  any
	    *M  methionine                      *  translation stop
	    *N  asparagine                      -  gap of indeterminate length
	*/
}

// Funtion that reads a matrix with amino acid letters (in uppercase and gaps as "#") and remove columns with more than a given percentage of gaps. It keeps the track of the original position numbers
void removing_gappy_columns(
	char** & letter_matrix,
	double percentage_of_gaps,
	vector<string> & position_vector,
	int & number_of_sequences, //number of sequences of the alignment
	int & number_of_positions, //number of positions of the alignment
	int & number_of_gappy_columns
)
	
{
	char** letter_matrix_without_gaps;
	vector<string> position_vector_without_gaps;
	number_of_gappy_columns=0;
	for (int j=0; j<(number_of_positions);j++){
		int number_of_gaps=0;
		for (int i=0; i<(number_of_sequences);i++ ){
			if (letter_matrix[i][j]=='#'){
				number_of_gaps++;
			}
		}
		double percentage_of_gappy_positions=(double (number_of_gaps)) /(double (number_of_sequences));
		if ( percentage_of_gappy_positions>percentage_of_gaps){
			number_of_gappy_columns++;
		}
	}
	
	letter_matrix_without_gaps=new char *[number_of_sequences];
	for(int i=0;i<number_of_sequences;i++){letter_matrix_without_gaps[i]=new char[number_of_positions-number_of_gappy_columns];}
	int position_vector_without_gaps_MARKER=-1;
	for (int j=0; j<(number_of_positions);j++){
		int number_of_gaps=0;
		for (int i=0; i<(number_of_sequences);i++ ){
			if (letter_matrix[i][j]=='#'){
				number_of_gaps++;
			}
		}
		if ((double (number_of_gaps)/double (number_of_sequences))<=percentage_of_gaps){
			position_vector_without_gaps_MARKER++;
			for (int i=0; i<(number_of_sequences);i++ ){
				letter_matrix_without_gaps[i][position_vector_without_gaps_MARKER]=letter_matrix[i][j];
			}
			position_vector_without_gaps.push_back(position_vector.at(j));
		}
	}
	for(int i=0;i<number_of_sequences;i++){delete letter_matrix[i];}
	delete [] letter_matrix;
	//NOTICE:
	letter_matrix=letter_matrix_without_gaps;
	position_vector=position_vector_without_gaps;
	number_of_positions=(number_of_positions-number_of_gappy_columns);
}

// Funtion that reads a matrix with amino acid letters (in uppercase and gaps as "#") and remove sequences with less than a given percentage of similarity to any of the rest of sequences
void removing_outliers(
	char** & letter_matrix,
	double percentage_of_similarity,
	int & number_of_sequences, //number of sequences of the alignment
	int & number_of_positions, //number of positions of the alignment
	int & number_of_outliers,
	vector<string> & protein_name_vector,
	vector<string> & outliers_name_vector
)	
{	
	double** similarity_matrix= new double *[number_of_sequences];
	for (int i = 0; i < number_of_sequences; i++){
		similarity_matrix[i] = new double [number_of_sequences];
	}
	for (int i=0; i<(number_of_sequences-1);i++ ){
		similarity_matrix[i][i]=1;
		for (int j=i+1; j<number_of_sequences;j++){
			int similarity=0;
			int number_of_non_empty_pairings=0;
			for(int k=0; k<(number_of_positions);k++ ){
				if((letter_matrix[i][k]!='#')and(letter_matrix[j][k]!='#')){
					number_of_non_empty_pairings++;
					if(letter_matrix[i][k]==letter_matrix[j][k]){
						similarity++;
					}
				}
			}
			double similarity_coef=double(similarity)/double(number_of_non_empty_pairings);
			similarity_matrix[i][j]=similarity_coef;
			similarity_matrix[j][i]=similarity_coef;
		}
	}
	similarity_matrix[number_of_sequences-1][number_of_sequences-1]=1;
	

	vector<int> outliers;
	number_of_outliers=0;
	for (int i=0; i<number_of_sequences;i++ ){
		int number_of_distant_seqs=0;
		for (int j=0; j<number_of_sequences;j++){
			if (similarity_matrix[i][j]<percentage_of_similarity){
				number_of_distant_seqs++;
			}
		}
		if(number_of_distant_seqs==(number_of_sequences-1)){
			number_of_outliers++;
			outliers.push_back(i);
		}
	}
	
	for (int i = 0; i <number_of_sequences ; i++) { delete [] similarity_matrix[i];}
	delete [] similarity_matrix;
	
	
	vector<string> protein_name_vector_without_outliers;
	
	char** letter_matrix_without_outliers=new char *[number_of_sequences-number_of_outliers];
	for(int i=0;i<number_of_sequences-number_of_outliers;i++){letter_matrix_without_outliers[i]=new char[number_of_positions];}
	
	int NOT_outlier_MARKER=-1;
	for (int i=0; i<number_of_sequences;i++ ){
		int flag=0;
		for(int k=0; k<outliers.size();k++){
			if(i==outliers.at(k)){
				flag=1;
			}
		}
		if (flag==0){
			NOT_outlier_MARKER++;
			protein_name_vector_without_outliers.push_back(protein_name_vector.at(i));
			for(int j=0; j<(number_of_positions);j++ ){
				letter_matrix_without_outliers[NOT_outlier_MARKER][j]=letter_matrix[i][j];
			}
		} else {
			outliers_name_vector.push_back(protein_name_vector.at(i));
		}
	}
	
	for(int i=0;i<number_of_sequences;i++){delete letter_matrix[i];}
	delete [] letter_matrix;
	
	//NOTICE
	protein_name_vector=protein_name_vector_without_outliers;
	number_of_sequences=number_of_sequences-number_of_outliers;
	letter_matrix=letter_matrix_without_outliers;
}	
	
//Dependencies: 	-L /home/arausell/Newmat10B/ -lnewmat -lm -I /home/arausell/Newmat10B/
void disjunctive_coding(
	Matrix & A,
	char** letter_matrix,
	vector<string> protein_name_vector,
	vector<string> position_vector,
	int number_of_sequences, //number of sequences of the alignment
	int number_of_positions, //number of positions of the alignment
	vector<string> & position_aminoacid_vector,//vector with aminoacid positions which has at least one non-null element
	vector<double> & sum_of_elements_by_column_vector,
	vector<double> & sum_of_elements_by_row_vector,
	double & total_sum,
	vector<int> & Conserved_positions,
	int conservation_threshold
)

{
	//Constructing the "position_aminoacid_vector" with elements the 21xL position-aminoacid flags, where L is the length of the alignment
//	char aminoacid[]={'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','#'};
//	char aminoacid[]={'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','#','B','J','O','U'};
//	char aminoacid[]={'A','R','D','C','E','H','I','L','F','T','V'};
        char aminoacid[]={'K','V','g','i','j','k','l','m','n','A','B','C','D','E','F','G','H','I','J','L','M','N','O','P','Q','R','S','T','U','W','X','Y','Z','a','b','c','d','e','f','h'};
//      char aminoacid[]={'A','B','D','C','E','F'};

      for (int p=0; p<(number_of_positions);p++){
//		for (int i=0;i<21;i++ ){
//		for (int i=0;i<25;i++ ){
//		for (int i=0;i<11;i++ ){	      
		for (int i=0;i<40;i++ ){
		      string position_aminoacid_string=position_vector.at(p)+"/"+aminoacid[i];//fusion "position" + "aminoacid" strings
			position_aminoacid_vector.push_back(position_aminoacid_string);
		}
	}

	//Constructing the disjunctive matrix A
	Matrix A_aux(number_of_sequences,1);
	A_aux=0;
	total_sum=0;	
	int position_aminoacid_vector_MARKER=-1;
	int columns_asigned =0;
	vector<double> sum_of_elements_by_row_vector_aux(number_of_sequences);
	vector<string> cleaned_position_aminoacid_vector;
	int limit_conserved= (number_of_sequences*conservation_threshold)/100;if(((number_of_sequences*conservation_threshold)%100)>0){limit_conserved++;}
	for (int j=0; j<(number_of_positions);j++){
//		for (int k=0; k<21; k++){
//		for (int k=0; k<25; k++){
//		for (int k=0; k<11; k++){			
                for (int k=0; k<40; k++){
			position_aminoacid_vector_MARKER += 1;
			Matrix C(number_of_sequences,1);//new column to add to A if has at least one non-null element
			C=0;
			int number_of_ones=0;
			
			for (int i=0;i<(number_of_sequences);i++ ){	
			
				if (letter_matrix[i][j]== aminoacid[k]){ //codifying for each aminoacid type k
					C[i][0]=1;
					sum_of_elements_by_row_vector_aux[i]+=1;
					number_of_ones +=1;
					
				} 
			}
			if (number_of_ones!=0){
				columns_asigned +=1;
				if (columns_asigned > 1){
					 //because the first time you need to overwrite the first default null column (honestly this statement is a shoddy piece of work)
					A_aux=A_aux|C;
				} else {	
					A_aux << C;
				}
				//constructing the vector with aminoacid positions which has at least one non-null element
				cleaned_position_aminoacid_vector.push_back(position_aminoacid_vector.at(position_aminoacid_vector_MARKER));
				sum_of_elements_by_column_vector.push_back(number_of_ones);
				total_sum+=number_of_ones;
				if(number_of_ones>limit_conserved){
					vector<string> tmp_aminoacid_position;
					Tokenize(position_aminoacid_vector.at(position_aminoacid_vector_MARKER),tmp_aminoacid_position,"/");
					Conserved_positions.push_back(atoi(tmp_aminoacid_position[0].c_str()));
				}
			}
		}
	}
	
	
	for(int i=0;i<number_of_sequences;i++){delete letter_matrix[i];}
	delete [] letter_matrix;
	
// 	NOTICE:
	position_aminoacid_vector=cleaned_position_aminoacid_vector;
	A=A_aux;
	ofstream fs2("disjunctive.txt"); //Create output file with the chosen name
	if (fs2.fail()){
        string out_file_string="disjunctive.txt";
        cout << "\t" << "\tError opening file " <<  out_file_string << endl << "\t\tFile not found or couldn't be opened" << endl;
        cout << "\t" << "\tExecution aborted" << endl;
        fs2.close();
}
	for (int j=0;j<A.Ncols(); j++){

		fs2 << "\t" << cleaned_position_aminoacid_vector[j];
	}
	fs2 << endl;
	for (int i=0; i<A.Nrows();i++){
  		fs2 << protein_name_vector[i]; 
                for (int j=0;j<A.Ncols(); j++){
                        fs2 << "\t"<< A[i][j];
                }
                fs2 << endl;
        }
     //   fs2 << endl;
	fs2.close();
	A_aux.Release();
	sum_of_elements_by_row_vector= sum_of_elements_by_row_vector_aux;
}

//Dependencies: 	-L /home/arausell/Newmat10B/ -lnewmat -lm -I /home/arausell/Newmat10B/
void eigen_descomposition(//NOTICE: Newmat descomposition leads to matrices where the first eigen* correspond to the last column
	Matrix & Z,
	DiagonalMatrix & D_Reduced,
	DiagonalMatrix & D_adjusted,//Only used in the MCA-Greenacre analysis (Analysis 2)
	Matrix & V_Reduced,
	DiagonalMatrix & Dc,
	Matrix & Column_Coordinates,//Matrix wich each file i has the coordinates of column i in factor j (columns arranging the coordinates for each factorj)
 	//Matrix & Row_Coordinates_standard,//Matrix wich each file i has the coordinates of row i in factor j (columns arranging the coordinates for each factorj)
	Matrix & Row_Coordinates_principal_indirect,
	Matrix & A,//disjunctive matrix
	int number_of_rows,
	int number_of_positions,
	int number_of_columns,
	vector<double> sum_of_elements_by_column_vector,
	vector<double> sum_of_elements_by_row_vector,
	double total_sum,
	int Analysis //Analysis type: 1.MCA 2.MCA-Greenacre 3.Non-centered-PCA (Notice: non identical to FASS) 4.PCA 
)
{
	for (int i=0;i<(number_of_rows);i++ ){
		for (int j=0;j<number_of_columns;j++){
			if(Analysis==1||Analysis==2){Z[i][j]=A[i][j]/(sqrt(sum_of_elements_by_row_vector[i])*sqrt(sum_of_elements_by_column_vector[j]));} //MCA
			if(Analysis==3){Z[i][j]=A[i][j];} // Centering of variables was also not done by FASS
			if(Analysis==4){Z[i][j]=A[i][j]-(sum_of_elements_by_column_vector[j]/number_of_rows);} //Centered PCA
		}
	}
	
	for (int i=1;i<(number_of_columns+1);i++ ){
			Dc(i,i)=1/(sqrt(sum_of_elements_by_column_vector[i-1]/total_sum));
	}
// 	A.Release();//The memory used by A is returned to the system, to improve efficiency
	
	if (number_of_rows <= number_of_columns){
		SymmetricMatrix X;
		DiagonalMatrix D;//Diagonal matrix containing the eigenvalues of X
		Matrix V;//Diagonal matrix containing the eigenvectors of X by columns
		
		X << Z*Z.t();
		if(Analysis==3||Analysis==4){ X << X/number_of_columns;}//NOTICE
		//NOTICE:Here FASS did: X << X/number_of_positions ; doing this is not possible anymore to perform an indirect assesment of V and D
		EigenValues(X,D,V);
		X.Release();
		
		int number_of_non_null_eigen_values=0;
		if(Analysis==1||Analysis==3||Analysis==4){
			for (int i=D.Nrows();((i>0) and (D(i,i)>0.0000000001));i--){
				number_of_non_null_eigen_values+=1;
			}
		}
		if(Analysis==2){
			for (int i=D.Nrows();((i>0) and (sqrt(D(i,i))>=(1/(double)number_of_positions)));i--){
				number_of_non_null_eigen_values+=1;
			}
		}
		// Just limitting the size of V and D 
		if(Analysis==1||Analysis==2){// In MCA the first eigenvalue (and thereby the first eigenvector) is removed because is the trivial "1"
			V_Reduced << V.SubMatrix(1,V.Nrows(),(V.Ncols()+1)-number_of_non_null_eigen_values,V.Ncols()-1); 
			D_Reduced << D.SubMatrix((D.Nrows()+1)-number_of_non_null_eigen_values,D.Nrows()-1,(D.Ncols()+1)-number_of_non_null_eigen_values,D.Ncols()-1);
		}
		if(Analysis==3||Analysis==4){
			V_Reduced << V.SubMatrix(1,V.Nrows(),(V.Ncols()+1)-number_of_non_null_eigen_values,V.Ncols());
			D_Reduced << D.SubMatrix((D.Nrows()+1)-number_of_non_null_eigen_values,D.Nrows(),(D.Ncols()+1)-number_of_non_null_eigen_values,D.Ncols());
		}
		V.Release();
		D.Release();
		
	} else {//Indirect assesment of V_Reduced and D_Reduced
		SymmetricMatrix Xt;
		DiagonalMatrix Dt;
		Matrix U,U_Reduced;
		
		Xt << Z.t()*Z;
		if(Analysis==3||Analysis==4){ Xt << Xt/number_of_rows;}//NOTICE
		EigenValues(Xt,Dt,U);
		Xt.Release();
		
		int number_of_non_null_eigen_values_transp=0;
		if(Analysis==1||Analysis==3||Analysis==4){
			for (int i=Dt.Ncols();((i>0) and (Dt(i,i)>0.0000000001));i--){
				number_of_non_null_eigen_values_transp+=1;
			}
		}
		if(Analysis==2){
			for (int i=Dt.Ncols();((i>0) and (sqrt(Dt(i,i))>=(1/(double)number_of_positions)));i--){
				number_of_non_null_eigen_values_transp+=1;
			}
		}
		
		if(Analysis==1||Analysis==2){
			U_Reduced << U.SubMatrix(1,U.Nrows(),(U.Ncols()+1)-number_of_non_null_eigen_values_transp,U.Ncols()-1);
			D_Reduced << Dt.SubMatrix((number_of_columns+1)-number_of_non_null_eigen_values_transp,number_of_columns-1,(number_of_columns+1)-number_of_non_null_eigen_values_transp,number_of_columns-1);
			U.Release();
			Dt.Release();
				
			DiagonalMatrix Inv_Sqrt_D_Reduced(D_Reduced.Ncols());
			Inv_Sqrt_D_Reduced=0;
			for (int i=1; i<=D_Reduced.Ncols(); i++){
				Inv_Sqrt_D_Reduced(i,i)=(1/sqrt(D_Reduced(i,i)));
			}
			V_Reduced=(Z*U_Reduced)*Inv_Sqrt_D_Reduced; //Indirect assesment of V_Reduced:
			
			//DiagonalMatrix Df(number_of_rows); //Constructing DiagonalMatrix Df
			//for (int i=1;i<(number_of_rows+1);i++ ){
			//		Df(i,i)=1/(sqrt(sum_of_elements_by_row_vector[i-1]/total_sum));
			//}
			//Row_Coordinates_principal_direct=Df*Z*U_Reduced;
			U_Reduced.Release();
		}
		if(Analysis==3||Analysis==4){
			DiagonalMatrix Dt_Reduced;
			U_Reduced << U.SubMatrix(1,U.Nrows(),(U.Ncols()+1)-number_of_non_null_eigen_values_transp,U.Ncols());
			Dt_Reduced << Dt.SubMatrix((number_of_columns+1)-number_of_non_null_eigen_values_transp,number_of_columns  ,(number_of_columns+1)-number_of_non_null_eigen_values_transp,number_of_columns);	
			U.Release();
			Dt.Release();
			
			DiagonalMatrix Inv_Sqrt_D_Reduced(Dt_Reduced.Ncols());
			Inv_Sqrt_D_Reduced=0;
			for (int i=1; i<=Dt_Reduced.Ncols(); i++){
				Inv_Sqrt_D_Reduced(i,i)=(1/sqrt(Dt_Reduced(i,i)));
			}
			
			double product=( (double) number_of_rows / (double) number_of_columns);
			D_Reduced << product*Dt_Reduced;
			V_Reduced=((Z*U_Reduced)*((1/sqrt((double)number_of_rows))*(Inv_Sqrt_D_Reduced)) ); //Indirect assesment of V_Reduced:
			Dt_Reduced.Release();
			//Row_Coordinates_principal_direct=Z*U_Reduced;
			U_Reduced.Release();
		}
	}
	DiagonalMatrix D_Reduced_sqrt(D_Reduced.Ncols());D_Reduced_sqrt=0;
	DiagonalMatrix Inv_Sqrt_D_Reduced(D_Reduced.Ncols());Inv_Sqrt_D_Reduced=0;
	for (int i=D_Reduced.Ncols();i>0; i--){
		D_Reduced_sqrt(i,i)=sqrt(D_Reduced(i,i));
		Inv_Sqrt_D_Reduced(i,i)=(1/sqrt(D_Reduced(i,i)));
	}
		
	if(Analysis==1){	//MCA
		Column_Coordinates=Dc*(Z.t()*V_Reduced);
 		//Row_Coordinates_standard=sqrt(number_of_rows)*V_Reduced;
		Row_Coordinates_principal_indirect=sqrt(number_of_rows)*(V_Reduced*D_Reduced_sqrt);
	}
	if(Analysis==2){	//MCA-Greenacre
		DiagonalMatrix D_adjusted_tmp(D_Reduced.Ncols());D_adjusted_tmp=0;
		D_adjusted=D_adjusted_tmp;
		DiagonalMatrix D_adjusted_sqrt(D_Reduced.Ncols());
		
		double total_inertia_burt=0;
		for (int i=1; i<=D_Reduced.Ncols(); i++){
			D_adjusted(i,i)=pow(((double) number_of_positions/(double)(number_of_positions-1)),2)*pow((D_Reduced(i,i)-(1/number_of_positions)),2);
			//total_inertia_burt+=pow(D_Reduced(i,i),2);
		}
		//double new_total_inertia=(number_of_positions/(number_of_positions-1))*(total_inertia_burt-((number_of_columns-number_of_positions)/(pow(number_of_positions,2))));
		//double new_total_inertia=(number_of_positions/(number_of_positions-1))*(total_inertia_burt-(((number_of_positions*21)-number_of_positions)/(pow(number_of_positions,2))));
		//D_adjusted_percentage=(1/new_total_inertia)*D_adjusted; NOTICE: This way is not used because the total sum could be greater than 100%
		//D_adjusted_percentage=(1/sum_eigenvalues_adjusted)*D_adjusted; NOTICE: This way is used in main
		
		for (int i=D_adjusted.Ncols();i>0; i--){
			D_adjusted_sqrt(i,i)=sqrt(D_adjusted(i,i));
		}
		Column_Coordinates=((Dc*(Z.t()*V_Reduced))*Inv_Sqrt_D_Reduced)*D_adjusted_sqrt;//Matrix wich each file i has the coordinates of position_aminoacid i in factor j (columns arranging the coordinates for each factorj)
		Row_Coordinates_principal_indirect=sqrt(number_of_rows)*(V_Reduced*D_adjusted_sqrt);
		//Row_Coordinates_standard=((sqrt(number_of_rows)*V_Reduced)*D_Reduced_sqrt.i())*D_adjusted_sqrt;
	}	
	if(Analysis==3){	// Non-Centered PCA
		Column_Coordinates=(Z.t()*V_Reduced);//Matrix wich each file i has the coordinates of position_aminoacid i in factor j (columns arranging the coordinates for each factorj)
		// Both next expressions are equivalent
		//Row_Coordinates_principal_indirect=(Z*Z.t())*V_Reduced*((1/sqrt((double)number_of_columns))*(Inv_Sqrt_D_Reduced));
		Row_Coordinates_principal_indirect=(Z*Column_Coordinates)*((1/sqrt((double)number_of_columns))*(Inv_Sqrt_D_Reduced));// THIS IS THE BEST WAY
		//Row_Coordinates_standard=V_Reduced*D_Reduced_sqrt; NOTICE: This way is how FASS assessed the coordinates of sequences
	}
	if(Analysis==4){	// PCA
		Row_Coordinates_principal_indirect=(Z*Z.t())*V_Reduced*((1/sqrt((double)number_of_columns))*(Inv_Sqrt_D_Reduced));
		Column_Coordinates=0;
		//cout << "NOTICE: PCA is not symmetrical. Only rows can be projected\n";
	}
}

//READING FUNCTION MATRIX
//Open the input file indicated as the second argument in the prompt which will contain the binary classification of the proteins of the alignment (by rows) corresponding a diferent functions (arranged by columns)
void loading_supervised_function_coding(
	char * function_file,
	vector<string> protein_name_vector,
	Matrix & Function_matrix,
	vector<string> & function_names_vector,
	int & number_of_functions,
	double & flag_for_interaction,
	int & error_couldnt_open,
	int & error_file_is_not_in_proper_format,
	int & error_sequence_not_found_in_function_file,
	int & flag_sequence_with_not_assigned_function,
	int & flag_sequence_with_sum_of_assigned_functions_different_from_1,
	string & protein_name_causing_error,
	int & flag_FUZZYClassification
)
{
	error_file_is_not_in_proper_format=0;
	error_sequence_not_found_in_function_file=0;
	flag_sequence_with_not_assigned_function=0;
	flag_sequence_with_sum_of_assigned_functions_different_from_1=0;
	protein_name_causing_error="";
	flag_FUZZYClassification=0;

	string line;
	ifstream infile2 (function_file);
	if(infile2.fail()){error_couldnt_open=1;return;}else{error_couldnt_open=0;}
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
		ifstream infile2 (function_file);
		int flag=0;
		getline(infile2,line,'\n');
		while ( (getline(infile2,line,'\n')) && (flag==0) ){
			if (boost::regex_match(line,e2)){
				vector<string> tokens;
				Tokenize(line,tokens,"\t");
				if(tokens.size()!=number_of_functions+1){
					error_file_is_not_in_proper_format=1;
					return;
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
								if((temp_double>0)and(temp_double<1)){
								  flag_FUZZYClassification=1;
								}
							}else{
								error_file_is_not_in_proper_format=1;
								return;
							}
						}else{
							error_file_is_not_in_proper_format=1;
							return;
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
			error_sequence_not_found_in_function_file=1;
			return;
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
}

//**************************************************************************************************************************************

double* Wilcoxon_pvalues;
void Wilcoxon_test(
	string order,
	Matrix Row_Coordinates_principal_indirect,
	string tmp_directory,
	string exec_directory,
	double* & Wilcoxon_pvalues,
	int Maximum_number_of_axes,
	int & error_couldnt_open
)
{
	Wilcoxon_pvalues=new double [Maximum_number_of_axes-1];
	string tmp_infile_string=tmp_directory;tmp_infile_string+="tmp_Wilcoxon_test_input.txt";
	const char* tmp_infile_const_char=tmp_infile_string.c_str();
	ofstream tmp1(tmp_infile_const_char);
	
	double squared_dist[Row_Coordinates_principal_indirect.Nrows()][Maximum_number_of_axes];
	int nc=Row_Coordinates_principal_indirect.Ncols()-1;
	for (int i=0;i<Row_Coordinates_principal_indirect.Nrows();i++){
		squared_dist[i][0]=pow(Row_Coordinates_principal_indirect[i][nc],2);
		tmp1 << setprecision(10) << squared_dist[i][0] ;
		for(int j=1;j<Maximum_number_of_axes;j++){
			squared_dist[i][j]=squared_dist[i][j-1]+pow(Row_Coordinates_principal_indirect[i][nc-j],2);
			tmp1 << "\t" << setprecision(10) << squared_dist[i][j];
		}
		tmp1 << endl;
	}
	tmp1.close();
	
	//i.e.:  system("/usr/local/lib64/R/bin/R --no-save </home/arausell/MCdet/Professional/C_trials_Wilcoxon_test.R");
	order+=" --no-save <";order+=exec_directory;order+="S3det_Wilcoxon_test.R >";
	order+=tmp_directory;order+="R_stdout ";
	order+=tmp_directory;
	const char* order_to_system=order.c_str();
	system(order_to_system);

	string tmp_outfile_string=tmp_directory;tmp_outfile_string+="tmp_Wilcoxon_test_output.txt";
	const char* tmp_outfile_const_char=tmp_outfile_string.c_str();
	
	ifstream tmp2(tmp_outfile_const_char);
	if(tmp2.fail()){error_couldnt_open=1;return;}else{error_couldnt_open=0;}
	string line;
	getline(tmp2,line,'\n');
	int index=-1;
	while (getline(tmp2,line,'\n')){	
		index++;
		vector<string> tokens;
		Tokenize(line,tokens,"\t");
		const char* caracter=tokens[1].c_str();
		double numero=atof(caracter);
		Wilcoxon_pvalues[index] = numero;
	}
	tmp2.close();
	
	order="rm ";order+=tmp_directory;order+="tmp_Wilcoxon_test_input.txt";
	order_to_system=order.c_str();
	system(order_to_system);
	
	order="rm ";order+=tmp_directory;order+="tmp_Wilcoxon_test_output.txt";
	order_to_system=order.c_str();
	system(order_to_system);

	order="rm ";order+=tmp_directory;order+="R_stdout";
	order_to_system=order.c_str();
	system(order_to_system);

}
//**************************************************************************************************************************************
void kmeans(
	int nrows,
	int ncols,
	double** data,
	int** mask,
	const int nclusters,
	int Analysis,
	string Index_type,
	double & Index_value,
	int* & clusterid,
	int npass,
	int & ifound,
	double & F_Fisher
)
{
	int i,j,k;
	const int transpose = 0;
	const char dist = 'e';
	const char method = 'a';
	ifound = 0;
	double error;
	int** index;
	int* count;
	double* weight = new double [ncols];
	double** cdata = new double *[nclusters];
	int** cmask = new int *[nclusters];
	for (i = 0; i < nclusters; i++){
		cdata[i] = new double [ncols];
		cmask[i] = new int [ncols];
	}
	for(i=0;i<ncols;i++) {weight[i]=1.0;}
	
	double* total_centroid = new double [ncols];
	if (Analysis==3||Analysis==4){
		for(i=0;i<ncols;i++){
			double sum=0;
			for (j = 0; j<nrows; j++){
				sum+=data[j][i];
			}
			total_centroid[i]=sum/nrows;
		}
	}
	kcluster(nclusters,nrows,ncols,data,mask,weight,transpose,npass,method,dist,clusterid,&error,&ifound);
	//printf ("Solution found %d times; ", ifound);
	//printf ("within-cluster sum of distances is %f\n", error); // (Unable to figure out what this measure corresponds to)
	//printf ("Cluster assignments:\n");
	//for (i = 0; i < nrows; i++){
	//	printf ("Gene %d: cluster %d\n", i, clusterid[i]);
	//}
	//printf ("\n");
	index = new int *[nclusters];
	count = new int [nclusters];
	for (i = 0; i < nclusters; i++) {count[i] = 0;}
	for (i = 0; i < nrows; i++) {count[clusterid[i]]++;} //"count" contains the number of elements in the cluster
	for (i = 0; i < nclusters; i++) {index[i] = new int [count[i]];}
	for (i = 0; i < nclusters; i++) {count[i] = 0;}
	for (i = 0; i < nrows; i++){
		int id = clusterid[i];
		index[id][count[id]] = i;
		count[id]++;
	}
	getclustercentroids(nclusters, nrows, ncols, data, mask, clusterid, cdata, cmask, 0, 'a');
	//_____________________________________________________________________________________	
	// Implementation of Calinski-Harabasz Index
		double CH_Index;
		double CH_Index_squared;
		double SSW=0; // Total Sum of Square distances Within clusters
		double SSB=0; // Total Sum of Square distances Between clusters
		double SDB=0; // Total Sum of simple Distances Within clusters
		double SDW=0; // Total Sum of simple Distances Between clusters (SSW=total_variation-SSB)
		//double total_variation=0;
		//for(i=0; i<nrows; i++){
		//	for (j=0; j<ncols; j++){
		//		total_variation+=pow((data[i][j]-total_centroid[j]),2);
		//	}
		//}
		for(i=0; i<nrows; i++){
			double SDW_temp=0;
			for (j=0; j<ncols; j++){
				SDW_temp+=pow((data[i][j]-cdata[clusterid[i]][j]),2);
			}
			SSW+=SDW_temp;
			SDW_temp=sqrt(SDW_temp);
			SDW+=SDW_temp;
		}
		// for(i=0; i<nclusters; i++){//An equivalent way to do assess SSW
		// 	double SSW_cluster=0;
		// 	for(int j=0;j<count[i];j++){
		// 		for(int k=0;k<ncols;k++){
		// 			SSW_cluster+=pow((data[index[i][j]][k]-cdata[i][k]),2);
		// 		}
		// 	}
		// 	SSW+=SSW_cluster;
		// }
		for(i=0; i<nclusters; i++){
			double SDB_temp=0;
			for (j=0; j<ncols; j++){
				if (Analysis==1||Analysis==2){
					SDB_temp+=pow((cdata[i][j]),2);
				}
				if (Analysis==3||Analysis==4){
					SDB_temp+=pow((cdata[i][j]-total_centroid[j]),2);
				}
			}
			SSB+=count[i]*SDB_temp;
			SDB_temp=sqrt(SDB_temp);
			SDB+=count[i]*SDB_temp;
		}
		CH_Index_squared=(SSB/(nclusters-1))/(SSW/(nrows-nclusters));
		CH_Index=(SDB/(nclusters-1))/(SDW/(nrows-nclusters));//This is the one used in FASS
	
	F_Fisher=CH_Index_squared;
	
	if(Index_type=="CH_Index"||Index_type=="CH_Index_squared"){
		if(Index_type=="CH_Index"        ){Index_value=CH_Index;}
		if(Index_type=="CH_Index_squared"){Index_value=CH_Index_squared;}
	}	
		//_____________________________________________________________________________________
	if(Index_type=="DB_Index"||Index_type=="DB_Index_squared"){// Implementation of Davies-Bouldin Index
		double DB_Index;
		double DB_Index_squared;
		double** R         = new double *[nclusters];
		double** R_squared = new double *[nclusters];
		for (i = 0; i < nclusters; i++){
			R[i]         = new double [nclusters];
			R_squared[i] = new double [nclusters];
		}
		double SDW_total         =0;
		double SDW_total_squared =0;
		double *SDW_each         = new double [nclusters];
		double *SDW_each_squared = new double [nclusters];
		
		for(i=0; i<nclusters; i++){
			double SDW_each_temp         =0;
			double SDW_each_temp_squared =0;
			
			for (j=0; j<count[i]; j++ ){
				double SDW_each_temp1=0;
				for (k=0; k<ncols; k++){
					SDW_each_temp1+=pow((data[index[i][j]][k]-cdata[i][k]),2);
				}
				SDW_each_temp_squared +=SDW_each_temp1;
				SDW_each_temp1         =sqrt(SDW_each_temp1);
				SDW_each_temp         +=SDW_each_temp1;
			}
			SDW_total            +=SDW_each_temp;
			SDW_each[i]           =SDW_each_temp/count[i] ;
			SDW_total_squared    +=SDW_each_temp_squared;
			SDW_each_squared[i]   =SDW_each_temp_squared/count[i] ;		
		}
		for(i=0; i<nclusters; i++){
			for (j=0; j<nclusters; j++){
				if(i==j){R[i][j]=0;R_squared[i][j]=0;}
				else{
					double dist=0;
					for (k=0; k<ncols; k++){
						dist+=pow(cdata[i][k]-cdata[j][k],2);
					}
					dist=sqrt(dist);
					R[i][j]        =(SDW_each[i]+SDW_each[j])/dist;
					R_squared[i][j]=(SDW_each_squared[i]+SDW_each_squared[j])/dist;
				}
			}
		}
		if(Index_type=="DB_Index"){
			DB_Index=0;
			for(i=0; i<nclusters; i++){
				double R_max=R[i][0];
				for (j=1; j<nclusters; j++){
					if (R[i][j]>R_max){
						R_max=R[i][j];
					}
				}
				DB_Index+=R_max;
			}
			DB_Index=DB_Index/nclusters;//This is the one used in FASS
			Index_value=DB_Index;
		}
		if(Index_type=="DB_Index_squared"){
			DB_Index_squared=0;
			for(i=0; i<nclusters; i++){
				double R_max_squared=R_squared[i][0];
				for (j=1; j<nclusters; j++){
					if (R_squared[i][j]>R_max_squared){
						R_max_squared=R_squared[i][j];
					}
				}
				DB_Index_squared+=R_max_squared;
			}
			DB_Index_squared=DB_Index_squared/nclusters;
			Index_value=DB_Index_squared;
		}
		for (i = 0; i < nclusters; i++){
			delete [] R[i];
			delete [] R_squared[i];
		}
		delete [] R;
		delete [] R_squared;
	}
		//_____________________________________________________________________________________
	if(Index_type=="C_Index"){
		// Implementation of C-Index
		double C_Index;
		double** distMatrix;
		distMatrix = distancematrix(nrows, ncols, data, mask, weight, 'e', 0);
		if (!distMatrix) { printf ("Insufficient memory to store the distance matrix\n");}
		
		int total_pairwise_distances=nrows*(nrows-1)/2;
		double *distances_distribution;
		distances_distribution= new double[total_pairwise_distances];
		int count_index=-1;
		double Sum_pairwise_distances_in_the_same_cluster=0;
		int number_of_pairwise_distances_in_the_same_cluster=0;
		for(i=1; i<nrows; i++){
			for(j=0; j<i; j++){
				count_index++;
				distances_distribution[count_index]=distMatrix[i][j];
				if(clusterid[i]==clusterid[j]){
					Sum_pairwise_distances_in_the_same_cluster+=distMatrix[i][j];
					number_of_pairwise_distances_in_the_same_cluster++;
				}
			}
		}
		combsort(distances_distribution,total_pairwise_distances);
		double Sum_of_smallest_p_pairwise_distances=0;
		double Sum_of_greatest_p_pairwise_distances=0;
		for (i=0; i<number_of_pairwise_distances_in_the_same_cluster; i++){
			Sum_of_smallest_p_pairwise_distances+=distances_distribution[i];
		}
		for (i=total_pairwise_distances-1; i>(total_pairwise_distances-1-number_of_pairwise_distances_in_the_same_cluster); i--){
			Sum_of_greatest_p_pairwise_distances+=distances_distribution[i];
		}
		C_Index=(Sum_pairwise_distances_in_the_same_cluster-Sum_of_smallest_p_pairwise_distances)/(Sum_of_greatest_p_pairwise_distances-Sum_of_smallest_p_pairwise_distances);
		Index_value=C_Index;
	}
	//_____________________________________________________________________________________
	
	for (i = 0; i < nclusters; i++) { delete [] index[i];}
	delete [] index;
	delete [] count;
	for (i = 0; i < nclusters; i++){
		delete [] cdata[i];
		delete [] cmask[i];
	}
	delete [] cdata;
	delete [] cmask;
	delete [] weight;
	return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//**************************************************************************************************************************************
void kmeans_v2_2(
	int nrows,
	int ncols,
	double** data,
	int** mask,
	const int nclusters,
	int Analysis,
	string Index_type,
	int npass,
	double npass_threshold,
	int kmeans_repeats,
	int & flag_unstable_clustering,
	int* & clusterid,
	int & ifound,
	double & Index_value,
	double & F_Fisher
)
{
        int i,j,k;
	const int transpose = 0;
	const char dist = 'e';
	const char method = 'a';
	ifound = 0;
	double error;
	int** index;
	int* count;
	double* weight = new double [ncols];
	double** cdata = new double *[nclusters];
	int** cmask = new int *[nclusters];
	for (i = 0; i < nclusters; i++){
		cdata[i] = new double [ncols];
		cmask[i] = new int [ncols];
	}
	for(i=0;i<ncols;i++) {weight[i]=1.0;}
	
	int* ifound_temp= new int [kmeans_repeats];
	int** clusterid_temp = new int *[kmeans_repeats];
	for(k=0; k<kmeans_repeats; k++){
	  clusterid_temp[k] = new int [nrows];
	}
	for(k=0; k<kmeans_repeats; k++){
	  kcluster(nclusters,nrows,ncols,data,mask,weight,transpose,npass,method,dist,clusterid_temp[k],&error,&ifound_temp[k]);
	  if(ifound_temp[k] < npass*npass_threshold){
	    flag_unstable_clustering=1;
	    //cout << "Kmeans_Internal_Messeage\tNum_of_groups: " << nclusters << ";\tNum_of_repeat: " << k+1 << ";\tStable_clustering: NO; ifound=" << ifound_temp[k] << " < threshold= " << npass*npass_threshold << endl;
	    for(int k2=0; k2<kmeans_repeats; k2++){ delete [] clusterid_temp[k2];}
	    delete [] clusterid_temp;
	    for (i = 0; i < nclusters; i++){
	      delete [] cdata[i];
	      delete [] cmask[i];
	    }
	    delete [] cdata;
	    delete [] cmask;
	    delete [] weight;
	    delete [] ifound_temp;
	    return;
	  }
	  if(k>0){
	    for(i=0;i<nrows-1;i++){
	      for(j=1+1;j<nrows;j++){
		if( ((clusterid_temp[k][i]== clusterid_temp[k][j]) and (clusterid_temp[k-1][i]!= clusterid_temp[k-1][j])) or
		    ((clusterid_temp[k][i]!= clusterid_temp[k][j]) and (clusterid_temp[k-1][i]== clusterid_temp[k-1][j])) )
		  {
		    flag_unstable_clustering=1;

		    for(int k2=0; k2<kmeans_repeats; k2++){ delete [] clusterid_temp[k2];}
		    delete [] clusterid_temp;
		    for (i = 0; i < nclusters; i++){
		      delete [] cdata[i];
		      delete [] cmask[i];
		    }
		    delete [] cdata;
		    delete [] cmask;
		    delete [] weight;
		    delete [] ifound_temp;
		    return;
		  }
	      }
	    }
	  }
	}
	flag_unstable_clustering=0;
	//Assessing ifound as the mean of ifounds from the kmeans_repeats resulting in clusterid_temp[] vector
	ifound=0;
	for(k=0; k<kmeans_repeats; k++){ifound=ifound+ifound_temp[k];}
	ifound=ifound/kmeans_repeats;
	//Relabeling first clusterid_temp as starting by 0
	int * equivalences_between_tags= new int [nclusters];
	for (i = 0; i < nclusters; i++) {equivalences_between_tags[i]=-1;}
	int equivalences_between_tags_filled=0;
	for (i = 0; i < nrows; i++) {
	  if(equivalences_between_tags[clusterid_temp[0][i]]==-1){equivalences_between_tags[clusterid_temp[0][i]]=equivalences_between_tags_filled;equivalences_between_tags_filled++;}
	}
	for (i = 0; i < nrows; i++) {
	  //cout << clusterid_temp[0][i] << "\t" << equivalences_between_tags[clusterid_temp[0][i]] << endl;
	  clusterid[i]=equivalences_between_tags[clusterid_temp[0][i]];
	}

	index = new int *[nclusters];
	count = new int [nclusters];
	for (i = 0; i < nclusters; i++) {count[i] = 0;}
	for (i = 0; i < nrows; i++) {count[clusterid[i]]++;} //"count" contains the number of elements in the cluster
	for (i = 0; i < nclusters; i++) {index[i] = new int [count[i]];}
	for (i = 0; i < nclusters; i++) {count[i] = 0;}
	for (i = 0; i < nrows; i++){
		int id = clusterid[i];
		index[id][count[id]] = i;
		count[id]++;
	}
	getclustercentroids(nclusters, nrows, ncols, data, mask, clusterid, cdata, cmask, 0, 'a');
	//_____________________________________________________________________________________	
	// Implementation of Calinski-Harabasz Index
		double CH_Index;
		double CH_Index_squared;
		double SSW=0; // Total Sum of Square distances Within clusters
		double SSB=0; // Total Sum of Square distances Between clusters
		double SDB=0; // Total Sum of simple Distances Within clusters
		double SDW=0; // Total Sum of simple Distances Between clusters (SSW=total_variation-SSB)
		//double total_variation=0;
		//for(i=0; i<nrows; i++){
		//	for (j=0; j<ncols; j++){
		//		total_variation+=pow((data[i][j]-total_centroid[j]),2);
		//	}
		//}
		for(i=0; i<nrows; i++){
			double SDW_temp=0;
			for (j=0; j<ncols; j++){
				SDW_temp+=pow((data[i][j]-cdata[clusterid[i]][j]),2);
			}
			SSW+=SDW_temp;
			SDW_temp=sqrt(SDW_temp);
			SDW+=SDW_temp;
		}
		// for(i=0; i<nclusters; i++){//An equivalent way to do assess SSW
		// 	double SSW_cluster=0;
		// 	for(int j=0;j<count[i];j++){
		// 		for(int k=0;k<ncols;k++){
		// 			SSW_cluster+=pow((data[index[i][j]][k]-cdata[i][k]),2);
		// 		}
		// 	}
		// 	SSW+=SSW_cluster;
		// }
		double* total_centroid = new double [ncols];
		if (Analysis==3||Analysis==4){
		  for(i=0;i<ncols;i++){
		    double sum=0;
		    for (j = 0; j<nrows; j++){
		      sum+=data[j][i];
		    }
		    total_centroid[i]=sum/nrows;
		  }
		}
		for(i=0; i<nclusters; i++){
			double SDB_temp=0;
			for (j=0; j<ncols; j++){
				if (Analysis==1||Analysis==2){
					SDB_temp+=pow((cdata[i][j]),2);
				}
				if (Analysis==3||Analysis==4){
					SDB_temp+=pow((cdata[i][j]-total_centroid[j]),2);
				}
			}
			SSB+=count[i]*SDB_temp;
			SDB_temp=sqrt(SDB_temp);
			SDB+=count[i]*SDB_temp;
		}
		CH_Index_squared=(SSB/(nclusters-1))/(SSW/(nrows-nclusters));
		CH_Index=(SDB/(nclusters-1))/(SDW/(nrows-nclusters));//This is the one used in FASS
	
	F_Fisher=CH_Index_squared;
	
	if(Index_type=="CH_Index"||Index_type=="CH_Index_squared"){
		if(Index_type=="CH_Index"        ){Index_value=CH_Index;}
		if(Index_type=="CH_Index_squared"){Index_value=CH_Index_squared;}
	}	
		//_____________________________________________________________________________________
	if(Index_type=="DB_Index"||Index_type=="DB_Index_squared"){// Implementation of Davies-Bouldin Index
		double DB_Index;
		double DB_Index_squared;
		double** R         = new double *[nclusters];
		double** R_squared = new double *[nclusters];
		for (i = 0; i < nclusters; i++){
			R[i]         = new double [nclusters];
			R_squared[i] = new double [nclusters];
		}
		double SDW_total         =0;
		double SDW_total_squared =0;
		double *SDW_each         = new double [nclusters];
		double *SDW_each_squared = new double [nclusters];
		
		for(i=0; i<nclusters; i++){
			double SDW_each_temp         =0;
			double SDW_each_temp_squared =0;
			
			for (j=0; j<count[i]; j++ ){
				double SDW_each_temp1=0;
				for (k=0; k<ncols; k++){
					SDW_each_temp1+=pow((data[index[i][j]][k]-cdata[i][k]),2);
				}
				SDW_each_temp_squared +=SDW_each_temp1;
				SDW_each_temp1         =sqrt(SDW_each_temp1);
				SDW_each_temp         +=SDW_each_temp1;
			}
			SDW_total            +=SDW_each_temp;
			SDW_each[i]           =SDW_each_temp/count[i] ;
			SDW_total_squared    +=SDW_each_temp_squared;
			SDW_each_squared[i]   =SDW_each_temp_squared/count[i] ;		
		}
		for(i=0; i<nclusters; i++){
			for (j=0; j<nclusters; j++){
				if(i==j){R[i][j]=0;R_squared[i][j]=0;}
				else{
					double dist=0;
					for (k=0; k<ncols; k++){
						dist+=pow(cdata[i][k]-cdata[j][k],2);
					}
					dist=sqrt(dist);
					R[i][j]        =(SDW_each[i]+SDW_each[j])/dist;
					R_squared[i][j]=(SDW_each_squared[i]+SDW_each_squared[j])/dist;
				}
			}
		}
		if(Index_type=="DB_Index"){
			DB_Index=0;
			for(i=0; i<nclusters; i++){
				double R_max=R[i][0];
				for (j=1; j<nclusters; j++){
					if (R[i][j]>R_max){
						R_max=R[i][j];
					}
				}
				DB_Index+=R_max;
			}
			DB_Index=DB_Index/nclusters;//This is the one used in FASS
			Index_value=DB_Index;
		}
		if(Index_type=="DB_Index_squared"){
			DB_Index_squared=0;
			for(i=0; i<nclusters; i++){
				double R_max_squared=R_squared[i][0];
				for (j=1; j<nclusters; j++){
					if (R_squared[i][j]>R_max_squared){
						R_max_squared=R_squared[i][j];
					}
				}
				DB_Index_squared+=R_max_squared;
			}
			DB_Index_squared=DB_Index_squared/nclusters;
			Index_value=DB_Index_squared;
		}
		for (i = 0; i < nclusters; i++){
			delete [] R[i];
			delete [] R_squared[i];
		}
		delete [] R;
		delete [] R_squared;
	}
		//_____________________________________________________________________________________
	if(Index_type=="C_Index"){
		// Implementation of C-Index
		double C_Index;
		double** distMatrix;
		distMatrix = distancematrix(nrows, ncols, data, mask, weight, 'e', 0);
		if (!distMatrix) { printf ("Insufficient memory to store the distance matrix\n");}
		
		int total_pairwise_distances=nrows*(nrows-1)/2;
		double *distances_distribution;
		distances_distribution= new double[total_pairwise_distances];
		int count_index=-1;
		double Sum_pairwise_distances_in_the_same_cluster=0;
		int number_of_pairwise_distances_in_the_same_cluster=0;
		for(i=1; i<nrows; i++){
			for(j=0; j<i; j++){
				count_index++;
				distances_distribution[count_index]=distMatrix[i][j];
				if(clusterid[i]==clusterid[j]){
					Sum_pairwise_distances_in_the_same_cluster+=distMatrix[i][j];
					number_of_pairwise_distances_in_the_same_cluster++;
				}
			}
		}
		combsort(distances_distribution,total_pairwise_distances);
		double Sum_of_smallest_p_pairwise_distances=0;
		double Sum_of_greatest_p_pairwise_distances=0;
		for (i=0; i<number_of_pairwise_distances_in_the_same_cluster; i++){
			Sum_of_smallest_p_pairwise_distances+=distances_distribution[i];
		}
		for (i=total_pairwise_distances-1; i>(total_pairwise_distances-1-number_of_pairwise_distances_in_the_same_cluster); i--){
			Sum_of_greatest_p_pairwise_distances+=distances_distribution[i];
		}
		C_Index=(Sum_pairwise_distances_in_the_same_cluster-Sum_of_smallest_p_pairwise_distances)/(Sum_of_greatest_p_pairwise_distances-Sum_of_smallest_p_pairwise_distances);
		Index_value=C_Index;
	}
	//_____________________________________________________________________________________
	
	for (i = 0; i < nclusters; i++) { delete [] index[i];}
	delete [] index;
	delete [] count;
	for (i = 0; i < nclusters; i++){
		delete [] cdata[i];
		delete [] cmask[i];
	}
	delete [] cdata;
	delete [] cmask;
	delete [] weight;
	delete [] total_centroid;
	return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void interaction_coding(
	Matrix A,
	vector<string> position_aminoacid_vector,
	vector<double> sum_of_elements_by_column_vector,
	double flag_for_interaction,
	Matrix & B,
	vector<int> & position_aminoacid_interaction_vector_first_element,
	vector<int> & position_aminoacid_interaction_vector_second_element,
	vector<double> & sum_of_elements_by_column_interaction_vector
)
{
	int maximum_size=(position_aminoacid_vector.size()*(position_aminoacid_vector.size()-1))/2;
	int number_of_sequences=A.Nrows();
	Matrix Auxiliar(number_of_sequences,maximum_size);
	position_aminoacid_interaction_vector_first_element.reserve(maximum_size);
	position_aminoacid_interaction_vector_second_element.reserve(maximum_size);
	
	Matrix Z=A.t()*A;
	int number_of_ones=0;
	
	double threshold=floor(flag_for_interaction*(2.0/3.0));
	int number_of_assignments=-1;
	for (int i=0; i<A.Ncols()-1;i++){
		for (int j=i+1;j<A.Ncols();j++ ){
			if(Z[i][j]>=threshold){
				number_of_ones=0;
				number_of_assignments++;
				for (int k=0;k<(number_of_sequences);k++ ){	
					Auxiliar[k][number_of_assignments]= A[k][i]*A[k][j];
					if (Auxiliar[k][number_of_assignments]==1){
						number_of_ones +=1;
					}
				}
				position_aminoacid_interaction_vector_first_element.push_back(i);
				position_aminoacid_interaction_vector_second_element.push_back(j);
				sum_of_elements_by_column_interaction_vector.push_back(number_of_ones);
			}
		}
	}
	B=Auxiliar.SubMatrix(1,A.Nrows(),1,number_of_assignments+1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void kmeans_index_evaluation(
	int nrows,
	int ncols,
	double** data,
	int** mask,
	const int nclusters,
	int Analysis,
	double & CH_Index_squared,
	double &  CH_Index,
	double &  DB_Index,
	double &  DB_Index_squared,
	double &  C_Index,
	int* & clusterid
)
{
	int i,j,k;
	const int transpose = 0;
	const char dist = 'e';
	const char method = 'a';
	int npass = 500;
	int ifound = 0;
	double error;
	int** index;
	int* count;
	double* weight = new double [ncols];
	double** cdata = new double *[nclusters];
	int** cmask = new int *[nclusters];
	for (i = 0; i < nclusters; i++){
		cdata[i] = new double [ncols];
		cmask[i] = new int [ncols];
	}
	for(i=0;i<ncols;i++) {weight[i]=1.0;}
	
	double* total_centroid = new double [ncols];
	if (Analysis==3||Analysis==4){
		for(i=0;i<ncols;i++){
			double sum=0;
			for (j = 0; j<nrows; j++){
				sum+=data[j][i];
			}
			total_centroid[i]=sum/nrows;
		}
	}
	kcluster(nclusters,nrows,ncols,data,mask,weight,transpose,npass,method,dist,clusterid,&error,&ifound);
	//printf ("Solution found %d times; ", ifound);
	//printf ("within-cluster sum of distances is %f\n", error); // (Unable to figure out what this measure corresponds to)
	//printf ("Cluster assignments:\n");
	//for (i = 0; i < nrows; i++){
	//	printf ("Gene %d: cluster %d\n", i, clusterid[i]);
	//}
	//printf ("\n");
	index = new int *[nclusters];
	count = new int [nclusters];
	for (i = 0; i < nclusters; i++) {count[i] = 0;}
	for (i = 0; i < nrows; i++) {count[clusterid[i]]++;} //"count" contains the number of elements in the cluster
	for (i = 0; i < nclusters; i++) {index[i] = new int [count[i]];}
	for (i = 0; i < nclusters; i++) {count[i] = 0;}
	for (i = 0; i < nrows; i++){
		int id = clusterid[i];
		index[id][count[id]] = i;
		count[id]++;
	}
	getclustercentroids(nclusters, nrows, ncols, data, mask, clusterid, cdata, cmask, 0, 'a');
	//_____________________________________________________________________________________	
	// Implementation of Calinski-Harabasz Index
		double SSW=0; // Total Sum of Square distances Within clusters
		double SSB=0; // Total Sum of Square distances Between clusters
		double SDB=0; // Total Sum of simple Distances Within clusters
		double SDW=0; // Total Sum of simple Distances Between clusters (SSW=total_variation-SSB)
		//double total_variation=0;
		//for(i=0; i<nrows; i++){
		//	for (j=0; j<ncols; j++){
		//		total_variation+=pow((data[i][j]-total_centroid[j]),2);
		//	}
		//}
		for(i=0; i<nrows; i++){
			double SDW_temp=0;
			for (j=0; j<ncols; j++){
				SDW_temp+=pow((data[i][j]-cdata[clusterid[i]][j]),2);
			}
			SSW+=SDW_temp;
			SDW_temp=sqrt(SDW_temp);
			SDW+=SDW_temp;
		}
		// for(i=0; i<nclusters; i++){//An equivalent way to do assess SSW
		// 	double SSW_cluster=0;
		// 	for(int j=0;j<count[i];j++){
		// 		for(int k=0;k<ncols;k++){
		// 			SSW_cluster+=pow((data[index[i][j]][k]-cdata[i][k]),2);
		// 		}
		// 	}
		// 	SSW+=SSW_cluster;
		// }
		for(i=0; i<nclusters; i++){
			double SDB_temp=0;
			for (j=0; j<ncols; j++){
				if (Analysis==1||Analysis==2){
					SDB_temp+=pow((cdata[i][j]),2);
				}
				if (Analysis==3||Analysis==4){
					SDB_temp+=pow((cdata[i][j]-total_centroid[j]),2);
				}
			}
			SSB+=count[i]*SDB_temp;
			SDB_temp=sqrt(SDB_temp);
			SDB+=count[i]*SDB_temp;
		}
		CH_Index_squared=(SSB/(nclusters-1))/(SSW/(nrows-nclusters));
		CH_Index=(SDB/(nclusters-1))/(SDW/(nrows-nclusters));//This is the one used in FASS	
	//_____________________________________________________________________________________
	// Implementation of Davies-Bouldin Index
		double** R         = new double *[nclusters];
		double** R_squared = new double *[nclusters];
		for (i = 0; i < nclusters; i++){
			R[i]         = new double [nclusters];
			R_squared[i] = new double [nclusters];
		}
		double SDW_total         =0;
		double SDW_total_squared =0;
		double *SDW_each         = new double [nclusters];
		double *SDW_each_squared = new double [nclusters];
		
		for(i=0; i<nclusters; i++){
			double SDW_each_temp         =0;
			double SDW_each_temp_squared =0;
			
			for (j=0; j<count[i]; j++ ){
				double SDW_each_temp1=0;
				for (k=0; k<ncols; k++){
					SDW_each_temp1+=pow((data[index[i][j]][k]-cdata[i][k]),2);
				}
				SDW_each_temp_squared +=SDW_each_temp1;
				SDW_each_temp1         =sqrt(SDW_each_temp1);
				SDW_each_temp         +=SDW_each_temp1;
			}
			SDW_total            +=SDW_each_temp;
			SDW_each[i]           =SDW_each_temp/count[i] ;
			SDW_total_squared    +=SDW_each_temp_squared;
			SDW_each_squared[i]   =SDW_each_temp_squared/count[i] ;		
		}
		for(i=0; i<nclusters; i++){
			for (j=0; j<nclusters; j++){
				if(i==j){R[i][j]=0;R_squared[i][j]=0;}
				else{
					double dist=0;
					for (k=0; k<ncols; k++){
						dist+=pow(cdata[i][k]-cdata[j][k],2);
					}
					dist=sqrt(dist);
					R[i][j]        =(SDW_each[i]+SDW_each[j])/dist;
					R_squared[i][j]=(SDW_each_squared[i]+SDW_each_squared[j])/dist;
				}
			}
		}

		DB_Index=0;
		for(i=0; i<nclusters; i++){
			double R_max=R[i][0];
			for (j=1; j<nclusters; j++){
				if (R[i][j]>R_max){
					R_max=R[i][j];
				}
			}
			DB_Index+=R_max;
		}
		DB_Index=DB_Index/nclusters;//This is the one used in FASS


		DB_Index_squared=0;
		for(i=0; i<nclusters; i++){
			double R_max_squared=R_squared[i][0];
			for (j=1; j<nclusters; j++){
				if (R_squared[i][j]>R_max_squared){
					R_max_squared=R_squared[i][j];
				}
			}
			DB_Index_squared+=R_max_squared;
		}
		DB_Index_squared=DB_Index_squared/nclusters;

		for (i = 0; i < nclusters; i++){
			delete [] R[i];
			delete [] R_squared[i];
		}
		delete [] R;
		delete [] R_squared;
	
	//_____________________________________________________________________________________
	// Implementation of C-Index
		double** distMatrix;
		distMatrix = distancematrix(nrows, ncols, data, mask, weight, 'e', 0);
		if (!distMatrix) { printf ("Insufficient memory to store the distance matrix\n");}
		
		int total_pairwise_distances=nrows*(nrows-1)/2;
		double *distances_distribution;
		distances_distribution= new double[total_pairwise_distances];
		int count_index=-1;
		double Sum_pairwise_distances_in_the_same_cluster=0;
		int number_of_pairwise_distances_in_the_same_cluster=0;
		for(i=1; i<nrows; i++){
			for(j=0; j<i; j++){
				count_index++;
				distances_distribution[count_index]=distMatrix[i][j];
				if(clusterid[i]==clusterid[j]){
					Sum_pairwise_distances_in_the_same_cluster+=distMatrix[i][j];
					number_of_pairwise_distances_in_the_same_cluster++;
				}
			}
		}
		combsort(distances_distribution,total_pairwise_distances);
		double Sum_of_smallest_p_pairwise_distances=0;
		double Sum_of_greatest_p_pairwise_distances=0;
		for (i=0; i<number_of_pairwise_distances_in_the_same_cluster; i++){
			Sum_of_smallest_p_pairwise_distances+=distances_distribution[i];
		}
		for (i=total_pairwise_distances-1; i>(total_pairwise_distances-1-number_of_pairwise_distances_in_the_same_cluster); i--){
			Sum_of_greatest_p_pairwise_distances+=distances_distribution[i];
		}
		C_Index=(Sum_pairwise_distances_in_the_same_cluster-Sum_of_smallest_p_pairwise_distances)/(Sum_of_greatest_p_pairwise_distances-Sum_of_smallest_p_pairwise_distances);
		if(C_Index<0.00000001){C_Index=0;}
	//_____________________________________________________________________________________
	
	for (i = 0; i < nclusters; i++) { delete [] index[i];}
	delete [] index;
	delete [] count;
	for (i = 0; i < nclusters; i++){
		delete [] cdata[i];
		delete [] cmask[i];
	}
	delete [] cdata;
	delete [] cmask;
	delete [] weight;
	return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// void ksub_next ( int n, int k, int a[], bool *more )
// //
// //  Purpose:
// //
// //    KSUB_NEXT generates the subsets of size K from a set of size N.
// //
// //  Modified:
// //
// //    29 May 2003
// //
// //  Reference:
// //
// //    Albert Nijenhuis, Herbert Wilf,
// //    Combinatorial Algorithms for Computers and Calculators,
// //    Second Edition,
// //    Academic Press, 1978,
// //    ISBN: 0-12-519260-6,
// //    LC: QA164.N54.
// //
// //  Parameters:
// //
// //    Input, int N, the size of the set from which subsets are drawn.
// //
// //    Input, int K, the desired size of the subsets.  K must
// //    be between 0 and N.
// //
// //    Output, int A[K].  A[I] is the I-th element of the
// //    subset.  Thus A[I] will be an integer between 1 and N.
// //    Note that the routine will return the values in A
// //    in sorted order: 1 <= A[0] < A[1] < ... < A[K-1] <= N
// //
// //    Input/output, bool *MORE.  Set MORE = FALSE before first call
// //    for a new sequence of subsets.  It then is set and remains
// //    TRUE as long as the subset computed on this call is not the
// //    final one.  When the final subset is computed, MORE is set to
// //    FALSE as a signal that the computation is done.
// //
// {
//   int j;
//   static int m = 0;
//   static int m2 = 0;
// 
//   if ( k < 0 || n < k )
//   {
//     cout << "\n";
//     cout << "KSUB_NEXT - Fatal error!\n";
//     cout << "N = " << n << "\n";
//     cout << "K = " << k << "\n";
//     cout << "but 0 <= K <= N is required!\n";
//     exit ( 1 );
//   }
// 
//   if ( !( *more ) )
//   {
//     m2 = 0;
//     m = k;
//   }
//   else
//   {
//     if ( m2 < n-m )
//     {
//       m = 0;
//     }
//     m = m + 1;
//     m2 = a[k-m];
//   }
// 
//   for ( j = 1; j <= m; j++ )
//   {
//     a[k+j-m-1] = m2 + j;
//   }
// 
//   *more = ( a[0] != (n-k+1) );
// 
//   return;
// }
// 
