//
//  main.cpp
//  dna
//
//  Created by Choo Xuan Wing on 26/10/2019.
//  Copyright © 2019 Choo Xuan Wing. All rights reserved.
//

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
using namespace std;

class tempDNA_DB{
public:
	string SEQ;
	string GID;
	string REF;
	string FileName;
	string Name;
	
};

class DNA_DB{
private:
//	string SEQ;
public:
	string SEQ;
	string GID;
	string REF;
	string FileName;
	string Name;

};

struct input{			// initialise struct to store function outputs
	int count=0;
	int check =0;
	string File;
	string initialSelect;
	string str;
	string tempstr;
	vector<string> fileNumber;
	vector<string> fileNames;
	vector<string> onlyNames;
};

struct initialOptions{
	vector<string> optionNum;
	vector<string> optionName;
	string optionSelect;
	
};

struct helpMenu{
	vector<string> helpNum;
	vector<string> helpName;
};

struct analyse{
	long regions=0, nRegions=0, cRegions=0, basePairs=0, G=0,A=0,T=0,C=0,R=0,Y=0,M=0,K=0,S=0,W=0,H=0,B=0,V=0,D=0,N=0,unknown=0;
	long nRegion, cRegion;
	vector<string> gapRegion;
	vector<string> codeRegion;
	vector<int> NindexStart;
	vector<int> NindexEnd;
	int nStart;
//	vector<int> CindexStart;
//	vector<int> CindexEnd;
};

class DNA_DB dna_db;			// make class global so values can be transferred throughout programme
class tempDNA_DB tempdna_db;
long indexG;			// gobal index value
struct input initial;			// make struct global so its easy to pass values around
struct initialOptions option;
struct helpMenu help;
struct analyse analyse;

void organiseFile(string File){			// function to collect and sort input file names

	stringstream file(File);			// ready string we entered fpr processing
		
	initial.fileNames.push_back("Summary statistics of the DNA database");
	initial.fileNumber.push_back("S");
	
	while (file.good()){			// check if file exists
		initial.count++;			// increment count for check
		initial.fileNumber.push_back(to_string(initial.count));			// numbers to oiption
		string subst;
		getline(file, subst, ',');			// takes input string, reads up till "," and saves to subst. repears untill threres none left
		initial.fileNames.push_back(subst);		// insert subst into vector
		initial.onlyNames.push_back(subst);		// insert subst into vector for summary processing
	}
	
	initial.fileNumber.push_back("Q");
	initial.fileNames.push_back("Quit");
	
}

void printOrganiseFile(){
	
	cout << "\nSelect one of the following options" << endl;
		for(int i = 0; i<initial.fileNames.size(); i++) {
			std::cout << "("<< initial.fileNumber.at(i)
			<< ")\t" << initial.fileNames.at(i) << endl;
		}
	cout<< ">";
}

int checkElement(){
	
	std::vector<string>::iterator check;
	check = std::find(initial.fileNumber.begin(), initial.fileNames.end(), initial.initialSelect);
	if (check == initial.fileNumber.end())
		return 1;
	else
		return 0;
}

void selectFile(string select){
	
	std::vector<string>::iterator command;			// make "command" as an iterator
	command = std::find(initial.fileNumber.begin(), initial.fileNumber.end(), select);			// finds specific element in vector
	if (command == initial.fileNumber.end()){
		cout << "Command not recognised" << endl;
	}
	else if (command != initial.fileNumber.end()){
		auto index = std::distance(initial.fileNumber.begin(), command );			// puts index of input and places them into filenames to print out filename
		indexG=index;
		cout << "Loading " << initial.fileNames.at(index) << "..."<< endl;
		std::ifstream input(initial.fileNames.at(index));		//streams selected file
		if (input.fail() ){
			cout << "Error, could not open file" << endl;
			initial.check =1;			// increment to tell int main() to halt
		}
	}
}

void assignClass(){
	
	string subst;
	std::ifstream input(initial.fileNames.at(indexG));		// streams selected file
	// Input file contents into classes
	string gid, seq, name, ref, temp, tempDNAseq;
	ostringstream ss;		//initialises ss as char array
	ss<<input.rdbuf();		// rdbuf points input(file conts) into ss
	initial.str = ss.str();		// ss as a string puts into struct initial str
	stringstream file(initial.str);
	
	while (!std::getline(file, subst, '|').eof()){		// run when not EOF
		getline(file, gid, '|');
		getline(file, temp, '|');		// temp used as junk variable to store unwanted characters
		getline(file, ref, '|');
		getline(file, name, '\n');
		while(!file.eof()){			// loop to read lines of dna Seq
			getline(file, tempDNAseq);
			seq.append(tempDNAseq);		// append adds new line to existing
		}
		dna_db.SEQ=seq;			// places elements into class
		dna_db.GID=gid;
		dna_db.Name=name;
		dna_db.REF=ref;
	}
	
};

void initialOptions(){
		
	option.optionNum.push_back("H");
	option.optionNum.push_back("S");
	option.optionNum.push_back("1");
	option.optionNum.push_back("2");
	option.optionNum.push_back("3");
	option.optionNum.push_back("4");
	option.optionNum.push_back("5");
	option.optionNum.push_back("R");
	option.optionNum.push_back("Q");
	
	option.optionName.push_back("Help");
	option.optionName.push_back("Summary statistics of the DNA sequence");
	option.optionName.push_back("Analyse gap region");
	option.optionName.push_back("Analyse coded region");
	option.optionName.push_back("Analyse base pair range");
	option.optionName.push_back("Find DNA sequence by manual input");
	option.optionName.push_back("Find DNA sequence by file input");
	option.optionName.push_back("Return to the previous menu");
	option.optionName.push_back("Quit");
	
}

void printInitialOptions(){
	
	cout << "\nSelect one of the following options" << endl;
		for(int i = 0; i<option.optionName.size(); i++) {
			std::cout << "("<< option.optionNum.at(i)
			<< ")\t" << option.optionName.at(i) << endl;
		}
	cout<< ">";
}

void firstSummary(){
	
	long size,i;
	string subst;
	string gid, seq, name, ref, temp, tempDNAseq;
	size= initial.onlyNames.size();
	
	cout << "\nThis DNA database holds " << size << " sequence(s)" << endl;
	for (i=0 ; i<size; i++){
		std::ifstream input(initial.onlyNames.at(i));				//read loaded files
		ostringstream ss;											// makes file into a string
		ss<< input.rdbuf();
		initial.tempstr = ss.str();						//loads string into temp dna class
		stringstream file(initial.tempstr);				// makes file a command for streaming in function
		
		while (!std::getline(file, subst, '|').eof()){
			getline(file,gid,'|');
			getline(file,temp,'|');
			getline(file,ref,'|');
			getline(file,name,'\n');
			
			while (!file.eof()){
				getline(file, tempDNAseq);
				seq.append(tempDNAseq);
			}
			tempdna_db.SEQ=seq;
			tempdna_db.GID=gid;
			tempdna_db.Name=name;
			tempdna_db.REF=ref;
		}
		cout << "\nSequence " << i+1 << ":" <<endl;
		cout << "Name:\t" << tempdna_db.Name << endl;
		cout << "GID:\t" << tempdna_db.GID << endl;
		cout << "REF:\t" << tempdna_db.REF << endl;
		cout << "# base pairs:\t" << tempdna_db.SEQ.length() <<"\n"<< endl;
	}
	
}


void helpMenu(){	// Help Menu
	
	help.helpNum.push_back("Code");
	help.helpNum.push_back("G");
	help.helpNum.push_back("A");
	help.helpNum.push_back("T");
	help.helpNum.push_back("C");
	help.helpNum.push_back("Y");
	help.helpNum.push_back("M");
	help.helpNum.push_back("K");
	help.helpNum.push_back("S");
	help.helpNum.push_back("W");
	help.helpNum.push_back("H");
	help.helpNum.push_back("B");
	help.helpNum.push_back("V");
	help.helpNum.push_back("D");
	help.helpNum.push_back("N");

	help.helpName.push_back("Base Description");
	help.helpName.push_back("Guanine");
	help.helpName.push_back("Adenine");
	help.helpName.push_back("Thymine (Uracil in RNA)");
	help.helpName.push_back("Cytosine");
	help.helpName.push_back("Purine (A or G)");
	help.helpName.push_back("Pyrimidine (C or T or U)");
	help.helpName.push_back("Amino (A or C)");
	help.helpName.push_back("Ketone (G or T)");
	help.helpName.push_back("Strong interaction (C or G)");
	help.helpName.push_back("Weak interaction (A or T)");
	help.helpName.push_back("Not-G (A or C or T) H follows G in the alphabet");
	help.helpName.push_back("Not-A (C or G or T) B follows A in the alphabet");
	help.helpName.push_back("Not-T (not-U) (A or C or G) V follows U in the alphabet");
	help.helpName.push_back("Not-C (A or G or T) D follows C in the alphabet");
	help.helpName.push_back("Any (A or C or G or T)");

	for(int i = 0; i<help.helpNum.size(); i++) {
		cout << help.helpNum.at(i)
		<< "\t\t" << help.helpName.at(i) << endl;
	}
}

void summary(){
	
	// Calculates base pair stuff
	cout << "\nLoading..."<< endl;
	analyse.basePairs= dna_db.SEQ.length();
	
	analyse.G= count(dna_db.SEQ.begin(), dna_db.SEQ.end(),'G');
	analyse.A= count(dna_db.SEQ.begin(), dna_db.SEQ.end(),'A');
	analyse.T= count(dna_db.SEQ.begin(), dna_db.SEQ.end(),'T');
	analyse.C= count(dna_db.SEQ.begin(), dna_db.SEQ.end(),'C');
	analyse.N= count(dna_db.SEQ.begin(), dna_db.SEQ.end(),'N');
	analyse.unknown= analyse.basePairs - analyse.G - analyse.C - analyse.T - analyse.A - analyse.N;
	
	cout << "Sequence identifiers:"<< endl;
	cout << "Name:\t" << dna_db.Name << endl;
	cout << "GID:\t" << dna_db.GID << endl;
	cout << "REF:\t" << dna_db.REF << endl;
	
	cout << "\nRegion characteristics:"<<endl;
	cout << "# regions:\t" << analyse.nRegion+analyse.cRegion << endl;
	cout << "# N regions:\t" << analyse.nRegion << endl;
	cout << "# C regions:\t" << analyse.cRegion << endl;
	
	cout << "\nBase pair characteristics:"<<endl;
	cout << "# base pairs\t" << dna_db.SEQ.length() << endl;
	cout << "G:\t" << analyse.G << endl;
	cout << "A:\t" << analyse.A  << endl;
	cout << "T:\t" << analyse.T  << endl;
	cout << "C:\t" << analyse.C  << endl;
	cout << "R:\t" << analyse.C + analyse.A  << endl;
	cout << "Y:\t" << analyse.C + analyse.T  << endl;
	cout << "M:\t" << analyse.A + analyse.C  << endl;
	cout << "K:\t" << analyse.G + analyse.T  << endl;
	cout << "S:\t" << analyse.C + analyse.G  << endl;
	cout << "W:\t" << analyse.A + analyse.T  << endl;
	cout << "H:\t" << analyse.A + analyse.C + analyse.T << endl;
	cout << "B:\t" << analyse.C + analyse.G + analyse.T << endl;
	cout << "V:\t" << analyse.C + analyse.G + analyse.A << endl;
	cout << "D:\t" << analyse.A + analyse.G + analyse.T << endl;
	cout << "N:\t" << analyse.N  << endl;
	cout << "Unknown:\t" << analyse.unknown  << endl;

}

void analyseRegions(){
	
	int countN=0,nRegion=0,countC=0,cRegion=0,j=0;
	if (dna_db.SEQ.at(j)=='N')			// if seq starts of with N
		analyse.nStart=1;				// makes that 1 so other functions know what to calculate
	else
		analyse.nStart=0;
	for(int i=0; i<dna_db.SEQ.size(); i++){
		if (dna_db.SEQ.at(i)=='N'){				// if char at i == N
			if(countN==0){					// if that is the 1st N encountered
				cRegion++;					// coded region +1
				analyse.NindexStart.push_back(i);		// and location added to vector
			}
			else;
			countN++;					// if not first N encounteres, keep counting number of N
		}
		else{							// if char at i is not  N
			if (countN !=0){			// and N counter is not, that means it has reached the end of non coding region
				analyse.NindexEnd.push_back(i);		// so location of end point into vector
				countN=0;				// reset n counter to 0
				nRegion++;				// increment the number of N regions found
			}
			countC++;		// if N is 0 and char at i is not N, keep incrementing coded char counter
		}
	}
	
	if (analyse.cRegion==0)			// if there is no no coding region
		analyse.cRegion=cRegion+1;
	else{
		analyse.nRegion=nRegion;
		analyse.cRegion=cRegion;
	}
}

void gapRegion(){
	
	int region,start,end,length;
	string sequence;
	cout << "Enter gap region number:\n>";
	cin >> region;
	
	if (analyse.nStart==1){				// if start is N
		start = analyse.NindexStart.at(region-1);		// start is 0th vector of start N index
		end = analyse.NindexEnd.at(region-1);			// End is 0th vector of end N index
		length = analyse.NindexEnd.at(region-1) - analyse.NindexStart.at(region-1);
	}
	else{									// if start is not N
		start = analyse.NindexEnd.at(region-1);			//	Start is 0th vector of end N index
		end = analyse.NindexStart.at(region);			// end is 1st vector of start N index
		length = analyse.NindexStart.at(region) - analyse.NindexEnd.at(region-1);
		
	}
	
	sequence = dna_db.SEQ.substr (start,length);
	
	cout << "Selected Sequence:\n" << "Base pair range: (" << start << "," << end << ")\n" << "Gap region number: " << region << "\n" <<endl;
	cout << "Sequence:\n" << sequence<< endl;
	
}

void codedRegion(){
	
	int region,start,end,length;
	string sequence;
	cout << "Enter coded region number:\n>";
	cin >> region;
	
	if (analyse.nStart==1){					// explaination same as codedregion function
		start = analyse.NindexEnd.at(region-1);
		end = analyse.NindexStart.at(region);
		length = analyse.NindexStart.at(region) - analyse.NindexEnd.at(region-1);
	}
	else{
		start = analyse.NindexStart.at(region-1);
		end = analyse.NindexEnd.at(region-1);
		length = analyse.NindexEnd.at(region-1) - analyse.NindexStart.at(region-1);
	}
	
	sequence = dna_db.SEQ.substr (start,length);

	cout << "Selected Sequence:\n" << "Base pair range: (" << start << "," << end << ")\n" << "Coded region number: " << region << "\n" <<endl;
	cout << "Sequence:\n" << sequence<< endl;
	
}

void findManual(){
	
	string findN;
	long count=0, length, found=0;
	cout << "\nSpecify the DNA sequence nucleotides you would like to find:\n>" ;
//	findN = "hello";
	cin >> findN;
//	dna_db.SEQ="hello";
	length = findN.size();
	for (int i=0 ; i < dna_db.SEQ.size(); i++){			// runs main seq char by char
		for (int j=0;j< findN.size();j++){				// runs 1 char of seq
			if(dna_db.SEQ.at(i) == findN.at(j)){		// if 1st char find == 1st char main seq
				count++;								// it will keep running incrementing char until theres no match
				i++;
				if (count == length){					// if all characters match, increment found and repeat
					found++;
					cout << "Base pair range(" << i-length << "," << i << ")" << endl;
				}
			}
			else{										// l'ets say if char 2 seq != char 2 main, it'll break and start again
				count=0;
				break;
			}
		}
	}
	
	cout << "Total Number of matches found: " << found << endl;
	
	
}

void findFile(){			// funtion to check file match with SEQ
	
	string file, tempDNAseq, seqCompare;
	long loc=0 ;
	stringstream seq(dna_db.SEQ);
	cout << "\n(Warning, do not key in large files as programme would be unresponsive trying to print out match)\n Specify the DNA sequence file you would like to find:\n>";
//	file= "gnFOXF1.fa";
	cin >> file;
	
	ifstream input(file);
	cout << "Loading " << file <<endl;
	if (!input.good()){
		cout << "File not found" << endl;
		findFile();
	}
	else
		cout << file << "..." << endl;
	
	input.ignore(5000,'\n');			// converts sample into string, ignores first line
	while(!input.eof()){			// loop to read lines of new dna Seq
		getline(input, tempDNAseq);
		seqCompare.append(tempDNAseq);		// append adds new line to existing
	}
	cout << "Successful loading of " << file << endl;

	loc = dna_db.SEQ.find(seqCompare);		// make loc the largest storable value initially, and equate to location of match
	if ( loc != dna_db.SEQ.length()){
		cout << "\nBase pair range: (" << loc << "," << dna_db.SEQ.length() << ")" << endl;
		cout << seqCompare;
	}
	else
		cout << "No Match" << endl;
	
	
}
void bpRange(){
	
	vector<int> vectRange;
	string ranges, subst, output;
	int length;
	cout << "Enter a comma ',' separating base pair ranges (NO SPACES):\n" << ">";
	cin >> ranges;
	
	stringstream range(ranges);			//converts range string into vector
	while (range.good()){
		getline (range, subst, ',');
		vectRange.push_back (stoi(subst));
	}
	length =vectRange.at(1)-vectRange.at(0);
	output = dna_db.SEQ.substr(vectRange.at(0),length);			// extracts string based on vector index, (start,length)
	cout << "Selected sequence: \nBase pair range: (" << vectRange.at(0) << "," << vectRange.at(1) << ")" << endl;
	cout << "\nSequence:\n" << output << endl;
									 
}


int main()
{
	cout << "DNA Sequence Database Software" << endl;
	cout << "Specify the name of DNA sequence file names you would like to load. For multiple files, add a ',' between each file name. (Add .fa extension after name, eg chr1.fa and no space in between files)\n>" ;
	
//	initial.File = "chr16.fa";
	cin >> initial.File;			//Enter string of diff files

	organiseFile(initial.File);		//string of files organised/separated
	initialOptions();				// initialise initial options

menu1:
	printOrganiseFile();			// prints which file to choose
//	initial.initialSelect ="1";
	cin >> initial.initialSelect;		//Select input file and put into struct
	
	if (initial.initialSelect=="Q"){			// reads and selects input options
		cout << "Programme ended." << endl;
		return 0;
	}
	else if(initial.initialSelect=="S"){
		firstSummary();
		goto menu1;
	}
	else{
		selectFile(initial.initialSelect);		//File being selected
		if (initial.check == 1){		// if cant open file, 1 would be output
			return 0;
		}
		else
			assignClass();					// selected file issued into class
		
	}
	
	analyseRegions();

option1:
	printInitialOptions();			// Prints options on what to do with file
	cin >> option.optionSelect;		// puts imputs for options into struct to process
	if (option.optionSelect=="Q"){
		cout << "Programme ended." << endl;
		return 0;
	}
	else if (option.optionSelect=="H"){
		helpMenu();
		goto option1;
	}
	else if (option.optionSelect=="S"){
		summary();
		goto option1;
	}
	else if (option.optionSelect=="R"){
		goto menu1;
	}
	else if (option.optionSelect=="1"){
		gapRegion();
		goto option1;

	}
	else if (option.optionSelect=="2"){
		codedRegion();
		goto option1;

	}
	else if (option.optionSelect=="3"){
		bpRange();
		goto option1;

	}
	else if (option.optionSelect == "4"){
		findManual();
		goto option1;

	}
	else if (option.optionSelect == "5"){
		findFile();
		goto option1;
	}
	else{
		cout << "Command not recognised, try again"<< endl;
		goto option1;
	}
}

