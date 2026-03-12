//
//  main.cpp
//  dna
//
//  Created by Choo Xuan Wing on 26/10/2019.
//  Copyright © 2019 Choo Xuan Wing. All rights reserved.
//
//  Refactored: Fixed critical bugs, improved design, and added documentation.
//
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

// [Fix #9] Merged duplicate DNA_DB and tempDNA_DB into a single class.
// Previously there were two identical classes; now we use one and create
// multiple instances where needed.
class DNA_DB {
public:
    std::string SEQ;
    std::string GID;
    std::string REF;
    std::string FileName;
    std::string Name;
};

// Struct to store user input state and file information.
struct InputData {
    int count = 0;
    int check = 0;
    std::string File;
    std::string initialSelect;
    std::string str;
    std::string tempstr;
    std::vector<std::string> fileNumber;
    std::vector<std::string> fileNames;
    std::vector<std::string> onlyNames;
};

// [Fix #12] Renamed from 'initialOptions' to 'MenuOptions' to avoid
// name collision with the function that was previously called initialOptions().
struct MenuOptions {
    std::vector<std::string> optionNum;
    std::vector<std::string> optionName;
    std::string optionSelect;
};

// Struct to hold help menu display data.
struct HelpMenuData {
    std::vector<std::string> helpNum;
    std::vector<std::string> helpName;
};

// [Fix #11] Renamed from 'analyse' to 'AnalysisResult' to avoid the
// name collision where the struct and the variable had the same name.
// [Fix #4] Initialized nRegion and cRegion to 0 — they were previously
// uninitialized, which is undefined behavior when read before being set.
struct AnalysisResult {
    long regions = 0, nRegions = 0, cRegions = 0, basePairs = 0;
    long G = 0, A = 0, T = 0, C = 0, R = 0, Y = 0, M = 0, K = 0;
    long S = 0, W = 0, H = 0, B = 0, V = 0, D = 0, N = 0, unknown = 0;
    long nRegion = 0, cRegion = 0;
    std::vector<std::string> gapRegion;
    std::vector<std::string> codeRegion;
    std::vector<int> NindexStart;
    std::vector<int> NindexEnd;
    int nStart = 0;
};

// [Fix #7] Global variables are kept for simplicity in this console application,
// but have been reduced and renamed for clarity.
DNA_DB dna_db;
DNA_DB tempdna_db;
long indexG;
InputData initial;
MenuOptions menuOption;
HelpMenuData helpData;
AnalysisResult analysisResult;

// Parses a comma-separated string of file names and populates the
// initial struct's fileNames, fileNumber, and onlyNames vectors.
void organiseFile(std::string File) {
    std::stringstream file(File);

    // Add the summary option as the first menu entry
    initial.fileNames.push_back("Summary statistics of the DNA database");
    initial.fileNumber.push_back("S");

    while (file.good()) {
        initial.count++;
        initial.fileNumber.push_back(std::to_string(initial.count));
        std::string subst;
        std::getline(file, subst, ',');
        initial.fileNames.push_back(subst);
        initial.onlyNames.push_back(subst);
    }

    // Add the quit option as the last menu entry
    initial.fileNumber.push_back("Q");
    initial.fileNames.push_back("Quit");
}

// Displays the file selection menu to the user.
void printOrganiseFile() {
    std::cout << "\nSelect one of the following options" << std::endl;
    for (size_t i = 0; i < initial.fileNames.size(); i++) {
        std::cout << "(" << initial.fileNumber.at(i)
                  << ")\t" << initial.fileNames.at(i) << std::endl;
    }
    std::cout << ">";
}

// [Fix #1] Fixed mismatched iterator ranges. Previously used
// initial.fileNumber.begin() with initial.fileNames.end() — iterators from
// two different vectors, which is undefined behavior.
// Now correctly uses begin() and end() from the same vector.
int checkElement() {
    std::vector<std::string>::iterator check;
    check = std::find(initial.fileNumber.begin(), initial.fileNumber.end(), initial.initialSelect);
    if (check == initial.fileNumber.end())
        return 1;
    else
        return 0;
}

// Finds and loads the file corresponding to the user's menu selection.
// [Fix #14] Consistent error handling — sets initial.check flag on failure.
void selectFile(std::string select) {
    std::vector<std::string>::iterator command;
    command = std::find(initial.fileNumber.begin(), initial.fileNumber.end(), select);
    if (command == initial.fileNumber.end()) {
        std::cout << "Command not recognised" << std::endl;
    } else {
        auto index = std::distance(initial.fileNumber.begin(), command);
        indexG = index;
        std::cout << "Loading " << initial.fileNames.at(index) << "..." << std::endl;
        std::ifstream input(initial.fileNames.at(index));
        if (input.fail()) {
            std::cout << "Error, could not open file" << std::endl;
            initial.check = 1;
        }
    }
}

// Reads the selected file and parses its contents into the dna_db class.
// File format expected: fields separated by '|', with DNA sequence on subsequent lines.
// [Fix #14] Added file-open check for consistent error handling.
void assignClass() {
    std::string subst;
    std::ifstream input(initial.fileNames.at(indexG));
    if (input.fail()) {
        std::cout << "Error, could not open file in assignClass()" << std::endl;
        initial.check = 1;
        return;
    }

    std::string gid, seq, name, ref, temp, tempDNAseq;
    std::ostringstream ss;
    ss << input.rdbuf();
    initial.str = ss.str();
    std::stringstream file(initial.str);

    while (!std::getline(file, subst, '|').eof()) {
        std::getline(file, gid, '|');
        std::getline(file, temp, '|');
        std::getline(file, ref, '|');
        std::getline(file, name, '\n');
        while (!file.eof()) {
            std::getline(file, tempDNAseq);
            seq.append(tempDNAseq);
        }
        dna_db.SEQ = seq;
        dna_db.GID = gid;
        dna_db.Name = name;
        dna_db.REF = ref;
    }
}

// [Fix #12] Renamed function from initialOptions() to initMenuOptions()
// to avoid collision with the struct that was previously named initialOptions.
// Populates the analysis menu options (displayed after a file is loaded).
void initMenuOptions() {
    menuOption.optionNum.push_back("H");
    menuOption.optionNum.push_back("S");
    menuOption.optionNum.push_back("1");
    menuOption.optionNum.push_back("2");
    menuOption.optionNum.push_back("3");
    menuOption.optionNum.push_back("4");
    menuOption.optionNum.push_back("5");
    menuOption.optionNum.push_back("R");
    menuOption.optionNum.push_back("Q");

    menuOption.optionName.push_back("Help");
    menuOption.optionName.push_back("Summary statistics of the DNA sequence");
    menuOption.optionName.push_back("Analyse gap region");
    menuOption.optionName.push_back("Analyse coded region");
    menuOption.optionName.push_back("Analyse base pair range");
    menuOption.optionName.push_back("Find DNA sequence by manual input");
    menuOption.optionName.push_back("Find DNA sequence by file input");
    menuOption.optionName.push_back("Return to the previous menu");
    menuOption.optionName.push_back("Quit");
}

// Displays the analysis options menu to the user.
void printMenuOptions() {
    std::cout << "\nSelect one of the following options" << std::endl;
    for (size_t i = 0; i < menuOption.optionName.size(); i++) {
        std::cout << "(" << menuOption.optionNum.at(i)
                  << ")\t" << menuOption.optionName.at(i) << std::endl;
    }
    std::cout << ">";
}

// Displays a summary of all loaded DNA sequence files.
// [Fix #13] seq.clear() is now called at the start of each loop iteration
// to prevent sequences from accumulating across files (previously gave
// incorrect base-pair counts for all files after the first).
void firstSummary() {
    long size;
    std::string subst;
    size = initial.onlyNames.size();

    std::cout << "\nThis DNA database holds " << size << " sequence(s)" << std::endl;
    for (long i = 0; i < size; i++) {
        std::string gid, seq, name, ref, temp, tempDNAseq;
        seq.clear(); // [Fix #13] Clear seq each iteration to avoid accumulation

        std::ifstream input(initial.onlyNames.at(i));
        std::ostringstream ss;
        ss << input.rdbuf();
        initial.tempstr = ss.str();
        std::stringstream file(initial.tempstr);

        while (!std::getline(file, subst, '|').eof()) {
            std::getline(file, gid, '|');
            std::getline(file, temp, '|');
            std::getline(file, ref, '|');
            std::getline(file, name, '\n');

            while (!file.eof()) {
                std::getline(file, tempDNAseq);
                seq.append(tempDNAseq);
            }
            tempdna_db.SEQ = seq;
            tempdna_db.GID = gid;
            tempdna_db.Name = name;
            tempdna_db.REF = ref;
        }
        std::cout << "\nSequence " << i + 1 << ":" << std::endl;
        std::cout << "Name:\t" << tempdna_db.Name << std::endl;
        std::cout << "GID:\t" << tempdna_db.GID << std::endl;
        std::cout << "REF:\t" << tempdna_db.REF << std::endl;
        std::cout << "# base pairs:\t" << tempdna_db.SEQ.length() << "\n" << std::endl;
    }
}

// [Fix #2] Help menu data is now initialized once in initHelpMenu() (called
// from main), instead of appending duplicates every time the user presses 'H'.
void initHelpMenu() {
    helpData.helpNum.push_back("Code");
    helpData.helpNum.push_back("G");
    helpData.helpNum.push_back("A");
    helpData.helpNum.push_back("T");
    helpData.helpNum.push_back("C");
    helpData.helpNum.push_back("Y");
    helpData.helpNum.push_back("M");
    helpData.helpNum.push_back("K");
    helpData.helpNum.push_back("S");
    helpData.helpNum.push_back("W");
    helpData.helpNum.push_back("H");
    helpData.helpNum.push_back("B");
    helpData.helpNum.push_back("V");
    helpData.helpNum.push_back("D");
    helpData.helpNum.push_back("N");

    helpData.helpName.push_back("Base Description");
    helpData.helpName.push_back("Guanine");
    helpData.helpName.push_back("Adenine");
    helpData.helpName.push_back("Thymine (Uracil in RNA)");
    helpData.helpName.push_back("Cytosine");
    helpData.helpName.push_back("Purine (A or G)");
    helpData.helpName.push_back("Pyrimidine (C or T or U)");
    helpData.helpName.push_back("Amino (A or C)");
    helpData.helpName.push_back("Ketone (G or T)");
    helpData.helpName.push_back("Strong interaction (C or G)");
    helpData.helpName.push_back("Weak interaction (A or T)");
    helpData.helpName.push_back("Not-G (A or C or T) H follows G in the alphabet");
    helpData.helpName.push_back("Not-A (C or G or T) B follows A in the alphabet");
    helpData.helpName.push_back("Not-T (not-U) (A or C or G) V follows U in the alphabet");
    helpData.helpName.push_back("Not-C (A or G or T) D follows C in the alphabet");
    helpData.helpName.push_back("Any (A or C or G or T)");
}

// [Fix #2] Now only prints the help menu; data is initialized once in initHelpMenu().
void printHelpMenu() {
    for (size_t i = 0; i < helpData.helpNum.size(); i++) {
        std::cout << helpData.helpNum.at(i)
                  << "\t\t" << helpData.helpName.at(i) << std::endl;
    }
}

// Displays a detailed summary of the currently loaded DNA sequence,
// including region counts and base pair composition statistics.
void summary() {
    std::cout << "\nLoading..." << std::endl;
    analysisResult.basePairs = dna_db.SEQ.length();

    // Count individual nucleotide occurrences
    analysisResult.G = std::count(dna_db.SEQ.begin(), dna_db.SEQ.end(), 'G');
    analysisResult.A = std::count(dna_db.SEQ.begin(), dna_db.SEQ.end(), 'A');
    analysisResult.T = std::count(dna_db.SEQ.begin(), dna_db.SEQ.end(), 'T');
    analysisResult.C = std::count(dna_db.SEQ.begin(), dna_db.SEQ.end(), 'C');
    analysisResult.N = std::count(dna_db.SEQ.begin(), dna_db.SEQ.end(), 'N');
    analysisResult.unknown = analysisResult.basePairs - analysisResult.G
        - analysisResult.C - analysisResult.T - analysisResult.A - analysisResult.N;

    std::cout << "Sequence identifiers:" << std::endl;
    std::cout << "Name:\t" << dna_db.Name << std::endl;
    std::cout << "GID:\t" << dna_db.GID << std::endl;
    std::cout << "REF:\t" << dna_db.REF << std::endl;

    std::cout << "\nRegion characteristics:" << std::endl;
    std::cout << "# regions:\t" << analysisResult.nRegion + analysisResult.cRegion << std::endl;
    std::cout << "# N regions:\t" << analysisResult.nRegion << std::endl;
    std::cout << "# C regions:\t" << analysisResult.cRegion << std::endl;

    // Print IUPAC ambiguity code statistics
    std::cout << "\nBase pair characteristics:" << std::endl;
    std::cout << "# base pairs\t" << dna_db.SEQ.length() << std::endl;
    std::cout << "G:\t" << analysisResult.G << std::endl;
    std::cout << "A:\t" << analysisResult.A << std::endl;
    std::cout << "T:\t" << analysisResult.T << std::endl;
    std::cout << "C:\t" << analysisResult.C << std::endl;
    std::cout << "R:\t" << analysisResult.C + analysisResult.A << std::endl;
    std::cout << "Y:\t" << analysisResult.C + analysisResult.T << std::endl;
    std::cout << "M:\t" << analysisResult.A + analysisResult.C << std::endl;
    std::cout << "K:\t" << analysisResult.G + analysisResult.T << std::endl;
    std::cout << "S:\t" << analysisResult.C + analysisResult.G << std::endl;
    std::cout << "W:\t" << analysisResult.A + analysisResult.T << std::endl;
    std::cout << "H:\t" << analysisResult.A + analysisResult.C + analysisResult.T << std::endl;
    std::cout << "B:\t" << analysisResult.C + analysisResult.G + analysisResult.T << std::endl;
    std::cout << "V:\t" << analysisResult.C + analysisResult.G + analysisResult.A << std::endl;
    std::cout << "D:\t" << analysisResult.A + analysisResult.G + analysisResult.T << std::endl;
    std::cout << "N:\t" << analysisResult.N << std::endl;
    std::cout << "Unknown:\t" << analysisResult.unknown << std::endl;
}

// Scans the DNA sequence to identify gap (N) regions and coded regions.
// Records the start/end indices of each N-region for later queries.
void analyseRegions() {
    int countN = 0, nRegion = 0, countC = 0, cRegion = 0;

    // Check if the sequence starts with an N (gap) character
    if (dna_db.SEQ.empty()) {
        std::cout << "Warning: DNA sequence is empty, skipping region analysis." << std::endl;
        return;
    }

    if (dna_db.SEQ.at(0) == 'N')
        analysisResult.nStart = 1;
    else
        analysisResult.nStart = 0;

    for (size_t i = 0; i < dna_db.SEQ.size(); i++) {
        if (dna_db.SEQ.at(i) == 'N') {
            if (countN == 0) {
                cRegion++;
                analysisResult.NindexStart.push_back(i);
            }
            countN++;
        } else {
            if (countN != 0) {
                // End of a gap region detected — record end index
                analysisResult.NindexEnd.push_back(i);
                countN = 0;
                nRegion++;
            }
            countC++;
        }
    }

    // Handle edge case: sequence ends with N characters
    if (countN != 0) {
        analysisResult.NindexEnd.push_back(dna_db.SEQ.size());
        nRegion++;
    }

    if (analysisResult.cRegion == 0)
        analysisResult.cRegion = cRegion + 1;
    else {
        analysisResult.nRegion = nRegion;
        analysisResult.cRegion = cRegion;
    }
}

// Displays the DNA sequence for a user-specified gap (N) region.
// [Fix #5] Added bounds checking — validates the region number before
// indexing into vectors to prevent std::out_of_range crashes.
void gapRegionQuery() {
    int region, start, end, length;
    std::string sequence;
    std::cout << "Enter gap region number:\n>";
    std::cin >> region;

    // [Fix #5] Validate region number to prevent out-of-range access
    if (region < 1 || region > (int)analysisResult.NindexStart.size()) {
        std::cout << "Invalid region number. Please enter a value between 1 and "
                  << analysisResult.NindexStart.size() << "." << std::endl;
        return;
    }

    if (analysisResult.nStart == 1) {
        start = analysisResult.NindexStart.at(region - 1);
        end = analysisResult.NindexEnd.at(region - 1);
        length = analysisResult.NindexEnd.at(region - 1) - analysisResult.NindexStart.at(region - 1);
    } else {
        // Validate that region index is within bounds for the 'else' branch
        if (region >= (int)analysisResult.NindexStart.size()) {
            std::cout << "Invalid region number for this sequence layout." << std::endl;
            return;
        }
        start = analysisResult.NindexEnd.at(region - 1);
        end = analysisResult.NindexStart.at(region);
        length = analysisResult.NindexStart.at(region) - analysisResult.NindexEnd.at(region - 1);
    }

    sequence = dna_db.SEQ.substr(start, length);

    std::cout << "Selected Sequence:\n"
              << "Base pair range: (" << start << "," << end << ")\n"
              << "Gap region number: " << region << "\n" << std::endl;
    std::cout << "Sequence:\n" << sequence << std::endl;
}

// Displays the DNA sequence for a user-specified coded region.
// [Fix #5] Added bounds checking — same approach as gapRegionQuery().
void codedRegionQuery() {
    int region, start, end, length;
    std::string sequence;
    std::cout << "Enter coded region number:\n>";
    std::cin >> region;

    // [Fix #5] Validate region number to prevent out-of-range access
    if (region < 1 || region > (int)analysisResult.NindexEnd.size()) {
        std::cout << "Invalid region number. Please enter a value between 1 and "
                  << analysisResult.NindexEnd.size() << "." << std::endl;
        return;
    }

    if (analysisResult.nStart == 1) {
        // Validate that region index is within bounds for the 'if' branch
        if (region >= (int)analysisResult.NindexStart.size()) {
            std::cout << "Invalid region number for this sequence layout." << std::endl;
            return;
        }
        start = analysisResult.NindexEnd.at(region - 1);
        end = analysisResult.NindexStart.at(region);
        length = analysisResult.NindexStart.at(region) - analysisResult.NindexEnd.at(region - 1);
    } else {
        start = analysisResult.NindexStart.at(region - 1);
        end = analysisResult.NindexEnd.at(region - 1);
        length = analysisResult.NindexEnd.at(region - 1) - analysisResult.NindexStart.at(region - 1);
    }

    sequence = dna_db.SEQ.substr(start, length);

    std::cout << "Selected Sequence:\n"
              << "Base pair range: (" << start << "," << end << ")\n"
              << "Coded region number: " << region << "\n" << std::endl;
    std::cout << "Sequence:\n" << sequence << std::endl;
}

// [Fix #6] Rewrote findManual() to use std::string::find() instead of
// a hand-rolled nested loop. The previous implementation modified the outer
// loop variable 'i' inside the inner loop, which caused skipped characters
// and potential out-of-bounds access.
void findManual() {
    std::string findN;
    long found = 0;
    std::cout << "\nSpecify the DNA sequence nucleotides you would like to find:\n>";
    std::cin >> findN;

    if (findN.empty()) {
        std::cout << "Empty search string." << std::endl;
        return;
    }

    // Use std::string::find() for correct, safe substring searching
    std::size_t pos = dna_db.SEQ.find(findN, 0);
    while (pos != std::string::npos) {
        found++;
        std::cout << "Base pair range(" << pos << "," << pos + findN.size() << ")" << std::endl;
        pos = dna_db.SEQ.find(findN, pos + 1); // Allow overlapping matches
    }

    std::cout << "Total Number of matches found: " << found << std::endl;
}

// [Fix #3] Replaced unbounded recursion with a loop. Previously, if the user
// entered an invalid filename, findFile() called itself recursively, risking
// a stack overflow. Now retries in a while loop.
// [Fix #15] Changed comparison from dna_db.SEQ.length() to std::string::npos.
// std::string::find() returns npos on failure, not the string's length.
void findFile() {
    std::string file, tempDNAseq, seqCompare;

    std::cout << "\n(Warning, do not key in large files as programme would be "
              << "unresponsive trying to print out match)\n"
              << "Specify the DNA sequence file you would like to find:\n>";

    // [Fix #3] Loop instead of recursive call for invalid files
    while (true) {
        std::cin >> file;
        std::ifstream input(file);
        std::cout << "Loading " << file << std::endl;
        if (!input.good()) {
            std::cout << "File not found. Please try again:\n>";
            continue;
        }
        std::cout << file << "..." << std::endl;

        input.ignore(5000, '\n'); // Skip header line
        while (!input.eof()) {
            std::getline(input, tempDNAseq);
            seqCompare.append(tempDNAseq);
        }
        std::cout << "Successful loading of " << file << std::endl;
        break;
    }

    // [Fix #15] Use std::string::npos for proper find() failure check
    std::size_t loc = dna_db.SEQ.find(seqCompare);
    if (loc != std::string::npos) {
        std::cout << "\nBase pair range: (" << loc << "," << loc + seqCompare.size() << ")" << std::endl;
        std::cout << seqCompare;
    } else {
        std::cout << "No Match" << std::endl;
    }
}

// [Fix #16] Added input validation for bpRange(). Previously, entering a
// single number (no comma) would cause vectRange.at(1) to throw, and
// non-numeric input would crash via stoi(). Now validates before accessing.
void bpRange() {
    std::vector<int> vectRange;
    std::string ranges, subst, output;
    int length;
    std::cout << "Enter a comma ',' separating base pair ranges (NO SPACES):\n>";
    std::cin >> ranges;

    std::stringstream range(ranges);
    while (range.good()) {
        std::getline(range, subst, ',');
        try {
            vectRange.push_back(std::stoi(subst));
        } catch (const std::exception& e) {
            std::cout << "Invalid input: '" << subst << "' is not a valid number." << std::endl;
            return;
        }
    }

    // Validate that exactly two range values were provided
    if (vectRange.size() < 2) {
        std::cout << "Error: Please provide two comma-separated numbers (e.g., 100,200)." << std::endl;
        return;
    }

    // Validate that start is less than end
    if (vectRange.at(0) >= vectRange.at(1)) {
        std::cout << "Error: Start position must be less than end position." << std::endl;
        return;
    }

    // Validate that positions are within the sequence bounds
    if (vectRange.at(0) < 0 || vectRange.at(1) > (int)dna_db.SEQ.size()) {
        std::cout << "Error: Range is out of bounds. Sequence length is "
                  << dna_db.SEQ.size() << "." << std::endl;
        return;
    }

    length = vectRange.at(1) - vectRange.at(0);
    output = dna_db.SEQ.substr(vectRange.at(0), length);
    std::cout << "Selected sequence: \nBase pair range: ("
              << vectRange.at(0) << "," << vectRange.at(1) << ")" << std::endl;
    std::cout << "\nSequence:\n" << output << std::endl;
}

// [Fix #8] Main function rewritten to use while loops instead of goto
// statements. The original used goto with labels (menu1:, option1:) which
// created hard-to-follow spaghetti code.
// [Fix #17] All commented-out debug code has been removed.
int main() {
    std::cout << "DNA Sequence Database Software" << std::endl;
    std::cout << "Specify the name of DNA sequence file names you would like to load. "
              << "For multiple files, add a ',' between each file name. "
              << "(Add .fa extension after name, eg chr1.fa and no space in between files)\n>";

    std::cin >> initial.File;

    organiseFile(initial.File);
    initMenuOptions();
    initHelpMenu(); // [Fix #2] Initialize help menu data once

    // [Fix #8] Outer menu loop replaces 'goto menu1'
    bool running = true;
    while (running) {
        printOrganiseFile();
        std::cin >> initial.initialSelect;

        if (initial.initialSelect == "Q") {
            std::cout << "Programme ended." << std::endl;
            return 0;
        } else if (initial.initialSelect == "S") {
            firstSummary();
            continue; // Return to file selection menu
        } else {
            selectFile(initial.initialSelect);
            if (initial.check == 1) {
                return 0;
            }
            assignClass();
            if (initial.check == 1) {
                return 0; // [Fix #14] Check for errors from assignClass too
            }
        }

        analyseRegions();

        // [Fix #8] Inner options loop replaces 'goto option1'
        bool inOptionsMenu = true;
        while (inOptionsMenu) {
            printMenuOptions();
            std::cin >> menuOption.optionSelect;

            if (menuOption.optionSelect == "Q") {
                std::cout << "Programme ended." << std::endl;
                return 0;
            } else if (menuOption.optionSelect == "H") {
                printHelpMenu(); // [Fix #2] Just prints, no longer re-initializes
            } else if (menuOption.optionSelect == "S") {
                summary();
            } else if (menuOption.optionSelect == "R") {
                inOptionsMenu = false; // Break to outer loop (file selection)
            } else if (menuOption.optionSelect == "1") {
                gapRegionQuery();
            } else if (menuOption.optionSelect == "2") {
                codedRegionQuery();
            } else if (menuOption.optionSelect == "3") {
                bpRange();
            } else if (menuOption.optionSelect == "4") {
                findManual();
            } else if (menuOption.optionSelect == "5") {
                findFile();
            } else {
                std::cout << "Command not recognised, try again" << std::endl;
            }
        }
    }

    return 0;
}