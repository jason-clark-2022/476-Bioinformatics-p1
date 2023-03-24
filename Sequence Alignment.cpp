// Jason Clark 
// CS476 Project 1: Sequence Alignment

#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <string>
#include <string.h>
#include <sys/stat.h>

using namespace std;

struct submatrix
{
    string description;
    vector<string> aaSequence;
    vector<vector<double>> scores;
};

struct userAASequence
{
    string description;
    vector<char> sequence;
};

struct processedAlignments
{
    vector<string> seq1;
    vector<string> seq2;
    double score;
    vector<vector<double>> opt;
};

struct coordinates
{
    int xcord;
    int ycord;
};

submatrix readSubMatrix(string inFile);
userAASequence readAASequence(string inFile);
processedAlignments globalAlignment(userAASequence seqX, userAASequence seqY, double gap, submatrix subs);
processedAlignments semiGlobalAlignment(userAASequence seqX, userAASequence seqY, double gap, submatrix subs);
processedAlignments localAlignment(userAASequence seqX, userAASequence seqY, double gap, submatrix subs);
void printVector(vector<auto> data);
void printMatrix(vector<vector<auto>> data);
void printOPT(vector<vector<auto>> matrix, userAASequence seq1, userAASequence seq2);
double findSubScore(char val1,  char val2, submatrix sub);
coordinates findSemiGlobalScore(processedAlignments alignments);
coordinates findLocalScore(processedAlignments alignment);
string submatrixFromInput(string input);
string sequenceFromInput(string input);
bool checkValidFile(string fileName);

int main()
{
    bool exit = false;
    bool firstLoop = false;
    bool validInput = false;
    bool global = false;
    bool semiglobal = false;
    bool local = false;

    string userinput;
    string submatrixFile;
    string seq1File;
    string seq2File;
    double gapPenalty;

    while(!exit)
    {
        if(!firstLoop)
        {   cout << endl << "***************************************************************" << endl;
            cout << "Welcome to the sequence alignment tool!" << endl;
            cout << "Any field that requests a filename will accept" << endl;
            cout << "the filename shortcuts or a valid filename as input" << endl;
            cout << "enter 'e' to exit at any time..." << endl;
            cout << "***************************************************************" << endl;
            firstLoop = true;
        }  
        
        cout << "Enter the substitution matrix filename (an/bs/hp/pm/tst)" << endl << ">";
        validInput = false;
        while(!validInput)
        {
            cin >> userinput;
            submatrixFile = submatrixFromInput(userinput);
            
            if(userinput == "e")
            {
                exit = true;
                break;
            }
            else if(checkValidFile(submatrixFile))
            {
                validInput = true;
            }
            else
            {
                cout << "Invalid filename" << endl << ">";
            }
        }
        if(exit)
            break;
    

        cout << "Enter the first sequence filename (a1/a2/b1/b2/c1/c2/cus1 - cus6)" << endl << ">";
        validInput = false;
        while(!validInput)
        {
            cin >> userinput;
            seq1File = sequenceFromInput(userinput);

            if(userinput == "e")
            {
                exit = true;
                break;
            }
            if(checkValidFile(seq1File))
            {
                validInput = true;
            }
            else
            {
                cout << "Invalid filename" << endl << ">";
            }
        }
        if(exit)
            break;


        cout << "Enter the second sequence filename (a1/a2/b1/b2/c1/c2/cus1 - cus6)" << endl << ">";
        validInput = false;
        while(!validInput)
        {
            cin >> userinput;
            seq2File = sequenceFromInput(userinput);
            if(userinput == "e")
            {
                exit = true;
                break;
            }
            if(checkValidFile(seq2File))
            {
                validInput = true;
            }
            else
            {
                cout << "Invalid filename" << endl << ">";
            }
        }
        if(exit)
            break;

        cout << "Enter the alignment type (gbl/sgb/loc)" << endl << ">";
        validInput = false;
        while(!validInput)
        {
            cin >> userinput;
            if(userinput == "e")
            {
                exit = true;
                break;
            }
            else if(userinput == "gbl")
            {
                global = true;
                validInput = true;
            }
            else if(userinput == "sgb")
            {
                semiglobal = true;
                validInput = true;
            }
            else if(userinput == "loc")
            {
                local = true;
                validInput = true;
            }
            else
            {
                cout << "invalid input" << endl << ">";
            }
        }
        if(exit)
            break;

        cout << "Enter the gap penalty (any number)" << endl << ">";
        validInput = false;
        while(!validInput)
        {
            cin >> gapPenalty;
            if(!cin)
            {
                cout << "Invalid input" << endl << ">";
            }
            else
                validInput = true;
        }
        if(exit)
            break;
        

        // THIS IS WHERE THE FUN BEGINS..........................................................................................................
        
        cout << "###############################################################" << endl;
        
        submatrix subscores = readSubMatrix(submatrixFile);
        userAASequence seq1 = readAASequence(seq1File);
        userAASequence seq2 = readAASequence(seq2File);
        processedAlignments alignment;

        if(global)
        {
            alignment = globalAlignment(seq1, seq2, gapPenalty, subscores);
            cout << ">Global Alignment" << endl;
            global = false;
        }
        else if(semiglobal)
        {
            alignment = semiGlobalAlignment(seq1, seq2, gapPenalty, subscores);
            cout << ">Semiglobal Alignment" << endl;
            semiglobal = false;
        }
        else if(local)
        {
            alignment = localAlignment(seq1, seq2, gapPenalty, subscores);
            cout << ">Local Alignment" << endl;
            local = false;
        }
        cout << ">" << submatrixFile << endl;
        cout << ">" << seq1File << endl;
        cout << ">" << seq2File << endl << endl;
        cout << ">Aligned Sequences" << endl;


        printVector(alignment.seq1);
        printVector(alignment.seq2);
        cout << "Score: " << alignment.score << endl << endl;
        printOPT(alignment.opt, seq1, seq2);
        cout << "###############################################################" << endl;
    
        cout << "Would you like to run a new sequence (y/n)" << endl << ">";
        
        validInput = false;
        while(!validInput)
        {
            cin >> userinput;
            if(userinput == "n")
            {
                validInput = true;
                exit = true;
            }
            else if(userinput == "y")
            {
                validInput = true;
            }
            else
            {
                cout << "invalid input (y/n)" << endl << ">";
            }
        }
    }

    cout << "Goodbye..." << endl;
    return 0;


    // Shortcuts
    //(an/bs/hp/pm/tst)
    //(a1/a2/b1/b2/c1/c2/cus1 - cus6)

    //submatrix subscorestest = readSubMatrix(submatrixFromInput("bs"));
    submatrix subscorestest = readSubMatrix("P1SubMatrices/PAM500scores.txt");
    userAASequence seq1test = readAASequence(sequenceFromInput("a1"));
    userAASequence seq2test = readAASequence(sequenceFromInput("b2"));


    processedAlignments globaltest = semiGlobalAlignment(seq1test, seq2test, -5, subscorestest);
    //processedAlignments globaltest = globalAlignment(seq1test, seq2test, -5, subscorestest);
    //processedAlignments globaltest = localAlignment(seq1test, seq2test, -5, subscorestest);
    
    printVector(globaltest.seq1);
    printVector(globaltest.seq2);
    cout << "score: " << globaltest.score << endl << endl;
    //printOPT(globaltest.opt, seq1test, seq2test);
    
    cout << endl << "DONE" << endl;
}

submatrix readSubMatrix(string inFile)
{
    submatrix submatrixValues;
    string line;
    string delim = ",";
    string token;
    size_t pos = 0;
    bool readDescription = false;
    bool readAASequence = false;
    int scoresLineCounter = 0;
    int scoresIndexCounter = 0;
    string::size_type sz;

    ifstream file;
    file.open(inFile);
    if(file.is_open())
    {
        while(getline(file, line))
        {
            if(!readDescription) //if description has not been read
            {
                submatrixValues.description = line;
                readDescription = true;
            }
            else if(!readAASequence) //if AASequence has not been read
            {
                while((pos = line.find(delim)) != string::npos)
                {
                    token = line.substr(0,pos);
                    submatrixValues.aaSequence.push_back(token);
                    line.erase(0, pos + delim.length());
                }
                submatrixValues.aaSequence.push_back(line);
                readAASequence = true;
            }
            else //Read Matrix Scores
            {
                scoresIndexCounter = 0;
                while((pos = line.find(delim)) != string::npos)
                {
                    
                    token = line.substr(0,pos);
                    if(scoresLineCounter == 0)                                
                    {
                        vector<double> score;
                        score.push_back(stod(token, &sz));
                        submatrixValues.scores.push_back(score);
                    }
                    else
                    {
                        submatrixValues.scores[scoresIndexCounter].push_back(stod(token, &sz));
                    }
                    
                    line.erase(0, pos + delim.length());
                    scoresIndexCounter++;
                }
                if(scoresLineCounter == 0)
                {
                        vector<double> score;
                        score.push_back(stod(line, &sz));
                        submatrixValues.scores.push_back(score);
                }
                else
                {
                        submatrixValues.scores[scoresIndexCounter].push_back(stod(line, &sz));
                }
                scoresLineCounter++;
            }
        }

        file.close();
    }
    else
    {
        cerr << "Error, unable to open file: " << inFile << endl;
    }
    return submatrixValues;
}

    else 
    {
        cerr << "Error, unable to open file: " << inFile << endl;
    }
    return sequence;
}

processedAlignments globalAlignment(userAASequence seqX, userAASequence seqY, double gap, submatrix subs)
{
    processedAlignments alignments;

    int xlen = seqX.sequence.size() + 1;
    int ylen = seqY.sequence.size() + 1;
    
    vector<double> initialize(ylen, 0);
    vector<vector<double>> alignmentTable(xlen, initialize);
    vector<vector<double>> seqdirect = alignmentTable;
    
    for(int i = 1; i < ylen; i++)
    {
        alignmentTable[0][i] = alignmentTable[0][i-1] + gap;
    }
    for(int i = 1; i < xlen; i++)
    {
        alignmentTable[i][0] = alignmentTable[i-1][0] + gap;
    }

    
    for(int i = 1; i < xlen; i++)
    {
        for(int j = 1; j < ylen; j++)
        {
            double subscore = findSubScore(seqX.sequence[i-1], seqY.sequence[j-1], subs);
            subscore = subscore+alignmentTable[i-1][j-1];
            double vertgap = gap + alignmentTable[i-1][j];
            double horzgap = gap + alignmentTable[i][j-1];

            if(subscore >= vertgap && subscore >= horzgap) 
            {
                alignmentTable[i][j] = subscore;
                seqdirect[i][j] = 1;    // set seqdirect to 1 to show match
            }
            else if(vertgap >= horzgap)
            {
                alignmentTable[i][j] = vertgap;
                seqdirect[i][j] = 0;    // set seqdirect to 0 to show vertical gap
            }
            else
            {
                alignmentTable[i][j] = horzgap;
                seqdirect[i][j] = -1;   // set seqdirect to -1 to show horizontal gap
            }
        }
    }
    alignments.opt = alignmentTable;
    alignments.score = alignmentTable[xlen-1][ylen-1];
    

    xlen = seqX.sequence.size();
    ylen = seqY.sequence.size();

    // Backtrace using sequence direct table to set both sequence alignments
    while(xlen > 0 && ylen > 0)
    {   
        if(seqdirect[xlen][ylen] == 1)
        {
            string v1(1, seqX.sequence[xlen-1]);
            string v2(1, seqY.sequence[ylen-1]);
            alignments.seq1.push_back(v1);
            alignments.seq2.push_back(v2);
            xlen--;
            ylen--;
        }
        else if(seqdirect[xlen][ylen] == 0)
        {
            string v1(1, seqX.sequence[xlen-1]);
            alignments.seq1.push_back(v1);
            alignments.seq2.push_back("_");
            xlen--;
        }
        else
        {
            string v1(1, seqY.sequence[ylen-1]);
            alignments.seq1.push_back("_");
            alignments.seq2.push_back(v1);
            ylen--;
        }
    }
    
    while(xlen > 0)
    {
        string v1(1, seqX.sequence[xlen-1]);
        alignments.seq1.push_back(v1);
        xlen--;
    }
    while(ylen > 0)
    {
        string v1(1, seqY.sequence[ylen-1]);
        alignments.seq2.push_back(v1);
        ylen--;
    }
    if(alignments.seq2.size() > alignments.seq1.size())
    {
        alignments.seq1.push_back("_");
    }
    else if(alignments.seq2.size() > alignments.seq1.size())
    {
        alignments.seq2.push_back("_");
    }

    reverse(alignments.seq1.begin(), alignments.seq1.end());
    reverse(alignments.seq2.begin(), alignments.seq2.end());
    return alignments;
}


processedAlignments semiGlobalAlignment(userAASequence seqX, userAASequence seqY, double gap, submatrix subs)
{
    processedAlignments alignments;

    int xlen = seqX.sequence.size() + 1;
    int ylen = seqY.sequence.size() + 1;
    
    vector<double> initialize(ylen, 0);
    vector<vector<double>> alignmentTable(xlen, initialize);
    vector<vector<double>> seqdirect = alignmentTable;


    for(int i = 1; i < xlen; i++)
    {
        for(int j = 1; j < ylen; j++)
        {
            double subscore = findSubScore(seqX.sequence[i-1], seqY.sequence[j-1], subs);
            subscore = subscore+alignmentTable[i-1][j-1];
            double vertgap = gap + alignmentTable[i-1][j];
            double horzgap = gap + alignmentTable[i][j-1];

            if(subscore >= vertgap && subscore >= horzgap) 
            {
                alignmentTable[i][j] = subscore;
                seqdirect[i][j] = 1;    // set seqdirect to 1 to show match
            }
            else if(vertgap >= horzgap)
            {
                alignmentTable[i][j] = vertgap;
                seqdirect[i][j] = 0;    // set seqdirect to 0 to show vertical gap
            }
            else
            {
                alignmentTable[i][j] = horzgap;
                seqdirect[i][j] = -1;   // set seqdirect to -1 to show horizontal gap
            }
        }
    }

    alignments.opt = alignmentTable;

    coordinates coord = findSemiGlobalScore(alignments);
    alignments.score = alignments.opt[coord.xcord][coord.ycord];


    xlen = coord.xcord;
    ylen = coord.ycord;


    // Backtrace using sequence direct table to set both sequence alignments
    while(xlen > 0 && ylen > 0)   
    {   
        if(seqdirect[xlen][ylen] == 1)
        {
            string v1(1, seqX.sequence[xlen-1]);
            string v2(1, seqY.sequence[ylen-1]);
            alignments.seq1.push_back(v1);
            alignments.seq2.push_back(v2);
            xlen--;
            ylen--;
        }
        else if(seqdirect[xlen][ylen] == 0)
        {
            string v1(1, seqX.sequence[xlen-1]);
            alignments.seq1.push_back(v1);
            alignments.seq2.push_back("_");
            xlen--;
        }
        else
        {
            string v1(1, seqY.sequence[ylen-1]);
            alignments.seq1.push_back("_");
            alignments.seq2.push_back(v1);
            ylen--;
        }
    }
    reverse(alignments.seq1.begin(), alignments.seq1.end());
    reverse(alignments.seq2.begin(), alignments.seq2.end());
    return alignments;  
}

processedAlignments localAlignment(userAASequence seqX, userAASequence seqY, double gap, submatrix subs)
{
    processedAlignments alignments;

    int xlen = seqX.sequence.size() + 1;
    int ylen = seqY.sequence.size() + 1;
    
    vector<double> initialize(ylen, 0);
    vector<vector<double>> alignmentTable(xlen, initialize);
    vector<vector<double>> seqdirect = alignmentTable;
    
    for(int i = 1; i < xlen; i++)
    {
        for(int j = 1; j < ylen; j++)
        {
            double subscore = findSubScore(seqX.sequence[i-1], seqY.sequence[j-1], subs);
            subscore = subscore+alignmentTable[i-1][j-1];
            double vertgap = gap + alignmentTable[i-1][j];
            double horzgap = gap + alignmentTable[i][j-1];

            if(subscore >= vertgap && subscore >= horzgap) 
            {
                if(subscore < 0)
                    subscore = 0;
                alignmentTable[i][j] = subscore;
                seqdirect[i][j] = 1;    // set seqdirect to 1 to show match
            }
            else if(vertgap >= horzgap)
            {
                if(vertgap < 0)
                    vertgap = 0;
                alignmentTable[i][j] = vertgap;
                seqdirect[i][j] = 0;    // set seqdirect to 0 to show vertical gap
            }
            else
            {
                if(horzgap < 0)
                    horzgap = 0;
                alignmentTable[i][j] = horzgap;
                seqdirect[i][j] = -1;   // set seqdirect to -1 to show horizontal gap
            }
        }
    }
    alignments.opt = alignmentTable;

    coordinates coord = findLocalScore(alignments);
    alignments.score = alignments.opt[coord.xcord][coord.ycord];


    xlen = coord.xcord;
    ylen = coord.ycord;

    // Backtrace using sequence direct table to set both sequence alignments
    while(alignments.opt[xlen][ylen] > 0)
    {   
        if(seqdirect[xlen][ylen] == 1)
        {
            string v1(1, seqX.sequence[xlen-1]);
            string v2(1, seqY.sequence[ylen-1]);
            alignments.seq1.push_back(v1);
            alignments.seq2.push_back(v2);
            xlen--;
            ylen--;
        }
        else if(seqdirect[xlen][ylen] == 0)
        {
            string v1(1, seqX.sequence[xlen-1]);
            alignments.seq1.push_back(v1);
            alignments.seq2.push_back("_");
            xlen--;
        }
        else
        {
            string v1(1, seqY.sequence[ylen-1]);
            alignments.seq1.push_back("_");
            alignments.seq2.push_back(v1);
            ylen--;
        }
    }
    
    reverse(alignments.seq1.begin(), alignments.seq1.end());
    reverse(alignments.seq2.begin(), alignments.seq2.end());
    return alignments;
}

void printVector(vector<auto> data)
{
    for(int i=0; i<data.size(); i++)
    {
        cout << data[i] << " ";
    }
    cout << endl;
}

void printMatrix(vector<vector<auto>> data)
{
    for(int i = 0; i < data.size(); i++)
    {
        for(int j = 0; j < data[i].size(); j++)
        {
            cout << data[i][j] << " ";
        }
        cout << endl;
    }
}

void printOPT(vector<vector<auto>> matrix, userAASequence seq1, userAASequence seq2)
{
    cout << "     ";
    for(int i = 0; i < seq2.sequence.size(); i ++)
    {
        cout << seq2.sequence[i] << " ";
    }
    cout << endl;

    for(int i = 0; i < matrix.size(); i++)
    {
        if(i == 0)
            cout << "   ";
        else
            cout << seq1.sequence[i-1] << ": ";
        for(int j = 0; j < matrix[i].size(); j++)
        {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

double findSubScore(char val1, char val2, submatrix sub)
{
    string v1(1, val1);
    string v2(1, val2);
    int i1;
    int i2;
    for(int i = 0; i < sub.aaSequence.size(); i++)
    {
        if(v1 == sub.aaSequence[i])
        {
            i1 = i;
        }
        if(v2 == sub.aaSequence[i])
        {
            i2 = i;
        }
    }
    return sub.scores[i1][i2];
}

coordinates findSemiGlobalScore(processedAlignments alignments)
{
    coordinates retCoords;
    retCoords.xcord = 0;
    retCoords.ycord = 0;
    int backIndex;
    double highScore = 0;

    backIndex = alignments.opt[0].size()-1;
    
    for(int i = 0; i < alignments.opt.size(); i++)
    {
        if(alignments.opt[i][backIndex] > highScore)
        {
            highScore = alignments.opt[i][backIndex];
                retCoords.xcord = i;
                retCoords.ycord = backIndex;
        }
    }

    backIndex = alignments.opt.size()-1;
    for(int i = 0; i < alignments.opt[0].size(); i++)
    {
        if(alignments.opt[backIndex][i] > highScore)
        {
            highScore = alignments.opt[backIndex][i];
            retCoords.xcord = backIndex;
            retCoords.ycord = i;
        }
    }
    return retCoords;
}

coordinates findLocalScore(processedAlignments alignment)
{
    coordinates retCoords;
    retCoords.xcord = 0;
    retCoords.ycord = 0;
    double highScore;
    for(int i = 0; i < alignment.opt.size(); i++)
    {
        for(int j = 0; j < alignment.opt[0].size(); j++)
        {
            if(highScore < alignment.opt[i][j])
            {
                highScore = alignment.opt[i][j];
                retCoords.xcord = i;
                retCoords.ycord = j;
            }
        }
    }
    return retCoords;
}

string submatrixFromInput(string input)
{
    string filename;
    string AAnucleopp = "P1SubMatrices/AAnucleoPP.txt";
    string blosum62 = "P1SubMatrices/BLOSUM62.txt";
    string hp = "P1SubMatrices/HP.txt";
    string pam250 = "P1SubMatrices/PAM250-scores.txt";
    string testdoc = "P1SubMatrices/test.txt";

    if(input == "an"){filename = AAnucleopp;}
    else if(input == "bs"){filename = blosum62;}
    else if(input == "hp"){filename = hp;}
    else if(input == "pm"){filename = pam250;}
    else if(input == "tst"){filename = testdoc;}
    else {filename = input;}
    return filename;
}

string sequenceFromInput(string input)
{
    string filename;
    string a1 = "P1AASeqs/sequenceA1.txt";
    string a2 = "P1AASeqs/sequenceA2.txt";
    string b1 = "P1AASeqs/sequenceB1.txt";
    string b2 = "P1AASeqs/sequenceB2.txt";
    string c1 = "P1AASeqs/sequenceC1.txt";
    string c2 = "P1AASeqs/sequenceC2.txt";
    string cus1 = "P1AASeqs/custom1.txt";
    string cus2 = "P1AASeqs/custom2.txt";
    string cus3 = "P1AASeqs/custom3.txt";
    string cus4 = "P1AASeqs/custom4.txt";
    string cus5 = "P1AASeqs/custom5.txt";
    string cus6 = "P1AASeqs/custom6.txt";

    if(input == "a1"){filename = a1;}
    else if(input == "a2"){filename = a2;}
    else if(input == "b1"){filename = b1;}
    else if(input == "b2"){filename = b2;}
    else if(input == "c1"){filename = c1;}
    else if(input == "c2"){filename = c2;}
    else if(input == "cus1"){filename = cus1;}
    else if(input == "cus2"){filename = cus2;}
    else if(input == "cus3"){filename = cus3;}
    else if(input == "cus4"){filename = cus4;}
    else if(input == "cus5"){filename = cus5;}
    else if(input == "cus6"){filename = cus6;}
    else {filename = input;}

    return filename;
}

bool checkValidFile(string fileName)
{
    ifstream inFile(fileName);
    return inFile.good();
}