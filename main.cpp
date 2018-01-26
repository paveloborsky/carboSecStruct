/* 
 * File:   main.cpp
 * Author: poborsky
 *
 * Created on March 11, 2014, 1:31 PM
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <assert.h>
#include <cstdlib>//stdlib.h>
#include <vector>
#include <utility>
#include <cstdio>
#include <unistd.h> //sleep
#include <sys/types.h> //to create directory
#include <sys/stat.h> //to create directory


using namespace std;

//
// CONSTANTS
//
const double PI =  3.14159265359;
const double kb =  1.380662e-26; // kJ.K^(-1)
const double T  =  298.15; //K
const double NA =  6.022045e23; // mol^(-1)
const double kbT=  kb*T; // kJ
const double NAkbT = kbT*NA; // (~2.47) kJ.mol^(-1)

//
// NAMES OF FILES
//
string cnfInPath; 
string cnfOutPath;
string angleInPath; 
string angleOutPath;
string tplPath;
string changePath;
string trcPath;
string outPath;

//
// STRUCTURE DEFINITION
//
struct vec {
    vec() : x(0.0), y(0.0), z(0.0) {}
    vec(double x, double y, double z) : x(x), y(y), z(z) {}
    double x;
    double y;
    double z;
};

struct Torsion_param {
    Torsion_param() : K(0.0), s(0.0), m(0.0) {}
    Torsion_param(double K, double s, double m) : K(K), s(s), m(m) {}
    double K;
    double s;
    double m;
};

struct LJparam {
    LJparam() : C12(0.0), C6(0.0), CS12(0.0), CS6(0.0) {}
    LJparam(double C12, double C6, double CS12, double CS6) : C12(C12), C6(C6), CS12(CS12), CS6(CS6) {}
    double C12;
    double C6;
    double CS12;
    double CS6;
};

struct atom {
    atom() : nrRes(0), res(" "), a(" "), nrA(0) {}
    int nrRes;
    string res;
    string a;
    int nrA;
    vec pos;
};

struct helix{
    helix() {}
    double phi;
    double psi;
    double omega;
    double R; //radius
    double p; //pitch
    double a; //rise per residue
    double g; //turn angle per residue
    double n; //residues per turn
    double maxH; //max posible Hbonds at same time
    double SteElecE; //Stereoelectronic energy
    double SteElecE_extended; //Stereoelectronic energy
    double SteElecE_just_TOR; //Stereoelectronic energy
    double SteElecE_just_LJ; //Stereoelectronic energy
    double SteElecE_just_LJ_extended; //Stereoelectronic energy
    vector<pair<string, int> > Hbond; //total number of posible Hbonds
    vector<pair<string, int> > HbondHelix; //total number of posible Hbonds
    vector<pair<string, double > > clash;
    vector<pair<string, double > > clashHelix;
};



//
// TORSIONAL AND LJ PARAMETERS OF HEXOPYRANOSES
//
//make torsion parameters
//56a6@CARBO-R
Torsion_param t34(5.920,      0.0,       3); //-CHn,SI-CHn-
Torsion_param t42(4.900,      0.0,       3); //generic CO
Torsion_param t44(7.000,    180.0,       1); //exo-anomeric
Torsion_param t46(2.500,     60.0,       1); //exocyclic methoxyl
Torsion_param t47(2.500,    -60.0,       1); //exocyclic methoxyl
Torsion_param t50(1.000,     60.0,       1); //exocyclic oxymethyl
Torsion_param t51(1.000,    -60.0,       1); //exocyclic oxymethyl
Torsion_param t52(4.500,    180.0,       1); //oxygen-oxygen gauche
    
map<int,Torsion_param> Torsions ={{34, t34},{42, t42}, {44, t44}, {46, t46}, 
{47, t47},{50, t50}, {51, t51}, {52, t52}};

//make map of C12, C6, CS12, CS6
//56a6@CARBO-R
//name type                 IAC
//OA   alcohol oxygen       3
//OE   ether oxygen         4
//CH2  CH2 - (C6)           15
//Or   ring oxygen          54
//Cr   ring carbon          56
//            C12            C6          CS12           CS6
LJparam   OO(1.505529e-06,  2.261954e-03,  1.265625e-06,  2.261954e-03); //OA-OA
LJparam  OOl(1.505529e-06,  2.261954e-03,  1.265625e-06,  2.261954e-03); //OA-OE
LJparam OlOl(1.210000e-06,  2.261954e-03,  1.265625e-06,  2.261954e-03); //OE-OE NOT USED
LJparam  OCe(6.410800e-06,  4.110135e-03,  2.449796e-06,  3.268799e-03); //OA-CH2
LJparam OlCe(6.410800e-06,  4.110135e-03,  2.449796e-06,  3.268799e-03); //OE-CH2
LJparam CeCe(3.396558e-05,  7.468416e-03,  4.741926e-06,  4.723813e-03); //CH2-CH2
LJparam  OOr(1.505529e-06,  2.261954e-03,  7.706250e-07,  2.261954e-03); //OA-Or
LJparam OlOr(1.210000e-06,  2.261954e-03,  7.706250e-07,  2.261954e-03); //OE-Or
LJparam CeOr(6.410800e-06,  4.110135e-03,  1.491654e-06,  3.268799e-03); //CH2-Or
LJparam OrOr(1.210000e-06,  2.261954e-03,  4.692250e-07,  2.261954e-03); //Or-Or
LJparam  OCr(1.083500e-05,  3.704924e-03,  1.530000e-06,  3.268799e-03); //OA-Cr
LJparam OlCr(1.083500e-05,  3.704924e-03,  1.530000e-06,  3.268799e-03); //OE-Cr
LJparam CeCr(5.740580e-05,  6.732118e-03,  2.961531e-06,  4.723813e-03); //CH2-Cr
LJparam OrCr(1.083500e-05,  3.704924e-03,  9.316000e-07,  3.268799e-03); //Or-Cr
LJparam CrCr(9.702250e-05,  6.068410e-03,  1.849600e-06,  4.723813e-03); //Cr-Cr

map<string,LJparam> LJpairs = {{"OO"  , OO   },
            {"OOl" , OOl  },
            {"OlOl", OlOl },
            {"OCe" , OCe  },
            {"OlCe", OlCe },
            {"CeCe", CeCe },
            {"OOr" , OOr  },
            {"OlOr", OlOr },
            {"CeOr", CeOr },
            {"OrOr", OrOr },
            {"OCr" , OCr  },
            {"OlCr", OlCr },
            {"CeCr", CeCr },
            {"OrCr", OrCr },
            {"CrCr", CrCr }};

//
// HBOND PARAMTERS
// 
const double Hdistmin = 0.0;
const double Hdistmax = 0.25;
const double Hangle = 135;
// LIST OF HYDROGENS AND H-ACCEPTORS
vector<string > Hlist {"HO1","HO2","HO3","HO4","HO6"};
vector<string > acceptList {"O1","O2","O3","O4","O5","O6"};

    
// bounds of Ramachandarn map
const float Rmin = -180;
const float Rmax =  180;



//
//DECLARATION OF ALL FUNCTIONS//
//


//basic geometric functions
vec getSubstrVec(const vec &a, const vec &b);
vec getSumVec(const vec &a, const vec &b);
double getAbs2Vec(const vec &v);
double getAbsVec(const vec &v);
double getDotofVec(const vec &a, const vec &b);
vec getCrossofVec(const vec &a, const vec &b);
vec getMultipVec(const vec &a, double x);
double calcDistance(const vec &a, const vec &b);
double calcAngle(const vec &a, const vec &b, const vec &c);
double calcDihedral(const vec &i, const vec &j, const vec &k, const vec &l);

//calculation of cartesians coordinate
double calcZ(const atom q,const atom c,const atom b,const atom a);
double calcY(const atom q,const atom c,const atom b,const atom a);
double calcX(const atom q,const atom c,const atom b,const atom a);

// functions for transformation
void translate(atom &a, atom &b, atom &c);
void rotate1(atom &a, atom &b, atom &c, bool inverse);
void rotate2(atom &a, atom &b, atom &c, bool inverse);
void transform(atom &a, atom &b, atom &c);
vec getMatMult(double rotM[3][3],const vec v);
void getMatrix1(atom const a, atom const b, atom const c, bool inverse, double rotM1[3][3]);
void getMatrix2(atom const a, atom const b, atom const c, bool inverse, double rotM2[3][3]);
void transformBack(atom const A, atom const B, atom const C, atom &d);

//functions just for testing used in -db option
void translatetest(atom &a, atom &b, atom &c, atom &d);
void rotate1test(atom &a, atom &b, atom &c,atom &d, bool inverse);
void rotate2test(atom &a, atom &b, atom &c,atom &d, bool inverse);
void transformtest(atom &a, atom &b, atom &c, atom &d);

//printing functions
void printInp(atom inp[17]);
void printTpl(vector<vector<string> > tpl);
void printCnf(vector<atom> cnf);
void printVec(vec v);
void printVec(vector<vector<string> > v);
void printVec(vector<vector<int> > v);
void printVec(vector<pair<string, int> > v);
void printVec(vector<vector<helix> > v, string par);
void printVec(vector<atom> v);
void printVec(vector<vector<atom> > v);
void printVec(vector<int> v);
void printVec(vector<string> v);
void printAtom(atom a);

//this functions are currently not used
void getDistance(vector<vector<string> > tpl, vector<atom> cnf);
void getAngle(vector<vector<string> > tpl, vector<atom> cnf);
void getDihedral(vector<vector<string> > tpl, vector<atom> cnf);
double chop(double d);
float round(float f,float prec);

//functions used in main//
//reading functions
void readInp(vector<atom > &angleIn);
void readChange(vector<vector<string> > &change);
void readTpl(vector<vector<string> > &tpl);
void readCnf(vector<atom> &cnf);
void readTrc(vector<atom> &cnf, vector<vector<atom> > &trc);

int changeAngle(vector<atom> &angle, vector<vector<string> > tpl, string flag, string instructions);
string merge(string str, char delim);
vector<string> split(const string &s, char delim);
vector<string> &split(const string &s, char delim, vector<string> &elems);
void setPsi(vector<atom> &angle, vector<vector<string> > tpl, double value, vector<pair <int,int> > con);
void setPhi(vector<atom> &angle, vector<vector<string> > tpl, double value, vector<pair <int,int> > con);

//creating cnf/angular matrix
int fillCnf(vector<atom> angleIn, vector<vector<string> > tpl, vector<atom> &cnfneu);
void fillCnf(vector<atom> angleIn, vector<vector<string> > tpl, vector<atom> &cnfneu, double phi, double psi);
void fillTpl(vector<vector<string> > tpl, vector<atom> cnf, vector<atom > &angle);
//helping function of creating
atom findAtom(string i, vector<atom> cnf, int Nres);
int findAtomPos(string name, vector<atom> cnf, int res);
void makeChain(int i, vector<vector<string> > tpl, atom &c, atom &b, atom &a, vector<atom> cnf, int res);
vector<atom> getRes(vector<atom> cnf, int Nres);
int findLink(vector<atom> cnf);
bool isIn(vector<atom> cnf, string str);
bool isIn(vector<vector<int> > skips, int i);
bool isIn(vector<string> vec, string str);
vector<vector<int> > findSkpis(vector<atom> cnf, int nrs);
//vector<int> findConect(vector<atom> cnf, int nrs);
vector<pair <int, int> > findConnect(vector<atom> cnf, int nrs);
vector<pair <int, int> > findConnect(vector<atom> cnf);
helix getHelix(vector<atom> cnf, int start, map<string,LJparam> LJpairs);
vector<pair<string, double > > getClash(vector<atom> res1, vector<atom> res2, vector<string > clashList, map<string,LJparam> LJpairs);
vector<pair<string, double > > getClashHelix(vector<atom> res1, vector<atom> res2, vector<string > clashList, map<string,LJparam> LJpairs);
vector<pair<string, double > > getClashTest(vector<atom> res1, vector<atom> res2, vector<string > clashList, map<string,LJparam> LJpairs);
double getLJ(atom i, atom j, double C12, double C6);
double getTorsionEnergy(atom i, atom j, atom k, atom l, double K, double s, double m);
double getSteElecE_just_LJ(vector<atom> res1, vector<atom> res2, long long int link,  map<string,LJparam> LJpairs);
double getSteElecE_just_LJ_extended(vector<atom> res1, vector<atom> res2, long long int link,  map<string,LJparam> LJpairs);
double getSteElecE_just_torsion(vector<atom> res1, vector<atom> res2, long long int ilink,  map<int,Torsion_param> Torsions);
vector<pair<string, int> > getHbond(vector<atom> angle, vector<vector<string> > tpl, double Hdistmin, double Hdistmax, 
        double Hangle, vector<string > Hlist,vector<string > acceptList);
vector<pair<string, int> > getHbond(vector<atom> angle, int firstR, int secondR, vector<vector<string> > tpl, 
        double Hdistmin, double Hdistmax, double Hangle, vector<string > Hlist,vector<string > acceptList);
void getHelixPar(vector<atom> cnfneu, helix & tmp);
int getDisacch(vector<atom > angleIn, vector<vector<string> > tpl, string link, float delta, vector<vector<helix> > &Hmatrix);
int getDisacch(vector<atom > angleIn, vector<vector<string> > tpl, string link, float delta, vector<vector<vector<helix> > > &Hmatrix);
int getHelixMap(vector<atom> angleIn, vector<vector<string> > tpl, string link, float delta, vector<vector<helix> > &Hmatrix);
int getHelixMap(vector<atom > angleIn, vector<vector<string> > tpl, string link, float delta, vector<vector<vector<helix> > > &Hmatrix);

//making files functions
vector<atom> orderCnf(vector<atom> res, vector<vector<string> > tpl, int nrAlast);
void makeOpt(vector<atom > angle);
void mkCnfOpt(vector<atom> cnf, vector<vector<string> > tpl);
void mkCnfOpt(vector<atom> cnf, vector<vector<string> > tpl, string path);
void makeHelix(vector<vector<helix> > h, string par, string name, string dir);
void makeHelix(vector<vector<vector<helix> > > h, string par, string name, string dir);

/*
 * The program has several options which are printed when program is called
 * with no argument or with --help
 */
int main(int argc, char** argv) {
    //if one wants to print arguments
//    for(int i = 1; i<argc; i++){
//        cout << argv[i] << endl;
//    }
    
    if (argc == 1 || string(argv[1]).compare("--help")==0){
        cout << "Usage: carbosecstruc [OPTION] FILES" << endl;
        cout << endl;
        cout << "  -t        transform cnf file to angular file" << endl;
        cout << "  -m        make cnf file from angular file" << endl;
        cout << "  -test     option to test conversion cartesian->angular->cartesian" << endl;
        cout << "  -test2    option to test conversion angular->cartesian->angular" << endl;
        cout << "  -help    to display info" << endl;
        cout << "  -helix  to calculate ramachandran plot of helix parameters" << endl;
        cout << "  -helmap to calculate ramachandran plot of helix parameters" << endl;
        cout << "  -disacch to calculate ramachandran plot of various parameters for disacharides" << endl;
        cout << "           variation of alpha and bete epimery possible" << endl;
        cout << endl;
        cout << "EXAMLPES:" << endl;
        cout << " ./carbosecstruc -t      source.cnf    template.tpl  anglefile.opt" << endl;
        cout << " ./carbosecstruc -m      anglefile.opt template.tpl  output.cnf [change file optional] [changed angle file optional]" << endl;
        cout << " ./carbosecstruc -test   source.cnf    template.tpl  anglefile.opt output.cnf" << endl;
        cout << " ./carbosecstruc -test2  source.ang    template.tpl  output.cnf output.ang" << endl;
        cout << " ./carbosecstruc -helix GLC.ang phex.tpl 5 GLC 4 beta" << endl;
        cout << " ./carbosecstruc -helmap  tetrasacch.ang phex.tpl delta_angle dissach_name linkNo" << endl;
        cout << " ./carbosecstruc -helmap  lib_tetrasacch/GlcA3.ang phex.tpl 5 GlcA3 3" << endl;
        cout << " ./carbosecstruc -disacch disacch.ang phex.tpl 5 dissach_name linkNo" << endl;
        cout << " ./carbosecstruc -disacch lib_disacch/GlcA3AllB.ang phex.tpl 30 GlcA3AllB 3" << endl;

        return 0;
    }
    
    ////////////////////////////////////////////////////////////////////////
    //******option to test conversion cartesian->angular->cartesian ******//
    ////////////////////////////////////////////////////////////////////////
    if (string(argv[1]).compare("-test")==0){
        //assign file path to program arguments
        cnfInPath = argv[2];
        tplPath = argv[3];
        angleOutPath = argv[4];
        angleInPath = argv[4];
        cnfOutPath = argv[5];
        cout << ">>>>>>>>>>testing<<<<<<<<<<<<<" << endl;
        //tpl file is read to tpl vector
        vector<string> one(4);
        vector<vector<string> > tpl;
        tpl.resize(17, one);
        readTpl(tpl);
        cout << "template has been read in-------------------" << endl;

        //cnf file is read to cnf[17] matrix
        vector<atom> cnf;
        readCnf(cnf);
        //printCnf(cnf);
        cout << "cnf file has been read in-------------------" << endl;

        //angular matrix is made and file is written
        vector<atom > angle;
        fillTpl(tpl, cnf, angle);
        cout << "angular coordinates have been calculated----" << endl;
        makeOpt(angle);
        cout << "angular file has been produced--------------" << endl;

        //reading in angular variable  into angleIN[17] matrix
        vector<atom > angleIn;
        readInp(angleIn);
        //printCnf(angleIn);
        cout << "angular coordinates has been read in--------" << endl;

        //coordinate matrix is created and file is written
        vector<atom> cnfneu;
        if (fillCnf(angleIn, tpl, cnfneu)==1){
            return 1;
        }
        cout << "cartesian coordinates have been calculated--" << endl;
        mkCnfOpt(cnfneu,tpl);
        cout << "output cnf file has beeen created-----------" << endl;
    } 
    else if  (string(argv[1]).compare("-test2")==0){
        //assign file path to program arguments
        cnfInPath = argv[4]; //not used
        tplPath = argv[3];
        angleOutPath = argv[5];
        angleInPath = argv[2];
        cnfOutPath = argv[4];
        cout << ">>>>>>>>>>testing<<<<<<<<<<<<<" << endl;
        //tpl file is read to tpl vector
        vector<string> one(4);
        vector<vector<string> > tpl;
        tpl.resize(17, one);
        readTpl(tpl);
        cout << "template has been read in-------------------" << endl;
        
        //reading in angular variable  into angleIN[17] matrix
        vector<atom > angleIn;
        readInp(angleIn);
        //printCnf(angleIn);
        cout << "angular coordinates has been read in--------" << endl;

        //coordinate matrix is created and file is written
        vector<atom> cnfneu;
        if (fillCnf(angleIn, tpl, cnfneu)==1){
            return 1;
        }
        cout << "cartesian coordinates have been calculated--" << endl;
        mkCnfOpt(cnfneu,tpl);
        cout << "output cnf file has beeen created-----------" << endl;
        
        //cnf file is read to cnf[17] matrix
        vector<atom> cnf;
        readCnf(cnf);
        //printCnf(cnf);
        cout << "cnf file has been read in-------------------" << endl;
        
        //angular matrix is made and file is written
        vector<atom > angle;
        fillTpl(tpl, cnf, angle);
        cout << "angular coordinates have been calculated----" << endl;
        makeOpt(angle);
        cout << "angular file has been produced--------------" << endl;

        
    /////////////////////////////////////////
    //****** option to make cnf file ******//
    /////////////////////////////////////////
    } 
    else if (string(argv[1]).compare("-m")==0){
        //assign file path to program arguments
        angleInPath  = argv[2];
        tplPath    = argv[3];
        cnfOutPath = argv[4];   
        if (argc >= 6) changePath = argv[5];
        if (argc == 7) angleOutPath = argv[6];
        cout << ">>>>>>>>>>make cnf from angleFile<<<<<<<<<<<<<" << endl;
        //tpl file is read to tpl vector
        vector<string> aone(4);
        vector<vector<string> > tpl;
        tpl.resize(17, aone);
        readTpl(tpl);
        cout << "template has been read in-------------------" << endl;
        
        //reading in angular variable  into angleIN[17] matrix
        vector<atom> angleIn;
        readInp(angleIn);
        //printCnf(inp);
        cout << "angular coordinates has been read in--------" << endl;
        
        if (argc >= 6) {
            vector<vector<string> > change;
            readChange(change);
            cout << "change has been read in---------------------" << endl;


            //printVec(change);
            for (vector<vector<string> >::iterator it = change.begin(); it < change.end(); it++){
                //printVec(*it);
                int state = changeAngle(angleIn, tpl, (*it)[0], (*it)[1]);
                //if state is >0 some error occured in changeAngle()
                if(state>0){
                    return 1;
                }
            }
            cout << "angular coordinates has been changed--------" << endl;
        }
        
        if (argc == 7) makeOpt(angleIn);
        //printCnf(angleIn);
        //coordinate matrix is created and file is written
        vector<atom> cnfneu;
        if (fillCnf(angleIn, tpl, cnfneu)==1){
            cerr << "error in fillCnf()" << endl;
            return 1;
        }
        
        //printCnf(cnfneu);
        cout << "cartesian coordinates have been calculated--" << endl;
        mkCnfOpt(cnfneu,tpl);
        cout << "output cnf file has beeen created-----------" << endl;
     
    /////////////////////////////////////////////
    //****** option to make angular file ******//
    /////////////////////////////////////////////
    } 
    else if (string(argv[1]).compare("-t")==0){
        //assign file path to program arguments
        cnfInPath = argv[2];
        tplPath = argv[3];
        angleOutPath = argv[4];
        cout << ">>>>>>>>>>transforming<<<<<<<<<<<<<" << endl;
        //tpl file is read to tpl[17][4] matrix
        vector<string> one(4);
        vector<vector<string> > tpl;
        tpl.resize(17, one);
        readTpl(tpl);
        //printTpl(tpl);
        cout << "template has been read in-------------------" << endl;

        //cnf file is read to cnf[17] matrix
        vector<atom> cnf;
        readCnf(cnf);
        //printCnf(cnf);
        cout << "cnf file has been read in-------------------" << endl;
        
        //angular matrix is made and file is written
        vector<atom > angle;
        fillTpl(tpl, cnf, angle);
        cout << "angular coordinates have been calculated----" << endl;
        makeOpt(angle);
        cout << "angular file has been produced--------------" << endl;
        //printCnf(angle);
        

    ///////////////////////////////////////////////////////////////
    //************************ HELIX ******************************
    ///////////////////////////////////////////////////////////////
    } 
    else if (string(argv[1]).compare("-helix")==0){
        ///////////////////////////////
        //reading in
        angleInPath = argv[2];
        tplPath = argv[3];
        float delta = stof(argv[4]); //resolution (deg)
        string sugar = argv[5];   //sugar name
        string link = argv[6];    //type of linking 1-?
        string anomer = argv[7];  //anomer
        cout << ">>>>>>>>>>Calculating ramachandran map<<<<<<<<<<<<<" << endl;
        //angular file is read into cnf vector
        vector<atom > angleIn;
        readInp(angleIn);
        cout << "angular file has been read in-------------------" << endl;
        
        //template is red into tpl
        vector<string> one(4);
        vector<vector<string> > tpl;
        tpl.resize(17, one);
        readTpl(tpl);
        cout << "template has been read in-------------------" << endl;
        ///////////////////////////////
        
        ///////////////////////////////
        //check anomer, made for D sugars
        string anflag; //anomer flag
        float DL = angleIn[findAtomPos("C6", angleIn, 1)].pos.z; //D or L sugar
        //cout << DL << endl;
        if(DL<180){
            if(anomer.compare("alpha")==0) anflag = "down";
            else if(anomer.compare("beta")==0) anflag = "up";
        }else if (DL>180){
            if(anomer.compare("alpha")==0) anflag = "up";
            else if(anomer.compare("beta")==0) anflag = "down";
        }
        
        string inst = "1 O1 " + anflag;
        int state = changeAngle(angleIn, tpl, "ena", inst);
        //if state is >0 some error occured in changeAngle()
        if(state>0)return 1;
        ///////////////////////////////

        
        ///////////////////////////////
        //check length of chain
        if (angleIn.back().nrRes!=1){
            cout << "-----------------------------" << endl;
            cout << "angle file consists of " << angleIn.back().nrRes << 
                    " residues. "  << angleIn.back().nrRes-1 << 
                    "residues are going to be removed" << endl;
            
            string inst = std::to_string(static_cast<long long>(angleIn.back().nrRes-1));
                    
            int state = changeAngle(angleIn, tpl, "del", inst);
            //if state is >0 some error occured in changeAngle()
            if(state>0) return 1;
        }
        
        if (link.compare("1")!=0){
            //make 4 residues long chain
            inst = "1 " + link + " 3";
            state = changeAngle(angleIn, tpl, "add", inst);
            //if state is >0 some error occured in changeAngle()
            if(state>0) return 1;
        } else {
            //make 2 residues long chain
            inst = "1 " + link + " 1";
            state = changeAngle(angleIn, tpl, "add", inst);
            //if state is >0 some error occured in changeAngle()
            if(state>0) return 1;
        }
        
        //////////////////////////////////////////////////////////////
        ///       create ramachandran plots   ////////////////////////
        //////////////////////////////////////////////////////////////
        
        //////////////////////////////////////////////////////////////
        // branch for 1-6 link
        //////////////////////////////////////////////////////////////       
        if (link.compare("6")==0){
            vector<vector<vector<helix> > > Hmatrix;
            ///////////////////////////////
            //loop over omega
            for (long double omega = -(180-delta/2); omega<180;omega += delta){
                //make psi string and change angle file
                string somega;
                if (omega<=0.0) somega = std::to_string(360+omega);//convert phi to positive angle
                else somega = std::to_string(omega);

                for (long long int r = 1; r<5; r++){
                    string par = std::to_string(r)+" o "+somega;
                    int state = changeAngle(angleIn, tpl, "rot", par);
                    //if state is >0 some error occured in changeAngle()
                    if(state>0){
                        return 1;
                    }
                }
                vector<vector<helix> > slice; //single row
                ///////////////////////////////
                //loop over psi
                for (long double psi = (180-delta/2); psi>-180;psi -= delta){
                    
                    //make phi string and change angle file
                    string spsi;
                    if (psi<=0.0) spsi = std::to_string(360+psi);//convert psi to positive angle
                    else spsi = std::to_string(psi);            

                    int state = changeAngle(angleIn, tpl, "allpsi", spsi);

                        //if state is >0 some error occured in changeAngle()
                    if(state>0){
                        return 1;
                    }

                    vector<helix> hrow; //single row
                    
                    ///////////////////////////////
                    //loop over phi
                    for (long double phi = -(180-delta/2); phi<180;phi += delta){

                        //make psi string and change angle file
                        string sphi;
                        if (phi<=0.0) sphi = std::to_string(360+phi);//convert phi to positive angle
                        else sphi = std::to_string(phi);

                        int state = changeAngle(angleIn, tpl, "allphi", sphi);

                        //if state is >0 some error occured in changeAngle()
                        if(state>0){
                            return 1;
                        }

                        vector<atom> cnfneu;
                        if (fillCnf(angleIn, tpl, cnfneu)==1){
                            cerr << "error in fillCnf()" << endl;
                            return 1;
                        }
                        //mkCnfOpt(cnfneu,tpl,"data/"+anomer+"_"+sugar+"_"+"1-"+link+"_phi"+sphi+"_psi"+spsi+"_omega"+somega+".cnf");
                        //printf("phi= %f ; psi = %f ; omega = %f \n", static_cast<double>(phi), static_cast<double>(psi), static_cast<double>(omega) );

                        ///////////////////////////////////////
                        //Helix calculation calculation
                        ///////////////////////////////////////
                        helix tmp;
                        tmp = getHelix(cnfneu, 2, LJpairs);
                        ///////////////////////////////////////
                        //Stereoelectronic energy calculation
                        ///////////////////////////////////////
                        vector<atom> res1 = getRes(cnfneu,1);
                        vector<atom> res2 = getRes(cnfneu,2);

                        long long int ilink = (long long) stoi(link);
                                       
                        tmp.SteElecE_just_TOR = getSteElecE_just_torsion(res1, res2, ilink, Torsions);
                        tmp.SteElecE_just_LJ_extended = getSteElecE_just_LJ_extended(res1, res2, ilink, LJpairs);
                        tmp.SteElecE_extended = tmp.SteElecE_just_LJ_extended + tmp.SteElecE_just_TOR;

                        ///////////////////////////////////////
                        //H bond calculation
                        ///////////////////////////////////////
                        double Hdistmin = 0.0;
                        double Hdistmax = 0.25;
                        double Hangle = 135;


                        tmp.Hbond = getHbond(angleIn, tpl, Hdistmin, Hdistmax, Hangle, Hlist, acceptList);

                        tmp.phi = phi;
                        //cout << tmp.phi << endl;
                        tmp.psi = psi;
                        tmp.omega = omega;
                        hrow.push_back(tmp);
                    }
                    slice.push_back(hrow);
                }
                Hmatrix.push_back(slice);
            }
            ///////////////////////////////
            //make outputs
            string name = sugar+anomer+"1-"+link;
            makeHelix(Hmatrix, "R", name, "data");
            makeHelix(Hmatrix, "p", name, "data");
            makeHelix(Hmatrix, "a", name, "data");
            makeHelix(Hmatrix, "g", name, "data");
            makeHelix(Hmatrix, "n", name, "data");
            makeHelix(Hmatrix, "clash", name, "data"); 
            makeHelix(Hmatrix, "Hbond", name, "data");
            makeHelix(Hmatrix, "SteElecE_extended", name, "data");
            makeHelix(Hmatrix, "SteElecE_just_LJ_extended", name, "data");
            makeHelix(Hmatrix, "SteElecE_just_TOR", name, "data");
            makeHelix(Hmatrix, "phi", name, "data");
            makeHelix(Hmatrix, "psi", name, "data");
            makeHelix(Hmatrix, "omega", name, "data");
            ///////////////////////////////
            
        //////////////////////////////////////////////////////////////
        // branch for OTHER THEN 1-6 link
        //////////////////////////////////////////////////////////////
        } else {   
            vector<vector<helix> > Hmatrix;
            for (long double psi = (180-delta/2); psi>-180;psi -= delta){
                //make phi string and change angle file
                string spsi;
                if (psi<=0.0) spsi = std::to_string(360+psi);//convert psi to positive angle
                else spsi = std::to_string(psi);            

                int state = changeAngle(angleIn, tpl, "allpsi", spsi);

                    //if state is >0 some error occured in changeAngle()
                if(state>0){
                    return 1;
                }

                vector<helix> hrow; //single row
                for (long double phi = -(180-delta/2); phi<180;phi += delta){

                    //make psi string and change angle file
                    string sphi;
                    if (phi<=0.0) sphi = std::to_string(360+phi);//convert phi to positive angle
                    else sphi = std::to_string(phi);

                    int state = changeAngle(angleIn, tpl, "allphi", sphi);

                    //if state is >0 some error occured in changeAngle()
                    if(state>0){
                        return 1;
                    }

                    vector<atom> cnfneu;
                    if (fillCnf(angleIn, tpl, cnfneu)==1){
                        cerr << "error in fillCnf()" << endl;
                        return 1;
                    }
                    //mkCnfOpt(cnfneu,tpl,"data/"+anomer+"_"+sugar+"_"+"1-"+link+"_psi"+spsi+"_phi"+sphi+".cnf");
                    //printf("phi= %f ; psi = %f \n", static_cast<double>(phi), static_cast<double>(psi) );

                    ///////////////////////////////////////
                    //Helix calculation calculation
                    ///////////////////////////////////////


                    helix tmp;
                    if (link.compare("1")==0){
                        ///////////////////////////////////////
                        //clash calculation
                        ///////////////////////////////////////
                        vector<atom> res1 = getRes(cnfneu,1);
                        vector<atom> res2 = getRes(cnfneu,2);

                        vector<string > clashList {"C1","C2","C3","C4","C5","C6","O1","O2","O3","O4","O5"};

                        tmp.clash = getClash(res1, res2, clashList, LJpairs);
                    } else {
                        tmp = getHelix(cnfneu, 2, LJpairs);                    
                    }
                    ///////////////////////////////////////
                    //Stereoelectronic energy calculation
                    ///////////////////////////////////////
                    vector<atom> res1 = getRes(cnfneu,1);
                    vector<atom> res2 = getRes(cnfneu,2);

                    long long int ilink = (long long) stoi(link);

//                    tmp.SteElecE_just_TOR = getSteElecE_just_torsion(res1, res2, ilink, Torsions);
//                    tmp.SteElecE_just_LJ_extended = getSteElecE_just_LJ_extended(res1, res2, ilink, LJpairs);
//                    tmp.SteElecE_extended = tmp.SteElecE_just_LJ_extended + tmp.SteElecE_just_TOR;

                    ///////////////////////////////////////
                    //H bond calculation
                    ///////////////////////////////////////
                    double Hdistmin = 0.0;
                    double Hdistmax = 0.25;
                    double Hangle = 135;



//                    tmp.Hbond = getHbond(angleIn, tpl, Hdistmin, Hdistmax, Hangle, Hlist, acceptList);

                    tmp.phi = phi;
                    tmp.psi = psi;
                    hrow.push_back(tmp);
                }
                Hmatrix.push_back(hrow);
            }
            ///////////////////////////////
            //make outputs
            string name = sugar+anomer+"1-"+link;
            if (link.compare("1")!=0){
                makeHelix(Hmatrix, "R", name, "data");
                makeHelix(Hmatrix, "p", name, "data");
                makeHelix(Hmatrix, "a", name, "data");
                makeHelix(Hmatrix, "g", name, "data");
                makeHelix(Hmatrix, "n", name, "data");
            }
//            makeHelix(Hmatrix, "clash", name, "data"); 
//            makeHelix(Hmatrix, "Hbond", name, "data");
//            makeHelix(Hmatrix, "SteElecE_extended", name, "data");
//            makeHelix(Hmatrix, "SteElecE_just_LJ_extended", name, "data");
//            makeHelix(Hmatrix, "SteElecE_just_TOR", name, "data");
//            makeHelix(Hmatrix, "phi", name, "data");
//            makeHelix(Hmatrix, "psi", name, "data");
            ///////////////////////////////
        
        }
            
    ///////////////////////////////////////////////////////////////
    //****** option to calculate disaccharide ramachandran maps ******
    // input is prepared disaccharide
    ///////////////////////////////////////////////////////////////
    }
    else if (string(argv[1]).compare("-disacch")==0){
        ///////////////////////////////
        //reading in
        angleInPath = argv[2];
        tplPath = argv[3];
        float delta = stof(argv[4]); //resolution (deg)
        string sugar = argv[5];   //sugar name
        string link = argv[6];    //type of linking 1-?

        cout << ">>>>>>>>>>Calculating ramachandran map<<<<<<<<<<<<<" << endl;
        //angular file is read into cnf vector
        vector<atom > angleIn;
        readInp(angleIn);
        cout << "angular file has been read in-------------------" << endl;
        
        //template is red into tpl
        vector<string> one(4);
        vector<vector<string> > tpl;
        tpl.resize(17, one);
        readTpl(tpl);
        cout << "template has been read in-------------------" << endl;
        ///////////////////////////////
        

        //////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////
        ///       create ramachandran plots   ////////////////////////
        //////////////////////////////////////////////////////////////
        
        //////////////////////////////////////////////////////////////
        // branch for 1-6 link
        if (link.compare("6")==0){
            vector<vector<vector<helix> > > Hmatrix;
            
            int state = getDisacch(angleIn, tpl, link, delta, Hmatrix);
            
            if (state != 0){
                printf("getDisacch() failed in 1-6 branch of -disacch2 option");
                return 1;
            }
            ///////////////////////////////
            //make outputs
            string name = sugar+"_";

            makeHelix(Hmatrix, "clash", name, "data_disacch"); 
            makeHelix(Hmatrix, "Hbond", name, "data_disacch");
            makeHelix(Hmatrix, "SteElecE", name, "data_disacch");
            makeHelix(Hmatrix, "SteElecE_extended", name, "data_disacch");
            makeHelix(Hmatrix, "SteElecE_just_TOR", name, "data_disacch");
            makeHelix(Hmatrix, "SteElecE_just_LJ", name, "data_disacch");
            makeHelix(Hmatrix, "SteElecE_just_LJ_extended", name, "data_disacch");
            makeHelix(Hmatrix, "phi", name, "data_disacch");
            makeHelix(Hmatrix, "psi", name, "data_disacch");
            makeHelix(Hmatrix, "omega", name, "data_disacch");

        //////////////////////////////////////////////////////////////
        // branch for OTHER THEN 1-6 link
        } else {
            vector<vector<helix> > Hmatrix;
            
            int state = getDisacch(angleIn, tpl, link, delta, Hmatrix);
            
            if (state != 0){
                printf("getDisacch() failed in 1-6 branch of -disacch2 option");
                return 1;
            }
            ///////////////////////////////
            //make outputs
            string name = sugar+"_";

            makeHelix(Hmatrix, "clash", name, "data_disacch"); 
            makeHelix(Hmatrix, "Hbond", name, "data_disacch");
            makeHelix(Hmatrix, "SteElecE", name, "data_disacch");
            makeHelix(Hmatrix, "SteElecE_extended", name, "data_disacch");
            makeHelix(Hmatrix, "SteElecE_just_TOR", name, "data_disacch");
            makeHelix(Hmatrix, "SteElecE_just_LJ", name, "data_disacch");
            makeHelix(Hmatrix, "SteElecE_just_LJ_extended", name, "data_disacch");
            //makeHelix(Hmatrix, "phi", name, "data_disacch");
            //makeHelix(Hmatrix, "psi", name, "data_disacch");
        }

        
        
    }
    else if (string(argv[1]).compare("-helmap")==0){
        ///////////////////////////////
        //reading in
        angleInPath = argv[2];
        tplPath = argv[3];
        float delta = stof(argv[4]); //resolution (deg)
        string sugar = argv[5];   //sugar name
        string link = argv[6];    //type of linking 1-?

        cout << ">>>>>>>>>>Calculating ramachandran map<<<<<<<<<<<<<" << endl;
        //angular file is read into cnf vector
        vector<atom > angleIn;
        readInp(angleIn);
        cout << "angular file has been read in-------------------" << endl;
        
        //template is red into tpl
        vector<string> one(4);
        vector<vector<string> > tpl;
        tpl.resize(17, one);
        readTpl(tpl);
        cout << "template has been read in-------------------" << endl;
        ///////////////////////////////
        
        if (angleIn.back().nrRes != 4){
            cout << "loaded sugar in -helmap should contain 4 residues" << endl;
            return 1;
        }

        //////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////
        ///       create ramachandran plots   ////////////////////////
        //////////////////////////////////////////////////////////////
        
        //////////////////////////////////////////////////////////////
        // branch for 1-6 link
        if (link.compare("6")==0){
            vector<vector<vector<helix> > > Hmatrix;
            
            int state = getHelixMap(angleIn, tpl, link, delta, Hmatrix);
            
            if (state != 0){
                printf("getHelixMap() failed in 1-6 branch of -helmap option");
                return 1;
            }
            
            //tmp = getHelix(cnfneu, 2, LJpairs);
            ///////////////////////////////
            //make outputs
            string name = sugar+"_";

            makeHelix(Hmatrix, "clash", name, "data_tetrasach"); 
            makeHelix(Hmatrix, "clashHelix", name, "data_tetrasach");
            makeHelix(Hmatrix, "Hbond", name, "data_tetrasach");
            makeHelix(Hmatrix, "HbondHelix", name, "data_tetrasach");
            makeHelix(Hmatrix, "SteElecE", name, "data_tetrasach");
            makeHelix(Hmatrix, "SteElecE_extended", name, "data_tetrasach");
            makeHelix(Hmatrix, "SteElecE_just_TOR", name, "data_tetrasach");
            makeHelix(Hmatrix, "SteElecE_just_LJ", name, "data_tetrasach");
            makeHelix(Hmatrix, "SteElecE_just_LJ_extended", name, "data_tetrasach");
            makeHelix(Hmatrix, "phi", name, "data_tetrasach");
            makeHelix(Hmatrix, "psi", name, "data_tetrasach");
            makeHelix(Hmatrix, "omega", name, "data_tetrasach");
            makeHelix(Hmatrix, "R", name, "data_tetrasach");
            makeHelix(Hmatrix, "p", name, "data_tetrasach");
            makeHelix(Hmatrix, "a", name, "data_tetrasach");
            makeHelix(Hmatrix, "g", name, "data_tetrasach");
            makeHelix(Hmatrix, "n", name, "data_tetrasach");

        //////////////////////////////////////////////////////////////
        // branch for OTHER THEN 1-6 link
        } 
        else {
            vector<vector<helix> > Hmatrix;
            
            
            int state = getHelixMap(angleIn, tpl, link, delta, Hmatrix);
            
            if (state != 0){
                printf("getHelixMap() failed in NON 1-6 branch of -helmap option");
                return 1;
            }
            
            ///////////////////////////////
            //make outputs
            string name = sugar+"_";

            makeHelix(Hmatrix, "clash", name, "data_tetrasach"); 
            makeHelix(Hmatrix, "clashHelix", name, "data_tetrasach"); 
            makeHelix(Hmatrix, "Hbond", name, "data_tetrasach");
            makeHelix(Hmatrix, "HbondHelix", name, "data_tetrasach");
            makeHelix(Hmatrix, "SteElecE", name, "data_tetrasach");
            makeHelix(Hmatrix, "SteElecE_extended", name, "data_tetrasach");
            makeHelix(Hmatrix, "SteElecE_just_TOR", name, "data_tetrasach");
            makeHelix(Hmatrix, "SteElecE_just_LJ", name, "data_tetrasach");
            makeHelix(Hmatrix, "SteElecE_just_LJ_extended", name, "data_tetrasach");
            makeHelix(Hmatrix, "phi", name, "data_tetrasach");
            makeHelix(Hmatrix, "psi", name, "data_tetrasach");
            makeHelix(Hmatrix, "R", name, "data_tetrasach");
            makeHelix(Hmatrix, "p", name, "data_tetrasach");
            makeHelix(Hmatrix, "a", name, "data_tetrasach");
            makeHelix(Hmatrix, "g", name, "data_tetrasach");
            makeHelix(Hmatrix, "n", name, "data_tetrasach");
        }

        
    }
    else if (string(argv[1]).compare("-fluct")==0){
        ///////////////////////////////
        //reading in
        cnfInPath = argv[2];
        tplPath = argv[3];
        trcPath = argv[4];
        outPath = argv[5];
        string link = argv[6];


        cout << ">>>>>>>>>>Calculating ramachandran map<<<<<<<<<<<<<" << endl;
        //angular file is read into cnf vector
        vector<atom > cnfIn;
        readCnf(cnfIn);
        cout << "angular file has been read in-------------------" << endl;
        
        //template is red into tpl
        vector<string> one(4);
        vector<vector<string> > tpl;
        tpl.resize(17, one);
        readTpl(tpl);
        cout << "template has been read in-------------------" << endl;
        ///////////////////////////////

        vector<vector<atom> > trc;
        
        readTrc(cnfIn, trc);
        cout << "trc has been read in-------------------" << endl;
        
        //printVec(trc);

        
        // opend output file
        fstream ofile(outPath.c_str(), ios::out);
        if ( ofile.fail() ){
            cout << "opt file cannot be opened.\n";
        }
        ofile << "#      phi       psi  Hbond  clash     SEE" <<
                "     SEE_extnd    SEE_j_LJ  SEE_j_LJ_extnd  SEE_j_TOR  \n";

            
        // calculate parameters
        for (vector<vector<atom> >::iterator it = trc.begin() ; it != trc.end(); ++it){
            
            vector<atom> res1 = getRes(*it,1);
            vector<atom> res2 = getRes(*it,2);

            vector<string > clashList {"C1","C2","C3","C4","C5","C6","O1","O2","O3","O4","O5"};

            vector<pair<string, double > > clash = getClash(res1, res2, clashList, LJpairs);
            ///////////////////////////////////////
            //Stereoelectronic energy calculation
            ///////////////////////////////////////
            res1 = getRes(*it,1);
            res2 = getRes(*it,2);

            long long int ilink = (long long) stoi(link);

            double SteElecE_just_TOR = getSteElecE_just_torsion(res1, res2, ilink, Torsions);

            double SteElecE_just_LJ = getSteElecE_just_LJ(res1, res2, ilink, LJpairs);
            double SteElecE = SteElecE_just_LJ + SteElecE_just_TOR;

            double SteElecE_just_LJ_extended = getSteElecE_just_LJ_extended(res1, res2, ilink, LJpairs);
            double SteElecE_extended = SteElecE_just_LJ_extended + SteElecE_just_TOR;
            ///////////////////////////////////////
            //H bond calculation
            ///////////////////////////////////////
            vector<atom> angle;
            fillTpl(tpl, *it, angle);
            
            vector<pair<string, int> > Hbond = getHbond(angle, tpl, Hdistmin, Hdistmax, Hangle, Hlist, acceptList);
            ///////////////////////////////////////
            // save phi, psi, and omega calculation
            ///////////////////////////////////////
            // this is made just for 1-4 link"
            atom O5 = findAtom("O5", *it, 1);
            atom C1 = findAtom("C1", *it, 1);
            atom O1 = findAtom("O1", *it, 1);
            atom C4 = findAtom("C4", *it, 2);
            atom C3 = findAtom("C3", *it, 2);
            
            double phi = calcDihedral(O5.pos, C1.pos, O1.pos, C4.pos);
            double psi = calcDihedral(C1.pos, O1.pos, C4.pos, C3.pos);
            ///////////////////////////////////////
            // Writing
            ///////////////////////////////////////
            //# phi    psi    Hbond  clash     SteElecE    SteElecE_extended    SteElecE_just_LJ    SteElecE_just_LJ_extended    SteElecE_just_TOR
            ofile.width(10);
            ofile.precision(4);
            ofile <<  phi ;
            ofile.width(10);
            ofile.precision(4);
            ofile <<  psi ;
            ofile.width(7);
            ofile.precision(3);
            ofile <<  Hbond.back().second  ;
            ofile.width(7);
            ofile.precision(3);
            ofile <<  clash.front().second ;
            ofile.width(11);
            ofile.precision(5);
            ofile <<  SteElecE ;
            ofile.width(11);
            ofile.precision(5);
            ofile <<  SteElecE_extended ;
            ofile.width(12);
            ofile.precision(5);
            ofile <<  SteElecE_just_LJ ;
            ofile.width(16);
            ofile.precision(5);
            ofile <<  SteElecE_just_LJ_extended ;
            ofile.width(11);
            ofile.precision(5);
            ofile <<  SteElecE_just_TOR << "\n";
        }
        ofile.close();
    }
    
//END of MAIN
} 




///////////////////////////////////////////////////////////////////////////////
////////////////////                               ////////////////////////////
////////////////////  DEFINITION OF ALL FUNCTIONS  ////////////////////////////
////////////////////                               ////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//basic geometric functions
void getVector(const vec &a){
    cout << " ( " << a.x << " , " <<a.y << " , "<< a.z << " ) \n";
}
vec getSubstrVec(const vec &a, const vec &b){
    vec v;
    v.x = b.x - a.x;
    v.y = b.y - a.y;
    v.z = b.z - a.z;
    return v;
}
vec getSumVec(const vec &a, const vec &b){
    vec v;
    v.x = b.x + a.x;
    v.y = b.y + a.y;
    v.z = b.z + a.z;
    return v;
}
double getAbs2Vec(const vec &v){
    const double abs2 = v.x * v.x + v.y * v.y + v.z * v.z;
    return abs2;
}
double getAbsVec(const vec &v){
    const double abs = sqrt(getAbs2Vec(v));
    return abs;
}
double getDotofVec(const vec &a, const vec &b){
    const double dot = a.x * b.x + a.y * b.y + a.z * b.z;
    return dot;
}
vec getCrossofVec(const vec &a, const vec &b){
    vec v;
    v.x = a.y * b.z - a.z * b.y;
    v.y = a.z * b.x - a.x * b.z;
    v.z = a.x * b.y - a.y * b.x;
    return v;
}
vec getMultipVec(const vec &a, double x){
    vec v;
    v.x = x*a.x;
    v.y = x*a.y;
    v.z = x*a.z;
    return v;
}
double calcDistance(const vec &a, const vec &b){
    const vec rab = getSubstrVec(a, b);
    const double distance = getAbsVec(rab);
    return distance;
    }
double calcAngle(const vec &a, const vec &b, const vec &c){
    const vec rab = getSubstrVec(a, b);
    const vec rcb = getSubstrVec(c, b);
    double angle = (getDotofVec(rab, rcb)/(getAbsVec(rab)*getAbsVec(rcb)));
    if (angle < -1.0) angle = -1.0;
    if (angle > 1.0) angle = 1.0;
    angle = acos(angle)* 180.0 / PI;
    return angle;
}
double calcDihedral(const vec &i, const vec &j, const vec &k, const vec &l){
    double dihedral = 0;
    vec rij, rkj, rkl, rim, rln, rnk;
    rij = getSubstrVec(j, i);
    rkj = getSubstrVec(j, k);
    rkl = getSubstrVec(l, k); 
    
    
    double dkj2 = getAbs2Vec(rkj);
    rnk = getCrossofVec(rkj, rkl);

    assert(dkj2 != 0.0);
    const double frim = getDotofVec(rij, rkj) / dkj2;
    const double frln = getDotofVec(rkl, rkj) / dkj2;

    rim = getSubstrVec(getMultipVec(rkj, frim), rij);
    rln = getSubstrVec(rkl, getMultipVec(rkj, frln));
    const double dim = getAbsVec(rim);
    const double dln = getAbsVec(rln);

    const double ip = getDotofVec(rim, rln);
    double cosphi = ip / (dim * dln);
    
    
    if (cosphi < -1.0) cosphi = -1.0;
    if (cosphi > 1.0) cosphi = 1.0;
    
    // get the value from 0-180
    dihedral = acos(cosphi)*180/PI;
    
    // get the value from 0-360
    const double sign = getDotofVec(rij, rnk);
    if (sign < 0) {
      dihedral = -dihedral + 2 * 180; // Put it into the range of [0; 360]
    }
    
    return dihedral;
}

void printTpl(vector<vector<string> > tpl){
    for(int i=0; i<17; i++){
        for (int j=0; j<4; j++){
            cout << tpl[i][j] << " ";
        }
        cout << endl;
    }
}
void printCnf(vector<atom> cnf){
    cout.setf( std::ios::fixed, std::ios::floatfield );
    cout.precision(9);
    for(int i=0; i<cnf.size(); i++){
        cout << cnf[i].nrRes  << " ";
        cout << cnf[i].res  << " ";
        cout << cnf[i].a   << " ";
        cout << cnf[i].nrA   << " ";
        cout << cnf[i].pos.x   << " ";
        cout << cnf[i].pos.y   << " ";
        cout << cnf[i].pos.z   << " ";
        cout << endl;
    }
}
void printInp(atom inp[17]){
    for(int i=0; i<17; i++){
        cout << inp[i].a   << " ";
        cout << inp[i].nrA   << " ";
        cout << inp[i].pos.x   << " ";
        cout << inp[i].pos.y   << " ";
        cout << inp[i].pos.z   << " ";
        cout << endl;
    }
}
void printVec(vec v){
    cout.setf( std::ios::fixed, std::ios::floatfield );
    cout.precision(9);
    cout << v.x << " " ;
    cout << v.y << " " ;
    cout << v.z << endl;
}
void printVec(vector<vector<string> > v){
    for (vector<vector<string> >::iterator it = v.begin() ; it != v.end(); ++it){
        cout << "subvector: ";
        for (vector<string>::iterator it2 = it->begin() ; it2 != it->end(); ++it2){
            cout << *it2 << ",";
        }
        cout << endl;
    }
}
void printVec(vector<vector<helix> > v, string par){
    for (vector<vector<helix> >::iterator it = v.begin() ; it != v.end(); ++it){
        cout << "subvector: ";
        for (vector<helix>::iterator it2 = it->begin() ; it2 != it->end(); ++it2){
            if     (par.compare("R")==0)   cout << (*it2).R  << ",";
            else if(par.compare("p")==0)   cout << (*it2).p  << ",";
            else if(par.compare("a")==0)   cout << (*it2).a  << ",";
            else if(par.compare("g")==0)   cout << (*it2).g  << ",";
            else if(par.compare("n")==0)   cout << (*it2).n  << ",";
            else if(par.compare("phi")==0) cout << (*it2).phi  << ",";
            else if(par.compare("psi")==0) cout << (*it2).psi  << ",";
            //else if(par.compare("clash")==0) cout << (*it2).clash  << ",";
        }
        cout << endl;
    }
}
void printVec(vector<vector<int> > v){
    for (vector<vector<int> >::iterator it = v.begin() ; it != v.end(); ++it){
        cout << "subvector: ";
        for (vector<int>::iterator it2 = it->begin() ; it2 != it->end(); ++it2){
            cout << *it2 << "  ";
        }
        cout << endl;
    }
}
void printVec(vector<pair<string, int> > v){
    for (vector<pair<string, int> >::iterator it = v.begin() ; it != v.end(); ++it){
        cout << "pair " << it->first << " " << it->second << endl;
    }
}
void printVec(vector<pair<int, int> > v){
    for (vector<pair<int, int> >::iterator it = v.begin() ; it != v.end(); ++it){
        cout << "pair " << it->first << " " << it->second << endl;
    }
}

void printVec(vector<atom> v){
    for (vector<atom>::iterator it = v.begin() ; it != v.end(); ++it){
        printAtom(*it);
    }
}

void printVec(vector<vector<atom> > v){
    for (vector<vector<atom> >::iterator it = v.begin() ; it != v.end(); ++it){
        printf("next vector \n");
        printVec(*it);
    }
}

void printVec(vector<int> v){
    cout << "vector: ";
    for (vector<int>::iterator it = v.begin() ; it != v.end(); ++it){
        cout << (*it) << " ";
    }
    cout << endl;
}

void printVec(vector<string> v){
    cout << "vector: ";
    for (vector<string>::iterator it = v.begin() ; it != v.end(); ++it){
        cout << (*it) << " ";
    }
    cout << endl;
}

void printAtom(atom a){
    printf("%d %s %s %d ", a.nrRes, a.res.c_str(), a.a.c_str(), a.nrA);
    //cout << a.nrRes << " " << a.res << " " << a.a << " " << a.nrA << " ";
    printVec(a.pos);
}

//functions used in main//
//reading functions
void readTpl(vector<vector<string> > &tpl ){
    fstream ifile(tplPath.c_str(), ios::in);
    if ( ifile.fail() ){
        cout << "tpl file cannot be opened.\n";
    }
    string input;
    ifile >> input;
    int i = 0, j =0;
    while (!ifile.fail()){
        if ( input.compare("#")==0){
                getline(ifile, input);
                ifile >> input;
                continue;
        }
        else{
            tpl[i][j] = input;
            j++;
            if (j == 4){
                j = 0;
                i++;
            }
            ifile >> input;
        }
    }
    ifile.close();
}

void readChange(vector<vector<string> > &change){
    fstream ifile(changePath.c_str(), ios::in);
    if ( ifile.fail() ){
        cout << "change file cannot be opened.\n";
    }
    string input;
    ifile >> input;
    int i = 0, j =0;
    while (!ifile.fail()){
        if ( input[0]=='#' ){
                getline(ifile, input);
                ifile >> input;
                continue;
        }
        else{
            vector<string> singleChange;    // make vector for single change
            singleChange.push_back(input);  // read flag
            getline(ifile, input);          
            singleChange.push_back(input);  //read instructions
            
            change.push_back(singleChange); //push back single change
            ifile >> input;
        }
    }
    ifile.close();
}

void readCnf(vector<atom> &cnf){
    fstream ifile(cnfInPath.c_str(), ios::in);
    if ( ifile.fail() ){
        cout << "cnf file cannot be opened.\n";
    }
    string input;
    ifile >> input;
    int i = 0;
    int j = 0;
    while (input != "POSITION"){
        ifile >> input;
        continue;
    }
    getline(ifile, input);
    while(getline(ifile, input)) { //can't be true here, I don't know why
        if (input[0] == '#'){
            continue;
        }
        atom a;
        istringstream ss(input);
        for (int j = 0; j < 7; j++){
            ss >> input;
            if (input.compare("END")==0 ) {
                a.a = input;
                break;
            }
            else if (j == 0) a.nrRes = atoi(input.c_str());
            else if (j == 1) a.res = input.c_str();
            else if (j == 2) a.a = input.c_str();
            else if (j == 3) a.nrA = atoi(input.c_str());
            else if (j == 4) a.pos.x = atof(input.c_str());
            else if (j == 5) a.pos.y = atof(input.c_str());
            else if (j == 6) a.pos.z = atof(input.c_str());
        }
        if (a.a.compare("END")==0 ) {
            break;
        }
        cnf.push_back(a);
    }
    ifile.close();
}

void readTrc(vector<atom> &cnf, vector<vector<atom> > &trc){
    readCnf(cnf);
    fstream ifile(trcPath.c_str(), ios::in);
    if ( ifile.fail() ){
        cout << "cnf file cannot be opened.\n";
    }
    string line;
    while ( getline(ifile,line) ){
        if (line.compare("POSITIONRED")!=0){
            continue;
        }
        vector<atom> tmp;
        int i = 0;
        getline(ifile, line);
        while(line.compare("END")!=0) {
            if (line[0] == '#'){
                getline(ifile, line);
                continue;
            }
            atom a = cnf[i++];
            
            istringstream ss(line);
            for (int j = 0; j < 3; j++){
                ss >> line;
                if      (j == 0) a.pos.x = atof(line.c_str());
                else if (j == 1) a.pos.y = atof(line.c_str());
                else if (j == 2) a.pos.z = atof(line.c_str());
            }

            tmp.push_back(a);
            getline(ifile, line);
        }
        trc.push_back(tmp);
    }
    ifile.close();
}

void readInp(vector<atom > &angleIn){
    fstream angleIfile(angleInPath.c_str(), ios::in);
    if ( angleIfile.fail() ){
        cout << "angleIfile file cannot be opened.\n";
    }
    string input;
    angleIfile >> input;
    while (input != "START"){
        angleIfile >> input;
        continue;
    }
    getline(angleIfile, input);
    while(getline(angleIfile, input)) { //can't be true here, I don't know why
        atom a;
        istringstream ss(input);
        for (int j = 0; j < 7; j++){
            ss >> input;
            if (input.compare("END")==0 ) {
                a.a = input; // a.a is assignet to END like a flag to jump out of second loop
                break;
            }
            else if (j == 0) a.nrRes = atoi(input.c_str());
            else if (j == 1) a.res = input.c_str();
            else if (j == 2) a.a = input.c_str();
            else if (j == 3) a.nrA = atoi(input.c_str());
            else if (j == 4) a.pos.x = atof(input.c_str());
            else if (j == 5) a.pos.y = atof(input.c_str());
            else if (j == 6) a.pos.z = atof(input.c_str());
        }
        if (a.a.compare("END")==0 ) { // using END flag to jump out of second loop
            break;
        }
        angleIn.push_back(a);
    }
    
    angleIfile.close();
}

vector<atom> getRes(vector<atom> cnf, int Nres){
    vector<atom> res;
    for(int i=0; i<cnf.size(); i++){
        if (Nres == cnf[i].nrRes){
            res.push_back(cnf[i]);
        }
    }
    return res;
}


vector<vector<int> > findSkpis(vector<atom> cnf, int nrs){
    vector<vector<int> > skips;
    //printVec(skips);
    for (int i=0; i<nrs; i++){              //loop over all residues
        vector<atom> res = getRes(cnf, i+1); //get a vector of atoms of a residue
        //printVec(res);
        vector<int> con;                   // vector of tpl positions of missing atoms
        
        // vector of possible atoms missing
        vector<pair<string, int> > miss;
        pair <string,int> h ("HO1", 16); //name of hydrogen and position in tpl
        miss.push_back(h);
        h = make_pair("O1", 15);
        miss.push_back(h);
        h = make_pair("HO2", 14);
        miss.push_back(h);
        h = make_pair("O2", 13);
        miss.push_back(h);
        h = make_pair("HO3", 12);
        miss.push_back(h);
        h = make_pair("O3", 11);
        miss.push_back(h);
        h = make_pair("HO4", 10);
        miss.push_back(h);
        h = make_pair("O4", 9);
        miss.push_back(h);
        h = make_pair("HO6",  8);
        miss.push_back(h);
        h = make_pair("O6", 7);
        miss.push_back(h);
        
        
        for (vector<atom>::iterator it = res.begin() ; it != res.end(); ++it){    //iterate over vector of residue atoms
            for (vector<std::pair<string, int> >::iterator itH = miss.begin() ; itH != miss.end(); itH++){ //iterate over vector of possible missing atoms
                if (itH->first.compare(it->a)==0){
                    miss.erase(itH);
                    itH--;
                }
            }
        }
        for (vector<std::pair<string, int> >::iterator itH = miss.begin() ; itH != miss.end(); ++itH){
                con.push_back((*itH).second + i*17);
        }
        skips.push_back(con);
    }
    return skips;
    
}

//gives back pair of residue number and atom number
vector<pair <int, int> > findConnect(vector<atom> cnf, int nrs){
    vector< pair <int, int> > con;
    vector<atom> res = getRes(cnf, nrs); //get a vector of atoms of a residue
    
    // vector of tpl positions of missing atoms
    // vector of possible atoms missing
    vector<pair<string, int> > miss;
    pair <string,int> h ("O2", 4); //name of oxygen and position of C atom in tpl
    miss.push_back(h);
    h = make_pair("O3", 3);
    miss.push_back(h);
    h = make_pair("O4", 2);
    miss.push_back(h);
    h = make_pair("O6", 6);
    miss.push_back(h);
    //O1 has to be in back, because it misses every time 
    //when it is not the last residue
    //it is kept here for 1-1 linking, then it is returned
    h = make_pair("O1", 5); 
    miss.push_back(h);


    for (vector<atom>::iterator it = res.begin() ; it != res.end(); ++it){    //iterate over vector of residue atoms
        for (vector<std::pair<string, int> >::iterator itH = miss.begin() ; itH != miss.end(); itH++){ //iterate over vector of possible missing atoms
            if (itH->first == it->a){
                miss.erase(itH);
                itH--;
            }
        }
    }
    
    pair <int, int> pcon;
    if (!miss.empty()){
        pcon.first = nrs;
        pcon.second = miss.front().second;
    }
    else {
        pcon.first = nrs;
        pcon.second = -1;
    }
    //pair <int, int> pcon(nrs,miss.front().second);
    
    con.push_back( pcon );
    return con;   
}
vector<pair <int, int> > findConnect(vector<atom> cnf){
    int nrs = cnf[cnf.size()-1].nrRes;
    vector<pair <int, int> > con;
    for (int i=1; i<nrs; i++){              //loop over all residues except first one
        vector<atom> res = getRes(cnf, i+1); //get a vector of atoms of a residue
        //printVec(res);
                           // vector of tpl positions of missing atoms
        
        // vector of possible atoms missing
        vector<pair<string, int> > miss;
        pair <string,int> h ("O2", 4); //name of oxygen and position of C atom in tpl
        miss.push_back(h);
        h = make_pair("O3", 3);
        miss.push_back(h);
        h = make_pair("O4", 2);
        miss.push_back(h);
        h = make_pair("O6", 6);
        miss.push_back(h);
        //O1 has to be in back, because it misses every time 
        //when it is not the last residue
        //it is kept here for 1-1 linking, then it is returned
        h = make_pair("O1", 5); 
        miss.push_back(h);
        
        
        for (vector<atom>::iterator it = res.begin() ; it != res.end(); ++it){    //iterate over vector of residue atoms
            for (vector<std::pair<string, int> >::iterator itH = miss.begin() ; itH != miss.end(); itH++){ //iterate over vector of possible missing atoms
                if (itH->first.compare(it->a)==0){
                    miss.erase(itH);
                    itH--;
                }
            }
        }
//        for (vector<std::pair<string, int> >::iterator itH = miss.begin() ; itH != miss.end(); ++itH){ //iterate over vector of missing atoms
//            pair <int, int> pcon(i+1,miss[0].second);
//            con.push_back( pcon );
//        }
        
        //just first missing "atom" is returned
        //see comment for ("O1", 5) above
        pair <int, int> pcon(i+1,miss.front().second);
        con.push_back( pcon );
    }
    return con;   
}


bool isIn(vector<atom> cnf, string str){
    for (vector<atom>::iterator it = cnf.begin() ; it != cnf.end(); ++it){
        if (str.compare((*it).a)==0) return true;
    }
    return false;
}

bool isIn(vector<vector<int> > skips, int i){
    for (vector<vector<int> >::iterator it = skips.begin() ; it != skips.end(); ++it){
        for (vector<int >::iterator it2 = (*it).begin() ; it2 != (*it).end(); ++it2){
            if (i==*it2) return true;
        }
    }
    return false;
}
bool isIn(vector<string> vec, string str){
    for (vector<string>::iterator it = vec.begin() ; it != vec.end(); ++it){
        if (str.compare(*it)==0) return true;
    }
    return false;
}

//creating cnf/angular matrix function
void fillTpl(vector<vector<string> > tpl, vector<atom> cnf, vector<atom > &angle){
    angle.resize(cnf.size());
    int nrs = cnf[cnf.size()-1].nrRes; //number of residues
    vector<vector<int> > skips = findSkpis(cnf, nrs); //vector for skips in each residue
    
    int j = 0; // position in angular file, keeps tract with right residue
    for(int i=0; i<(17*nrs); ++i, ++j){
        if (isIn(skips,i)){
            j--;
            continue;
        }
        vec a, b, c, d;
        atom temp;

            if (i%17 == 0 && i != 0){ //in case of new residue
                int ar = i/17; //actual residue-1
                
                int r; //ring position
                if (isIn(skips,(ar)*17+7)) { //TO DO this is not fully implemented!!!
                    r = 6; //if 06 in skips build ring from 7rd atom in tpl - C6 - but it is not in the ring!
                }
                else if (isIn(skips,(ar)*17+9)) {
                    r = 2; //if 04 in skips build ring from 3rd atom in tpl - C4
                }
                else if (isIn(skips,(ar)*17+11)) {
                    r = 3; //if 03 in skips build ring from 4th atom in tpl - C3
                }
                else if (isIn(skips,(ar)*17+13)) {
                    r = 4; //if 02 in skips build ring from 5th atom in tpl - C2
                }
                else if (isIn(skips,(ar)*17+15)) {
                    r = 5; //if 01 in skips build ring from 6th atom in tpl - C1
                }
                else cout << "conection to residue " << ar << " not found" << endl;
                
                c = findAtom(tpl[16][1], cnf, ar ).pos;
                b = findAtom(tpl[16][2], cnf, ar ).pos;
                a = findAtom(tpl[16][3], cnf, ar ).pos;
                
                for (int g = 0; g < 6; g++, r++){ // loop over 6 three ring atoms anticlockwise
                    // i+r is position of solving atom in angleIn
                    
                    temp = findAtom(tpl[r%6][0], cnf, ar+1 ); // i/17);
                    
                    angle[j+r%6].nrRes = temp.nrRes;
                    angle[j+r%6].res = temp.res;
                    angle[j+r%6].a = temp.a;
                    angle[j+r%6].nrA = temp.nrA;
                    d = temp.pos;
                    
                    angle[j+r%6].pos.x = calcDistance(d,c);
                    angle[j+r%6].pos.y = calcAngle(d,c,b);
                    angle[j+r%6].pos.z = calcDihedral(d,c,b,a);
                    
                    a=b;
                    b=c;
                    c=d;
                }
                //skips no included because we stay in the ring and there are no missing atoms
                j=j+5;
                i=i+5;
                continue;
            }
            else {
                temp = findAtom(tpl[i%17][0], cnf, i/17+1 ); 
                angle[j].nrRes = temp.nrRes;
                angle[j].res = temp.res;
                angle[j].a = temp.a;
                angle[j].nrA = temp.nrA;
                d = temp.pos;
                c = findAtom(tpl[i%17][1], cnf, i/17+1 ).pos;
                b = findAtom(tpl[i%17][2], cnf, i/17+1 ).pos;
                a = findAtom(tpl[i%17][3], cnf, i/17+1 ).pos;
            }
            angle[j].pos.x = calcDistance(d,c);
            angle[j].pos.y = calcAngle(d,c,b);
            angle[j].pos.z = calcDihedral(d,c,b,a);
        }
    
}

// Peter C. Kahn 1989 - lars:KA89.4
// 
helix getHelix(vector<atom> cnf, int start, map<string,LJparam> LJpairs){
    helix h;
    atom Pa1 = findAtom("O1", cnf, start-1);
    atom Pa2 = findAtom("O1", cnf, start);
    atom Pa3 = findAtom("O1", cnf, start+1);
    vec Aa = getSubstrVec(Pa2.pos, Pa1.pos);
    vec Ba = getSubstrVec(Pa2.pos, Pa3.pos);
    vec Va = getSumVec(Aa, Ba);
    Va = getMultipVec(Va,  1/getAbsVec(Va));
    //printAtom(Pa1);
//    printAtom(Pa2);
//    printAtom(Pa3);
//    printVec(Aa);
//    printVec(Ba);
//    cout << "Aa " << getAbsVec(Aa) << " Ba " << getAbsVec(Ba) << endl;
//    cout << "Va " << getAbsVec(Va) <<  endl;
//    printVec(Va);

    atom Pb1 = findAtom("O1", cnf, start);
    atom Pb2 = findAtom("O1", cnf, start+1);
    atom Pb3 = findAtom("O1", cnf, start+2);
    vec Ab = getSubstrVec(Pb2.pos, Pb1.pos);
    vec Bb = getSubstrVec(Pb2.pos, Pb3.pos);
    vec Vb = getSumVec(Ab, Bb);
    Vb = getMultipVec(Vb,  1/getAbsVec(Vb));
    //printAtom(Pb1);
    //printAtom(Pb2);
    //printAtom(Pb3);
    //printVec(Ab);
    //printVec(Bb);
    //cout << "Ab " << getAbsVec(Ab) << " Bb " << getAbsVec(Bb) << endl;
    //cout << "Vb " << getAbsVec(Vb) <<  endl;
    //printVec(Vb);
    
    //axis
    vec H = getCrossofVec(Va, Vb);
    H = getMultipVec(H, 1/getAbsVec(H));
    //printVec(H);

    vec SUB = getSubstrVec(Pa2.pos,Pb2.pos);
    
    // rise per residue
    double a = getDotofVec(SUB,H);
    h.a = a;
    if (a<0.0) { // residues are making over 180deg therefore H is miss oriented
        vec H2 = getCrossofVec(Vb, Va);
        H2 = getMultipVec(H2, 1/getAbsVec(H2));
        a = getDotofVec(SUB,H2);
    }
    h.a = a;
    //printf("rise per res %f \n", a);
        
    // radius
    vec SUB2 = getSubstrVec(Pb2.pos,Pa2.pos);
    //the minus is added empiricaly, It is not according to the paper
    double r = -( (a*a*getAbs2Vec(H)) - getAbs2Vec(SUB) ) / (2*abs(getDotofVec(SUB2,Vb))); 
    h.R = r;
    //cout << "r " << r <<  endl;
    
    // phase angle angular rotation per residue APR
    // can be used just for successive residues 0-360deg;
    // http://math.stackexchange.com/questions/878785/how-to-find-an-angle-in-range0-360-between-2-vectors
    double dot = getDotofVec(Va,Vb);
    //H is changing orientation with Va and Vb (left/right handed) -> det is always positive
    double det = getDotofVec(H,getCrossofVec(Va, Vb));
    double gamma = atan2(det,dot)*180/PI;
    
    //check left and right handing
    //g- without consequences
//    if(getDotofVec(SUB,H)<0) {
//        h.g = -gamma;
//    } else h.g = gamma;
    
    //g- consequences
    if(getDotofVec(SUB,H)<0) { 
        gamma = -gamma;
    } 
    h.g = gamma;
    //printf("angle per res %f\n", gamma);
    
    //number of residues per turn
    double RPT = 360/gamma;
    h.n = RPT;
    //printf("residues per turn %f\n", RPT);
    
    // pitch
    double p = RPT*a;
    h.p = p;
    //printf("pitch %f\n ----- \n", p);
    
    
    ///////////////////////////////////////
    //clash calculation
    ///////////////////////////////////////
    vector<atom> res1 = getRes(cnf,1);
    vector<atom> res2 = getRes(cnf,2);
    
    vector<string > clashList {"C1","C2","C3","C4","C5","C6","O1","O2","O3","O4","O5"};
        
    h.clash = getClash(res1, res2, clashList, LJpairs);
    
    return h;
}

bool isHbond(atom accept, atom H, atom donor, double Hdistmin, double Hdistmax, double Hangle){
    bool bond = false;
    
    double distance = calcDistance(accept.pos,H.pos);
    double angle = calcAngle(accept.pos, H.pos, donor.pos);
    

    if ( Hdistmin < distance && distance < Hdistmax && Hangle < angle) bond = true;
    
    return bond;
}

vector<pair<string, int> > getHbond(vector<atom> angle, vector<vector<string> > tpl, 
        double Hdistmin, double Hdistmax, double Hangle, vector<string > Hlist,vector<string > acceptList){
    //get just first two residues, better for later angle-> cnf conversion
    vector<atom> one = getRes(angle, 1);
    vector<atom> two = getRes(angle, 2);
    angle.clear();
    angle.reserve( one.size() + two.size());
    angle.insert(angle.end(),one.begin(),one.end() );
    angle.insert(angle.end(),two.begin(),two.end() );
    
    vector<pair<string, int> > Hbond;
    
    pair<string, int> HBcumul;
    HBcumul.first= "cumulative";
    HBcumul.second=0;
    
    //3 cases when combination of H bond is not possible
    //1. case H can access 2 acceptors
    //2. case two acceptor/donor groups can interchange (flip-flop) 
    //3. case (?) can three acceptor/donor groups form three Hbonds?
    
    vector<string> Hdonor1; //holds H which are already in hbonds - to prevent 1.case
    vector<string> Hdonor2; //holds Hbond group pairs - to prevent 2.case
    
    pair<string, int> HBcompat1; //omits just 1.case
    HBcompat1.first= "compatible1";
    HBcompat1.second=0;

    pair<string, int> HBcompat2; //omits just 2.case
    HBcompat2.first= "compatible2";
    HBcompat2.second=0;
    
    pair<string, int> HBcompat; //omits 1. and 2. case
    HBcompat.first= "compatible";
    HBcompat.second=0;
    
    vector<atom> cnf;
    if (fillCnf(angle, tpl, cnf)==1){
        cerr << "function getHbond(), error in fillCnf()" << endl;
        return Hbond;
    }
    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    //for res1 as acceptor and res2 as donor 
    //////////////////////////////////////////////////////////////////////////////////////
    int ac_res = 1;
    int do_res = 2;
    
    for (vector<string>::iterator ita = acceptList.begin(); ita != acceptList.end(); ++ita){
        //acceptor atom in Cartesian coordinates
        atom accept = findAtom(*ita, cnf, ac_res); 
        if(accept.nrA==0) continue; //if the atom doesn't exist go next
        for (vector<string>::iterator ith = Hlist.begin(); ith != Hlist.end(); ++ith){
            //H atom in Cartesian coordinates
            atom H = findAtom(*ith, cnf, do_res);
            if(H.nrA==0) continue; //if the atom doesn't exist go next
            
            ///////////////////////////
            //H bond pair defining name and presence of Hbond
            pair<string, int> Hpair;
            Hpair.first= "1" + *ita + "-2" + *ith;
            Hpair.second=0;
            Hbond.push_back(Hpair);
            ///////////////////////////
            //check if H bonding is even possible
            //calculate distance of acceptor to H distance
            double dis0 = calcDistance(accept.pos,H.pos);          
            if(dis0>(Hdistmax+0.19) && H.a.compare("HO6")!=0) continue; // 0.1885 is diameter of circle made by OH rotation
            
            ///////////////////////////
            //Lets search for Hbond
            ///////////////////////////
            //set some variables for H
            //donor atom of H in Cartesian coordinates
            atom donor = findAtom((*ith).substr(1,2), cnf, do_res);
            //flag for changing right chi angle
            string ch = to_string(static_cast<long long int>(do_res))+" ch"+(*ith).substr(2,1)+" "; 
            //initial angle step
            double angleStep = 120; 
            //chi angle value
            long double chAngle = findAtom(*ith, angle, do_res).pos.z;
            string s_chAngle;
            ////////////////////////////
            
            ////////////////////////////
            //if H is HO6, optimize omega tilde angle
            if(H.a.compare("HO6")==0){
                ////////////////////
                //set some variables
                //C6 atom
                atom C6 = findAtom("C6", cnf, do_res);
                //calculate angle of C6-donor-acceptor - this we want to make close to 109.5
                double ang0 = calcAngle(C6.pos,donor.pos,accept.pos);
                //difference between current and last ang0
                double angDiff = abs(109.5 - calcAngle(C6.pos,donor.pos,accept.pos));
                //flag for changing omega tilde angle
                string ot = to_string(static_cast<long long int>(do_res))+" ot "; 
                //initial ot angle step
                double otStep = 120; 
                //ot angle value
                long double otAngle = findAtom(donor.a, angle, do_res).pos.z;
                string s_otAngle;
                
                ////////////////////
                //now lets optimize
                while(otStep < -10 || otStep > 10 ){
                    //if we are getting further -> change direction and make smaller step
                    if (angDiff < abs(109.5 - calcAngle(C6.pos,donor.pos,accept.pos ) ) )
                        otStep = -otStep/2;
                    angDiff = abs( 109.5 - calcAngle(C6.pos,donor.pos,accept.pos) );

                    //get new angle
                    otAngle += angleStep;
                    s_otAngle = ot+std::to_string(otAngle);

                    //change angle
                    int state = changeAngle(angle, tpl, "rot", s_otAngle);
                    if(state>0) {
                        printf("problem in H bond analysis changing ot angle");
                        return Hbond;
                    }

                    //calculate new cnf
                    cnf.clear();
                    if (fillCnf(angle, tpl, cnf)==1){
                        cerr << "function getHbond(), error in fillCnf()" << endl;
                        return Hbond;
                    }
                    // get new donor
                    donor = findAtom("O6", cnf, do_res);
                    
                }
                H = findAtom(*ith, cnf, do_res);
                //cout << "angle: " << calcAngle(C6.pos,donor.pos,accept.pos) << endl;
            }
            ////////////////////////////
            
            //now continue searching for H bond changing chi angle
            while(angleStep < -10 || angleStep > 10 ){
                //do we have H bond -> increase Hbond and start with new pair
                if(isHbond(accept, H, donor, Hdistmin, Hdistmax, Hangle)){
                    Hbond.back().second += 1; //increase current H bond type
                    HBcumul.second += 1;  //increase cumulative H bond
                    
                    //create flag for this Hbond - 2nd case
                    stringstream ss;
                    string pairCode;//number of group on first residue + number of group on second residue
                    char tmp = (*ita)[(*ita).length()-1];
                    ss << tmp;
                    tmp = (*ith)[(*ith).length()-1];
                    ss << tmp;
                    ss >> pairCode; 
                    //omits 1. and 2. case (see above the definitions)
                    if (!isIn(Hdonor2, pairCode ) && !isIn(Hdonor1, "2"+*ith ) ) {
                        HBcompat.second += 1;
                    }
                    //omits just 1. case (see above the definitions)
                    if (!isIn(Hdonor1, "2"+*ith ) ) {
                        Hdonor1.push_back("2"+*ith);
                        HBcompat1.second += 1;
                    }
                    //omits just 2. case (see above the definitions)
                    if (!isIn(Hdonor2, pairCode ) ) {
                        Hdonor2.push_back(pairCode);
                        HBcompat2.second += 1;
                    }
                                       
                    break;
                } else {
                    //if we are getting further -> change direction and make smaller step
                    if (dis0 < calcDistance(accept.pos,H.pos)) angleStep = -angleStep/2;
                    dis0 = calcDistance(accept.pos,H.pos);

                    //get new angle
                    chAngle +=  angleStep;
                    s_chAngle = ch+std::to_string(chAngle);

                    //change angle
                    int state = changeAngle(angle, tpl, "rot", s_chAngle);
                    if(state>0) {
                        printf("problem in H bond analysis changing chi angle");
                        return Hbond;
                    }

                    //calculate new cnf
                    cnf.clear();
                    if (fillCnf(angle, tpl, cnf)==1){
                        cerr << "function getHbond(), error in fillCnf()" << endl;
                        return Hbond;
                    }
                    // get new H
                    H = findAtom(*ith, cnf, do_res);
                    //to get right orientation for the first step
                    if (angleStep == 120 && dis0 < calcDistance(accept.pos,H.pos)) angleStep *= 2;
                    
                }
            }
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    //for res2 as acceptor and res1 as donor 
    //////////////////////////////////////////////////////////////////////////////////////
    ac_res = 2;
    do_res = 1;
    
    
    for (vector<string>::iterator ita = acceptList.begin(); ita != acceptList.end(); ++ita){
        atom accept = findAtom(*ita, cnf, ac_res);
        if(accept.nrA==0) continue; //if the atom doesn't exist go next
        for (vector<string>::iterator ith = Hlist.begin(); ith != Hlist.end(); ++ith){
            atom H = findAtom(*ith, cnf, do_res);
            if(H.nrA==0) continue; //if the atom doesn't exist go next
            
            ///////////////////////////
            pair<string, int> Hpair;
            Hpair.first= "2" + *ita + "-1" + *ith;
            Hpair.second=0;
            Hbond.push_back(Hpair);
            ///////////////////////////
            //check if H bonding is even possible
            //calculate distance of acceptor to H distance
            double dis0 = calcDistance(accept.pos,H.pos);          
            if(dis0>(Hdistmax+0.19) && H.a.compare("HO6")!=0) continue; // 0.1885 is diameter of circle made by OH rotation
            
            ///////////////////////////
            //Lets search for Hbond
            ///////////////////////////
            //set some variables for H
            atom donor = findAtom((*ith).substr(1,2), cnf, do_res);
            string ch = to_string(static_cast<long long int>(do_res))+" ch"+(*ith).substr(2,1)+" "; //flag for changing right angle
            double angleStep = 120; //initial angle step
            long double chAngle = findAtom(*ith, angle, do_res).pos.z;
            string s_chAngle;
            int count = 0;
            ////////////////////////////
            
            ////////////////////////////
            //if H is HO6, optimize omega tilde angle
            if(H.a.compare("HO6")==0){
                ////////////////////
                //set some variables
                //C6 atom
                atom C6 = findAtom("C6", cnf, do_res);
                //calculate angle of C6-donor-acceptor - this we want to make close to 109.5
                double ang0 = calcAngle(C6.pos,donor.pos,accept.pos);
                //difference between current and last ang0
                double angDiff = abs(109.5 - calcAngle(C6.pos,donor.pos,accept.pos));
                //flag for changing omega tilde angle
                string ot = to_string(static_cast<long long int>(do_res))+" ot "; 
                //initial ot angle step
                double otStep = 120; 
                //ot angle value
                long double otAngle = findAtom(donor.a, angle, do_res).pos.z;
                string s_otAngle;
                
                ////////////////////
                //now lets optimize
                while(otStep < -10 || otStep > 10 ){
                    //if we are getting further -> change direction and make smaller step
                    if (angDiff < abs(109.5 - calcAngle(C6.pos,donor.pos,accept.pos ) ) )
                        otStep = -otStep/2;
                    angDiff = abs( 109.5 - calcAngle(C6.pos,donor.pos,accept.pos) );

                    //get new angle
                    otAngle += angleStep;
                    s_otAngle = ot+std::to_string(otAngle);

                    //change angle
                    int state = changeAngle(angle, tpl, "rot", s_otAngle);
                    if(state>0) {
                        printf("problem in H bond analysis changing ot angle");
                        return Hbond;
                    }

                    //calculate new cnf
                    cnf.clear();
                    if (fillCnf(angle, tpl, cnf)==1){
                        cerr << "function getHbond(), error in fillCnf()" << endl;
                        return Hbond;
                    }
                    // get new donor
                    donor = findAtom("O6", cnf, do_res);
                    
                }
                H = findAtom(*ith, cnf, do_res);
            }
            ////////////////////////////
            
            //now continue searching for H bond changing chi angle
            while(angleStep < -10 || angleStep > 10 ){
                count += 1;
                //do we have H bond -> increase Hbond and start with new pair
                if(isHbond(accept, H, donor, Hdistmin, Hdistmax, Hangle)){
                    Hbond.back().second += 1;
                    HBcumul.second += 1;
                    stringstream ss;
                    string pairCode; //number of group on first residue + number of group on second residue
                    char tmp = (*ith)[(*ith).length()-1];
                    ss << tmp;
                    tmp = (*ita)[(*ita).length()-1];
                    ss << tmp;
                    ss >> pairCode;
                    if (!isIn(Hdonor2, pairCode ) && !isIn(Hdonor1, "1"+*ith ) ) {
                        HBcompat.second += 1;
                    }
                    if (!isIn(Hdonor1, "1"+*ith ) ) {
                        Hdonor1.push_back("1"+*ith);
                        HBcompat1.second += 1;
                    }
                    if (!isIn(Hdonor2, pairCode ) ) {
                        Hdonor2.push_back(pairCode);
                        HBcompat2.second += 1;
                    }
                    break;
                } else {
                    //if we are getting further -> change direction and make smaller step
                    if (dis0 < calcDistance(accept.pos,H.pos)) angleStep = -angleStep/2;
                    dis0 = calcDistance(accept.pos,H.pos);

                    //get new angle
                    chAngle +=  angleStep;
                    s_chAngle = ch+std::to_string(chAngle);

                    //change angle
                    int state = changeAngle(angle, tpl, "rot", s_chAngle);
                    if(state>0) {
                        printf("problem in H bond analysis changing chi angle");
                        return Hbond;
                    }

                    //calculate new cnf
                    cnf.clear();
                    if (fillCnf(angle, tpl, cnf)==1){
                        cerr << "function getHbond(), error in fillCnf()" << endl;
                        return Hbond;
                    }
                    // get new H
                    H = findAtom(*ith, cnf, do_res);
                    //to get right orientation for the first step
                    if (angleStep == 120 && dis0 < calcDistance(accept.pos,H.pos)) angleStep *= 2;
                    
                }
            }
        }
    }
    
    Hbond.push_back(HBcumul);
    Hbond.push_back(HBcompat1);
    Hbond.push_back(HBcompat2);
    Hbond.push_back(HBcompat);
    
    return Hbond;
}

vector<pair<string, int> > getHbond(vector<atom> angle, int firstR, int secondR, vector<vector<string> > tpl, 
        double Hdistmin, double Hdistmax, double Hangle, vector<string > Hlist,vector<string > acceptList){

    //printVec(angle);
    vector<pair<string, int> > Hbond;
    
    pair<string, int> HBcumul;
    HBcumul.first= "cumulative";
    HBcumul.second=0;
    
    //3 cases when combination of H bond is not possible
    //1. case H can access 2 acceptors
    //2. case two acceptor/donor groups can interchange (flip-flop) 
    //3. case (?) can three acceptor/donor groups form three Hbonds?
    
    vector<string> Hdonor1; //holds H which are already in hbonds - to prevent 1.case
    vector<string> Hdonor2; //holds Hbond group pairs - to prevent 2.case
    
    pair<string, int> HBcompat1; //omits just 1.case
    HBcompat1.first= "compatible1";
    HBcompat1.second=0;

    pair<string, int> HBcompat2; //omits just 2.case
    HBcompat2.first= "compatible2";
    HBcompat2.second=0;
    
    pair<string, int> HBcompat; //omits 1. and 2. case
    HBcompat.first= "compatible";
    HBcompat.second=0;
    
    vector<atom> cnf;
    if (fillCnf(angle, tpl, cnf)==1){
        cerr << "function getHbond(), error in fillCnf()" << endl;
        return Hbond;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    //for res1 as acceptor and res2 as donor 
    //////////////////////////////////////////////////////////////////////////////////////
    int ac_res = firstR;
    int do_res = secondR;
    
    for (vector<string>::iterator ita = acceptList.begin(); ita != acceptList.end(); ++ita){
        //acceptor atom in Cartesian coordinates
        atom accept = findAtom(*ita, cnf, ac_res); 
        if(accept.nrA==0) continue; //if the atom doesn't exist go next
        for (vector<string>::iterator ith = Hlist.begin(); ith != Hlist.end(); ++ith){
            //H atom in Cartesian coordinates
            atom H = findAtom(*ith, cnf, do_res);
            if(H.nrA==0) continue; //if the atom doesn't exist go next
            
            ///////////////////////////
            //H bond pair defining name and presence of Hbond
            pair<string, int> Hpair;
            Hpair.first= "1" + *ita + "-2" + *ith;
            Hpair.second=0;
            Hbond.push_back(Hpair);
            ///////////////////////////
            // if secondR==2 then HbondHelix is same as in Hbond
            if (secondR==2 ) continue;
            ///////////////////////////
            //check if H bonding is even possible
            //calculate distance of acceptor to H distance
            double dis0 = calcDistance(accept.pos,H.pos);          
            if(dis0>(Hdistmax+0.19) && H.a.compare("HO6")!=0) continue; // 0.1885 is diameter of circle made by OH rotation
            
            ///////////////////////////
            //Lets search for Hbond
            ///////////////////////////
            //set some variables for H
            //donor atom of H in Cartesian coordinates
            atom donor = findAtom((*ith).substr(1,2), cnf, do_res);
            //flag for changing right chi angle
            string ch = to_string(static_cast<long long int>(do_res))+" ch"+(*ith).substr(2,1)+" "; 
            //initial angle step
            double angleStep = 120; 
            //chi angle value
            long double chAngle = findAtom(*ith, angle, do_res).pos.z;
            string s_chAngle;
            ////////////////////////////
            
            ////////////////////////////
            //if H is HO6, optimize omega tilde angle
            if(H.a.compare("HO6")==0){
                ////////////////////
                //set some variables
                //C6 atom
                atom C6 = findAtom("C6", cnf, do_res);
                //calculate angle of C6-donor-acceptor - this we want to make close to 109.5
                double ang0 = calcAngle(C6.pos,donor.pos,accept.pos);
                //difference between current and last ang0
                double angDiff = abs(109.5 - calcAngle(C6.pos,donor.pos,accept.pos));
                //flag for changing omega tilde angle
                string ot = to_string(static_cast<long long int>(do_res))+" ot "; 
                //initial ot angle step
                double otStep = 120; 
                //ot angle value
                long double otAngle = findAtom(donor.a, angle, do_res).pos.z;
                string s_otAngle;
                
                ////////////////////
                //now lets optimize
                while(otStep < -10 || otStep > 10 ){
                    //if we are getting further -> change direction and make smaller step
                    if (angDiff < abs(109.5 - calcAngle(C6.pos,donor.pos,accept.pos ) ) )
                        otStep = -otStep/2;
                    angDiff = abs( 109.5 - calcAngle(C6.pos,donor.pos,accept.pos) );

                    //get new angle
                    otAngle += angleStep;
                    s_otAngle = ot+std::to_string(otAngle);

                    //change angle
                    int state = changeAngle(angle, tpl, "rot", s_otAngle);
                    if(state>0) {
                        printf("problem in H bond analysis changing ot angle");
                        return Hbond;
                    }

                    //calculate new cnf
                    cnf.clear();
                    if (fillCnf(angle, tpl, cnf)==1){
                        cerr << "function getHbond(), error in fillCnf()" << endl;
                        return Hbond;
                    }
                    // get new donor
                    donor = findAtom("O6", cnf, do_res);
                    
                }
                H = findAtom(*ith, cnf, do_res);
                //cout << "angle: " << calcAngle(C6.pos,donor.pos,accept.pos) << endl;
            }
            ////////////////////////////
            
            //now continue searching for H bond changing chi angle
            while(angleStep < -10 || angleStep > 10 ){
                //do we have H bond -> increase Hbond and start with new pair
                if(isHbond(accept, H, donor, Hdistmin, Hdistmax, Hangle)){
                    Hbond.back().second += 1; //increase current H bond type
                    HBcumul.second += 1;  //increase cumulative H bond
                    
                    //create flag for this Hbond - 2nd case
                    stringstream ss;
                    string pairCode;//number of group on first residue + number of group on second residue
                    char tmp = (*ita)[(*ita).length()-1];
                    ss << tmp;
                    tmp = (*ith)[(*ith).length()-1];
                    ss << tmp;
                    ss >> pairCode; 
                    //omits 1. and 2. case (see above the definitions)
                    if (!isIn(Hdonor2, pairCode ) && !isIn(Hdonor1, "2"+*ith ) ) {
                        HBcompat.second += 1;
                    }
                    //omits just 1. case (see above the definitions)
                    if (!isIn(Hdonor1, "2"+*ith ) ) {
                        Hdonor1.push_back("2"+*ith);
                        HBcompat1.second += 1;
                    }
                    //omits just 2. case (see above the definitions)
                    if (!isIn(Hdonor2, pairCode ) ) {
                        Hdonor2.push_back(pairCode);
                        HBcompat2.second += 1;
                    }
                                       
                    break;
                } else {
                    //if we are getting further -> change direction and make smaller step
                    if (dis0 < calcDistance(accept.pos,H.pos)) angleStep = -angleStep/2;
                    dis0 = calcDistance(accept.pos,H.pos);

                    //get new angle
                    chAngle +=  angleStep;
                    s_chAngle = ch+std::to_string(chAngle);

                    //change angle
                    int state = changeAngle(angle, tpl, "rot", s_chAngle);
                    if(state>0) {
                        printf("problem in H bond analysis changing chi angle");
                        return Hbond;
                    }

                    //calculate new cnf
                    cnf.clear();
                    if (fillCnf(angle, tpl, cnf)==1){
                        cerr << "function getHbond(), error in fillCnf()" << endl;
                        return Hbond;
                    }
                    // get new H
                    H = findAtom(*ith, cnf, do_res);
                    //to get right orientation for the first step
                    if (angleStep == 120 && dis0 < calcDistance(accept.pos,H.pos)) angleStep *= 2;
                    
                }
            }
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    //for res2 as acceptor and res1 as donor 
    //////////////////////////////////////////////////////////////////////////////////////
    ac_res = secondR;
    do_res = firstR;
    
    
    for (vector<string>::iterator ita = acceptList.begin(); ita != acceptList.end(); ++ita){
        atom accept = findAtom(*ita, cnf, ac_res);
        if(accept.nrA==0) continue; //if the atom doesn't exist go next
        for (vector<string>::iterator ith = Hlist.begin(); ith != Hlist.end(); ++ith){
            atom H = findAtom(*ith, cnf, do_res);
            if(H.nrA==0) continue; //if the atom doesn't exist go next
            
            ///////////////////////////
            pair<string, int> Hpair;
            Hpair.first= "2" + *ita + "-1" + *ith;
            Hpair.second=0;
            Hbond.push_back(Hpair);
            ///////////////////////////
            // if secondR==2 then HbondHelix is same as in Hbond
            if (secondR==2 ) continue;
            ///////////////////////////
            //check if H bonding is even possible
            //calculate distance of acceptor to H distance
            double dis0 = calcDistance(accept.pos,H.pos);          
            if(dis0>(Hdistmax+0.19) && H.a.compare("HO6")!=0) continue; // 0.1885 is diameter of circle made by OH rotation
            
            ///////////////////////////
            //Lets search for Hbond
            ///////////////////////////
            //set some variables for H
            atom donor = findAtom((*ith).substr(1,2), cnf, do_res);
            string ch = to_string(static_cast<long long int>(do_res))+" ch"+(*ith).substr(2,1)+" "; //flag for changing right angle
            double angleStep = 120; //initial angle step
            long double chAngle = findAtom(*ith, angle, do_res).pos.z;
            string s_chAngle;
            int count = 0;
            ////////////////////////////
            
            ////////////////////////////
            //if H is HO6, optimize omega tilde angle
            if(H.a.compare("HO6")==0){
                ////////////////////
                //set some variables
                //C6 atom
                atom C6 = findAtom("C6", cnf, do_res);
                //calculate angle of C6-donor-acceptor - this we want to make close to 109.5
                double ang0 = calcAngle(C6.pos,donor.pos,accept.pos);
                //difference between current and last ang0
                double angDiff = abs(109.5 - calcAngle(C6.pos,donor.pos,accept.pos));
                //flag for changing omega tilde angle
                string ot = to_string(static_cast<long long int>(do_res))+" ot "; 
                //initial ot angle step
                double otStep = 120; 
                //ot angle value
                long double otAngle = findAtom(donor.a, angle, do_res).pos.z;
                string s_otAngle;
                
                ////////////////////
                //now lets optimize
                while(otStep < -10 || otStep > 10 ){
                    //if we are getting further -> change direction and make smaller step
                    if (angDiff < abs(109.5 - calcAngle(C6.pos,donor.pos,accept.pos ) ) )
                        otStep = -otStep/2;
                    angDiff = abs( 109.5 - calcAngle(C6.pos,donor.pos,accept.pos) );

                    //get new angle
                    otAngle += angleStep;
                    s_otAngle = ot+std::to_string(otAngle);

                    //change angle
                    int state = changeAngle(angle, tpl, "rot", s_otAngle);
                    if(state>0) {
                        printf("problem in H bond analysis changing ot angle");
                        return Hbond;
                    }

                    //calculate new cnf
                    cnf.clear();
                    if (fillCnf(angle, tpl, cnf)==1){
                        cerr << "function getHbond(), error in fillCnf()" << endl;
                        return Hbond;
                    }
                    // get new donor
                    donor = findAtom("O6", cnf, do_res);
                    
                }
                H = findAtom(*ith, cnf, do_res);
            }
            ////////////////////////////
            
            //now continue searching for H bond changing chi angle
            while(angleStep < -10 || angleStep > 10 ){
                count += 1;
                //do we have H bond -> increase Hbond and start with new pair
                if(isHbond(accept, H, donor, Hdistmin, Hdistmax, Hangle)){
                    Hbond.back().second += 1;
                    HBcumul.second += 1;
                    stringstream ss;
                    string pairCode; //number of group on first residue + number of group on second residue
                    char tmp = (*ith)[(*ith).length()-1];
                    ss << tmp;
                    tmp = (*ita)[(*ita).length()-1];
                    ss << tmp;
                    ss >> pairCode;
                    if (!isIn(Hdonor2, pairCode ) && !isIn(Hdonor1, "1"+*ith ) ) {
                        HBcompat.second += 1;
                    }
                    if (!isIn(Hdonor1, "1"+*ith ) ) {
                        Hdonor1.push_back("1"+*ith);
                        HBcompat1.second += 1;
                    }
                    if (!isIn(Hdonor2, pairCode ) ) {
                        Hdonor2.push_back(pairCode);
                        HBcompat2.second += 1;
                    }
                    break;
                } else {
                    //if we are getting further -> change direction and make smaller step
                    if (dis0 < calcDistance(accept.pos,H.pos)) angleStep = -angleStep/2;
                    dis0 = calcDistance(accept.pos,H.pos);

                    //get new angle
                    chAngle +=  angleStep;
                    s_chAngle = ch+std::to_string(chAngle);

                    //change angle
                    int state = changeAngle(angle, tpl, "rot", s_chAngle);
                    if(state>0) {
                        printf("problem in H bond analysis changing chi angle");
                        return Hbond;
                    }

                    //calculate new cnf
                    cnf.clear();
                    if (fillCnf(angle, tpl, cnf)==1){
                        cerr << "function getHbond(), error in fillCnf()" << endl;
                        return Hbond;
                    }
                    // get new H
                    H = findAtom(*ith, cnf, do_res);
                    //to get right orientation for the first step
                    if (angleStep == 120 && dis0 < calcDistance(accept.pos,H.pos)) angleStep *= 2;
                    
                }
            }
        }
    }
    
    Hbond.push_back(HBcumul);
    Hbond.push_back(HBcompat1);
    Hbond.push_back(HBcompat2);
    Hbond.push_back(HBcompat);
    
    if (secondR==2 ) {
        for (vector< pair<string, int> >::iterator it = Hbond.begin(); it != Hbond.end(); ++it) {
            it->second = -1;
        }
    }
    
    return Hbond;
}


double getLJ(atom i, atom j, double C12, double C6){
    double LJpot;
    double r = calcDistance(i.pos,j.pos);
    double r6 = r*r*r;
    r6 *= r6; //r6
    double r12 = r6*r6; //r12
    LJpot= C12/r12 - C6/r6;
    return LJpot;
}

int findLink(vector<atom> cnf){
    if (!isIn(cnf, "HO6") ){   
        return 6;
    }else if (!isIn(cnf, "HO4") ){   
        return 4;
    }else if (!isIn(cnf, "HO3") ){   
        return 3;
    }else if (!isIn(cnf, "HO2") ){   
        return 2;
    }else if (!isIn(cnf, "HO1") ){   
        return 1;
    }
    return 0;
}

double getTorsionEnergy(atom i, atom j, atom k, atom l, double K, double s, double m){
    double a = calcDihedral(i.pos,j.pos,k.pos,l.pos) / 180 * PI;
    return K * (1 + cos(m * a - s));
}


double getSteElecE_just_LJ(vector<atom> res1, vector<atom> res2, long long int link,  map<string,LJparam> LJpairs){
    double E = 0.0;
    //names are given like it would be for 1-4 linking
    atom O5, C2, C1, O1; //atoms of residue 1
    atom C4, C3, C5; //atoms of residue 2
    O5 = findAtom("O5", res1, 1);
    C2 = findAtom("C2", res1, 1);
    C1 = findAtom("C1", res1, 1);
    O1 = findAtom("O1", res1, 1);
    
    C4 = findAtom("C"+std::to_string(link), res2, 2);
    
    if(link==1) C3 =findAtom("O5", res2, 2);
    else C3 = findAtom("C"+std::to_string(link-1), res2, 2);
    
    C5 = findAtom("C"+std::to_string(link+1), res2, 2); //not used for 1-6
    
    
    //////////////////////////
    //CALCULATE LJ
    //////////////////////////
    if (link==1){
        E += getLJ(O5, C3, LJpairs["OrOr"].C12, LJpairs["OrOr"].C6);
        E += getLJ(C2, C3, LJpairs["OrCr"].C12, LJpairs["OrCr"].C6);
        E += getLJ(C1, C3, LJpairs["OrCr"].CS12, LJpairs["OrCr"].CS6);
    } else {
        E += getLJ(O5, C3, LJpairs["OrCr"].C12, LJpairs["OrCr"].C6);
        E += getLJ(C2, C3, LJpairs["CrCr"].C12, LJpairs["CrCr"].C6);
        E += getLJ(C1, C3, LJpairs["CrCr"].CS12, LJpairs["CrCr"].CS6);
    }
    
    if (link!=6){
        E += getLJ(O5, C4, LJpairs["OrCr"].CS12, LJpairs["OrCr"].CS6);
        E += getLJ(C2, C4, LJpairs["CrCr"].CS12, LJpairs["CrCr"].CS6);
    
        E += getLJ(O5, C5, LJpairs["OrCr"].C12, LJpairs["OrCr"].C6);
        E += getLJ(C2, C5, LJpairs["CrCr"].C12, LJpairs["CrCr"].C6);
        E += getLJ(C1, C5, LJpairs["CrCr"].CS12, LJpairs["CrCr"].CS6);
    } else {
        E += getLJ(O5, C4, LJpairs["OrCe"].CS12, LJpairs["OrCe"].CS6);
        E += getLJ(C2, C4, LJpairs["CrCe"].CS12, LJpairs["CrCe"].CS6);
    }
    
    return E;
}

double getSteElecE_just_LJ_extended(vector<atom> res1, vector<atom> res2, long long int link,  map<string,LJparam> LJpairs){
    double E = 0.0;
    //names are given like it would be for 1-4 linking !!!!
    atom O5, C2, O2, C1, O1; //atoms of residue 1
    atom C4, C3, O3, O3b, C5, C6; //atoms of residue 2 (O3b is used for 1-6 linking only)
    O5 = findAtom("O5", res1, 1);
    C2 = findAtom("C2", res1, 1);
    O2 = findAtom("O2", res1, 1); //extend
    C1 = findAtom("C1", res1, 1);
    O1 = findAtom("O1", res1, 1);
    
    //linking atom
    C4 = findAtom("C"+std::to_string(link), res2, 2);
    
    //atoms with link-1
    if(link==1) C3 =findAtom("O5", res2, 2);
    else {
        C3 = findAtom("C"+std::to_string(link-1), res2, 2);
        O3 = findAtom("O"+std::to_string(link-1), res2, 2);
    }
    if(link==6) O3b = findAtom("C4", res2, 2); //special, just for 1-6 linking
    
    //atoms with link+1
    if(link!=6 ) C5 = findAtom("C"+std::to_string(link+1), res2, 2); //  1-6 link doesn't have "C5"   
    
    if(link==4 ) C6 = findAtom("C6", res2, 2);
    else if(link!=6 ) C6 = findAtom("O"+std::to_string(link+1), res2, 2); // 1-6 link doesn't have "C6"
    
    
    //////////////////////////
    //CALCULATE LJ
    //////////////////////////
    if (link==4){
        //O5
        E += getLJ(O5, C4, LJpairs["OrCr"].CS12, LJpairs["OrCr"].CS6);
        E += getLJ(O5, C3, LJpairs["OrCr"].C12, LJpairs["OrCr"].C6);
        E += getLJ(O5, C5, LJpairs["OrCr"].C12, LJpairs["OrCr"].C6);
        E += getLJ(O5, C6, LJpairs["CeOr"].C12, LJpairs["CeOr"].C6);
        E += getLJ(O5, O3, LJpairs["OOr"].C12, LJpairs["OOr"].C6);
        //O2
        E += getLJ(O2, C4, LJpairs["OCr"].C12, LJpairs["OCr"].C6);
        E += getLJ(O2, C3, LJpairs["OCr"].C12, LJpairs["OCr"].C6);
        E += getLJ(O2, C5, LJpairs["OCr"].C12, LJpairs["OCr"].C6);
        E += getLJ(O2, C6, LJpairs["OCe"].C12, LJpairs["OCe"].C6);
        E += getLJ(O2, O3, LJpairs["OO"].C12, LJpairs["OO"].C6);
        //C2
        E += getLJ(C2, C4, LJpairs["CrCr"].CS12, LJpairs["CrCr"].CS6);
        E += getLJ(C2, C3, LJpairs["CrCr"].C12, LJpairs["CrCr"].C6);
        E += getLJ(C2, C5, LJpairs["CrCr"].C12, LJpairs["CrCr"].C6);
        E += getLJ(C2, C6, LJpairs["CeCr"].C12, LJpairs["CeCr"].C6);
        E += getLJ(C2, O3, LJpairs["OCr"].C12, LJpairs["OCr"].C6);
        //C1
        E += getLJ(C1, C3, LJpairs["CrCr"].CS12, LJpairs["CrCr"].CS6);
        E += getLJ(C1, C5, LJpairs["CrCr"].CS12, LJpairs["CrCr"].CS6);
        E += getLJ(C1, C6, LJpairs["CeCr"].C12, LJpairs["CeCr"].C6);
        E += getLJ(C1, O3, LJpairs["OCr"].C12, LJpairs["OCr"].C6);
    } else if (link==2 || link==3){
        //O5
        E += getLJ(O5, C4, LJpairs["OrCr"].CS12, LJpairs["OrCr"].CS6);
        E += getLJ(O5, C3, LJpairs["OrCr"].C12, LJpairs["OrCr"].C6);
        E += getLJ(O5, C5, LJpairs["OrCr"].C12, LJpairs["OrCr"].C6);
        E += getLJ(O5, C6, LJpairs["OOr"].C12, LJpairs["OOr"].C6); // C6 is here alcohol oxygen
        E += getLJ(O5, O3, LJpairs["OOr"].C12, LJpairs["OOr"].C6);
        //O2
        E += getLJ(O2, C4, LJpairs["OCr"].C12, LJpairs["OCr"].C6);
        E += getLJ(O2, C3, LJpairs["OCr"].C12, LJpairs["OCr"].C6);
        E += getLJ(O2, C5, LJpairs["OCr"].C12, LJpairs["OCr"].C6);
        E += getLJ(O2, C6, LJpairs["OO"].C12, LJpairs["OO"].C6); // C6 is here alcohol oxygen
        E += getLJ(O2, O3, LJpairs["OO"].C12, LJpairs["OO"].C6);
        //C2
        E += getLJ(C2, C4, LJpairs["CrCr"].CS12, LJpairs["CrCr"].CS6);
        E += getLJ(C2, C3, LJpairs["CrCr"].C12, LJpairs["CrCr"].C6);
        E += getLJ(C2, C5, LJpairs["CrCr"].C12, LJpairs["CrCr"].C6);
        E += getLJ(C2, C6, LJpairs["OCr"].C12, LJpairs["OCr"].C6); // C6 is here alcohol oxygen
        E += getLJ(C2, O3, LJpairs["OCr"].C12, LJpairs["OCr"].C6);
        //C1
        E += getLJ(C1, C3, LJpairs["CrCr"].CS12, LJpairs["CrCr"].CS6);
        E += getLJ(C1, C5, LJpairs["CrCr"].CS12, LJpairs["CrCr"].CS6);
        E += getLJ(C1, C6, LJpairs["OCr"].C12, LJpairs["OCr"].C6); // C6 is here alcohol oxygen
        E += getLJ(C1, O3, LJpairs["OCr"].C12, LJpairs["OCr"].C6);
    } else if (link==1){
        //O5
        E += getLJ(O5, C4, LJpairs["OrCr"].CS12, LJpairs["OrCr"].CS6); // C3 is here ring oxygen
        E += getLJ(O5, C3, LJpairs["OrOr"].C12, LJpairs["OrOr"].C6);
        E += getLJ(O5, C5, LJpairs["OrCr"].C12, LJpairs["OrCr"].C6);
        E += getLJ(O5, C6, LJpairs["OOr"].C12, LJpairs["OOr"].C6); // C6 is here alcohol oxygen
        // O3 doesn't exist here
        //O2
        E += getLJ(O2, C4, LJpairs["OCr"].C12, LJpairs["OCr"].C6);
        E += getLJ(O2, C3, LJpairs["OOr"].C12, LJpairs["OOr"].C6); // C3 is here ring oxygen
        E += getLJ(O2, C5, LJpairs["OCr"].C12, LJpairs["OCr"].C6);
        E += getLJ(O2, C6, LJpairs["OO"].C12, LJpairs["OO"].C6); // C6 is here alcohol oxygen
        // O3 doesn't exist here
        //C2
        E += getLJ(C2, C4, LJpairs["CrCr"].CS12, LJpairs["CrCr"].CS6);
        E += getLJ(C2, C3, LJpairs["OrCr"].C12, LJpairs["OrCr"].C6); // C3 is here ring oxygen
        E += getLJ(C2, C5, LJpairs["CrCr"].C12, LJpairs["CrCr"].C6);
        E += getLJ(C2, C6, LJpairs["OCr"].C12, LJpairs["OCr"].C6); // C6 is here alcohol oxygen
        // O3 doesn't exist here
        //C1
        E += getLJ(C1, C3, LJpairs["OrCr"].CS12, LJpairs["OrCr"].CS6); // C3 is here ring oxygen
        E += getLJ(C1, C5, LJpairs["CrCr"].CS12, LJpairs["CrCr"].CS6);
        E += getLJ(C1, C6, LJpairs["OCr"].C12, LJpairs["OCr"].C6); // C6 is here alcohol oxygen
        // O3 doesn't exist here
    } else if (link==6){
        //O5------------------------------------------------------------------
        E += getLJ(O5, C4, LJpairs["CeOr"].CS12, LJpairs["CeOr"].CS6); //C4 is here C6 -(CH2)
        E += getLJ(O5, C3, LJpairs["OrCr"].C12, LJpairs["OrCr"].C6);
        // C5 and C6 doesn't exist here
        E += getLJ(O5, O3, LJpairs["OrOr"].C12, LJpairs["OrOr"].C6); // O3 is here ring oxygen
        E += getLJ(O5, O3b, LJpairs["OrCr"].C12, LJpairs["OrCr"].C6); // O3b is C4 for 1-6 link
        //O2------------------------------------------------------------------
        E += getLJ(O2, C4, LJpairs["OCe"].C12, LJpairs["OCe"].C6); //C4 is here C6 -(CH2)
        E += getLJ(O2, C3, LJpairs["OCr"].C12, LJpairs["OCr"].C6);
        // C5 and C6 doesn't exist here
        E += getLJ(O2, O3, LJpairs["OOr"].C12, LJpairs["OOr"].C6); // O3 is here ring oxygen
        E += getLJ(O2, O3b, LJpairs["OCr"].C12, LJpairs["OCr"].C6); // O3b is C4 for 1-6 link
        //C2------------------------------------------------------------------
        E += getLJ(C2, C4, LJpairs["CeCr"].CS12, LJpairs["CeCr"].CS6); //C4 is here C6 -(CH2)
        E += getLJ(C2, C3, LJpairs["CrCr"].C12, LJpairs["CrCr"].C6);
        // C5 and C6 doesn't exist here
        E += getLJ(C2, O3, LJpairs["OrCr"].C12, LJpairs["OrCr"].C6); // O3 is here ring oxygen
        E += getLJ(C2, O3b, LJpairs["CrCr"].C12, LJpairs["CrCr"].C6); // O3b is C4 for 1-6 link
        //C1------------------------------------------------------------------
        E += getLJ(C1, C3, LJpairs["CrCr"].CS12, LJpairs["CrCr"].CS6);
        // C5 and C6 doesn't exist here
        E += getLJ(C1, O3, LJpairs["OrCr"].C12, LJpairs["OrCr"].C6); // O3 is here ring oxygen
        E += getLJ(C1, O3b, LJpairs["CrCr"].C12, LJpairs["CrCr"].C6); // O3b is C4 for 1-6 link
    }
    
    return E;
}

double getSteElecE_just_torsion(vector<atom> res1, vector<atom> res2, long long int link,  map<int,Torsion_param> Torsions){
    double E = 0.0;
    //names are given like it would be for 1-4 linking
    atom O5, C2, C1, O1; //atoms of residue 1
    atom C4, C3, C5, C2b; //atoms of residue 2
    O5 = findAtom("O5", res1, 1);
    C2 = findAtom("C2", res1, 1);
    C1 = findAtom("C1", res1, 1);
    O1 = findAtom("O1", res1, 1);
    
    C4 = findAtom("C"+std::to_string(link), res2, 2);
    
    if(link==1) C3 =findAtom("O5", res2, 2);
    else C3 = findAtom("C"+std::to_string(link-1), res2, 2);
    
    if(link!=6) {
        C5 = findAtom("C"+std::to_string(link+1), res2, 2); 
    } else {
        C5 = findAtom("O5", res2, 2); 
        C2b = findAtom("C4", res2, 2); //used just for 1-6
    }
    
    //////////////////////////
    //CALCULATE TORSIONS
    //////////////////////////
    //PHI
    E += getTorsionEnergy(C2,C1,O1,C4,Torsions[46].K,Torsions[46].s,Torsions[46].m);
    E += getTorsionEnergy(O5,C1,O1,C4,Torsions[42].K,Torsions[42].s,Torsions[42].m);
    E += getTorsionEnergy(O5,C1,O1,C4,Torsions[44].K,Torsions[44].s,Torsions[44].m);
    E += getTorsionEnergy(O5,C1,O1,C4,Torsions[47].K,Torsions[47].s,Torsions[47].m);
    
    //PSI
    E += getTorsionEnergy(C1,O1,C4,C3,Torsions[42].K,Torsions[42].s,Torsions[42].m);
    if (link!=6){
        E += getTorsionEnergy(C1,O1,C4,C3,Torsions[47].K,Torsions[47].s,Torsions[47].m);
        E += getTorsionEnergy(C1,O1,C4,C5,Torsions[46].K,Torsions[46].s,Torsions[46].m);
    } else { //omega
        E += getTorsionEnergy(O1,C4,C3,C2b,Torsions[34].K,Torsions[34].s,Torsions[34].m);
        E += getTorsionEnergy(O1,C4,C3,C2b,Torsions[50].K,Torsions[50].s,Torsions[50].m);
        
        E += getTorsionEnergy(O1,C4,C3,C5,Torsions[51].K,Torsions[51].s,Torsions[51].m);
        E += getTorsionEnergy(O1,C4,C3,C5,Torsions[52].K,Torsions[52].s,Torsions[52].m);
    }
    if (link==1){
        E += getTorsionEnergy(C1,O1,C4,C3,Torsions[44].K,Torsions[44].s,Torsions[44].m);
    }

    return E;
}

vector<pair<string, double > > getClash(vector<atom> res1, vector<atom> res2, vector<string > clashList, map<string,LJparam> LJpairs){
    vector<pair<string, double > > Clashvector;
    long long int link = findLink(res2);

    pair< string, double > NrClash16("NrClash16",0);
    
    //loop over first residue
    for (vector<atom>::iterator it1 = res1.begin(); it1 != res1.end(); ++it1){
        if ((*it1).a.compare("01")==0){ //skip O1 for first residue
            continue;
        }
        if (isIn(clashList,(*it1).a)){ //check if atom in 1st residue is in clash List
            //loop over second residue
            for (vector<atom>::iterator it2 = res2.begin(); it2 != res2.end(); ++it2){
                if (isIn(clashList,(*it2).a)){ //check if atom in 2st residue is in clash List
                    
                    bool third = false; //third neighbor
                    if ((*it1).a.compare("C1")==0 && (*it2).a.compare("C"+std::to_string(link))==0){ //second neighbor
                        continue;
                    }
                    if ((*it1).a.compare("C1")==0 && ( (*it2).a.compare("C"+std::to_string(link-1))==0 || 
                            (*it2).a.compare("C"+std::to_string(link+1))==0) ){ //third neighbor
                        third = true;
                    } 
                    else if ( ( (*it1).a.compare("C2")==0 || (*it1).a.compare("O5")==0 ) && ( (*it2).a.compare("C"+std::to_string(link))==0  ) ){ //third neighbor
                        third = true;
                    }
                    string a,b;
                    //make code for first residue
                    if((*it1).a.compare("O5")==0){
                        a="Or";
                    }
                    else if((*it1).a.compare("C6")==0){
                        a="Ce";
                    } 
                    else if((*it1).a.compare("O1")==0){
                        a="Ol";
                    } 
                    else if ((*it1).a[0]=='C'){
                        a="Cr";
                    } 
                    else {
                        a="O";
                    }
                    //make code for second residue
                    if((*it2).a.compare("O5")==0){
                        b="Or";
                    }
                    else if((*it2).a.compare("C6")==0){
                        b="Ce";
                    }
                    else if ((*it1).a[0]=='C'){
                        b="Cr";
                    } 
                    else {
                        b="O";
                    }
                    string code = a+b;
                    LJparam tmp = LJpairs[code];
                    double C12, C6;
                    if (third){
                        C12 = tmp.CS12;
                        C6  = tmp.CS6;
                    } 
                    else {
                        C12 = tmp.C12;
                        C6  = tmp.C6;
                    }
                    
                    double LJpot = getLJ((*it1),(*it2), C12, C6);
                    
                    if(LJpot > 16*NAkbT ){
                        NrClash16.second += 1.0;
                    }
                }

            }
        }
    }

    Clashvector.push_back(NrClash16);

    return Clashvector;
}

vector<pair<string, double > > getClashHelix(vector<atom> res1, vector<atom> res2, vector<string > clashList, map<string,LJparam> LJpairs){
    vector<pair<string, double > > Clashvector;
    long long int link = findLink(res2);

    pair< string, double > NrClash16("NrClash16",0);
    
    //loop over first residue
    for (vector<atom>::iterator it1 = res1.begin(); it1 != res1.end(); ++it1){
        if ((*it1).a.compare("01")==0){ //skip O1 for first residue
            continue;
        }
        if (isIn(clashList,(*it1).a)){ //check if atom in 1st residue is in clash List
            //loop over second residue
            for (vector<atom>::iterator it2 = res2.begin(); it2 != res2.end(); ++it2){
                if (isIn(clashList,(*it2).a)){ //check if atom in 2st residue is in clash List
                    
                    string a,b;
                    //make code for first residue
                    if((*it1).a.compare("O5")==0){
                        a="Or";
                    }
                    else if((*it1).a.compare("C6")==0){
                        a="Ce";
                    } 
                    else if((*it1).a.compare("O1")==0){
                        a="Ol";
                    } 
                    else if ((*it1).a[0]=='C'){
                        a="Cr";
                    } 
                    else {
                        a="O";
                    }
                    //make code for second residue
                    if((*it2).a.compare("O5")==0){
                        b="Or";
                    }
                    else if((*it2).a.compare("C6")==0){
                        b="Ce";
                    }
                    else if ((*it1).a[0]=='C'){
                        b="Cr";
                    } 
                    else {
                        b="O";
                    }
                    string code = a+b;
                    LJparam tmp = LJpairs[code];
                    double C12, C6;
                    C12 = tmp.C12;
                    C6  = tmp.C6;
                    
                    double LJpot = getLJ((*it1),(*it2), C12, C6);
                    
                    if(LJpot > 16*NAkbT ){
                        NrClash16.second += 1.0;
                    }
                }

            }
        }
    }

    Clashvector.push_back(NrClash16);

    return Clashvector;
}

vector<pair<string, double > > getClashTest(vector<atom> res1, vector<atom> res2, vector<string > clashList, map<string,LJparam> LJpairs){
    vector<pair<string, double > > Clashvector;
    long long int link = findLink(res2);

    pair< string, double > NrClash("NrClash",0);
    pair< string, double > NrClash2("NrClash2",0);
    pair< string, double > NrClash4("NrClash4",0);
    pair< string, double > NrClash8("NrClash8",0);
    pair< string, double > NrClash16("NrClash16",0);
    pair< string, double > NrClash32("NrClash32",0);
    pair< string, double > NrClash64("NrClash64",0);
    pair< string, double > PotClash("PotClash",0);
    
    //loop over first residue
    for (vector<atom>::iterator it1 = res1.begin(); it1 != res1.end(); ++it1){
        if ((*it1).a.compare("01")==0){ //skip O1 for first residue
            continue;
        }
        if (isIn(clashList,(*it1).a)){ //check if atom in 1st residue is in clash List
            //loop over second residue
            for (vector<atom>::iterator it2 = res2.begin(); it2 != res2.end(); ++it2){
                if (isIn(clashList,(*it2).a)){ //check if atom in 2st residue is in clash List
                    //
                    //printAtom(*it1);
                    //printAtom(*it2);
                    //
                    bool third = false; //third neighbor
                    if ((*it1).a.compare("C1")==0 && (*it2).a.compare("C"+std::to_string(link))==0){ //second neighbor
                        continue;
                    }
                    if ((*it1).a.compare("C1")==0 && ( (*it2).a.compare("C"+std::to_string(link-1))==0 || 
                            (*it2).a.compare("C"+std::to_string(link+1))==0) ){ //third neighbor
                        third = true;
                    } else if ( ( (*it1).a.compare("C2")==0 || (*it1).a.compare("O5")==0 ) && ( (*it2).a.compare("C"+std::to_string(link))==0  ) ){ //third neighbor
                        third = true;
                    }
                    string a,b;
                    //make code for first residue
                    if((*it1).a.compare("O5")==0){
                        a="Or";
                    }else if((*it1).a.compare("C6")==0){
                        a="Ce";
                    } else if((*it1).a.compare("O1")==0){
                        a="Ol";
                    } else if ((*it1).a[0]=='C'){
                        a="Cr";
                    } else {
                        a="O";
                    }
                    //make code for second residue
                    if((*it2).a.compare("O5")==0){
                        b="Or";
                    }else if((*it2).a.compare("C6")==0){
                        b="Ce";
                    }else if ((*it1).a[0]=='C'){
                        b="Cr";
                    } else {
                        b="O";
                    }
                    string code = a+b;
                    LJparam tmp = LJpairs[code];
                    double C12, C6;
                    if (third){
                        C12 = tmp.CS12;
                        C6  = tmp.CS6;
                    } else {
                        C12 = tmp.C12;
                        C6  = tmp.C6;
                    }
                    
                    double LJpot = getLJ((*it1),(*it2), C12, C6);
                    PotClash.second += LJpot;
                    
                    if(LJpot > NAkbT ){
                        NrClash.second += 1.0;
                    }
                    if(LJpot > 2*NAkbT ){
                        NrClash2.second += 1.0;
                    }
                    if(LJpot > 4*NAkbT ){
                        NrClash4.second += 1.0;
                    }
                    if(LJpot > 8*NAkbT ){
                        NrClash8.second += 1.0;
                    }
                    if(LJpot > 16*NAkbT ){
                        NrClash16.second += 1.0;
                    }
                    if(LJpot > 32*NAkbT ){
                        NrClash32.second += 1.0;
                    }
                    if(LJpot > 64*NAkbT ){
                        NrClash64.second += 1.0;
                    }
                }

            }
        }
    }
    Clashvector.push_back(PotClash);
    Clashvector.push_back(NrClash);
    Clashvector.push_back(NrClash2);
    Clashvector.push_back(NrClash4);
    Clashvector.push_back(NrClash8);
    Clashvector.push_back(NrClash16);
    Clashvector.push_back(NrClash32);
    Clashvector.push_back(NrClash64);
    return Clashvector;
}

int changeAngle(vector<atom> &angle, vector<vector<string> > tpl, string flag, string instruct){
    //prepare input
    instruct = merge(instruct, ' '); //merge multispaces
    if (instruct[0]==' '){ //if first space is present, remove
        instruct.assign(instruct,1,instruct.length());
    }
    vector<string> vecInstr = split(instruct, ' ');
////////////////////////////////////////////////////////
//                        PHI                         //
////////////////////////////////////////////////////////
    if (flag.compare("allphi")==0){
        // set input variables
        double value= atof(vecInstr[0].c_str());
        
        vector<pair <int,int> > con = findConnect(angle); //get connection
        setPhi(angle, tpl, value, con);
        //set last phi
        if(con.back().second!=5){
            angle[findAtomPos("HO1", angle, angle.back().nrRes)].pos.z = value;
        }
        

    } else if (flag.compare("phi")==0){
        // set input variables
        int res = atoi(vecInstr[0].c_str())+1; //phi is defined in the residue after the linkage
        double value= atof(vecInstr[1].c_str());
  
        if (res == angle.back().nrRes+1) { //if trying to access one too much res make phi for HO1
            angle[findAtomPos("HO1", angle, angle.back().nrRes)].pos.z = value;
        } else {
            vector<pair <int,int> > con = findConnect(angle, res); //get connection
            setPhi(angle, tpl, value, con);
        }

        
////////////////////////////////////////////////////////
//                        PSI                         //
////////////////////////////////////////////////////////
    } else if (flag.compare("allpsi")==0){
        // set input variables
        double value= atof(vecInstr[0].c_str());
        vector<pair <int,int> > con = findConnect(angle); //get connection
        setPsi(angle, tpl, value, con);
        
    } else if (flag.compare("psi")==0){
        // set input variables
        int res = atoi(vecInstr[0].c_str())+1;//phi is defined in the residue after the linkage
        double value= atof(vecInstr[1].c_str());
        
        vector<pair <int,int> > con = findConnect(angle, res); //get connection
        setPsi(angle, tpl, value, con);
        
////////////////////////////////////////////////////////
//                        ADD                         //
////////////////////////////////////////////////////////
    } else if (flag.compare("add")==0){
        //Defining some variables
        vector<atom> res = getRes(angle, atoi(vecInstr[0].c_str())); //get residue to multiply
        string link = vecInstr[1];              // linkage (1-?)
        string  Cn, Cnn, Cnnn;                  // carbons after linking carbon:"C next", "C next next" and for 6link "C next next next"
        atom O;                                 // O link to be removed used to check enantiomer
//        vector<string > ring {"C4","C3","C2","C1","O5","C5"};
//        for(vector<string>::iterator it = ring.begin();it != ring.end();it++){
//            if((*it).compare("C"+link)==0){
//                Cn  = (*++it);
//                Cnn = (*++it);
//                break;  
//            }
//        }
        switch(stoi(link)){
            case 1:
                Cn  = "O5";
                Cnn = "C5"; 
                break;
            case 2:
                Cn  = "C1";
                Cnn = "O5"; 
                break;
            case 3:
                Cn  = "C2";
                Cnn = "C1"; 
                break;
            case 4:
                Cn  = "C3";
                Cnn = "C2"; 
                break;
            case 6:
                Cn  = "C5";
                Cnn = "C4"; 
                Cnnn= "C3";
                break;
        }
        
        //get rid of last atom form res - HO1
        if (res.back().a.compare("HO1")==0){
            res.pop_back();
        }
        //get rid of connectin HO and O; save O
        for(vector<atom>::iterator it = res.begin();it != res.end();it++){
            if((*it).a.compare("HO"+link)==0 ){
                res.erase(it);
                it--;
            } else if ((*it).a.compare("O"+link)==0){
                O = (*it);
                res.erase(it);
                it--;
            } 
        }
        
        if (link.compare("6")==0){
            // make default phi and psi and translate torsion for Cnn depending on enantiomer at the linkage
            for(vector<atom>::iterator it = res.begin();it != res.end();it++){
                double rac = 0.0;
                if ((*it).a.compare("C"+link)==0){
                    rac = (*it).pos.z;
                    (*it).pos.z = 60.0;
                } else if ((*it).a.compare(Cn)==0){
                    (*it).pos.z = 180.0;
                } else if ((*it).a.compare(Cnn)==0){
                    (*it).pos.z = 300.0;
                } else if ((*it).a.compare(Cnnn)==0){
                    if (rac > 180.0){ // check for enantiomer in linkage
                        (*it).pos.z += 120.0;
                    } else {
                        (*it).pos.z -= 120.0;
                    }
                }
            }
        } else {
            // make default phi and psi and translate torsion for Cnn depending on enantiomer at the linkage
            for(vector<atom>::iterator it = res.begin();it != res.end();it++){
                if ((*it).a.compare("C"+link)==0){
                    (*it).pos.z = 60.0;
                } else if ((*it).a.compare(Cn)==0){
                    (*it).pos.z = 180.0;
                } else if ((*it).a.compare(Cnn)==0){
                    if (O.pos.z > 180.0){ // check for enantiomer in linkage
                        (*it).pos.z += 120.0;
                    } else {
                        (*it).pos.z -= 120.0;
                    }
                }
            }
        }
        
        //take last atom from angle - HO1
        atom last = angle.back();
        angle.pop_back();
        
        for(int i = 0; i < atoi(vecInstr[2].c_str()); i++){ //loop over number of added residue
            for(int j = 0; j < res.size(); j++){ //loop over size of added residue
                res[j].nrA = angle.back().nrA+1; //increase atom number by one
                res[j].nrRes = last.nrRes+1+i; //increase residue number
                angle.push_back(res[j]); 
            }
        }
        
        if (link.compare("1")!=0){
            last.nrA = angle.back().nrA+1;  //fix correct nrA for last hydrogen
            last.nrRes = last.nrRes+atoi(vecInstr[2].c_str()); //fix correct nrRes for last hydrogen
            angle.push_back(last);//add last hydrogen
        }
        
        //printCnf(angle);
        
////////////////////////////////////////////////////////
//                    DELETION                        //
////////////////////////////////////////////////////////    
    } else if (flag.compare("del")==0){
        atom last = angle.back(); // get HO1
        angle.pop_back();
        int res = last.nrRes; //holding number of last residue
        
        //check number or residues to be removed
        if (atoi(vecInstr[0].c_str())>=res){
                cerr << "-----------------------------" << endl;
                cerr << "Error:" << endl;
                cerr << "In procedure:" << endl;
                cerr << flag << " " << instruct << endl;
                cerr << "ERROR: Trying to remove all residues or more residues than are pressent in the molecule" << endl;
                cerr << "-----------------------------" << endl;
                return 1;
        }
        for(int i = 0; i < atoi(vecInstr[0].c_str()); i++, res--){ //loop over number of residues to be removed
            while (res==angle.back().nrRes){ // if last atom is in removing residue
                angle.pop_back();     //remove atom
            }
        }
        //printCnf(angle);
        last.nrA = angle.back().nrA+1;  //fix correct nrA for last hydrogen
        last.nrRes = angle.back().nrRes; //fix correct nrRes for last hydrogen
        angle.push_back(last);  //add last hydrogen
        
////////////////////////////////////////////////////////
//                    ENANTIOMER                      //
////////////////////////////////////////////////////////
    }  else if (flag.compare("ena")==0){
        int res = atoi(vecInstr[0].c_str());
        string name = vecInstr[1];

        //check if name is possible
        const char *aname[] = {"O4","O3","O2","O1","C6"};
        vector<string > vname (aname, aname+5);
        for(vector<string>::iterator it = vname.begin();it != vname.end();it++){
            if((*it).compare(name)==0){
                break;  
            } else if (it==--vname.end()){
                cerr << "-----------------------------" << endl;
                cerr << "Errorr:" << endl;
                cerr << "Ivalid procedure:" << endl;
                cerr << "-----------------------------" << endl;
                cerr << flag << " " << instruct << endl;
                cerr << "-----------------------------" << endl;
                cerr << "Flag: \""  << vecInstr[1] << "\" unknown"<< endl;
                return 1;
            }
        }
        
        
        //check if enantiomer site is connection site
        //---
        //get connection
        vector<pair <int,int> > con = findConnect(angle, res); 
        
        bool connect=false; //variable to hold if we are at enantiomer site
        string D,C,B,A,H; //atom names for calculating orientation of site and H atom to make change
        
        //check if "name" in connection
        //if yes make "connect" true and give correct names to D,C,B,A,H atoms
        if (con[0].second == 5 && name.compare("O1")==0) {
            connect = true;
            D = "O1"; //of preceding residue
            C = "C1";
            B = "C2";
            A = "O5";
            H = "C5";
        }
        else if (con[0].second == 4 && name.compare("O2")==0) {
            connect = true;
            D = "O1"; //of preceding residue
            C = "C2";
            B = "C3";
            A = "C1";
            H = "O5";
        }
        else if (con[0].second == 3 && name.compare("O3")==0) {
            connect = true;
            D = "O1"; //of preceding residue
            C = "C3";
            B = "C4";
            A = "C2";
            H = "C1";
        }
        else if (con[0].second == 2 && name.compare("O4")==0){
            connect = true;
            D = "O1"; //of preceding residue
            C = "C4";
            B = "C5";
            A = "C3";
            H = "C2";
        }
        else if (con[0].second == 6 && name.compare("C6")==0){
            connect = true;
            D = "C6";
            C = "C5";
            B = "O5";
            A = "C4";
            H = "C3";
        }
        //if we are at connection site
        //we have to check orientation at the site
        if(connect){
            //create coordinate matrix
            vector<atom> cnftmp;
            if (fillCnf(angle, tpl, cnftmp)==1){
                return 1;
            }
            //create vectors to calculate dih - angle at enantimeric center
            vec d = cnftmp[findAtomPos(D, cnftmp, res-1)].pos;
            vec c = cnftmp[findAtomPos(C, cnftmp, res)].pos;
            vec b = cnftmp[findAtomPos(B, cnftmp, res)].pos;
            vec a = cnftmp[findAtomPos(A, cnftmp, res)].pos;
            //calculate the angle
            double dih = calcDihedral(d, c, b, a);
            //check if client asks for orientation which is same or different than present
            if (dih<180 && vecInstr[2].compare("up")==0){
                //nothing to do
                cout << "-----------------------------" << endl;
                cout << "Warning" << endl;
                cout << "Procedure: " << flag << " " << instruct << endl;
                cout << "System is already in demanded conformation." << endl;
                cout << "Nothing done." << endl;
                cout << "-----------------------------" << endl;
            }else if (dih>180 && vecInstr[2].compare("down")==0){
                //nothing to do
                cout << "-----------------------------" << endl;
                cout << "Warning" << endl;
                cout << "Procedure: " << flag << " " << instruct << endl;
                cout << "System is already in demanded conformation." << endl;
                cout << "Nothing done." << endl;
                cout << "-----------------------------" << endl;
            }else if (dih<180 && vecInstr[2].compare("down")==0){
                //system is up and client wants it down
                //angle[findAtomPos(A, angle, res)].pos.z += 60.0;
                angle[findAtomPos(H, angle, res)].pos.z -= 120.0;
            }else if (dih>180 && vecInstr[2].compare("up")==0){
                //system is down and client wants it up
                //angle[findAtomPos(A, angle, res)].pos.z -= 60.0;
                angle[findAtomPos(H, angle, res)].pos.z += 120.0;
            }
            return 0;
        }
        
        //case we have correct name and it is not at connection site
        //make the change
        double value;
        if(vecInstr[2].compare("up")==0){
            value = 120.0;
        } else if(vecInstr[2].compare("down")==0){
            value = 240.0;
        } else {
            cerr << "-----------------------------" << endl;
            cerr << "Error:" << endl;
            cerr << "Ivalid procedure:" << endl;
            cerr << "-----------------------------" << endl;
            cerr << flag << " " << instruct << endl;
            cerr << "-----------------------------" << endl;
            cerr << "Flag: \""  << vecInstr[2] << "\" unknown"<< endl;
            return 1;
        }
        angle[findAtomPos(name, angle, res)].pos.z = value;
        
////////////////////////////////////////////////////////
//                    ROTATION                        //
////////////////////////////////////////////////////////
    } else if (flag.compare("rot")==0){
        int res = atoi(vecInstr[0].c_str());
        if (res > angle.back().nrRes){
            cerr << "Ivalid procedure:" << endl;
            cerr << "-----------------------------" << endl;
            cerr << flag << " " << instruct << endl;
            cerr << "residue number is higher than total number of residues" << endl;
            cerr << "-----------------------------" << endl;
        }
        string name = vecInstr[1];
        double value = atof(vecInstr[2].c_str());
        double value2 = value; //used just in case "ot" and 1-6 linking see below
        string atomn;
        if(name.compare("ch1")==0){
            atomn = "HO1";
        } else if(name.compare("ch2")==0){
            atomn = "HO2";
        } else if(name.compare("ch3")==0){
            atomn = "HO3";
        } else if(name.compare("ch4")==0){
            atomn = "HO4";
        } else if(name.compare("ch6")==0){
            atomn = "HO6";
        //omega tilda - used in template (O6-C6-C5-O5)
        } else if(name.compare("ot")==0){
            atomn = "O6";
        //omega (IUPAC) - (O6-C6-C5-C4)
        //need to change value to fit omega tilda
        //it depends on L or D serie (orientation of C6)
        } else if(name.compare("o")==0){
            // construct saccharide in cartesian coordinates
            vector<atom> cnftmp;
            if (fillCnf(angle, tpl, cnftmp)==1){
                cerr << "error in fillCnf()" << endl;
                return 1;
            }
            // get atoms to to calculate chirality
            atom tmpC6 = findAtom("C6", cnftmp, res);
            atom tmpC5 = findAtom("C5", cnftmp, res);
            atom tmpO5 = findAtom("O5", cnftmp, res);
            atom tmpC4 = findAtom("C4", cnftmp, res);
            // calculate dihedral angle at chiral center
            double dih = calcDihedral(tmpC6.pos, tmpC5.pos, tmpO5.pos, tmpC4.pos);
            if(dih > 180.0){
                value += 120.0;
            } else{
                value -= 120.0;
            }
            atomn = "O6";
        }
        //check for 1-6 linking
        // if we have 1-6 linking, O6 doesn't exist 
        // then we have to change angle on C4
        // because thats fillCnf() procedure for successive residue
        // then "o" should get original value and "ot" has to be changed based on chirality on C6
        vector<pair <int,int> > con = findConnect(angle, res);
        if (con[0].second==6){ 
            if(name.compare("o")==0){
                atomn = "C4";
                value = value2; 
            } else if (name.compare("ot")==0){
                atomn = "C4";
                // construct saccharide in cartesian coordinates
                vector<atom> cnftmp;
                if (fillCnf(angle, tpl, cnftmp)==1){
                    cerr << "error in fillCnf()" << endl;
                    return 1;
                }
                // get atoms to to calculate chirality
                atom tmpC6 = findAtom("C6", cnftmp, res);
                atom tmpC5 = findAtom("C5", cnftmp, res);
                atom tmpO5 = findAtom("O5", cnftmp, res);
                atom tmpC4 = findAtom("C4", cnftmp, res);
                // calculate dihedral angle at chiral center
                double dih = calcDihedral(tmpC6.pos, tmpC5.pos, tmpO5.pos, tmpC4.pos);
                if(dih > 180.0){
                    value += 120.0;
                } else{
                    value -= 120.0;
                }
            }
        }
        
        //check if we can find the atom
        if (findAtomPos(atomn, angle, res)==-1){
            cerr << "Ivalid procedure:" << endl;
            cerr << "-----------------------------" << endl;
            cerr << flag << " " << instruct << endl;
            //cerr << "-----------------------------" << endl;
            //cerr << "Trying to acces nonexisting atom through " <<
            //"flag: \"" << name << "\"" << endl;
            cerr << "Instead is called:" << endl;
            stringstream sstream;
            sstream << (res-1);
            string sres = sstream.str();
            cerr << "psi " << sres+" "+vecInstr[2] << endl;
            cerr << "-----------------------------" << endl;
            changeAngle(angle, tpl, "psi", sres+" "+vecInstr[2]);
            return 0;
        }
        angle[findAtomPos(atomn, angle, res)].pos.z = value;
////////////////////////////////////////////////////////
//                      RING                          //
////////////////////////////////////////////////////////
    } else if (flag.compare("ring")==0){
        int res = atoi(vecInstr[0].c_str());
        string conf = vecInstr[1];
        double value;
        if(conf.compare("4C1")==0){
            angle[findAtomPos("O5", angle, res)].pos.z =  60;
            angle[findAtomPos("C5", angle, res)].pos.z = 300;
            angle[findAtomPos("C4", angle, res)].pos.z =  60;
            angle[findAtomPos("C3", angle, res)].pos.z = 300;
            angle[findAtomPos("C2", angle, res)].pos.z =  60;
            angle[findAtomPos("C1", angle, res)].pos.z = 300;
        } else if(conf.compare("1C4")==0){
            angle[findAtomPos("O5", angle, res)].pos.z = 300;
            angle[findAtomPos("C5", angle, res)].pos.z =  60;
            angle[findAtomPos("C4", angle, res)].pos.z = 300;
            angle[findAtomPos("C3", angle, res)].pos.z =  60;
            angle[findAtomPos("C2", angle, res)].pos.z = 300;
            angle[findAtomPos("C1", angle, res)].pos.z =  60;
        } else {
            cerr << "-----------------------------" << endl;
            cerr << "Error:" << endl;
            cerr << "Ivalid procedure:" << endl;
            cerr << "-----------------------------" << endl;
            cerr << flag << " " << instruct << endl;
            cerr << "-----------------------------" << endl;
            cerr << "Flag: \""  << conf << "\" unknown"<< endl;
            return 1;
        }
        
    }else{
        cerr << "-----------------------------" << endl;
        cerr << "Error:" << endl;
        cerr << "Ivalid procedure:" << endl;
        cerr << "-----------------------------" << endl;
        cerr << flag << " " << instruct << endl;
        cerr << "-----------------------------" << endl;
        cerr << "flag: \"" << flag << "\" unknown" << endl;
        return 1;
    }
    return 0;
}

void setPhi(vector<atom> &angle, vector<vector<string> > tpl, double value, vector<pair <int,int> > con){
    // set input variables
    
    for(vector<pair <int, int> >::iterator it = con.begin(); it < con.end();it++){
        string name = tpl[it->second][0]; //get name of connected C atom
        angle[findAtomPos(name, angle, it->first)].pos.z = value; // set the value
    }
}

void setPsi(vector<atom> &angle, vector<vector<string> > tpl, double value, vector<pair <int,int> > con){
    // set input variables
    
    for(vector<pair <int, int> >::iterator it = con.begin(); it < con.end();it++){
        string name;
        if (it->second > 4){
            name = tpl[it->second%5][0];
        } else { 
            name = tpl[it->second+1][0]; //get name of connected C atom
        }
        angle[findAtomPos(name, angle, it->first)].pos.z = value; // set the value
        
    }
}

string merge(string str, char delim){
    for (string::iterator it = str.begin(); it < str.end(); it++){
        if (*it!=delim) continue;
        string::iterator it2 = it;
        it2++;
        while (*it==*it2){
            //cout << "," << *it2 << "," << endl;
            str.erase(it2);
            //cout << "," << *it2 << "," << endl;
        }
    }
    return str;
}

vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

int findAtomPos(string name, vector<atom> cnf, int res){
    for(int i=0; i<cnf.size(); i++){
        if(cnf[i].a.compare(name)==0  && cnf[i].nrRes == res){
            atom a = cnf[i];
            return i;
        }
    }
    cerr << "-----------------------------" << endl;
    cerr << "Error: atom " << name << " in residue nr. " << res << " doesn't exist" << endl;
    return -1;
}

int fillCnf( vector<atom> angleIn, vector<vector<string> > tpl, vector<atom> &cnfneu){
    // d - atom to calculate coordinate of
    // c,b,a - atoms to which respect is "d" defined
    int nrs = angleIn.back().nrRes; //number of residues
    vector<vector<int> > skips = findSkpis(angleIn, nrs); //vector for skips in each residue
    atom d,c,b,a;
    // make right size of new coordinate vector
    cnfneu.resize(angleIn.size());
    // loop over angular file
    int i = 0; // for access angleIn and cnfneu
    for(int j = 0; j < (17*nrs); j++, i++){ // j for access tpl
        //printVec(cnfneu);
        
        if (isIn(skips,j)){
            i--;
            continue;
        }
        // copy atom info from angular file to right position (from original cnf) to new cart. coord vector
        cnfneu[angleIn[i].nrA-1].a = angleIn[i].a;
        cnfneu[angleIn[i].nrA-1].nrA = angleIn[i].nrA;
        cnfneu[angleIn[i].nrA-1].nrRes = angleIn[i].nrRes;
        cnfneu[angleIn[i].nrA-1].res = angleIn[i].res;
        if(i==0){ //for first anchoring atom
            //copy positions from angular file
            //cnfneu[angleIn[i].nrA-1].pos = angleIn[i].pos;
            cnfneu[angleIn[i].nrA-1].pos.x = 0.0;
            cnfneu[angleIn[i].nrA-1].pos.y = 0.0;
            cnfneu[angleIn[i].nrA-1].pos.z = 0.0;
            // c used as reference to second anchoring atom
            c = cnfneu[angleIn[i].nrA-1];
        }
        else if(i==1){ //for second anchoring atom
            // anchoring reference positions same as in fillTpl()
            a.pos.x = 0.000000000;
            a.pos.y = -1.000000000;
            a.pos.z = -1.000000000;
            b.pos.x = 0.000000000;
            b.pos.y = -1.000000000;
            b.pos.z = 0.000000000;
            //c is set from first anchoring atom
            
            //const atoms for backward conversion
            const atom A = a;
            const atom B = b;
            const atom C = c;
            
            //transform reference atoms
            transform(a,b,c);
            //calculate coordinates
            cnfneu[angleIn[i].nrA-1].pos.z = calcZ(angleIn[i],c ,b ,a);
            cnfneu[angleIn[i].nrA-1].pos.y = calcY(angleIn[i],c ,b ,a);
            cnfneu[angleIn[i].nrA-1].pos.x = calcX(angleIn[i],c ,b ,a);
            //transform back the calculated atom
            transformBack(A,B,C,cnfneu[angleIn[i].nrA-1]);
            
            //shift reference atom to be used by third anchoring atom
            a = B;
            b = C;
            c = cnfneu[angleIn[i].nrA-1];
        }
        else if(i==2){ //for third anchoring atom
            //const atoms for backward conversion
            const atom A = a;
            const atom B = b;
            const atom C = c;
            
            //transform reference atoms
            transform(a,b,c);
            //calculate coordinates
            cnfneu[angleIn[i].nrA-1].pos.z = calcZ(angleIn[i],c ,b ,a);
            cnfneu[angleIn[i].nrA-1].pos.y = calcY(angleIn[i],c ,b ,a);
            cnfneu[angleIn[i].nrA-1].pos.x = calcX(angleIn[i],c ,b ,a);
            //transform back the calculated atom
            transformBack(A,B,C,cnfneu[angleIn[i].nrA-1]);
        }
        else{ //for regular atom
            //get reference atoms "c,b,a" acording to chain rule specified in "tpl"
            //coordinates for ref. atoms are taken from "cnfneu"

            if (j%17 == 0){ //in case of new residue
                int ar = angleIn[i].nrRes; //actual residue

                int r; //ring position
                if (isIn(skips,(ar-1)*17+7)) { //TO DO this is not fully implemented!!!
                    r = 6; //if 06 in skips build ring from 7rd atom in tpl - C6 - but it is not in the ring!
                }
                else if (isIn(skips,(ar-1)*17+9)) {
                    r = 2; //if 04 in skips build ring from 3rd atom in tpl - C4
                }
                else if (isIn(skips,(ar-1)*17+11)) {
                    r = 3; //if 03 in skips build ring from 4th atom in tpl - C3
                }
                else if (isIn(skips,(ar-1)*17+13)) {
                    r = 4; //if 02 in skips build ring from 5th atom in tpl - C2
                }
                else if (isIn(skips,(ar-1)*17+15)) {
                    r = 5; //if 01 in skips build ring from 6th atom in tpl - C1
                }
                else {
                    cerr << "-------------------------------------------" << endl;
                    cerr << "Error: " << endl;
                    cerr << "conection to residue " << ar << " not found" << endl;
                    cerr << "-------------------------------------------" << endl;
                    return 1;
                }
                
                makeChain(16, tpl, c, b, a, cnfneu, angleIn[i].nrRes-1 ); //make a chain for HO1 of preceding residue
                bool six = false;
                if (r==6){ //option for 6 link
                    //calculate C6
                    six = true;
                    cnfneu[angleIn[i+r].nrA-1].a     = angleIn[i+r].a;
                    cnfneu[angleIn[i+r].nrA-1].nrA   = angleIn[i+r].nrA;
                    cnfneu[angleIn[i+r].nrA-1].nrRes = angleIn[i+r].nrRes;
                    cnfneu[angleIn[i+r].nrA-1].res   = angleIn[i+r].res;
                    const atom A = a;
                    const atom B = b;
                    const atom C = c;
                    //printAtom(angleIn[i+r%6]);
                    //printAtom(c);
                    //printAtom(b);
                    //printAtom(a);
                    transform(a,b,c);

                    cnfneu[angleIn[i+r].nrA-1].pos.z = calcZ(angleIn[i+r],c ,b ,a);
                    cnfneu[angleIn[i+r].nrA-1].pos.y = calcY(angleIn[i+r],c ,b ,a);
                    cnfneu[angleIn[i+r].nrA-1].pos.x = calcX(angleIn[i+r],c ,b ,a);
                    transformBack(A,B,C,cnfneu[angleIn[i+r].nrA-1]);

                    a=B;
                    b=C;
                    c=cnfneu[angleIn[i+r].nrA-1];
                    r=1;
                }
                //calculate ring
                for (int g = 0; g < 6; g++, r++){ // loop over 6 ring atoms anticlockwise
                    // i+r%6 is position of solving atom in angleIn
                    cnfneu[angleIn[i+r%6].nrA-1].a     = angleIn[i+r%6].a;
                    cnfneu[angleIn[i+r%6].nrA-1].nrA   = angleIn[i+r%6].nrA;
                    cnfneu[angleIn[i+r%6].nrA-1].nrRes = angleIn[i+r%6].nrRes;
                    cnfneu[angleIn[i+r%6].nrA-1].res   = angleIn[i+r%6].res;
                    const atom A = a;
                    const atom B = b;
                    const atom C = c;
                    //printAtom(angleIn[i+r%6]);
                    //printAtom(c);
                    //printAtom(b);
                    //printAtom(a);
                    transform(a,b,c);

                    cnfneu[angleIn[i+r%6].nrA-1].pos.z = calcZ(angleIn[i+r%6],c ,b ,a);
                    cnfneu[angleIn[i+r%6].nrA-1].pos.y = calcY(angleIn[i+r%6],c ,b ,a);
                    cnfneu[angleIn[i+r%6].nrA-1].pos.x = calcX(angleIn[i+r%6],c ,b ,a);
                    transformBack(A,B,C,cnfneu[angleIn[i+r%6].nrA-1]);

                    a=B;
                    b=C;
                    c=cnfneu[angleIn[i+r%6].nrA-1];
                }
                //skips no included because we stay in the ring and there are no missing atoms
                j=j+5;
                i=i+5;
                if(six){
                    j=j+1;
                    i=i+1;
                }
            } 
            else {
                makeChain(j, tpl, c, b, a, cnfneu, angleIn[i].nrRes );
                //cout << angleIn[i].nrRes << " " << angleIn[i].a << endl;
                const atom A = a;
                const atom B = b;
                const atom C = c;             
                transform(a,b,c);

                cnfneu[angleIn[i].nrA-1].pos.z = calcZ(angleIn[i],c ,b ,a);
                cnfneu[angleIn[i].nrA-1].pos.y = calcY(angleIn[i],c ,b ,a);
                cnfneu[angleIn[i].nrA-1].pos.x = calcX(angleIn[i],c ,b ,a);

                transformBack(A,B,C,cnfneu[angleIn[i].nrA-1]);
            }
        }
    }
}

void fillCnf( vector<atom> angleIn, vector<vector<string> > tpl, vector<atom> &cnfneu, double phi, double psi){
    // d - atom to calculate coordinate of
    // c,b,a - atoms to which respect id "d" defined
    int nrs = angleIn[angleIn.size()-1].nrRes; //number of residues
    vector<vector<int> > skips = findSkpis(angleIn, nrs); //vector for skips in each residue
    atom d,c,b,a;
    // make right size of new coordinate vector
    cnfneu.resize(angleIn.size());
    // loop over angular file
    int i = 0; // for access angleIn and cnfneu
    for(int j = 0; j < (17*nrs); j++, i++){ // j for access tpl
        if (isIn(skips,j)){
            i--;
            continue;
        }
        // copy atom info from angular file to right position (from original cnf) to new cart. coord vector
        cnfneu[angleIn[i].nrA-1].a = angleIn[i].a;
        cnfneu[angleIn[i].nrA-1].nrA = angleIn[i].nrA;
        cnfneu[angleIn[i].nrA-1].nrRes = angleIn[i].nrRes;
        cnfneu[angleIn[i].nrA-1].res = angleIn[i].res;
        if(i==0){ //for first anchoring atom
            //copy positions from angular file
            cnfneu[angleIn[i].nrA-1].pos = angleIn[i].pos;
            // c used as reference to second anchoring atom
            c = cnfneu[angleIn[i].nrA-1];
        }
        else if(i==1){ //for second anchoring atom
            // anchoring reference positions same as in fillTpl()
            a.pos.x = 0.000000000;
            a.pos.y = -1.000000000;
            a.pos.z = -1.000000000;
            b.pos.x = 0.000000000;
            b.pos.y = -1.000000000;
            b.pos.z = 0.000000000;
            //c is set from first anchoring atom
            
            //const atoms for backward conversion
            const atom A = a;
            const atom B = b;
            const atom C = c;
            
            //transform reference atoms
            transform(a,b,c);
            //calculate coordinates
            cnfneu[angleIn[i].nrA-1].pos.z = calcZ(angleIn[i],c ,b ,a);
            cnfneu[angleIn[i].nrA-1].pos.y = calcY(angleIn[i],c ,b ,a);
            cnfneu[angleIn[i].nrA-1].pos.x = calcX(angleIn[i],c ,b ,a);
            //transform back the calculated atom
            transformBack(A,B,C,cnfneu[angleIn[i].nrA-1]);
            
            //shift reference atom to be used by third anchoring atom
            a = B;
            b = C;
            c = cnfneu[angleIn[i].nrA-1];
        }
        else if(i==2){ //for third anchoring atom
            //const atoms for backward conversion
            const atom A = a;
            const atom B = b;
            const atom C = c;
            
            //transform reference atoms
            transform(a,b,c);
            //calculate coordinates
            cnfneu[angleIn[i].nrA-1].pos.z = calcZ(angleIn[i],c ,b ,a);
            cnfneu[angleIn[i].nrA-1].pos.y = calcY(angleIn[i],c ,b ,a);
            cnfneu[angleIn[i].nrA-1].pos.x = calcX(angleIn[i],c ,b ,a);
            //transform back the calculated atom
            transformBack(A,B,C,cnfneu[angleIn[i].nrA-1]);
        }
        else{ //for regular atom
            //get reference atoms "c,b,a" acording to chain rule specified in "tpl"
            //coordinates for ref. atoms are taken from "cnfneu"

            if (j%17 == 0){ //in case of new residue
                int ar = angleIn[i].nrRes; //actual residue

                int r; //ring position
                if (isIn(skips,(ar-1)*17+7)) { //TO DO this is not fully implemented!!!
                    r = 6; //if 06 in skips build ring from 7rd atom in tpl - C6 - but it is not in the ring!
                }
                else if (isIn(skips,(ar-1)*17+9)) {
                    r = 2; //if 04 in skips build ring from 3rd atom in tpl - C4
                }
                else if (isIn(skips,(ar-1)*17+11)) {
                    r = 3; //if 03 in skips build ring from 4th atom in tpl - C3
                }
                else if (isIn(skips,(ar-1)*17+13)) {
                    r = 4; //if 02 in skips build ring from 5th atom in tpl - C2
                }
                else if (isIn(skips,(ar-1)*17+15)) {
                    r = 5; //if 01 in skips build ring from 6th atom in tpl - C1
                }
                else cout << "conection to residue " << ar << " not found" << endl;
                
                angleIn[i+r%6].pos.z = phi;     //set phi
                angleIn[i+r%6+1].pos.z = psi;   //set psi
                makeChain(16, tpl, c, b, a, cnfneu, angleIn[i].nrRes-1 ); //make a chain for HO1 of preceding residue
                
                for (int g = 0; g < 6; g++, r++){ // loop over 6 ring atoms anticlockwise
                    // i+r%6 is position of solving atom in angleIn
                    cnfneu[angleIn[i+r%6].nrA-1].a     = angleIn[i+r%6].a;
                    cnfneu[angleIn[i+r%6].nrA-1].nrA   = angleIn[i+r%6].nrA;
                    cnfneu[angleIn[i+r%6].nrA-1].nrRes = angleIn[i+r%6].nrRes;
                    cnfneu[angleIn[i+r%6].nrA-1].res   = angleIn[i+r%6].res;
                    const atom A = a;
                    const atom B = b;
                    const atom C = c;
                
                    transform(a,b,c);

                    cnfneu[angleIn[i+r%6].nrA-1].pos.z = calcZ(angleIn[i+r%6],c ,b ,a);
                    cnfneu[angleIn[i+r%6].nrA-1].pos.y = calcY(angleIn[i+r%6],c ,b ,a);
                    cnfneu[angleIn[i+r%6].nrA-1].pos.x = calcX(angleIn[i+r%6],c ,b ,a);
                    transformBack(A,B,C,cnfneu[angleIn[i+r%6].nrA-1]);
                    
                    a=B;
                    b=C;
                    c=cnfneu[angleIn[i+r%6].nrA-1];
                }
                //skips no included because we stay in the ring and there are no missing atoms
                j=j+5;
                i=i+5;
            } 
            else {
                makeChain(j, tpl, c, b, a, cnfneu, angleIn[i].nrRes );
                const atom A = a;
                const atom B = b;
                const atom C = c;             
                transform(a,b,c);

                cnfneu[angleIn[i].nrA-1].pos.z = calcZ(angleIn[i],c ,b ,a);
                cnfneu[angleIn[i].nrA-1].pos.y = calcY(angleIn[i],c ,b ,a);
                cnfneu[angleIn[i].nrA-1].pos.x = calcX(angleIn[i],c ,b ,a);

                transformBack(A,B,C,cnfneu[angleIn[i].nrA-1]);
            }
        }
    }
}
//helping function of creating function
atom findAtom(string name, vector<atom> cnf, int res){
    for(int i=0; i<cnf.size(); i++){
        if(cnf[i].a == name && cnf[i].nrRes == res){
            atom a = cnf[i];
            return a;
        }
    }
    atom a;
    return a;
}

void makeChain(int j, vector<vector<string> > tpl, atom &c, atom &b, atom &a, vector<atom> cnfneu, int res){
    c = findAtom(tpl[j%17][1], cnfneu, res );
    b = findAtom(tpl[j%17][2], cnfneu, res );
    a = findAtom(tpl[j%17][3], cnfneu, res );
 }

//making files functions
void makeOpt(vector<atom > angle){
    fstream ofile(angleOutPath.c_str(), ios::out);
    if ( ofile.fail() ){
        cout << "opt file cannot be opened.\n";
    }
    string output;
//    ofile << "# template for beta-D-Glucopyranose \n";
//    ofile << "# d  - distance to d_atom \n";
//    ofile << "# a  - angle calculated acording \"distance path\" \n";
//    ofile << "# tn - torsion angle calculated acording \"distance path\" \n";
//    ofile << "# te - torsion exocyclic (Cn-1_Cn_On_Hn)-in disacharide one is psi ,and phi (O5_C1_O1_H1) \n";
//    ofile << "# ch - improper diherdal for chiral center \n";
//    ofile << "# o  - out of plane dihedral - ring puckering Pickett and Strauss \n";
//    ofile << "# C  - one or two or tree Cartesian coordinate on place of torsion angle(1), angle(2) and distance(3) \n";
//    ofile << "\n";
    ofile << "# Nres res atom    Natom         distance          angle               torsion\n";
    ofile << "START\n";

    for (vector<atom>::iterator it = angle.begin(); it!=angle.end(); it++){
        ofile.flags(std::ios::right);
        ofile.fill(' ');
        ofile.width(5);
        ofile << (*it).nrRes;
        ofile << ' ';
        ofile.flags(std::ios::left);
        ofile.width(6);
        ofile << (*it).res;
        ofile.width(6);
        ofile << (*it).a;
        ofile.width(6);
        ofile.flags(std::ios::right);
        ofile << (*it).nrA;
        ofile.width(22);
        ofile.precision(15);
        ofile.flags(std::ios::fixed);
        ofile << (*it).pos.x;
        ofile.width(22);
        ofile.precision(15);
        ofile << (*it).pos.y;
        ofile.width(22);
        ofile.precision(15);
        ofile << (*it).pos.z << "\n";
    }
    ofile << "END";
    ofile.close();
}

vector<atom> orderCnf(vector<atom> res, vector<vector<string> > tpl, int nrAlast){
    vector<atom> newres; // variable to be returned
    //order in topology
    vector<string> top = {"C2","O2","HO2","C3","O3","HO3","C4","O4","HO4","C6","O6","HO6","C5","O5","C1","O1","HO1"};
         
    //start to order with first group
    if (res.front().nrRes==1){ //case we have first residue
        newres.push_back(findAtom("HO4", res, res.front().nrRes));
        top.erase(top.begin()+8);
        newres.push_back(findAtom("O4", res, res.front().nrRes));
        top.erase(top.begin()+7);
        newres.push_back(findAtom("C4", res, res.front().nrRes));
        top.erase(top.begin()+6);
    } else{                   //case we don't have first residue
        vector<pair<int, int > > con = findConnect(res, res.front().nrRes); // first is nres (not used!) second is position of C in tpl
        int i = con.front().second;   
        string c = tpl[i][0]; //get string of Carbon atom
        newres.push_back(findAtom(c, res, res.front().nrRes));
        switch(i){  //get carbon position in top vector (link between tpl and top)
            case 2: //C4
                i = 6;
                break;
            case 3: //C3
                i = 3;
                break;
            case 4: //C2
                i = 0;
                break;
            case 5: //C1
                i = 14;
                break;
            case 6: //C6
                i = 9;
                break;
        }
        if (i==14){
            //there is different order in topology for 1->1 link
            top = {"O5","C5","C4","O4","HO4","C3","O3","HO3","C2","O2","HO2","C6","O6","HO6"};
        }
        else 
            top.erase(top.begin()+i,top.begin()+i+3); // remove not existing linking atoms
    }
    
    //loop over the rest of top vector
    for (vector<string>::iterator it = top.begin(); it != top.end(); it++){
        atom a = findAtom(*it, res, res.front().nrRes);
        if (a.nrRes!=0) newres.push_back(a); //check if atom exists
    }
    
    //order atom numbers
    for (vector<atom>::iterator it = newres.begin(); it != newres.end(); it++){
        if (it==newres.begin()){
            (*it).nrA = nrAlast+1;
        } else {
            (*it).nrA = (*(it-1)).nrA+1;
        }
    }
    return newres;
}

void mkCnfOpt(vector<atom> cnf, vector<vector<string> > tpl){
    fstream ofile(cnfOutPath.c_str(), ios::out);
    if ( ofile.fail() ){
        cout << "opt file cannot be opened.\n";
    }
    
    vector<vector<atom > > cnfsplit; //vector to hold split and ordered residues
    int r = 1;                       //residue number 
    vector<atom > tmp;               //tmp vector to hold single residue
    //loop over whole cnf
    for (vector<atom >::iterator it = cnf.begin(); it!=cnf.end(); it++){
        if ((*it).nrRes==r){ //loading tmp vector
            tmp.push_back(*it);
        }
        if ((*it).nrRes>r || (it+1)==cnf.end()){ //end of single residue
            //get last atom number
            int nrAlast;
            if(cnfsplit.empty()){ //case it is first residue
                nrAlast = 0;
            }else{
                nrAlast = cnfsplit.back().back().nrA;
            }
            //order tmp vector and load it in cnfsplit
            tmp = orderCnf(tmp,tpl,nrAlast);
            cnfsplit.push_back(tmp);
            //prepare for the next residue
            tmp.clear();                  //clear tmp for the next residue
            tmp.push_back(*it);           //add first atom of not first residue
            r++;
        }
    }
    
    string output;
    ofile << "TITLE\n";
    ofile << "# CONVERTED  CNF FILE\n";
    ofile << "END\n";
    ofile << "POSITION\n";
    //ofile << "# first 24 chars ignored\n";
    for (vector<vector<atom > >::iterator it = cnfsplit.begin(); it!=cnfsplit.end(); it++){
        vector<atom > res = *it;
        for (int i = 0; i<res.size(); i++){
            ofile.flags(std::ios::right);
            ofile.fill(' ');
            ofile.width(5);
            ofile << res[i].nrRes;
            ofile << ' ';
            ofile.flags(std::ios::left);
            ofile.width(6);
            ofile << res[i].res;
            ofile.width(6);
            ofile << res[i].a;
            ofile.width(6);
            ofile.flags(std::ios::right);
            ofile << res[i].nrA;
            ofile.width(27);
            ofile.precision(20);
            ofile.flags(std::ios::fixed);
            ofile << res[i].pos.x;
            ofile.width(27);
            ofile.precision(20);
            ofile << res[i].pos.y;
            ofile.width(27);
            ofile.precision(20);
            ofile << res[i].pos.z << "\n";
        }
    }
    ofile << "END\n";
    ofile.close();
}

void mkCnfOpt(vector<atom> cnf, vector<vector<string> > tpl, string path){
    //create directory if it doesn't exist
    string dir = path.substr(0, path.rfind("/"));
    struct stat st = {0};
    if (stat(dir.c_str(), &st) == -1) {
        mkdir(dir.c_str(), 0755);
    }
    
    fstream ofile(path.c_str(), ios::out);
    if ( ofile.fail() ){
        cout << "opt file cannot be opened.\n";
    }
    
    vector<vector<atom > > cnfsplit; //vector to hold split and ordered residues
    int r = 1;                       //residue number 
    vector<atom > tmp;               //tmp vector to hold single residue
    //loop over whole cnf
    for (vector<atom >::iterator it = cnf.begin(); it!=cnf.end(); it++){
        if ((*it).nrRes==r){ //loading tmp vector
            tmp.push_back(*it);
        }
        if ((*it).nrRes>r || (it+1)==cnf.end()){ //end of single residue
            //get last atom number
            int nrAlast;
            if(cnfsplit.empty()){ //case it is first residue
                nrAlast = 0;
            }else{
                nrAlast = cnfsplit.back().back().nrA;
            }
            //order tmp vector and load it in cnfsplit
            tmp = orderCnf(tmp,tpl,nrAlast);
            cnfsplit.push_back(tmp);
            //prepare for the next residue
            tmp.clear();                  //clear tmp for the next residue
            tmp.push_back(*it);           //add first atom of not first residue
            r++;
        }
    }
    
    string output;
    ofile << "TITLE\n";
    ofile << "# CONVERTED  CNF FILE\n";
    ofile << "END\n";
    ofile << "POSITION\n";
    //ofile << "# first 24 chars ignored\n";
    for (vector<vector<atom > >::iterator it = cnfsplit.begin(); it!=cnfsplit.end(); it++){
        vector<atom > res = *it;
        for (int i = 0; i<res.size(); i++){
            ofile.flags(std::ios::right);
            ofile.fill(' ');
            ofile.width(5);
            ofile << res[i].nrRes;
            ofile << ' ';
            ofile.flags(std::ios::left);
            ofile.width(6);
            ofile << res[i].res;
            ofile.width(6);
            ofile << res[i].a;
            ofile.width(6);
            ofile.flags(std::ios::right);
            ofile << res[i].nrA;
            ofile.width(27);
            ofile.precision(20);
            ofile.flags(std::ios::fixed);
            ofile << res[i].pos.x;
            ofile.width(27);
            ofile.precision(20);
            ofile << res[i].pos.y;
            ofile.width(27);
            ofile.precision(20);
            ofile << res[i].pos.z << "\n";
        }
    }
    ofile << "END\n";
    ofile.close();
}
 
//making files functions
void makeHelix(vector<vector<helix> > h, string par, string name, string dir){
    //create directory data if it doesn't exist
    struct stat st = {0};
    if (stat(dir.c_str(), &st) == -1) {
        mkdir(dir.c_str(), 0755);
    }

    
    //output file is different for  Hbond (Hbond is a vector)
    if (par.compare("Hbond")==0){
        //create directory data/Hbond if it doesn't exist
        struct stat st = {0};
        if (stat( (dir+"/Hbond").c_str(), &st) == -1) {
            mkdir( (dir+"/Hbond").c_str(), 0755);
        }
        //delete Hbond files
        for (vector<pair<string, int> >::iterator HBit = h[1][1].Hbond.begin(); HBit!= h[1][1].Hbond.end(); HBit++){
            fstream ofile(dir+"/Hbond/"+name+par+"_"+(*HBit).first+".dat", std::fstream::out);
            if ( ofile.fail() ){
                cout << par+".dat file cannot be opened.\n";
            }
            ofile.close();
        }
        
        for (vector<vector<helix> >::iterator it = h.begin(); it!= h.end(); it++){
            for (vector<helix> ::iterator it2 = it->begin(); it2!= it->end(); it2++){
                for (vector<pair<string, int> >::iterator HBit=(*it2).Hbond.begin(); HBit!= (*it2).Hbond.end(); HBit++){
                    fstream ofile(dir+"/Hbond/"+name+par+"_"+(*HBit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                    if ( ofile.fail() ){
                        cout << par+".dat file cannot be opened.\n";
                    }
                    ofile.flags(std::ios::left);
                    ofile.fill(' ');
                    ofile.width(12);
                    ofile.precision(5);

                    ofile << (*HBit).second ;
                    ofile.close();
                }
            }
            //end of row, add "\n" to each file
            
            for (vector<pair<string, int> >::iterator HBit = h[1][1].Hbond.begin(); HBit!= h[1][1].Hbond.end(); HBit++){
                fstream ofile(dir+"/Hbond/"+name+par+"_"+(*HBit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                if ( ofile.fail() ){
                    cout << par+".dat file cannot be opened.\n";
                }

                ofile << "\n";
                ofile.close();
            }
            
        }
    } 
    else if (par.compare("HbondHelix")==0){
        //create directory data/HbondHelix if it doesn't exist
        struct stat st = {0};
        if (stat( (dir+"/HbondHelix").c_str(), &st) == -1) {
            mkdir( (dir+"/HbondHelix").c_str(), 0755);
        }
        //delete Hbond files
        for (vector<pair<string, int> >::iterator HBit = h[1][1].HbondHelix.begin(); HBit!= h[1][1].HbondHelix.end(); HBit++){
            fstream ofile(dir+"/HbondHelix/"+name+par+"_"+(*HBit).first+".dat", std::fstream::out);
            if ( ofile.fail() ){
                cout << par+".dat file cannot be opened.\n";
            }
            ofile.close();
        }
        
        for (vector<vector<helix> >::iterator it = h.begin(); it!= h.end(); it++){
            for (vector<helix> ::iterator it2 = it->begin(); it2!= it->end(); it2++){
                for (vector<pair<string, int> >::iterator HBit=(*it2).HbondHelix.begin(); HBit!= (*it2).HbondHelix.end(); HBit++){
                    fstream ofile(dir+"/HbondHelix/"+name+par+"_"+(*HBit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                    if ( ofile.fail() ){
                        cout << par+".dat file cannot be opened.\n";
                    }
                    ofile.flags(std::ios::left);
                    ofile.fill(' ');
                    ofile.width(12);
                    ofile.precision(5);

                    ofile << (*HBit).second ;
                    ofile.close();
                }
            }
            //end of row, add "\n" to each file
            
            for (vector<pair<string, int> >::iterator HBit = h[1][1].HbondHelix.begin(); HBit!= h[1][1].HbondHelix.end(); HBit++){
                fstream ofile(dir+"/HbondHelix/"+name+par+"_"+(*HBit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                if ( ofile.fail() ){
                    cout << par+".dat file cannot be opened.\n";
                }

                ofile << "\n";
                ofile.close();
            }
            
        }
    //make standard output for clash
    } 
    else if (par.compare("clash")==0){
        //delete clash files
        for (vector<pair<string, double> >::iterator Cit = h[1][1].clash.begin(); Cit!= h[1][1].clash.end(); Cit++){
            fstream ofile(dir+"/"+name+par+"_"+(*Cit).first+".dat", std::fstream::out);
            if ( ofile.fail() ){
                cout << par+".dat file cannot be opened.\n";
            }
            ofile.close();
        }
        
        for (vector<vector<helix> >::iterator it = h.begin(); it!= h.end(); it++){
            for (vector<helix> ::iterator it2 = it->begin(); it2!= it->end(); it2++){
                for (vector<pair<string, double> >::iterator Cit=(*it2).clash.begin(); Cit!= (*it2).clash.end(); Cit++){
                    fstream ofile(dir+"/"+name+par+"_"+(*Cit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                    if ( ofile.fail() ){
                        cout << par+".dat file cannot be opened.\n";
                    }
                    ofile.flags(std::ios::left);
                    ofile.fill(' ');
                    ofile.width(12);
                    ofile.precision(5);

                    ofile << (*Cit).second ;
                    ofile.close();
                }
            }
            //end of row, add "\n" to each file
            
            for (vector<pair<string, double> >::iterator Cit = h[1][1].clash.begin(); Cit!= h[1][1].clash.end(); Cit++){
                fstream ofile(dir+"/"+name+par+"_"+(*Cit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                if ( ofile.fail() ){
                    cout << par+".dat file cannot be opened.\n";
                }

                ofile << "\n";
                ofile.close();
            }
            
        }
    //make standard output for clashHelix
    } 
    else if (par.compare("clashHelix")==0){
        //delete clash files
        for (vector<pair<string, double> >::iterator Cit = h[1][1].clashHelix.begin(); Cit!= h[1][1].clashHelix.end(); Cit++){
            fstream ofile(dir+"/"+name+par+"_"+(*Cit).first+".dat", std::fstream::out);
            if ( ofile.fail() ){
                cout << par+".dat file cannot be opened.\n";
            }
            ofile.close();
        }
        
        for (vector<vector<helix> >::iterator it = h.begin(); it!= h.end(); it++){
            for (vector<helix> ::iterator it2 = it->begin(); it2!= it->end(); it2++){
                for (vector<pair<string, double> >::iterator Cit=(*it2).clashHelix.begin(); Cit!= (*it2).clashHelix.end(); Cit++){
                    fstream ofile(dir+"/"+name+par+"_"+(*Cit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                    if ( ofile.fail() ){
                        cout << par+".dat file cannot be opened.\n";
                    }
                    ofile.flags(std::ios::left);
                    ofile.fill(' ');
                    ofile.width(12);
                    ofile.precision(5);

                    ofile << (*Cit).second ;
                    ofile.close();
                }
            }
            //end of row, add "\n" to each file
            
            for (vector<pair<string, double> >::iterator Cit = h[1][1].clashHelix.begin(); Cit!= h[1][1].clashHelix.end(); Cit++){
                fstream ofile(dir+"/"+name+par+"_"+(*Cit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                if ( ofile.fail() ){
                    cout << par+".dat file cannot be opened.\n";
                }

                ofile << "\n";
                ofile.close();
            }
            
        }
    //make standard output for all other parameters
    } 
    else {
        fstream ofile(dir+"/"+name+par+".dat", ios::out);
        if ( ofile.fail() ){
            cout << par+".dat file cannot be opened.\n";
        }
        for (vector<vector<helix> >::iterator it = h.begin(); it!= h.end(); it++){
            for (vector<helix> ::iterator it2 = it->begin(); it2!= it->end(); it2++){
                ofile.flags(std::ios::left);
                ofile.fill(' ');
                ofile.width(12);
                ofile.precision(5);
                if     (par.compare("R")==0)   ofile << (*it2).R ;
                else if(par.compare("p")==0)   ofile << (*it2).p ;
                else if(par.compare("a")==0)   ofile << (*it2).a ;
                else if(par.compare("g")==0)   ofile << (*it2).g ;
                else if(par.compare("n")==0)   ofile << (*it2).n ;
                else if(par.compare("SteElecE")==0) ofile << (*it2).SteElecE ;
                else if(par.compare("SteElecE_extended")==0) ofile << (*it2).SteElecE_extended ;
                else if(par.compare("SteElecE_just_TOR")==0) ofile << (*it2).SteElecE_just_TOR ;
                else if(par.compare("SteElecE_just_LJ")==0) ofile << (*it2).SteElecE_just_LJ ;
                else if(par.compare("SteElecE_just_LJ_extended")==0) ofile << (*it2).SteElecE_just_LJ_extended ;
                else if(par.compare("phi")==0) ofile << (*it2).phi ;
                else if(par.compare("psi")==0) ofile << (*it2).psi ;
            }
            ofile << "\n";
        }
        ofile.close();
    }
    
}

//making files functions
void makeHelix(vector<vector<vector<helix> > > h, string par, string name, string dir){
    //create directory data if it doesn't exist
    struct stat st = {0};
    if (stat(dir.c_str(), &st) == -1) {
        mkdir(dir.c_str(), 0755);
    }

    
    //output file is different for  Hbond (Hbond is a vector)
    if (par.compare("Hbond")==0){
        //create directory data/Hbond if it doesn't exist
        struct stat st = {0};
        if (stat( (dir+"/Hbond").c_str(), &st) == -1) {
            mkdir( (dir+"/Hbond").c_str(), 0755);
        }
        //delete Hbond files
        for (vector<pair<string, int> >::iterator HBit = h[1][1][1].Hbond.begin(); HBit!= h[1][1][1].Hbond.end(); HBit++){
            fstream ofile(dir+"/Hbond/"+name+par+"_"+(*HBit).first+".dat", std::fstream::out);
            if ( ofile.fail() ){
                cout << par+".dat file cannot be opened.\n";
            }
            ofile.close();
        }
        
        for (vector<vector<vector<helix> > >::iterator it0 = h.begin(); it0!= h.end(); it0++){
            for (vector<vector<helix> >::iterator it = it0->begin(); it!= it0->end(); it++){
                for (vector<helix> ::iterator it2 = it->begin(); it2!= it->end(); it2++){
                    for (vector<pair<string, int> >::iterator HBit=(*it2).Hbond.begin(); HBit!= (*it2).Hbond.end(); HBit++){
                        fstream ofile(dir+"/Hbond/"+name+par+"_"+(*HBit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                        if ( ofile.fail() ){
                            cout << par+".dat file cannot be opened.\n";
                        }
                        ofile.flags(std::ios::left);
                        ofile.fill(' ');
                        ofile.width(12);
                        ofile.precision(5);

                        ofile << (*HBit).second ;
                        ofile.close();
                    }
                }
                //end of row, add "\n" to each file
                for (vector<pair<string, int> >::iterator HBit = h[1][1][1].Hbond.begin(); HBit!= h[1][1][1].Hbond.end(); HBit++){
                    fstream ofile(dir+"/Hbond/"+name+par+"_"+(*HBit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                    if ( ofile.fail() ){
                        cout << par+".dat file cannot be opened.\n";
                    }

                    ofile << "\n";
                    ofile.close();
                }
            }
            //end of slice, add "# next slice\n" to each file
            for (vector<pair<string, int> >::iterator HBit = h[1][1][1].Hbond.begin(); HBit!= h[1][1][1].Hbond.end(); HBit++){
                fstream ofile(dir+"/Hbond/"+name+par+"_"+(*HBit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                if ( ofile.fail() ){
                    cout << par+".dat file cannot be opened.\n";
                }

                ofile << "# next slice\n";
                ofile.close();
            }
        }
    //make standard output for clash
    } else if (par.compare("HbondHelix")==0){
        //create directory data/Hbond if it doesn't exist
        struct stat st = {0};
        if (stat( (dir+"/HbondHelix").c_str(), &st) == -1) {
            mkdir( (dir+"/HbondHelix").c_str(), 0755);
        }
        //delete Hbond files
        for (vector<pair<string, int> >::iterator HBit = h[1][1][1].HbondHelix.begin(); HBit!= h[1][1][1].HbondHelix.end(); HBit++){
            fstream ofile(dir+"/HbondHelix/"+name+par+"_"+(*HBit).first+".dat", std::fstream::out);
            if ( ofile.fail() ){
                cout << par+".dat file cannot be opened.\n";
            }
            ofile.close();
        }
        
        for (vector<vector<vector<helix> > >::iterator it0 = h.begin(); it0!= h.end(); it0++){
            for (vector<vector<helix> >::iterator it = it0->begin(); it!= it0->end(); it++){
                for (vector<helix> ::iterator it2 = it->begin(); it2!= it->end(); it2++){
                    for (vector<pair<string, int> >::iterator HBit=(*it2).HbondHelix.begin(); HBit!= (*it2).HbondHelix.end(); HBit++){
                        fstream ofile(dir+"/HbondHelix/"+name+par+"_"+(*HBit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                        if ( ofile.fail() ){
                            cout << par+".dat file cannot be opened.\n";
                        }
                        ofile.flags(std::ios::left);
                        ofile.fill(' ');
                        ofile.width(12);
                        ofile.precision(5);

                        ofile << (*HBit).second ;
                        ofile.close();
                    }
                }
                //end of row, add "\n" to each file
                for (vector<pair<string, int> >::iterator HBit = h[1][1][1].HbondHelix.begin(); HBit!= h[1][1][1].HbondHelix.end(); HBit++){
                    fstream ofile(dir+"/HbondHelix/"+name+par+"_"+(*HBit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                    if ( ofile.fail() ){
                        cout << par+".dat file cannot be opened.\n";
                    }

                    ofile << "\n";
                    ofile.close();
                }
            }
            //end of slice, add "# next slice\n" to each file
            for (vector<pair<string, int> >::iterator HBit = h[1][1][1].HbondHelix.begin(); HBit!= h[1][1][1].HbondHelix.end(); HBit++){
                fstream ofile(dir+"/HbondHelix/"+name+par+"_"+(*HBit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                if ( ofile.fail() ){
                    cout << par+".dat file cannot be opened.\n";
                }

                ofile << "# next slice\n";
                ofile.close();
            }
        }
    //make standard output for clash
    } else if (par.compare("clash")==0){
        //delete clash files
        for (vector<pair<string, double> >::iterator Cit = h[1][1][1].clash.begin(); Cit!= h[1][1][1].clash.end(); Cit++){
            fstream ofile(dir+"/"+name+par+"_"+(*Cit).first+".dat", std::fstream::out);
            if ( ofile.fail() ){
                cout << par+".dat file cannot be opened.\n";
            }
            ofile.close();
        }
        
        for (vector<vector<vector<helix> > >::iterator it0 = h.begin(); it0!= h.end(); it0++){
            for (vector<vector<helix> >::iterator it = it0->begin(); it!= it0->end(); it++){
                for (vector<helix> ::iterator it2 = it->begin(); it2!= it->end(); it2++){
                    for (vector<pair<string, double> >::iterator Cit=(*it2).clash.begin(); Cit!= (*it2).clash.end(); Cit++){
                        fstream ofile(dir+"/"+name+par+"_"+(*Cit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                        if ( ofile.fail() ){
                            cout << par+".dat file cannot be opened.\n";
                        }
                        ofile.flags(std::ios::left);
                        ofile.fill(' ');
                        ofile.width(12);
                        ofile.precision(5);

                        ofile << (*Cit).second ;
                        ofile.close();
                    }
                }
                //end of row, add "\n" to each file
                for (vector<pair<string, double> >::iterator Cit = h[1][1][1].clash.begin(); Cit!= h[1][1][1].clash.end(); Cit++){
                    fstream ofile(dir+"/"+name+par+"_"+(*Cit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                    if ( ofile.fail() ){
                        cout << par+".dat file cannot be opened.\n";
                    }

                    ofile << "\n";
                    ofile.close();
                }
            }
            //end of row, add "# next slice\n" to each file
            for (vector<pair<string, double> >::iterator Cit = h[1][1][1].clash.begin(); Cit!= h[1][1][1].clash.end(); Cit++){
                fstream ofile(dir+"/"+name+par+"_"+(*Cit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                if ( ofile.fail() ){
                    cout << par+".dat file cannot be opened.\n";
                }

                ofile << "# next slice\n";
                ofile.close();
            }
        }
        //make clasHelix output
    } else if (par.compare("clashHelix")==0){
        //delete clash files
        for (vector<pair<string, double> >::iterator Cit = h[1][1][1].clashHelix.begin(); Cit!= h[1][1][1].clashHelix.end(); Cit++){
            fstream ofile(dir+"/"+name+par+"_"+(*Cit).first+".dat", std::fstream::out);
            if ( ofile.fail() ){
                cout << par+".dat file cannot be opened.\n";
            }
            ofile.close();
        }
        
        for (vector<vector<vector<helix> > >::iterator it0 = h.begin(); it0!= h.end(); it0++){
            for (vector<vector<helix> >::iterator it = it0->begin(); it!= it0->end(); it++){
                for (vector<helix> ::iterator it2 = it->begin(); it2!= it->end(); it2++){
                    for (vector<pair<string, double> >::iterator Cit=(*it2).clashHelix.begin(); Cit!= (*it2).clashHelix.end(); Cit++){
                        fstream ofile(dir+"/"+name+par+"_"+(*Cit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                        if ( ofile.fail() ){
                            cout << par+".dat file cannot be opened.\n";
                        }
                        ofile.flags(std::ios::left);
                        ofile.fill(' ');
                        ofile.width(12);
                        ofile.precision(5);

                        ofile << (*Cit).second ;
                        ofile.close();
                    }
                }
                //end of row, add "\n" to each file
                for (vector<pair<string, double> >::iterator Cit = h[1][1][1].clashHelix.begin(); Cit!= h[1][1][1].clashHelix.end(); Cit++){
                    fstream ofile(dir+"/"+name+par+"_"+(*Cit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                    if ( ofile.fail() ){
                        cout << par+".dat file cannot be opened.\n";
                    }

                    ofile << "\n";
                    ofile.close();
                }
            }
            //end of row, add "# next slice\n" to each file
            for (vector<pair<string, double> >::iterator Cit = h[1][1][1].clashHelix.begin(); Cit!= h[1][1][1].clashHelix.end(); Cit++){
                fstream ofile(dir+"/"+name+par+"_"+(*Cit).first+".dat", std::fstream::in | std::fstream::out | std::fstream::app);
                if ( ofile.fail() ){
                    cout << par+".dat file cannot be opened.\n";
                }

                ofile << "# next slice\n";
                ofile.close();
            }
        }
    //make standard output for all other parameters
    } else {
        fstream ofile(dir+"/"+name+par+".dat", ios::out);
        if ( ofile.fail() ){
            cout << par+".dat file cannot be opened.\n";
        }
        for (vector<vector<vector<helix> > >::iterator it0 = h.begin(); it0!= h.end(); it0++){
            for (vector<vector<helix> >::iterator it = it0->begin(); it!= it0->end(); it++){
                for (vector<helix> ::iterator it2 = it->begin(); it2!= it->end(); it2++){
                    ofile.flags(std::ios::left);
                    ofile.fill(' ');
                    ofile.width(12);
                    ofile.precision(5);
                    if     (par.compare("R")==0)   ofile << (*it2).R ;
                    else if(par.compare("p")==0)   ofile << (*it2).p ;
                    else if(par.compare("a")==0)   ofile << (*it2).a ;
                    else if(par.compare("g")==0)   ofile << (*it2).g ;
                    else if(par.compare("n")==0)   ofile << (*it2).n ;
                    else if(par.compare("SteElecE")==0) ofile << (*it2).SteElecE ;
                    else if(par.compare("SteElecE_extended")==0) ofile << (*it2).SteElecE_extended ;
                    else if(par.compare("SteElecE_just_TOR")==0) ofile << (*it2).SteElecE_just_TOR ;
                    else if(par.compare("SteElecE_just_LJ")==0) ofile << (*it2).SteElecE_just_LJ ;
                    else if(par.compare("SteElecE_just_LJ_extended")==0) ofile << (*it2).SteElecE_just_LJ_extended ;
                    else if(par.compare("phi")==0) ofile << (*it2).phi ;
                    else if(par.compare("psi")==0) ofile << (*it2).psi ;
                    else if(par.compare("omega")==0) ofile << (*it2).omega ;
                }
                ofile << "\n";
            }
            ofile << "# next slice\n";
        }
        ofile.close();
    }
}

void getHelixPar(vector<atom> cnfneu, helix & tmp){
    ///////////////////////////////////////
    // calculate helix parameters
    ///////////////////////////////////////
    atom Pa1 = findAtom("O1", cnfneu, 1);
    atom Pa2 = findAtom("O1", cnfneu, 2);
    atom Pa3 = findAtom("O1", cnfneu, 3);
    vec Aa = getSubstrVec(Pa2.pos, Pa1.pos);
    vec Ba = getSubstrVec(Pa2.pos, Pa3.pos);
    vec Va = getSumVec(Aa, Ba);
    Va = getMultipVec(Va,  1/getAbsVec(Va));

    atom Pb1 = findAtom("O1", cnfneu, 2);
    atom Pb2 = findAtom("O1", cnfneu, 3);
    atom Pb3 = findAtom("O1", cnfneu, 4);
    vec Ab = getSubstrVec(Pb2.pos, Pb1.pos);
    vec Bb = getSubstrVec(Pb2.pos, Pb3.pos);
    vec Vb = getSumVec(Ab, Bb);
    Vb = getMultipVec(Vb,  1/getAbsVec(Vb));

    //axis
    vec H = getCrossofVec(Va, Vb);
    H = getMultipVec(H, 1/getAbsVec(H));
    //printVec(H);

    vec SUB = getSubstrVec(Pa2.pos,Pb2.pos);

    // rise per residue
    double a = getDotofVec(SUB,H);
    tmp.a = a;
    if (a<0.0) { // residues are making over 180deg therefore H is miss oriented
        vec H2 = getCrossofVec(Vb, Va);
        H2 = getMultipVec(H2, 1/getAbsVec(H2));
        a = getDotofVec(SUB,H2);
    }
    tmp.a = a;
    //printf("rise per res %f \n", a);

    // radius
    vec SUB2 = getSubstrVec(Pb2.pos,Pa2.pos);
    //the minus is added empiricaly, It is not according to the paper
    double r = -( (a*a*getAbs2Vec(H)) - getAbs2Vec(SUB) ) / (2*abs(getDotofVec(SUB2,Vb))); 
    tmp.R = r;
    //cout << "r " << r <<  endl;

    // phase angle angular rotation per residue APR
    // can be used just for successive residues 0-360deg;
    // http://math.stackexchange.com/questions/878785/how-to-find-an-angle-in-range0-360-between-2-vectors
    double dot = getDotofVec(Va,Vb);
    //H is changing orientation with Va and Vb (left/right handed) -> det is always positive
    double det = getDotofVec(H,getCrossofVec(Va, Vb));
    double gamma = atan2(det,dot)*180/PI;

    //check left and right handing
    //g- without consequences
//    if(getDotofVec(SUB,H)<0) {
//        h.g = -gamma;
//    } else h.g = gamma;

    //g- consequences
    if(getDotofVec(SUB,H)<0) { 
        gamma = -gamma;
    } 
    tmp.g = gamma;
    //printf("angle per res %f\n", gamma);

    //number of residues per turn
    double RPT = 360/gamma;
    tmp.n = RPT;
    //tmp.n = floor(abs(RPT)+0.5);
    //printf("residues per turn %f\n", RPT);

    // pitch
    double p = RPT*a;
    tmp.p = p;
    //printf("pitch %f\n ----- \n", p);
}

int getHelixMap(vector<atom> angleIn, vector<vector<string> > tpl, string link, float delta, vector<vector<helix> > &Hmatrix){
    int maxLength = 20;
    for (long double psi = (Rmax-delta/2); psi>Rmin;psi -= delta){
        //make phi string and change angle file
        string spsi;
        if (psi<=0.0) spsi = std::to_string(360+psi);//convert psi to positive angle
        else spsi = std::to_string(psi);            

        int state = changeAngle(angleIn, tpl, "allpsi", spsi);

        //if state is >0 some error occured in changeAngle()
        if(state>0){
            return 1;
        }

        vector<helix> hrow; //single row
        for (long double phi = (Rmin+delta/2); phi<Rmax;phi += delta){
            //make psi string and change angle file
            string sphi;
            if (phi<=0.0) sphi = std::to_string(360+phi);//convert phi to positive angle
            else sphi = std::to_string(phi);
            
            int state = changeAngle(angleIn, tpl, "allphi", sphi);

            //if state is >0 some error occured in changeAngle()
            if(state>0){
                return 1;
            }

            vector<atom> cnfneu;
            if (fillCnf(angleIn, tpl, cnfneu)==1){
                cerr << "error in fillCnf()" << endl;
                return 1;
            }
            //printf("phi= %f ; psi = %f \n", static_cast<double>(phi), static_cast<double>(psi) );

            ///////////////////////////////////////
            //Helix calculation calculation
            ///////////////////////////////////////
            helix tmp;
            
            vector<atom> res1 = getRes(cnfneu,1);
            vector<atom> res2 = getRes(cnfneu,2);
            
            ///////////////////////////////////////
            // calculate helix parameters
            ///////////////////////////////////////
            getHelixPar(cnfneu, tmp);
            
            ///////////////////////////////////////
            //clash 
            ///////////////////////////////////////
            vector<string > clashList {"C1","C2","C3","C4","C5","C6","O1","O2","O3","O4","O5"};

            tmp.clash = getClash(res1, res2, clashList, LJpairs);
            
            ///////////////////////////////////////
            //clashHelix 
            ///////////////////////////////////////
            int Nturn = (int)floor(abs(tmp.n)+0.5); //number of residues in one turn
            
            vector<atom> cnftmp; //temporal cnf for clashHelix
            if (Nturn < maxLength && Nturn > 4){
                long long int nr = Nturn-4;
                
                vector<atom> angletmp = angleIn;
                if(changeAngle(angletmp, tpl, "add", "1 " + link + " " + to_string(nr)) == 1) return 1;
                if(changeAngle(angletmp, tpl, "allpsi", spsi)>0) return 1;
                if(changeAngle(angletmp, tpl, "allphi", sphi)>0)return 1;
                
                if (fillCnf(angletmp, tpl, cnftmp) == 1) return 1;
            } 
            else {
                if (fillCnf(angleIn, tpl, cnftmp) == 1) return 1;
            }
            
            tmp.clashHelix = tmp.clash;
            for (int i = 3; i <= (Nturn+1); i++){
                vector<atom> resN =  getRes(cnftmp,i);
                tmp.clashHelix.back().second = tmp.clashHelix.back().second + 
                        getClashHelix(res1, resN, clashList, LJpairs).back().second;
            }
            

            ///////////////////////////////////////
            //Stereoelectronic energy calculation
            ///////////////////////////////////////
            res1 = getRes(cnfneu,1);
            res2 = getRes(cnfneu,2);

            long long int ilink = (long long) stoi(link);

            tmp.SteElecE_just_TOR = getSteElecE_just_torsion(res1, res2, ilink, Torsions);
            
            tmp.SteElecE_just_LJ = getSteElecE_just_LJ(res1, res2, ilink, LJpairs);
            tmp.SteElecE = tmp.SteElecE_just_LJ + tmp.SteElecE_just_TOR;

            tmp.SteElecE_just_LJ_extended = getSteElecE_just_LJ_extended(res1, res2, ilink, LJpairs);
            tmp.SteElecE_extended = tmp.SteElecE_just_LJ_extended + tmp.SteElecE_just_TOR;


            ///////////////////////////////////////
            //H bond calculation
            ///////////////////////////////////////
            tmp.Hbond = getHbond(angleIn, tpl, Hdistmin, Hdistmax, Hangle, Hlist, acceptList);
            
            //Different behavior depending on Nturn
            if (Nturn<4){ // calculate Hbonds (function knows that if Ntrun==2 it's not helical Hbond)
                tmp.HbondHelix = getHbond(angleIn, 1, Nturn, tpl, Hdistmin, Hdistmax, Hangle, Hlist, acceptList);
            } 
            else if (Nturn <= maxLength && Nturn > 3){ // if the second residue is more than Nturn, make angle longer
                vector<atom> angleHelix = angleIn;
                long long int nr = Nturn-3;
                if(changeAngle(angleHelix, tpl, "add", "1 " + link + " " + to_string(nr)) == 1) return 1;
                if(changeAngle(angleHelix, tpl, "allpsi", spsi)>0) return 1;
                if(changeAngle(angleHelix, tpl, "allphi", sphi)>0)return 1;
                tmp.HbondHelix = getHbond(angleHelix, 1, Nturn, tpl, Hdistmin, Hdistmax, Hangle, Hlist, acceptList);
            }
            else if (Nturn > maxLength){ //if Nturn is big make all Hbonds 0
                vector<atom> angleHelix = angleIn;
                long long int nr = maxLength-3;
                if(changeAngle(angleHelix, tpl, "add", "1 " + link + " " + to_string(nr)) == 1) return 1;
                if(changeAngle(angleHelix, tpl, "allpsi", spsi)>0) return 1;
                if(changeAngle(angleHelix, tpl, "allphi", sphi)>0)return 1;
                tmp.HbondHelix = getHbond(angleHelix, 1, maxLength, tpl, Hdistmin, Hdistmax, Hangle, Hlist, acceptList);
            }
            else {
                printf("Some problem while calculating HbondHelix \n");
                printf("Nturn= %f", static_cast<double>(Nturn));
            }
            
            
            ///////////////////////////////////////
            // save phi, psi, and omega calculation
            ///////////////////////////////////////
            tmp.phi = phi;
            tmp.psi = psi;

            hrow.push_back(tmp);
        }
        Hmatrix.push_back(hrow);
    }
    return 0;
}
int getHelixMap(vector<atom > angleIn, vector<vector<string> > tpl, string link, float delta, vector<vector<vector<helix> > > &Hmatrix){
    int maxLength = 20;  
    ///////////////////////////////
    //loop over omega
    for (long double omega = (Rmin+delta/2); omega<Rmax;omega += delta){
        //make psi string and change angle file
        string somega;
        if (omega<=0.0) somega = std::to_string(360+omega);//convert phi to positive angle
        else somega = std::to_string(omega);

        for (long long int r = 1; r<3; r++){
            string par = std::to_string(r)+" o "+somega;
            int state = changeAngle(angleIn, tpl, "rot", par);
            //if state is >0 some error occured in changeAngle()
            if(state>0){
                return 1;
            }
        }
        vector<vector<helix> > slice; //single row
        ///////////////////////////////
        //loop over psi
        for (long double psi = (Rmax-delta/2); psi>Rmin; psi -= delta){

            //make phi string and change angle file
            string spsi;
            if (psi<=0.0) spsi = std::to_string(360+psi);//convert psi to positive angle
            else spsi = std::to_string(psi);            

            int state = changeAngle(angleIn, tpl, "allpsi", spsi);

                //if state is >0 some error occured in changeAngle()
            if(state>0){
                return 1;
            }

            vector<helix> hrow; //single row

            ///////////////////////////////
            //loop over phi
            for (long double phi = (Rmin+delta/2); phi<Rmax;phi += delta){

                //make psi string and change angle file
                string sphi;
                if (phi<=0.0) sphi = std::to_string(360+phi);//convert phi to positive angle
                else sphi = std::to_string(phi);

                int state = changeAngle(angleIn, tpl, "allphi", sphi);

                //if state is >0 some error occured in changeAngle()
                if(state>0){
                    return 1;
                }

                vector<atom> cnfneu;
                if (fillCnf(angleIn, tpl, cnfneu)==1){
                    cerr << "error in fillCnf()" << endl;
                    return 1;
                }

                ///////////////////////////////////////
                //Helix variable
                ///////////////////////////////////////
                helix tmp;

                vector<atom> res1 = getRes(cnfneu,1);
                vector<atom> res2 = getRes(cnfneu,2);
                
                ///////////////////////////////////////
                // calculate helix parameters
                ///////////////////////////////////////
                getHelixPar(cnfneu, tmp);
                
                
                ///////////////////////////////////////
                //clash 
                ///////////////////////////////////////
                vector<string > clashList {"C1","C2","C3","C4","C5","C6","O1","O2","O3","O4","O5"};

                tmp.clash = getClash(res1, res2, clashList, LJpairs);

                ///////////////////////////////////////
                //clashHelix 
                ///////////////////////////////////////
                int Nturn = (int)floor(abs(tmp.n)+0.5); //number of residues in one turn

                vector<atom> cnftmp; //temporal cnf for clashHelix
                if (Nturn < maxLength && Nturn > 4){
                    long long int nr = Nturn-4;

                    vector<atom> angletmp = angleIn;
                    if(changeAngle(angletmp, tpl, "add", "1 " + link + " " + to_string(nr)) == 1) return 1;
                    if(changeAngle(angletmp, tpl, "allpsi", spsi)>0) return 1;
                    if(changeAngle(angletmp, tpl, "allphi", sphi)>0)return 1;

                    if (fillCnf(angletmp, tpl, cnftmp) == 1) return 1;
                } 
                else {
                    if (fillCnf(angleIn, tpl, cnftmp) == 1) return 1;
                }

                tmp.clashHelix = tmp.clash;
                for (int i = 3; i <= (Nturn+1); i++){
                    vector<atom> resN =  getRes(cnftmp,i);
                    tmp.clashHelix.back().second = tmp.clashHelix.back().second + 
                            getClashHelix(res1, resN, clashList, LJpairs).back().second;
                }
            

                ///////////////////////////////////////
                //Stereoelectronic energy calculation
                ///////////////////////////////////////
                long long int ilink = (long long) stoi(link);

                tmp.SteElecE_just_TOR = getSteElecE_just_torsion(res1, res2, ilink, Torsions);
                
                tmp.SteElecE_just_LJ = getSteElecE_just_LJ(res1, res2, ilink, LJpairs);
                tmp.SteElecE = tmp.SteElecE_just_LJ + tmp.SteElecE_just_TOR;
                
                tmp.SteElecE_just_LJ_extended = getSteElecE_just_LJ_extended(res1, res2, ilink, LJpairs);
                tmp.SteElecE_extended = tmp.SteElecE_just_LJ_extended + tmp.SteElecE_just_TOR;

                ///////////////////////////////////////
                //H bond calculation
                ///////////////////////////////////////
                tmp.Hbond = getHbond(angleIn, tpl, Hdistmin, Hdistmax, Hangle, Hlist, acceptList);
                
                
                //Different behavior depending on Nturn
                if (Nturn<4){ // calculate Hbonds (function knows that if Ntrun==2 it's not helical Hbond)
                    tmp.HbondHelix = getHbond(angleIn, 1, Nturn, tpl, Hdistmin, Hdistmax, Hangle, Hlist, acceptList);
                } 
                else if (Nturn < maxLength && Nturn > 3){ // if the second residue is more than Nturn, make angle longer
                    vector<atom> angleHelix = angleIn;
                    long long int nr = Nturn-3;
                    if(changeAngle(angleHelix, tpl, "add", "1 " + link + " " + to_string(nr)) == 1) return 1;
                    if(changeAngle(angleHelix, tpl, "allpsi", spsi)>0) return 1;
                    if(changeAngle(angleHelix, tpl, "allphi", sphi)>0)return 1;
                    tmp.HbondHelix = getHbond(angleHelix, 1, Nturn, tpl, Hdistmin, Hdistmax, Hangle, Hlist, acceptList);
                }
                else if (Nturn > maxLength){ //if Nturn is big make all Hbonds 0
                    vector<atom> angleHelix = angleIn;
                    long long int nr = maxLength-3;
                    if(changeAngle(angleHelix, tpl, "add", "1 " + link + " " + to_string(nr)) == 1) return 1;
                    if(changeAngle(angleHelix, tpl, "allpsi", spsi)>0) return 1;
                    if(changeAngle(angleHelix, tpl, "allphi", sphi)>0)return 1;
                    tmp.HbondHelix = getHbond(angleHelix, 1, maxLength, tpl, Hdistmin, Hdistmax, Hangle, Hlist, acceptList);
                } 
                else {
                    printf("Some problem while calculating HbondHelix");
                }
                

                ///////////////////////////////////////
                // save phi, psi, and omega calculation
                ///////////////////////////////////////
                tmp.phi = phi;
                tmp.psi = psi;
                tmp.omega = omega;


                hrow.push_back(tmp);
            }
            slice.push_back(hrow);
        }
        Hmatrix.push_back(slice);
    }

    return 0;
}

///////////////////////////////

    //////////////////////////////////////////////////////////////
    // branch for OTHER THEN 1-6 link
    //////////////////////////////////////////////////////////////
int getDisacch(vector<atom > angleIn, vector<vector<string> > tpl, string link, float delta, vector<vector<helix> > &Hmatrix){

    for (long double psi = (Rmax-delta/2); psi>Rmin;psi -= delta){
        //make phi string and change angle file
        string spsi;
        if (psi<=0.0) spsi = std::to_string(360+psi);//convert psi to positive angle
        else spsi = std::to_string(psi);            

        int state = changeAngle(angleIn, tpl, "allpsi", spsi);

        //if state is >0 some error occured in changeAngle()
        if(state>0){
            return 1;
        }

        vector<helix> hrow; //single row
        for (long double phi = (Rmin+delta/2); phi<Rmax;phi += delta){
            //make psi string and change angle file
            string sphi;
            if (phi<=0.0) sphi = std::to_string(360+phi);//convert phi to positive angle
            else sphi = std::to_string(phi);

            int state = changeAngle(angleIn, tpl, "allphi", sphi);

            //if state is >0 some error occured in changeAngle()
            if(state>0){
                return 1;
            }

            vector<atom> cnfneu;
            if (fillCnf(angleIn, tpl, cnfneu)==1){
                cerr << "error in fillCnf()" << endl;
                return 1;
            }
            //printf("phi= %f ; psi = %f \n", static_cast<double>(phi), static_cast<double>(psi) );

            ///////////////////////////////////////
            //Helix calculation calculation
            ///////////////////////////////////////
            helix tmp;

            ///////////////////////////////////////
            //clash calculation
            ///////////////////////////////////////
            vector<atom> res1 = getRes(cnfneu,1);
            vector<atom> res2 = getRes(cnfneu,2);

            vector<string > clashList {"C1","C2","C3","C4","C5","C6","O1","O2","O3","O4","O5"};

            tmp.clash = getClash(res1, res2, clashList, LJpairs);
            //tmp.clash = getClashTest(res1, res2, clashList, LJpairs);

            ///////////////////////////////////////
            //Stereoelectronic energy calculation
            ///////////////////////////////////////
            res1 = getRes(cnfneu,1);
            res2 = getRes(cnfneu,2);

            long long int ilink = (long long) stoi(link);

            tmp.SteElecE_just_TOR = getSteElecE_just_torsion(res1, res2, ilink, Torsions);
            
            tmp.SteElecE_just_LJ = getSteElecE_just_LJ(res1, res2, ilink, LJpairs);
            tmp.SteElecE = tmp.SteElecE_just_LJ + tmp.SteElecE_just_TOR;

            tmp.SteElecE_just_LJ_extended = getSteElecE_just_LJ_extended(res1, res2, ilink, LJpairs);
            tmp.SteElecE_extended = tmp.SteElecE_just_LJ_extended + tmp.SteElecE_just_TOR;


            ///////////////////////////////////////
            //H bond calculation
            ///////////////////////////////////////
            tmp.Hbond = getHbond(angleIn, tpl, Hdistmin, Hdistmax, Hangle, Hlist, acceptList);

            ///////////////////////////////////////
            // save phi, psi, and omega calculation
            ///////////////////////////////////////
            tmp.phi = phi;
            tmp.psi = psi;

            hrow.push_back(tmp);
        }
        Hmatrix.push_back(hrow);
    }
    return 0;
}
//////////////////////////////////////////////////////////////
///       create ramachandran plots   ////////////////////////
///       case for 1-6 link           ////////////////////////
//////////////////////////////////////////////////////////////
int getDisacch(vector<atom > angleIn, vector<vector<string> > tpl, string link, float delta, vector<vector<vector<helix> > > &Hmatrix){

           
    ///////////////////////////////
    //loop over omega
    for (long double omega = (Rmin+delta/2); omega<Rmax;omega += delta){
        //make psi string and change angle file
        string somega;
        if (omega<=0.0) somega = std::to_string(360+omega);//convert phi to positive angle
        else somega = std::to_string(omega);

        for (long long int r = 1; r<3; r++){
            string par = std::to_string(r)+" o "+somega;
            int state = changeAngle(angleIn, tpl, "rot", par);
            //if state is >0 some error occured in changeAngle()
            if(state>0){
                return 1;
            }
        }
        vector<vector<helix> > slice; //single row
        ///////////////////////////////
        //loop over psi
        for (long double psi = (Rmax-delta/2); psi>Rmin; psi -= delta){

            //make phi string and change angle file
            string spsi;
            if (psi<=0.0) spsi = std::to_string(360+psi);//convert psi to positive angle
            else spsi = std::to_string(psi);            

            int state = changeAngle(angleIn, tpl, "allpsi", spsi);

                //if state is >0 some error occured in changeAngle()
            if(state>0){
                return 1;
            }

            vector<helix> hrow; //single row

            ///////////////////////////////
            //loop over phi
            for (long double phi = (Rmin+delta/2); phi<Rmax;phi += delta){

                //make psi string and change angle file
                string sphi;
                if (phi<=0.0) sphi = std::to_string(360+phi);//convert phi to positive angle
                else sphi = std::to_string(phi);

                int state = changeAngle(angleIn, tpl, "allphi", sphi);

                //if state is >0 some error occured in changeAngle()
                if(state>0){
                    return 1;
                }

                vector<atom> cnfneu;
                if (fillCnf(angleIn, tpl, cnfneu)==1){
                    cerr << "error in fillCnf()" << endl;
                    return 1;
                }
                //mkCnfOpt(cnfneu,tpl,"data/"+anomer+"_"+sugar+"_"+"1-"+link+"_phi"+sphi+"_psi"+spsi+"_omega"+somega+".cnf");
                //printf("phi= %f ; psi = %f ; omega = %f \n", static_cast<double>(phi), static_cast<double>(psi), static_cast<double>(omega) );

                ///////////////////////////////////////
                //Helix variable
                ///////////////////////////////////////
                helix tmp;

                vector<atom> res1 = getRes(cnfneu,1);
                vector<atom> res2 = getRes(cnfneu,2);

                ///////////////////////////////////////
                //clash calculation
                ///////////////////////////////////////
                vector<string > clashList {"C1","C2","C3","C4","C5","C6","O1","O2","O3","O4","O5"};

                tmp.clash = getClash(res1, res2, clashList, LJpairs);

                ///////////////////////////////////////
                //Stereoelectronic energy calculation
                ///////////////////////////////////////
                long long int ilink = (long long) stoi(link);

                tmp.SteElecE_just_TOR = getSteElecE_just_torsion(res1, res2, ilink, Torsions);
                
                tmp.SteElecE_just_LJ = getSteElecE_just_LJ(res1, res2, ilink, LJpairs);
                tmp.SteElecE = tmp.SteElecE_just_LJ + tmp.SteElecE_just_TOR;
                
                tmp.SteElecE_just_LJ_extended = getSteElecE_just_LJ_extended(res1, res2, ilink, LJpairs);
                tmp.SteElecE_extended = tmp.SteElecE_just_LJ_extended + tmp.SteElecE_just_TOR;

                ///////////////////////////////////////
                //H bond calculation
                ///////////////////////////////////////
                tmp.Hbond = getHbond(angleIn, tpl, Hdistmin, Hdistmax, Hangle, Hlist, acceptList);

                ///////////////////////////////////////
                // save phi, psi, and omega calculation
                ///////////////////////////////////////
                tmp.phi = phi;
                tmp.psi = psi;
                tmp.omega = omega;


                hrow.push_back(tmp);
            }
            slice.push_back(hrow);
        }
        Hmatrix.push_back(slice);
    }

    return 0;
}


//calculation of cartesians coordinate
//vec q is vector of distance(x), angle(y) and torsion(z)
double calcZ(const atom q,const atom c,const atom b,const atom a){
    double Z, nom,  den; 
    int s = (a.pos.y > 0) - (a.pos.y < 0); //signum
    
    nom = (-pow(b.pos.x,2) + pow(calcDistance(b.pos,c.pos),2)*pow(cos(q.pos.y*PI/180),2));
    nom = nom*(-1+pow(1/cos(q.pos.z*PI/180),2)*pow(s,2));

    nom = sqrt(-nom);    // Martin's "i" taken into account
    nom = q.pos.x* cos(q.pos.z*PI/180)*nom; //So I remove minus here. RIGHT?
    den = b.pos.x * s;
    Z = nom/den;
    
    double z = fmod(q.pos.z,360.0);
    //empirical fixing, don't know where in the calculation is a mistake
    if(z > 0.0 && z < 90.0 || z > 180.0 && z < 270.0){
        Z = -Z;
    }
    return Z;
}
double calcY(const atom q,const atom c,const atom b,const atom a){
    double Y;
    
    //calculated with z expression
    Y = -pow(b.pos.x,2)+pow(calcDistance(c.pos,b.pos),2)*pow(cos(q.pos.y*PI/180),2);
    Y = Y*pow(cos(q.pos.z*PI/180),2);

    Y = sqrt(-Y); // Martin's "i" taken into account
    
    Y = Y*q.pos.x/b.pos.x; 
    
    double z = fmod(q.pos.z,360.0);
    //HERE IS PROBLEM WITH SIGN OF "Y" REGARDING TORSION ANGLE (IS IT JUST IN CASE OF 0 AND 180?)
    if (z > 179.9){
        Y = -Y;
    }
    //empirical fixing, don't know where in the calculation is a mistake
    if(z > 90.0 && z < 180.0 || z > 270.0 && z < 360.0){
        Y = -Y;
    }
    return Y;
}
double calcX(const atom q,const atom c,const atom b,const atom a){
    double X;
    X = calcDistance(c.pos,b.pos)*q.pos.x*sqrt(pow(cos(q.pos.y*PI/180),2))/sqrt(pow(b.pos.x,2));
    
    double y = fmod(q.pos.y,180.0);
    if (y < 90.0){
        X = -X;
    }
    
    return X;
}

void transform(atom &a, atom &b, atom &c){
    translate(a,b,c);
    rotate1(a,b,c, false);
    rotate2(a,b,c, false);
}
void translate(atom &a, atom &b, atom &c){
    a.pos.x -= c.pos.x;
    a.pos.y -= c.pos.y;
    a.pos.z -= c.pos.z;
    b.pos.x -= c.pos.x;
    b.pos.y -= c.pos.y;
    b.pos.z -= c.pos.z;
    c.pos.x -= c.pos.x;
    c.pos.y -= c.pos.y;
    c.pos.z -= c.pos.z;
}
//a is the third neighbor, b is second neighbor, c is first neighbor
//http://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
void rotate1(atom &a, atom &b, atom &c, bool inverse){
    double rotM1[3][3];
    getMatrix1(a, b, c, inverse, rotM1);
    a.pos = getMatMult(rotM1,a.pos);
    b.pos = getMatMult(rotM1,b.pos);
    c.pos = getMatMult(rotM1,c.pos);
}
void rotate2(atom &a, atom &b, atom &c, bool inverse){
    double rotM2[3][3];
    getMatrix2(a, b, c, inverse, rotM2);
    a.pos = getMatMult(rotM2,a.pos);
    b.pos = getMatMult(rotM2,b.pos);
    c.pos = getMatMult(rotM2,c.pos);
}
void getMatrix1(atom const a, atom const b, atom const c, bool inverse, double rotM1[3][3]){
    vec m = getCrossofVec(getSubstrVec(b.pos,a.pos),getSubstrVec(b.pos,c.pos));
    vec n = {0.0, 0.0, 1.0};
    vec RotAxis = getMultipVec(getCrossofVec(m,n),1/getAbsVec(getCrossofVec(m,n)));
    double theta = acos(getDotofVec(m,n)/(getAbsVec(m)*getAbsVec(n)));
    if(inverse) theta=-theta;
    double co = cos(theta);
    double s = sin(theta);
    double CO = 1-co;
    rotM1[0][0] = pow(RotAxis.x,2)*CO+co;
    rotM1[0][1] = RotAxis.x*RotAxis.y*CO-RotAxis.z*s;
    rotM1[0][2] = RotAxis.x*RotAxis.z*CO+RotAxis.y*s;
    rotM1[1][0] = RotAxis.x*RotAxis.y*CO+RotAxis.z*s;
    rotM1[1][1] = pow(RotAxis.y,2)*CO+co;
    rotM1[1][2] = RotAxis.y*RotAxis.z*CO-RotAxis.x*s;
    rotM1[2][0] = RotAxis.x*RotAxis.z*CO-RotAxis.y*s;
    rotM1[2][1] = RotAxis.y*RotAxis.z*CO+RotAxis.x*s;
    rotM1[2][2] = pow(RotAxis.z,2)*CO+co;
}
void getMatrix2(atom const a, atom const b, atom const c, bool inverse, double rotM2[3][3]){
    vec z = {1.0 , 0.0, 0.0};
    vec fr;
    double alpha;
    fr = getSubstrVec(c.pos,b.pos);
    alpha = getDotofVec(fr,z);
    alpha = alpha/(getAbsVec(getSubstrVec(c.pos,b.pos))*getAbsVec(z));
    alpha = acos(alpha);
    alpha = alpha-PI;
    if (b.pos.y > 0) alpha = -alpha;
    //cout << "alpha " << alpha << endl;
    if(inverse) alpha=-alpha;
    rotM2[0][0] = cos(alpha);
    rotM2[0][1] = -sin(alpha);
    rotM2[0][2] = 0.0;
    rotM2[1][0] = sin(alpha);
    rotM2[1][1] = cos(alpha);
    rotM2[1][2] = 0.0;
    rotM2[2][0] = 0.0;
    rotM2[2][1] = 0.0;
    rotM2[2][2] = 1.0;
}
vec getMatMult(double rotM[3][3],const vec v){
    vec a;
    a.x = (rotM[0][0]*v.x+rotM[0][1]*v.y+rotM[0][2]*v.z);
    a.y = (rotM[1][0]*v.x+rotM[1][1]*v.y+rotM[1][2]*v.z);
    a.z = (rotM[2][0]*v.x+rotM[2][1]*v.y+rotM[2][2]*v.z);
    return a;
}
void transformBack(atom const A, atom const B, atom const C, atom &d){
    atom a,b,c;
    a=A;
    b=B;
    c=C;
    double rotM1[3][3];
    double rotM2[3][3];
    translate(a,b,c);
    getMatrix1(a,b,c, true, rotM1);
    rotate1(a,b,c, false);
    getMatrix2(a,b,c, true, rotM2);
    //now transform d
    d.pos = getMatMult(rotM2,d.pos);
    d.pos = getMatMult(rotM1,d.pos);
    d.pos.x += C.pos.x;
    d.pos.y += C.pos.y;
    d.pos.z += C.pos.z;
}

//this function is just used under -db option
void transformtest(atom &a, atom &b, atom &c, atom &d){
    translatetest(a,b,c,d);
    rotate1test(a,b,c,d, false);
    rotate2test(a,b,c,d, false);
}
void translatetest(atom &a, atom &b, atom &c, atom &d){
    a.pos.x -= c.pos.x;
    a.pos.y -= c.pos.y;
    a.pos.z -= c.pos.z;
    b.pos.x -= c.pos.x;
    b.pos.y -= c.pos.y;
    b.pos.z -= c.pos.z;
    d.pos.x -= c.pos.x;
    d.pos.y -= c.pos.y;
    d.pos.z -= c.pos.z;
    c.pos.x -= c.pos.x;
    c.pos.y -= c.pos.y;
    c.pos.z -= c.pos.z;
    
}
void rotate2test(atom &a, atom &b, atom &c,atom &d, bool inverse){
    double rotM2[3][3];
    getMatrix2(a, b, c, inverse, rotM2);
    a.pos = getMatMult(rotM2,a.pos);
    b.pos = getMatMult(rotM2,b.pos);
    c.pos = getMatMult(rotM2,c.pos);
    d.pos = getMatMult(rotM2,d.pos);
}
void rotate1test(atom &a, atom &b, atom &c,atom &d, bool inverse){
    double rotM1[3][3];
    getMatrix1(a, b, c, inverse, rotM1);
    a.pos = getMatMult(rotM1,a.pos);
    b.pos = getMatMult(rotM1,b.pos);
    c.pos = getMatMult(rotM1,c.pos);
    d.pos = getMatMult(rotM1,d.pos);
}


//this functions are not used
void getDistance(vector<vector<string> > tpl, vector<atom> cnf){
    for(int i=0; i<17; i++){
        if (tpl[i][2].compare("d")==0){
            int a1 = atoi(tpl[i][1].c_str())-1;
            vec a = {cnf[i].pos.x, cnf[i].pos.y, cnf[i].pos.z};
            vec b = {cnf[a1].pos.x, cnf[a1].pos.y, cnf[a1].pos.z};
            char buffer [50];
            sprintf (buffer, "%f", calcDistance(a,b));
            tpl[i][2] = buffer;
        }
        else {
            char buffer [50];
            sprintf (buffer, "%f", cnf[i].pos.x);
            tpl[i][2] = buffer;
        }
    }
}
void getAngle(vector<vector<string> > tpl, vector<atom> cnf){
    for(int i=0; i<17; i++){
        if (tpl[i][3].compare("a")==0){
            int a1 = atoi(tpl[i][1].c_str())-1;
            int a2 = atoi(tpl[a1][1].c_str())-1;
            //cout << i+1 << " " << a1+1 << " " << a2+1 << endl;
            vec a = {cnf[i].pos.x, cnf[i].pos.y, cnf[i].pos.z};
            vec b = {cnf[a1].pos.x, cnf[a1].pos.y, cnf[a1].pos.z};
            vec c = {cnf[a2].pos.x, cnf[a2].pos.y, cnf[a2].pos.z};
            char buffer [50];
            sprintf (buffer, "%f", calcAngle(a,b,c));
            tpl[i][3] = buffer;
        }
        else {
            char buffer [50];
            sprintf (buffer, "%f", cnf[i].pos.y);
            tpl[i][3] = buffer;
        }
    }
}
void getDihedral(vector<vector<string> > tpl, vector<atom> cnf){
    for(int i=0; i<17; i++){
        if (tpl[i][4].compare("t")==0){
            int a1 = atoi(tpl[i][1].c_str())-1;
            int a2 = atoi(tpl[a1][1].c_str())-1;
            int a3 = atoi(tpl[a2][1].c_str())-1;
            //cout << i+1 << " " << a1+1 << " " << a2+1 << " " << a3+1 << endl;
            vec a = {cnf[i].pos.x, cnf[i].pos.y, cnf[i].pos.z};
            vec b = {cnf[a1].pos.x, cnf[a1].pos.y, cnf[a1].pos.z};
            vec c = {cnf[a2].pos.x, cnf[a2].pos.y, cnf[a2].pos.z};
            vec d = {cnf[a3].pos.x, cnf[a3].pos.y, cnf[a3].pos.z};
            char buffer [50];
            sprintf (buffer, "%f", calcDihedral(a,b,c,d));
            tpl[i][4] = buffer;
        }
        else {
            char buffer [50];
            sprintf (buffer, "%f", cnf[i].pos.z);
            tpl[i][4] = buffer;
        }
    }
}
double chop(double d){
    double back;
    if(d<0.00000000001 && d>-0.00000000001) back = 0;
    else back = d;
    return back;
}
float round(float f,float prec){
    return (float) (floor(f*(1.0/prec) + 0.5)/(1.0/prec));
}

