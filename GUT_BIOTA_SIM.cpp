


//
//  main.cpp
//  plasmid
//
//  Created by Burcu Tepekule on 25/10/17.
//  Copyright © 2017 Burcu Tepekule. All rights reserved.
//
# pragma warning(push, 0)
# include <iostream>
# include <iomanip>
# include <cmath>
# include <stdio.h>
# include <random>
# include <fstream>
# include <vector>
# include <tgmath.h>
# include <cstdlib>
# include <string>
# include <ctime>
# include <iomanip>
# include <stack>
# include <math.h>
# include <stdlib.h>
# include <algorithm>
# include <vector>
# include <map>
# include <set>
# include <sstream>
# include <chrono>
# include <sys/types.h>
# include <sys/stat.h>
# include <sys/time.h>
# include <boost/array.hpp>
# include <boost/numeric/odeint.hpp>
# include <boost/numeric/ublas/vector.hpp>
# include <boost/numeric/ublas/matrix.hpp>
# include <boost/numeric/ublas/matrix_proxy.hpp>
# include <boost/numeric/ublas/io.hpp>
# include <boost/numeric/ublas/storage.hpp>
# include <boost/numeric/ublas/assignment.hpp>
# include <boost/fusion/algorithm/transformation/push_back.hpp>
# include <boost/fusion/include/push_back.hpp>
# include <gsl/gsl_math.h>
# include <gsl/gsl_rng.h>
# include <gsl/gsl_randist.h>

using namespace std;
namespace odeint = boost::numeric::odeint;
namespace ublas  = boost::numeric::ublas;
//namespace plt    = matplotlibcpp;
typedef ublas::vector<double> popStates;
int numOfPhyla   = 4;
int numOfSpecies = 2*numOfPhyla;
popStates x(numOfSpecies);
popStates xInt(2);
popStates extinct(2);
popStates c0counter(2);
popStates g0(1);
ublas::vector<double> compMult(numOfSpecies);
ublas::vector<double> conjMult (numOfSpecies);
ublas::vector<double> xt(numOfSpecies+1);
ublas::vector<double> xICSum(numOfPhyla);
ublas::vector<double> ds(numOfPhyla);
ublas::vector<double> dr(numOfPhyla);
ublas::vector<double> is(numOfPhyla);
ublas::vector<double> ir(numOfPhyla);
ublas::vector<double> drugKillRates (numOfSpecies);
ublas::vector<double> infAbundVec (numOfSpecies);
ublas::vector<double> rVec (numOfSpecies);
ublas::vector<double> aVec (numOfPhyla*numOfPhyla);
ublas::matrix<double> hMat (numOfPhyla,numOfPhyla);
ublas::vector<double> sampleTimes(72);
ublas::vector<double> samplePops(72); // [time, C_{0}, C_{0}^{+}, C_{1}, C_{1}^{+}]
ublas::vector<double> hVec (numOfPhyla*numOfPhyla);
ublas::vector<double> h (4);
ublas::vector<double> resFracs(numOfPhyla);
ublas::vector<double> resAbs(numOfSpecies);
ublas::vector<double> otherParams(4);
ublas::vector<double> parameterVector(5);
ublas::vector<double> c(numOfSpecies);
ublas::vector<double> gama(numOfSpecies);
ublas::vector<int> trtLengthVector(4);
ublas::vector<int> trtCountVector(50);
ublas::vector<int> infCountVector(50);
ublas::vector<int> trtInit(50);
ublas::vector<int> trtDura(50);
ublas::vector<int> infInit(50);
ublas::matrix<int> tspanTrt(100,3);
ublas::matrix<int> tspanInf(100,3);
ublas::matrix<int> epsMat(88,8);

const double popScale = 1e13;
const double dt  = 0.01; // dt=1 is 1.12 days (can take as 1)
double t   = 0.0;
double mu;
double kappa;
double h0,h1;
double Ca;
double Cx;
double sol01,sol02,sol03,sol04;
int drugPressure;
int inf;
int save;
int trtTime;
int obsTime;
int trtLength;
int trtPeriod;
int totalTrts;
int totalInfs;
int numReactions_c;
int numReactions_d;
int updateExt;
string folderChar;
////// FOR HYBRID GILLESPIE /////
ublas::vector<double> Md (88);
ublas::vector<double> Mc (88);
ublas::vector<double> propensities(88);
ublas::vector<double> deltaNormVec(88);

class RandomDevice{
private:
    unsigned long rand_seed;
    default_random_engine engine;
public:
    RandomDevice(unsigned long n);
    int randBinomInt(double p);
    int randPoissInt(double lam);
    double randGauss0 (double mu);
    double randUniform01 (double a,double b);
};

RandomDevice::RandomDevice(unsigned long n){
    rand_seed = n;
    engine.seed(n);
}

int RandomDevice::randBinomInt(double p){
    binomial_distribution<int> distribution(12,p);
    return distribution(engine);
}

int RandomDevice::randPoissInt(double lam){
    poisson_distribution<int> distribution(lam);
    return distribution(engine);
}

double RandomDevice::randGauss0 (double mu){
    normal_distribution<double> distribution(mu,0);
    return distribution(engine);
}

double RandomDevice::randUniform01 (double a, double b){
    uniform_real_distribution<double> distribution(a,b);
    return distribution(engine);
}

const string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d, %X", &tstruct);
    return buf;
}

popStates addStatVecs(popStates s1, popStates s2){
    int sizeStateVec = (int)s1.size();
    popStates res(sizeStateVec);
    for (int i=0; i<sizeStateVec; i++){
        res(i) = s1(i)+s2(i);
    }
    return res;
}

popStates subStatVecs(popStates s1, popStates s2){
    int sizeStateVec = (int)s1.size();
    popStates res(sizeStateVec);
    for (int i=0; i<sizeStateVec; i++){
        res(i) = s1(i)-s2(i);
    }
    return res;
}

int delNormCalc(popStates deltaVec){
    int ret = 0;
    for (int i=0;i<deltaVec.size();i++){
        ret = ret + abs(deltaVec(i));
    }
    return ret;
}

class Reaction {
public:
    double rate;
    popStates in, out, delta;
    int deltaNorm;
    
    // Constructors:
    Reaction(popStates p_in, popStates p_out, double p_rate){
        rate  = p_rate;
        in    = p_in;
        out   = p_out;
        delta = subStatVecs(out,in);
        deltaNorm = delNormCalc(delta);
        
    }
    Reaction() {};
    
    // Calculate reaction propensity for given state:
    double propensity(popStates x)
    {
        double a = rate;
        
        for (int i=0; i<x.size(); i++) {
            for (int m=0; m<in(i); m++)
                a *= x(i)-(double)m;
        }
        
        if (a<0)
            return 0;
        
        return a;
    }
    
    // Implement reaction on given state:
    popStates implement(popStates x)
    {
        return addStatVecs(x,delta);
    }
};
//////////////////////////

template <typename T>
bool my_isnan(const T x)
{
#if __cplusplus >= 201103L
    using std::isnan;
#endif
    return isnan(x);
}

int factorial(int n){
    int nFac = 1;
    for(unsigned i=0; i<n-1; ++i){
        nFac = nFac*(n-i);
    }
    return nFac;
}

int choose(int n, int k){
    if(n==k){
        return 1;
    }else{
        return factorial(n)/(factorial(k)*factorial(n-k));
    }
}

double unifRand(double a, double b){
    return (b-a)*((double) rand() / (RAND_MAX))+a;
}

void setRNG(){
    const gsl_rng_type * T2;
    gsl_rng * r2;
    srand( static_cast<unsigned int>(time(NULL)));
    time_t Seed2 = rand();
    gsl_rng_env_setup();
    T2 = gsl_rng_default;
    r2 = gsl_rng_alloc (T2);
    gsl_rng_set (r2, Seed2);
}

template<typename T>
void resizeVec(ublas::vector<T> &m, int len){
    m.resize(len,1);
}

template<typename T>
void resizeMat(ublas::matrix<T> &m, int len){
    m.resize(2*len+1,3);
}

string setDirectoriesFullRandomScheduling(string saveDirectory, int trtTime, int totalTrts, int totalInfs, int save){
    stringstream folderName;
    string folderChar;
    if(totalInfs==0){
        folderName << saveDirectory<< "TL_"<<trtTime<<"_N_"<<totalTrts; //
    }else{
        folderName << saveDirectory<< "TL_"<<trtTime<<"_N_"<<totalTrts<<"_RC_"<<totalInfs; // In case recolonization events are added
    }
    folderChar = folderName.str();
    if(save==1){mkdir(folderChar.c_str(),0777);}
    return folderChar;
}

template <typename T>
ublas::vector<T> mixVecsAbs (const ublas::vector<T> &vs, const ublas::vector<T> &vr){
    double k=0;
    ublas::vector<double> vReturn(numOfSpecies);
    for(unsigned i=0; i<numOfPhyla; ++i){
        vReturn(k)=vs(i)-vr(i); k++;
        vReturn(k)=vr(i); k++;
    }
    return vReturn;
}

template <typename T>
ublas::vector<T> mixVecsFrac (const ublas::vector<T> &vs, const ublas::vector<T> &vr){
    double k=0;
    ublas::vector<double> vReturn(numOfSpecies);
    for(unsigned i=0; i<numOfPhyla; ++i){
        vReturn(k)=vs(i)*(1-vr(i)); k++;
        vReturn(k)=vs(i)*vr(i); k++;
    }
    return vReturn;
}

template <typename T>
ublas::vector<T> mixVecs (const ublas::vector<T> &vs, const ublas::vector<T> &vr){
    double k=0;
    ublas::vector<double> vReturn(numOfSpecies);
    for(unsigned i=0; i<numOfPhyla; ++i){
        vReturn(k)=vs(i); k++;
        vReturn(k)=vr(i); k++;
    }
    return vReturn;
}

template <typename T>
void saveUblasVecAppend(const ublas::vector<T> &v, string folderName, string fileName, const char * precision, int save){
    if(save==1){
        double *tempVec = new double[v.size()]; // TO SAVE THE DATA
        for(unsigned k=0;k<v.size();++k){tempVec[k]=v(k);}
        
        FILE *result;
        stringstream stream;
        stringstream myRoot;
        myRoot << folderName<< "/";
        stream.str(string());
        stream << myRoot.str() <<fileName<<".txt";
        fileName = stream.str();
        result = fopen(fileName.c_str(),"a+");
        for(unsigned i=0;i<v.size();++i){
            fprintf(result, precision, tempVec[i]);
        }
        fprintf(result,"\n ");
        fclose(result);
        delete[] tempVec;
    }
}


void sysNP(popStates &x, popStates &dxdt , double t )
{
    for (unsigned i=0;i<x.size();++i){
        if(i%2==0){// SENSITIVE VERSION OF THE SAME SPECIES
            dxdt(i) =  rVec(i)*x(i)+rVec(i)*(1-c(i))*gama(i)*(1-kappa)*x(i+1)+compMult(i)*x(i)-x(i)*conjMult(i)-drugPressure*drugKillRates(i)*x(i);
        }else{ // RESISTANT VERSION OF THE SAME SPECIES
            dxdt(i) =  rVec(i)*(1-c(i))*(1-gama(i))*x(i)+compMult(i)*x(i)+x(i-1)*conjMult(i)-drugPressure*drugKillRates(i)*x(i);
        }
    }
}

template <typename T>
popStates implementReaction(ublas::matrix<T> &tspanTrt, int reactIdx){
    ////////////////////////////////////////////// REACTIONS ////////////////////////////////////////////////
    
    // in(0) -> C0
    // in(1) -> C0+
    // in(2) -> C1
    // in(3) -> C1+
    // in(4) -> C2
    // in(5) -> C2+ (=0)
    // in(6) -> C3
    // in(7) -> C3+ (=0)
    
    ///////// REACTTIONS ///////////
    Reaction reactions[88];
    popStates in(numOfSpecies), out(numOfSpecies);
    /////// BIRTH PROCESSES ////////
    // REACTION 0 -> BIRTH OF C0
    in  <<= 1,0,0,0,0,0,0,0;
    out <<= 2,0,0,0,0,0,0,0;
    reactions[0] = Reaction(in, out, rVec(0));
    // REACTION 1 -> BIRTH OF C0+
    in  <<= 0,1,0,0,0,0,0,0;
    out <<= 0,2,0,0,0,0,0,0;
    reactions[1] = Reaction(in, out, rVec(1)*(1-c(1))*(1-gama(1)));
    // REACTION 2 -> BIRTH OF C1
    in  <<= 0,0,1,0,0,0,0,0;
    out <<= 0,0,2,0,0,0,0,0;
    reactions[2] = Reaction(in, out, rVec(2));
    // REACTION 3 -> BIRTH OF C1+
    in  <<= 0,0,0,1,0,0,0,0;
    out <<= 0,0,0,2,0,0,0,0;
    reactions[3] = Reaction(in, out, rVec(3)*(1-c(3))*(1-gama(3)));
    // REACTION 4 -> BIRTH OF C2
    in  <<= 0,0,0,0,1,0,0,0;
    out <<= 0,0,0,0,2,0,0,0;
    reactions[4] = Reaction(in, out, rVec(4));
    // REACTION 5 -> BIRTH OF C2+
    in  <<= 0,0,0,0,0,1,0,0;
    out <<= 0,0,0,0,0,2,0,0;
    reactions[5] = Reaction(in, out, rVec(5)*(1-c(5))*(1-gama(5)));
    // REACTION 6 -> BIRTH OF C3
    in  <<= 0,0,0,0,0,0,1,0;
    out <<= 0,0,0,0,0,0,2,0;
    reactions[6] = Reaction(in, out, rVec(6));
    // REACTION 7 -> BIRTH OF C3+
    in  <<= 0,0,0,0,0,0,0,1;
    out <<= 0,0,0,0,0,0,0,2;
    reactions[7] = Reaction(in, out, rVec(7)*(1-c(7))*(1-gama(7)));
    /////// INTERACTION PROCESSES //////
    // TAKE THE ABS VALUE OF AVEC //
    for (int i=0;i<numOfPhyla*numOfPhyla;i++){
        if(aVec(i)<0){
            aVec(i)=-1*aVec(i);
        }
    }
    //// INTRA COMPETITION AMONG C0 ////////////////////////////////////////////////////////////////////
    // REACTION 8 -> COMPETITION C0 + C0 -> C0
    in  <<= 2,0,0,0,0,0,0,0;
    out <<= 1,0,0,0,0,0,0,0;
    reactions[8] = Reaction(in, out, aVec(0));
    // REACTION 9 -> COMPETITION C0 + C0+ -> C0
    in  <<= 1,1,0,0,0,0,0,0;
    out <<= 1,0,0,0,0,0,0,0;
    reactions[9] = Reaction(in, out, aVec(0));
    // REACTION 10 -> COMPETITION C0 + C0+ -> C0+
    in  <<= 1,1,0,0,0,0,0,0;
    out <<= 0,1,0,0,0,0,0,0;
    reactions[10] = Reaction(in, out, aVec(0));
    // REACTION 11 -> COMPETITION C0+ + C0+ -> C0+
    in  <<= 0,2,0,0,0,0,0,0;
    out <<= 0,1,0,0,0,0,0,0;
    reactions[11] = Reaction(in, out, aVec(0));
    //// INTRA COMPETITION AMONG C1 ////////////////////////////////////////////////////////////////////
    // REACTION 12 -> COMPETITION C1 + C1 -> C1
    in  <<= 0,0,2,0,0,0,0,0;
    out <<= 0,0,1,0,0,0,0,0;
    reactions[12] = Reaction(in, out, aVec(5));
    // REACTION 13 -> COMPETITION C1 + C1+ -> C1
    in  <<= 0,0,1,1,0,0,0,0;
    out <<= 0,0,1,0,0,0,0,0;
    reactions[13] = Reaction(in, out, aVec(5));
    // REACTION 14 -> COMPETITION C1 + C1+ -> C1+
    in  <<= 0,0,1,1,0,0,0,0;
    out <<= 0,0,0,1,0,0,0,0;
    reactions[14] = Reaction(in, out, aVec(5));
    // REACTION 15 -> COMPETITION C1+ + C1+ -> C1+
    in  <<= 0,0,0,2,0,0,0,0;
    out <<= 0,0,0,1,0,0,0,0;
    reactions[15] = Reaction(in, out, aVec(5));
    //// INTRA COMPETITION AMONG C2 ////////////////////////////////////////////////////////////////////
    // REACTION 16 -> COMPETITION C2 + C2 -> C2
    in  <<= 0,0,0,0,2,0,0,0;
    out <<= 0,0,0,0,1,0,0,0;
    reactions[16] = Reaction(in, out, aVec(10));
    // REACTION 17 -> COMPETITION C2 + C2+ -> C2
    in  <<= 0,0,0,0,1,1,0,0;
    out <<= 0,0,0,0,1,0,0,0;
    reactions[17] = Reaction(in, out, aVec(10));
    // REACTION 18 -> COMPETITION C2 + C2+ -> C2+
    in  <<= 0,0,0,0,1,1,0,0;
    out <<= 0,0,0,0,0,1,0,0;
    reactions[18] = Reaction(in, out, aVec(10));
    // REACTION 19 -> COMPETITION C2+ + C2+ -> C2+
    in  <<= 0,0,0,0,0,2,0,0;
    out <<= 0,0,0,0,0,1,0,0;
    reactions[19] = Reaction(in, out, aVec(10));
    //// INTRA COMPETITION AMONG C3 ////////////////////////////////////////////////////////////////////
    // REACTION 20 -> COMPETITION C3 + C3 -> C3
    in  <<= 0,0,0,0,0,0,2,0;
    out <<= 0,0,0,0,0,0,1,0;
    reactions[20] = Reaction(in, out, aVec(15));
    // REACTION 21 -> COMPETITION C3 + C3+ -> C3
    in  <<= 0,0,0,0,0,0,1,1;
    out <<= 0,0,0,0,0,0,1,0;
    reactions[21] = Reaction(in, out, aVec(15));
    // REACTION 22 -> COMPETITION C3 + C3+ -> C3+
    in  <<= 0,0,0,0,0,0,1,1;
    out <<= 0,0,0,0,0,0,0,1;
    reactions[22] = Reaction(in, out, aVec(15));
    // REACTION 23 -> COMPETITION C3+ + C3+ -> C3+
    in  <<= 0,0,0,0,0,0,0,2;
    out <<= 0,0,0,0,0,0,0,1;
    reactions[23] = Reaction(in, out, aVec(15));
    //// INTER COMPETITION C1->C0 ////////////////////////////////////////////////////////////////////
    // REACTION 24 -> COMPETITION C0 + C1 -> C1
    in  <<= 1,0,1,0,0,0,0,0;
    out <<= 0,0,1,0,0,0,0,0;
    reactions[24] = Reaction(in, out, aVec(1));
    // REACTION 25 -> COMPETITION C0 + C1+ -> C1+
    in  <<= 1,0,0,1,0,0,0,0;
    out <<= 0,0,0,1,0,0,0,0;
    reactions[25] = Reaction(in, out, aVec(1));
    // REACTION 26 -> COMPETITION C0+ + C1 -> C1
    in  <<= 0,1,1,0,0,0,0,0;
    out <<= 0,0,1,0,0,0,0,0;
    reactions[26] = Reaction(in, out, aVec(1));
    // REACTION 27 -> COMPETITION C0+ + C1+ -> C1+
    in  <<= 0,1,0,1,0,0,0,0;
    out <<= 0,0,0,1,0,0,0,0;
    reactions[27] = Reaction(in, out, aVec(1));
    //// INTER COMPETITION C2->C0 ////////////////////////////////////////////////////////////////////
    // REACTION 28 -> COMPETITION C0 + C2 -> C2
    in  <<= 1,0,0,0,1,0,0,0;
    out <<= 0,0,0,0,1,0,0,0;
    reactions[28] = Reaction(in, out, aVec(2));
    // REACTION 29 -> COMPETITION C0 + C2+ -> C2+
    in  <<= 1,0,0,0,0,1,0,0;
    out <<= 0,0,0,0,0,1,0,0;
    reactions[29] = Reaction(in, out, aVec(2));
    // REACTION 30 -> COMPETITION C0+ + C2 -> C2
    in  <<= 0,1,0,0,1,0,0,0;
    out <<= 0,0,0,0,1,0,0,0;
    reactions[30] = Reaction(in, out, aVec(2));
    // REACTION 31 -> COMPETITION C0+ + C2+ -> C2+
    in  <<= 0,1,0,0,0,1,0,0;
    out <<= 0,0,0,0,0,1,0,0;
    reactions[31] = Reaction(in, out, aVec(2));
    //// INTER COMPETITION C3->C0 ////////////////////////////////////////////////////////////////////
    // REACTION 32 -> COMPETITION C0 + C3 -> C3
    in  <<= 1,0,0,0,0,0,1,0;
    out <<= 0,0,0,0,0,0,1,0;
    reactions[32] = Reaction(in, out, aVec(3));
    // REACTION 33 -> COMPETITION C0 + C3+ -> C3+
    in  <<= 1,0,0,0,0,0,0,1;
    out <<= 0,0,0,0,0,0,0,1;
    reactions[33] = Reaction(in, out, aVec(3));
    // REACTION 34 -> COMPETITION C0+ + C3 -> C3
    in  <<= 0,1,0,0,0,0,1,0;
    out <<= 0,0,0,0,0,0,1,0;
    reactions[34] = Reaction(in, out, aVec(3));
    // REACTION 35 -> COMPETITION C0+ + C3+ -> C3+
    in  <<= 0,1,0,0,0,0,0,1;
    out <<= 0,0,0,0,0,0,0,1;
    reactions[35] = Reaction(in, out, aVec(3));
    //// INTER COMPETITION C0->C1 ////////////////////////////////////////////////////////////////////
    // REACTION 36 -> COMPETITION C0 + C1 -> C0
    in  <<= 1,0,1,0,0,0,0,0;
    out <<= 1,0,0,0,0,0,0,0;
    reactions[36] = Reaction(in, out, aVec(4));
    // REACTION 37 -> COMPETITION C0 + C1+ -> C0
    in  <<= 1,0,0,1,0,0,0,0;
    out <<= 1,0,0,0,0,0,0,0;
    reactions[37] = Reaction(in, out, aVec(4));
    // REACTION 38 -> COMPETITION C0+ + C1 -> C0+
    in  <<= 0,1,1,0,0,0,0,0;
    out <<= 0,1,0,0,0,0,0,0;
    reactions[38] = Reaction(in, out, aVec(4));
    // REACTION 39 -> COMPETITION C0+ + C1+ -> C0+
    in  <<= 0,1,0,1,0,0,0,0;
    out <<= 0,1,0,0,0,0,0,0;
    reactions[39] = Reaction(in, out, aVec(4));
    //// INTER COMPETITION C2->C1 ////////////////////////////////////////////////////////////////////
    // REACTION 40 -> COMPETITION C2 + C1 -> C2
    in  <<= 0,0,1,0,1,0,0,0;
    out <<= 0,0,0,0,1,0,0,0;
    reactions[40] = Reaction(in, out, aVec(6));
    // REACTION 41 -> COMPETITION C2 + C1+ -> C2
    in  <<= 0,0,0,1,1,0,0,0;
    out <<= 0,0,0,0,1,0,0,0;
    reactions[41] = Reaction(in, out, aVec(6));
    // REACTION 42 -> COMPETITION C2+ + C1 -> C2+
    in  <<= 0,0,1,0,0,1,0,0;
    out <<= 0,0,0,0,0,1,0,0;
    reactions[42] = Reaction(in, out, aVec(6));
    // REACTION 43 -> COMPETITION C2+ + C1+ -> C2+
    in  <<= 0,0,0,1,0,1,0,0;
    out <<= 0,0,0,0,0,1,0,0;
    reactions[43] = Reaction(in, out, aVec(6));
    //// INTER COMPETITION C3->C1 ////////////////////////////////////////////////////////////////////
    // REACTION 44 -> COMPETITION C3 + C1 -> C3
    in  <<= 0,0,1,0,0,0,1,0;
    out <<= 0,0,0,0,0,0,1,0;
    reactions[44] = Reaction(in, out, aVec(7));
    // REACTION 45 -> COMPETITION C3 + C1+ -> C3
    in  <<= 0,0,0,1,0,0,1,0;
    out <<= 0,0,0,0,0,0,1,0;
    reactions[45] = Reaction(in, out, aVec(7));
    // REACTION 46 -> COMPETITION C3+ + C1 -> C3+
    in  <<= 0,0,1,0,0,0,0,1;
    out <<= 0,0,0,0,0,0,0,1;
    reactions[46] = Reaction(in, out, aVec(7));
    // REACTION 47 -> COMPETITION C3+ + C1+ -> C3+
    in  <<= 0,0,0,1,0,0,0,1;
    out <<= 0,0,0,0,0,0,0,1;
    reactions[47] = Reaction(in, out, aVec(7));
    //// INTER COMPETITION C0->C2 ////////////////////////////////////////////////////////////////////
    // REACTION 48 -> COMPETITION C0 + C2 -> C0
    in  <<= 1,0,0,0,1,0,0,0;
    out <<= 1,0,0,0,0,0,0,0;
    reactions[48] = Reaction(in, out, aVec(8));
    // REACTION 49 -> COMPETITION C0 + C2+ -> C0
    in  <<= 1,0,0,0,0,1,0,0;
    out <<= 1,0,0,0,0,0,0,0;
    reactions[49] = Reaction(in, out, aVec(8));
    // REACTION 50 -> COMPETITION C0+ + C2 -> C0+
    in  <<= 0,1,0,0,1,0,0,0;
    out <<= 0,1,0,0,0,0,0,0;
    reactions[50] = Reaction(in, out, aVec(8));
    // REACTION 51 -> COMPETITION C0+ + C2+ -> C0+
    in  <<= 0,1,0,0,0,1,0,0;
    out <<= 0,1,0,0,0,0,0,0;
    reactions[51] = Reaction(in, out, aVec(8));
    //// INTER COMPETITION C1->C2 ////////////////////////////////////////////////////////////////////
    // REACTION 52 -> COMPETITION C2 + C1 -> C1
    in  <<= 0,0,1,0,1,0,0,0;
    out <<= 0,0,1,0,0,0,0,0;
    reactions[52] = Reaction(in, out, aVec(9));
    // REACTION 53 -> COMPETITION C2 + C1+ -> C1+
    in  <<= 0,0,0,1,1,0,0,0;
    out <<= 0,0,0,1,0,0,0,0;
    reactions[53] = Reaction(in, out, aVec(9));
    // REACTION 54 -> COMPETITION C2+ + C1 -> C1
    in  <<= 0,0,1,0,0,1,0,0;
    out <<= 0,0,1,0,0,0,0,0;
    reactions[54] = Reaction(in, out, aVec(9));
    // REACTION 55 -> COMPETITION C2+ + C1+ -> C1+
    in  <<= 0,0,0,1,0,1,0,0;
    out <<= 0,0,0,1,0,0,0,0;
    reactions[55] = Reaction(in, out, aVec(9));
    //// INTER COMPETITION C3->C2 ////////////////////////////////////////////////////////////////////
    // REACTION 56 -> COMPETITION C3 + C2 -> C3
    in  <<= 0,0,0,0,1,0,1,0;
    out <<= 0,0,0,0,0,0,1,0;
    reactions[56] = Reaction(in, out, aVec(11));
    // REACTION 57 -> COMPETITION C3 + C2+ -> C3
    in  <<= 0,0,0,0,0,1,1,0;
    out <<= 0,0,0,0,0,0,1,0;
    reactions[57] = Reaction(in, out, aVec(11));
    // REACTION 58 -> COMPETITION C3+ + C2 -> C3+
    in  <<= 0,0,0,0,1,0,0,1;
    out <<= 0,0,0,0,0,0,0,1;
    reactions[58] = Reaction(in, out, aVec(11));
    // REACTION 59 -> COMPETITION C3+ + C2+ -> C3+
    in  <<= 0,0,0,0,0,1,0,1;
    out <<= 0,0,0,0,0,0,0,1;
    reactions[59] = Reaction(in, out, aVec(11));
    //// INTER COMPETITION C0->C3 ////////////////////////////////////////////////////////////////////
    // REACTION 60 -> COMPETITION C0 + C3 -> C0
    in  <<= 1,0,0,0,0,0,1,0;
    out <<= 1,0,0,0,0,0,0,0;
    reactions[60] = Reaction(in, out, aVec(12));
    // REACTION 61 -> COMPETITION C0 + C3+ -> C0
    in  <<= 1,0,0,0,0,0,0,1;
    out <<= 1,0,0,0,0,0,0,0;
    reactions[61] = Reaction(in, out, aVec(12));
    // REACTION 62 -> COMPETITION C0+ + C3 -> C0+
    in  <<= 0,1,0,0,0,0,1,0;
    out <<= 0,1,0,0,0,0,0,0;
    reactions[62] = Reaction(in, out, aVec(12));
    // REACTION 63 -> COMPETITION C0+ + C3+ -> C0+
    in  <<= 0,1,0,0,0,0,0,1;
    out <<= 0,1,0,0,0,0,0,0;
    reactions[63] = Reaction(in, out, aVec(12));
    //// INTER COMPETITION C1->C3 (ONLY POSITIVE VALUE) ////////////////////////////////////////////////
    // REACTION 64 -> COMPETITION C3 + C1 -> 2C3 + C1
    in  <<= 0,0,1,0,0,0,1,0;
    out <<= 0,0,1,0,0,0,2,0;
    reactions[64] = Reaction(in, out, aVec(13));
    // REACTION 65 -> COMPETITION C3 + C1+ -> 2C3 + C1+
    in  <<= 0,0,0,1,0,0,1,0;
    out <<= 0,0,0,1,0,0,2,0;
    reactions[65] = Reaction(in, out, aVec(13));
    // REACTION 66 -> COMPETITION C3+ + C1 -> 2C3+ C1
    in  <<= 0,0,1,0,0,0,0,1;
    out <<= 0,0,1,0,0,0,0,2;
    reactions[66] = Reaction(in, out, aVec(13));
    // REACTION 67 -> COMPETITION C3+ + C1+ -> 2C3+ C1
    in  <<= 0,0,0,1,0,0,0,1;
    out <<= 0,0,0,1,0,0,0,2;
    reactions[67] = Reaction(in, out, aVec(13));
    //// INTER COMPETITION C2->C3 ////////////////////////////////////////////////////////////////////
    // REACTION 68 -> COMPETITION C2 + C3 -> C2
    in  <<= 0,0,0,0,1,0,1,0;
    out <<= 0,0,0,0,1,0,0,0;
    reactions[68] = Reaction(in, out, aVec(14));
    // REACTION 69 -> COMPETITION C2 + C3+ -> C2
    in  <<= 0,0,0,0,1,0,0,1;
    out <<= 0,0,0,0,1,0,0,0;
    reactions[69] = Reaction(in, out, aVec(14));
    // REACTION 70 -> COMPETITION C2+ + C3 -> C2+
    in  <<= 0,0,0,0,0,1,1,0;
    out <<= 0,0,0,0,0,1,0,0;
    reactions[70] = Reaction(in, out, aVec(14));
    // REACTION 71 -> COMPETITION C2+ + C3+ -> C2+
    in  <<= 0,0,0,0,0,1,0,1;
    out <<= 0,0,0,0,0,1,0,0;
    reactions[71] = Reaction(in, out, aVec(14));
    //// CONJUGATION ////////////////////////////////////////////////////////////////////
    // REACTION 72 -> CONJUGATION C0 + C0+ -> 2C0+
    in  <<= 1,1,0,0,0,0,0,0;
    out <<= 0,2,0,0,0,0,0,0;
    reactions[72] = Reaction(in, out, hVec(0));
    // REACTION 73 -> CONJUGATION C0 + C1+ -> C0+ +C1+
    in  <<= 1,0,0,1,0,0,0,0;
    out <<= 0,1,0,1,0,0,0,0;
    reactions[73] = Reaction(in, out, hVec(1));
    // REACTION 74 -> CONJUGATION C1 + C0+ -> C1+ + C0+
    in  <<= 0,1,1,0,0,0,0,0;
    out <<= 0,1,0,1,0,0,0,0;
    reactions[74] = Reaction(in, out, hVec(4));
    // REACTION 75 -> CONJUGATION C1 + C1+ -> 2C1+
    in  <<= 0,0,1,1,0,0,0,0;
    out <<= 0,0,0,2,0,0,0,0;
    reactions[75] = Reaction(in, out, hVec(5));
    //// SEGREGATION LOSS (ALSO BIRTH PROCESS) ///////////////////////////////////
    // REACTION 76 -> C0+ -> C0+ + C0
    in  <<= 0,1,0,0,0,0,0,0;
    out <<= 1,1,0,0,0,0,0,0;
    reactions[76] = Reaction(in, out, rVec(1)*(1-c(1))*gama(1));
    // REACTION 77 -> C1+ -> C1+ + C1
    in  <<= 0,0,0,1,0,0,0,0;
    out <<= 0,0,1,1,0,0,0,0;
    reactions[77] = Reaction(in, out, rVec(3)*(1-c(3))*gama(3));
    // REACTION 78 -> C2+ -> C2+ + C2
    in  <<= 0,0,0,0,0,1,0,0;
    out <<= 0,0,0,0,1,1,0,0;
    reactions[78] = Reaction(in, out, rVec(5)*(1-c(5))*gama(5));
    // REACTION 79 -> C3+ -> C3+ + C3
    in  <<= 0,0,0,0,0,0,0,1;
    out <<= 0,0,0,0,0,0,1,1;
    reactions[79] = Reaction(in, out, rVec(7)*(1-c(7))*gama(7));
    
    //// DEATH DUE TO DRUG PRESSURE ///////////////////////////////////
    // REACTION 80 -> C0 ->
    in  <<= 1,0,0,0,0,0,0,0;
    out <<= 0,0,0,0,0,0,0,0;
    reactions[80] = Reaction(in, out, drugPressure*drugKillRates(0));
    // REACTION 81 -> C0+ ->
    in  <<= 0,1,0,0,0,0,0,0;
    out <<= 0,0,0,0,0,0,0,0;
    reactions[81] = Reaction(in, out, drugPressure*drugKillRates(1));
    // REACTION 82 -> C1 ->
    in  <<= 0,0,1,0,0,0,0,0;
    out <<= 0,0,0,0,0,0,0,0;
    reactions[82] = Reaction(in, out, drugPressure*drugKillRates(2));
    // REACTION 83 -> C1+ ->
    in  <<= 0,0,0,1,0,0,0,0;
    out <<= 0,0,0,0,0,0,0,0;
    reactions[83] = Reaction(in, out, drugPressure*drugKillRates(3));
    // REACTION 84 -> C2 ->
    in  <<= 0,0,0,0,1,0,0,0;
    out <<= 0,0,0,0,0,0,0,0;
    reactions[84] = Reaction(in, out, drugPressure*drugKillRates(4));
    // REACTION 85 -> C2+ ->
    in  <<= 0,0,0,0,0,1,0,0;
    out <<= 0,0,0,0,0,0,0,0;
    reactions[85] = Reaction(in, out, drugPressure*drugKillRates(5));
    // REACTION 86 -> C3 ->
    in  <<= 0,0,0,0,0,0,1,0;
    out <<= 0,0,0,0,0,0,0,0;
    reactions[86] = Reaction(in, out, drugPressure*drugKillRates(6));
    // REACTION 87 -> C3+ ->
    in  <<= 0,0,0,0,0,0,0,1;
    out <<= 0,0,0,0,0,0,0,0;
    reactions[87] = Reaction(in, out, drugPressure*drugKillRates(7));
    //////////////////////////////////////////////////////////////////////////////////////////////
    x = reactions[reactIdx].implement(x);
    return x;
}

template <typename T>
void calculatePropensities(ublas::matrix<T> &tspanTrt, ublas::matrix<T> &tspanInf){
    
    
    // in(0) -> C0
    // in(1) -> C0+
    // in(2) -> C1
    // in(3) -> C1+
    // in(4) -> C2
    // in(5) -> C2+ (=0)
    // in(6) -> C3
    // in(7) -> C3+ (=0)
    
    ////////////////////////////////////////////// REACTIONS ////////////////////////////////////////////////
    Reaction reactions[88];
    popStates in(numOfSpecies), out(numOfSpecies);
    /////// BIRTH PROCESSES ////////
    // REACTION 0 -> BIRTH OF C0
    in  <<= 1,0,0,0,0,0,0,0;
    out <<= 2,0,0,0,0,0,0,0;
    reactions[0] = Reaction(in, out, rVec(0));
    // REACTION 1 -> BIRTH OF C0+
    in  <<= 0,1,0,0,0,0,0,0;
    out <<= 0,2,0,0,0,0,0,0;
    reactions[1] = Reaction(in, out, rVec(1)*(1-c(1))*(1-gama(1)));
    // REACTION 2 -> BIRTH OF C1
    in  <<= 0,0,1,0,0,0,0,0;
    out <<= 0,0,2,0,0,0,0,0;
    reactions[2] = Reaction(in, out, rVec(2));
    // REACTION 3 -> BIRTH OF C1+
    in  <<= 0,0,0,1,0,0,0,0;
    out <<= 0,0,0,2,0,0,0,0;
    reactions[3] = Reaction(in, out, rVec(3)*(1-c(3))*(1-gama(3)));
    // REACTION 4 -> BIRTH OF C2
    in  <<= 0,0,0,0,1,0,0,0;
    out <<= 0,0,0,0,2,0,0,0;
    reactions[4] = Reaction(in, out, rVec(4));
    // REACTION 5 -> BIRTH OF C2+
    in  <<= 0,0,0,0,0,1,0,0;
    out <<= 0,0,0,0,0,2,0,0;
    reactions[5] = Reaction(in, out, rVec(5)*(1-c(5))*(1-gama(5)));
    // REACTION 6 -> BIRTH OF C3
    in  <<= 0,0,0,0,0,0,1,0;
    out <<= 0,0,0,0,0,0,2,0;
    reactions[6] = Reaction(in, out, rVec(6));
    // REACTION 7 -> BIRTH OF C3+
    in  <<= 0,0,0,0,0,0,0,1;
    out <<= 0,0,0,0,0,0,0,2;
    reactions[7] = Reaction(in, out, rVec(7)*(1-c(7))*(1-gama(7)));
    /////// INTERACTION PROCESSES //////
    // TAKE THE ABS VALUE OF AVEC //
    for (int i=0;i<numOfPhyla*numOfPhyla;i++){
        if(aVec(i)<0){
            aVec(i)=-1*aVec(i);
        }
    }
    //// INTRA COMPETITION AMONG C0 ////////////////////////////////////////////////////////////////////
    // REACTION 8 -> COMPETITION C0 + C0 -> C0
    in  <<= 2,0,0,0,0,0,0,0;
    out <<= 1,0,0,0,0,0,0,0;
    reactions[8] = Reaction(in, out, aVec(0));
    // REACTION 9 -> COMPETITION C0 + C0+ -> C0
    in  <<= 1,1,0,0,0,0,0,0;
    out <<= 1,0,0,0,0,0,0,0;
    reactions[9] = Reaction(in, out, aVec(0));
    // REACTION 10 -> COMPETITION C0 + C0+ -> C0+
    in  <<= 1,1,0,0,0,0,0,0;
    out <<= 0,1,0,0,0,0,0,0;
    reactions[10] = Reaction(in, out, aVec(0));
    // REACTION 11 -> COMPETITION C0+ + C0+ -> C0+
    in  <<= 0,2,0,0,0,0,0,0;
    out <<= 0,1,0,0,0,0,0,0;
    reactions[11] = Reaction(in, out, aVec(0));
    //// INTRA COMPETITION AMONG C1 ////////////////////////////////////////////////////////////////////
    // REACTION 12 -> COMPETITION C1 + C1 -> C1
    in  <<= 0,0,2,0,0,0,0,0;
    out <<= 0,0,1,0,0,0,0,0;
    reactions[12] = Reaction(in, out, aVec(5));
    // REACTION 13 -> COMPETITION C1 + C1+ -> C1
    in  <<= 0,0,1,1,0,0,0,0;
    out <<= 0,0,1,0,0,0,0,0;
    reactions[13] = Reaction(in, out, aVec(5));
    // REACTION 14 -> COMPETITION C1 + C1+ -> C1+
    in  <<= 0,0,1,1,0,0,0,0;
    out <<= 0,0,0,1,0,0,0,0;
    reactions[14] = Reaction(in, out, aVec(5));
    // REACTION 15 -> COMPETITION C1+ + C1+ -> C1+
    in  <<= 0,0,0,2,0,0,0,0;
    out <<= 0,0,0,1,0,0,0,0;
    reactions[15] = Reaction(in, out, aVec(5));
    //// INTRA COMPETITION AMONG C2 ////////////////////////////////////////////////////////////////////
    // REACTION 16 -> COMPETITION C2 + C2 -> C2
    in  <<= 0,0,0,0,2,0,0,0;
    out <<= 0,0,0,0,1,0,0,0;
    reactions[16] = Reaction(in, out, aVec(10));
    // REACTION 17 -> COMPETITION C2 + C2+ -> C2
    in  <<= 0,0,0,0,1,1,0,0;
    out <<= 0,0,0,0,1,0,0,0;
    reactions[17] = Reaction(in, out, aVec(10));
    // REACTION 18 -> COMPETITION C2 + C2+ -> C2+
    in  <<= 0,0,0,0,1,1,0,0;
    out <<= 0,0,0,0,0,1,0,0;
    reactions[18] = Reaction(in, out, aVec(10));
    // REACTION 19 -> COMPETITION C2+ + C2+ -> C2+
    in  <<= 0,0,0,0,0,2,0,0;
    out <<= 0,0,0,0,0,1,0,0;
    reactions[19] = Reaction(in, out, aVec(10));
    //// INTRA COMPETITION AMONG C3 ////////////////////////////////////////////////////////////////////
    // REACTION 20 -> COMPETITION C3 + C3 -> C3
    in  <<= 0,0,0,0,0,0,2,0;
    out <<= 0,0,0,0,0,0,1,0;
    reactions[20] = Reaction(in, out, aVec(15));
    // REACTION 21 -> COMPETITION C3 + C3+ -> C3
    in  <<= 0,0,0,0,0,0,1,1;
    out <<= 0,0,0,0,0,0,1,0;
    reactions[21] = Reaction(in, out, aVec(15));
    // REACTION 22 -> COMPETITION C3 + C3+ -> C3+
    in  <<= 0,0,0,0,0,0,1,1;
    out <<= 0,0,0,0,0,0,0,1;
    reactions[22] = Reaction(in, out, aVec(15));
    // REACTION 23 -> COMPETITION C3+ + C3+ -> C3+
    in  <<= 0,0,0,0,0,0,0,2;
    out <<= 0,0,0,0,0,0,0,1;
    reactions[23] = Reaction(in, out, aVec(15));
    //// INTER COMPETITION C1->C0 ////////////////////////////////////////////////////////////////////
    // REACTION 24 -> COMPETITION C0 + C1 -> C1
    in  <<= 1,0,1,0,0,0,0,0;
    out <<= 0,0,1,0,0,0,0,0;
    reactions[24] = Reaction(in, out, aVec(1));
    // REACTION 25 -> COMPETITION C0 + C1+ -> C1+
    in  <<= 1,0,0,1,0,0,0,0;
    out <<= 0,0,0,1,0,0,0,0;
    reactions[25] = Reaction(in, out, aVec(1));
    // REACTION 26 -> COMPETITION C0+ + C1 -> C1
    in  <<= 0,1,1,0,0,0,0,0;
    out <<= 0,0,1,0,0,0,0,0;
    reactions[26] = Reaction(in, out, aVec(1));
    // REACTION 27 -> COMPETITION C0+ + C1+ -> C1+
    in  <<= 0,1,0,1,0,0,0,0;
    out <<= 0,0,0,1,0,0,0,0;
    reactions[27] = Reaction(in, out, aVec(1));
    //// INTER COMPETITION C2->C0 ////////////////////////////////////////////////////////////////////
    // REACTION 28 -> COMPETITION C0 + C2 -> C2
    in  <<= 1,0,0,0,1,0,0,0;
    out <<= 0,0,0,0,1,0,0,0;
    reactions[28] = Reaction(in, out, aVec(2));
    // REACTION 29 -> COMPETITION C0 + C2+ -> C2+
    in  <<= 1,0,0,0,0,1,0,0;
    out <<= 0,0,0,0,0,1,0,0;
    reactions[29] = Reaction(in, out, aVec(2));
    // REACTION 30 -> COMPETITION C0+ + C2 -> C2
    in  <<= 0,1,0,0,1,0,0,0;
    out <<= 0,0,0,0,1,0,0,0;
    reactions[30] = Reaction(in, out, aVec(2));
    // REACTION 31 -> COMPETITION C0+ + C2+ -> C2+
    in  <<= 0,1,0,0,0,1,0,0;
    out <<= 0,0,0,0,0,1,0,0;
    reactions[31] = Reaction(in, out, aVec(2));
    //// INTER COMPETITION C3->C0 ////////////////////////////////////////////////////////////////////
    // REACTION 32 -> COMPETITION C0 + C3 -> C3
    in  <<= 1,0,0,0,0,0,1,0;
    out <<= 0,0,0,0,0,0,1,0;
    reactions[32] = Reaction(in, out, aVec(3));
    // REACTION 33 -> COMPETITION C0 + C3+ -> C3+
    in  <<= 1,0,0,0,0,0,0,1;
    out <<= 0,0,0,0,0,0,0,1;
    reactions[33] = Reaction(in, out, aVec(3));
    // REACTION 34 -> COMPETITION C0+ + C3 -> C3
    in  <<= 0,1,0,0,0,0,1,0;
    out <<= 0,0,0,0,0,0,1,0;
    reactions[34] = Reaction(in, out, aVec(3));
    // REACTION 35 -> COMPETITION C0+ + C3+ -> C3+
    in  <<= 0,1,0,0,0,0,0,1;
    out <<= 0,0,0,0,0,0,0,1;
    reactions[35] = Reaction(in, out, aVec(3));
    //// INTER COMPETITION C0->C1 ////////////////////////////////////////////////////////////////////
    // REACTION 36 -> COMPETITION C0 + C1 -> C0
    in  <<= 1,0,1,0,0,0,0,0;
    out <<= 1,0,0,0,0,0,0,0;
    reactions[36] = Reaction(in, out, aVec(4));
    // REACTION 37 -> COMPETITION C0 + C1+ -> C0
    in  <<= 1,0,0,1,0,0,0,0;
    out <<= 1,0,0,0,0,0,0,0;
    reactions[37] = Reaction(in, out, aVec(4));
    // REACTION 38 -> COMPETITION C0+ + C1 -> C0+
    in  <<= 0,1,1,0,0,0,0,0;
    out <<= 0,1,0,0,0,0,0,0;
    reactions[38] = Reaction(in, out, aVec(4));
    // REACTION 39 -> COMPETITION C0+ + C1+ -> C0+
    in  <<= 0,1,0,1,0,0,0,0;
    out <<= 0,1,0,0,0,0,0,0;
    reactions[39] = Reaction(in, out, aVec(4));
    //// INTER COMPETITION C2->C1 ////////////////////////////////////////////////////////////////////
    // REACTION 40 -> COMPETITION C2 + C1 -> C2
    in  <<= 0,0,1,0,1,0,0,0;
    out <<= 0,0,0,0,1,0,0,0;
    reactions[40] = Reaction(in, out, aVec(6));
    // REACTION 41 -> COMPETITION C2 + C1+ -> C2
    in  <<= 0,0,0,1,1,0,0,0;
    out <<= 0,0,0,0,1,0,0,0;
    reactions[41] = Reaction(in, out, aVec(6));
    // REACTION 42 -> COMPETITION C2+ + C1 -> C2+
    in  <<= 0,0,1,0,0,1,0,0;
    out <<= 0,0,0,0,0,1,0,0;
    reactions[42] = Reaction(in, out, aVec(6));
    // REACTION 43 -> COMPETITION C2+ + C1+ -> C2+
    in  <<= 0,0,0,1,0,1,0,0;
    out <<= 0,0,0,0,0,1,0,0;
    reactions[43] = Reaction(in, out, aVec(6));
    //// INTER COMPETITION C3->C1 ////////////////////////////////////////////////////////////////////
    // REACTION 44 -> COMPETITION C3 + C1 -> C3
    in  <<= 0,0,1,0,0,0,1,0;
    out <<= 0,0,0,0,0,0,1,0;
    reactions[44] = Reaction(in, out, aVec(7));
    // REACTION 45 -> COMPETITION C3 + C1+ -> C3
    in  <<= 0,0,0,1,0,0,1,0;
    out <<= 0,0,0,0,0,0,1,0;
    reactions[45] = Reaction(in, out, aVec(7));
    // REACTION 46 -> COMPETITION C3+ + C1 -> C3+
    in  <<= 0,0,1,0,0,0,0,1;
    out <<= 0,0,0,0,0,0,0,1;
    reactions[46] = Reaction(in, out, aVec(7));
    // REACTION 47 -> COMPETITION C3+ + C1+ -> C3+
    in  <<= 0,0,0,1,0,0,0,1;
    out <<= 0,0,0,0,0,0,0,1;
    reactions[47] = Reaction(in, out, aVec(7));
    //// INTER COMPETITION C0->C2 ////////////////////////////////////////////////////////////////////
    // REACTION 48 -> COMPETITION C0 + C2 -> C0
    in  <<= 1,0,0,0,1,0,0,0;
    out <<= 1,0,0,0,0,0,0,0;
    reactions[48] = Reaction(in, out, aVec(8));
    // REACTION 49 -> COMPETITION C0 + C2+ -> C0
    in  <<= 1,0,0,0,0,1,0,0;
    out <<= 1,0,0,0,0,0,0,0;
    reactions[49] = Reaction(in, out, aVec(8));
    // REACTION 50 -> COMPETITION C0+ + C2 -> C0+
    in  <<= 0,1,0,0,1,0,0,0;
    out <<= 0,1,0,0,0,0,0,0;
    reactions[50] = Reaction(in, out, aVec(8));
    // REACTION 51 -> COMPETITION C0+ + C2+ -> C0+
    in  <<= 0,1,0,0,0,1,0,0;
    out <<= 0,1,0,0,0,0,0,0;
    reactions[51] = Reaction(in, out, aVec(8));
    //// INTER COMPETITION C1->C2 ////////////////////////////////////////////////////////////////////
    // REACTION 52 -> COMPETITION C2 + C1 -> C1
    in  <<= 0,0,1,0,1,0,0,0;
    out <<= 0,0,1,0,0,0,0,0;
    reactions[52] = Reaction(in, out, aVec(9));
    // REACTION 53 -> COMPETITION C2 + C1+ -> C1+
    in  <<= 0,0,0,1,1,0,0,0;
    out <<= 0,0,0,1,0,0,0,0;
    reactions[53] = Reaction(in, out, aVec(9));
    // REACTION 54 -> COMPETITION C2+ + C1 -> C1
    in  <<= 0,0,1,0,0,1,0,0;
    out <<= 0,0,1,0,0,0,0,0;
    reactions[54] = Reaction(in, out, aVec(9));
    // REACTION 55 -> COMPETITION C2+ + C1+ -> C1+
    in  <<= 0,0,0,1,0,1,0,0;
    out <<= 0,0,0,1,0,0,0,0;
    reactions[55] = Reaction(in, out, aVec(9));
    //// INTER COMPETITION C3->C2 ////////////////////////////////////////////////////////////////////
    // REACTION 56 -> COMPETITION C3 + C2 -> C3
    in  <<= 0,0,0,0,1,0,1,0;
    out <<= 0,0,0,0,0,0,1,0;
    reactions[56] = Reaction(in, out, aVec(11));
    // REACTION 57 -> COMPETITION C3 + C2+ -> C3
    in  <<= 0,0,0,0,0,1,1,0;
    out <<= 0,0,0,0,0,0,1,0;
    reactions[57] = Reaction(in, out, aVec(11));
    // REACTION 58 -> COMPETITION C3+ + C2 -> C3+
    in  <<= 0,0,0,0,1,0,0,1;
    out <<= 0,0,0,0,0,0,0,1;
    reactions[58] = Reaction(in, out, aVec(11));
    // REACTION 59 -> COMPETITION C3+ + C2+ -> C3+
    in  <<= 0,0,0,0,0,1,0,1;
    out <<= 0,0,0,0,0,0,0,1;
    reactions[59] = Reaction(in, out, aVec(11));
    //// INTER COMPETITION C0->C3 ////////////////////////////////////////////////////////////////////
    // REACTION 60 -> COMPETITION C0 + C3 -> C0
    in  <<= 1,0,0,0,0,0,1,0;
    out <<= 1,0,0,0,0,0,0,0;
    reactions[60] = Reaction(in, out, aVec(12));
    // REACTION 61 -> COMPETITION C0 + C3+ -> C0
    in  <<= 1,0,0,0,0,0,0,1;
    out <<= 1,0,0,0,0,0,0,0;
    reactions[61] = Reaction(in, out, aVec(12));
    // REACTION 62 -> COMPETITION C0+ + C3 -> C0+
    in  <<= 0,1,0,0,0,0,1,0;
    out <<= 0,1,0,0,0,0,0,0;
    reactions[62] = Reaction(in, out, aVec(12));
    // REACTION 63 -> COMPETITION C0+ + C3+ -> C0+
    in  <<= 0,1,0,0,0,0,0,1;
    out <<= 0,1,0,0,0,0,0,0;
    reactions[63] = Reaction(in, out, aVec(12));
    //// INTER COMPETITION C1->C3 (ONLY POSITIVE VALUE) ////////////////////////////////////////////////
    // REACTION 64 -> COMPETITION C3 + C1 -> 2C3 + C1
    in  <<= 0,0,1,0,0,0,1,0;
    out <<= 0,0,1,0,0,0,2,0;
    reactions[64] = Reaction(in, out, aVec(13));
    // REACTION 65 -> COMPETITION C3 + C1+ -> 2C3 + C1+
    in  <<= 0,0,0,1,0,0,1,0;
    out <<= 0,0,0,1,0,0,2,0;
    reactions[65] = Reaction(in, out, aVec(13));
    // REACTION 66 -> COMPETITION C3+ + C1 -> 2C3+ C1
    in  <<= 0,0,1,0,0,0,0,1;
    out <<= 0,0,1,0,0,0,0,2;
    reactions[66] = Reaction(in, out, aVec(13));
    // REACTION 67 -> COMPETITION C3+ + C1+ -> 2C3+ C1
    in  <<= 0,0,0,1,0,0,0,1;
    out <<= 0,0,0,1,0,0,0,2;
    reactions[67] = Reaction(in, out, aVec(13));
    //// INTER COMPETITION C2->C3 ////////////////////////////////////////////////////////////////////
    // REACTION 68 -> COMPETITION C2 + C3 -> C2
    in  <<= 0,0,0,0,1,0,1,0;
    out <<= 0,0,0,0,1,0,0,0;
    reactions[68] = Reaction(in, out, aVec(14));
    // REACTION 69 -> COMPETITION C2 + C3+ -> C2
    in  <<= 0,0,0,0,1,0,0,1;
    out <<= 0,0,0,0,1,0,0,0;
    reactions[69] = Reaction(in, out, aVec(14));
    // REACTION 70 -> COMPETITION C2+ + C3 -> C2+
    in  <<= 0,0,0,0,0,1,1,0;
    out <<= 0,0,0,0,0,1,0,0;
    reactions[70] = Reaction(in, out, aVec(14));
    // REACTION 71 -> COMPETITION C2+ + C3+ -> C2+
    in  <<= 0,0,0,0,0,1,0,1;
    out <<= 0,0,0,0,0,1,0,0;
    reactions[71] = Reaction(in, out, aVec(14));
    //// CONJUGATION ////////////////////////////////////////////////////////////////////
    // REACTION 72 -> CONJUGATION C0 + C0+ -> 2C0+
    in  <<= 1,1,0,0,0,0,0,0;
    out <<= 0,2,0,0,0,0,0,0;
    reactions[72] = Reaction(in, out, hVec(0));
    // REACTION 73 -> CONJUGATION C0 + C1+ -> C0+ +C1+
    in  <<= 1,0,0,1,0,0,0,0;
    out <<= 0,1,0,1,0,0,0,0;
    reactions[73] = Reaction(in, out, hVec(1));
    // REACTION 74 -> CONJUGATION C1 + C0+ -> C1+ + C0+
    in  <<= 0,1,1,0,0,0,0,0;
    out <<= 0,1,0,1,0,0,0,0;
    reactions[74] = Reaction(in, out, hVec(4));
    // REACTION 75 -> CONJUGATION C1 + C1+ -> 2C1+
    in  <<= 0,0,1,1,0,0,0,0;
    out <<= 0,0,0,2,0,0,0,0;
    reactions[75] = Reaction(in, out, hVec(5));
    //// SEGREGATION LOSS (ALSO BIRTH PROCESS) ///////////////////////////////////
    // REACTION 76 -> C0+ -> C0+ + C0
    in  <<= 0,1,0,0,0,0,0,0;
    out <<= 1,1,0,0,0,0,0,0;
    reactions[76] = Reaction(in, out, rVec(1)*(1-c(1))*gama(1));
    // REACTION 77 -> C1+ -> C1+ + C1
    in  <<= 0,0,0,1,0,0,0,0;
    out <<= 0,0,1,1,0,0,0,0;
    reactions[77] = Reaction(in, out, rVec(3)*(1-c(3))*gama(3));
    // REACTION 78 -> C2+ -> C2+ + C2
    in  <<= 0,0,0,0,0,1,0,0;
    out <<= 0,0,0,0,1,1,0,0;
    reactions[78] = Reaction(in, out, rVec(5)*(1-c(5))*gama(5));
    // REACTION 79 -> C3+ -> C3+ + C3
    in  <<= 0,0,0,0,0,0,0,1;
    out <<= 0,0,0,0,0,0,1,1;
    reactions[79] = Reaction(in, out, rVec(7)*(1-c(7))*gama(7));
    
    //// DEATH DUE TO DRUG PRESSURE ///////////////////////////////////
    // REACTION 80 -> C0 ->
    in  <<= 1,0,0,0,0,0,0,0;
    out <<= 0,0,0,0,0,0,0,0;
    reactions[80] = Reaction(in, out, drugPressure*drugKillRates(0));
    // REACTION 81 -> C0+ ->
    in  <<= 0,1,0,0,0,0,0,0;
    out <<= 0,0,0,0,0,0,0,0;
    reactions[81] = Reaction(in, out, drugPressure*drugKillRates(1));
    // REACTION 82 -> C1 ->
    in  <<= 0,0,1,0,0,0,0,0;
    out <<= 0,0,0,0,0,0,0,0;
    reactions[82] = Reaction(in, out, drugPressure*drugKillRates(2));
    // REACTION 83 -> C1+ ->
    in  <<= 0,0,0,1,0,0,0,0;
    out <<= 0,0,0,0,0,0,0,0;
    reactions[83] = Reaction(in, out, drugPressure*drugKillRates(3));
    // REACTION 84 -> C2 ->
    in  <<= 0,0,0,0,1,0,0,0;
    out <<= 0,0,0,0,0,0,0,0;
    reactions[84] = Reaction(in, out, drugPressure*drugKillRates(4));
    // REACTION 85 -> C2+ ->
    in  <<= 0,0,0,0,0,1,0,0;
    out <<= 0,0,0,0,0,0,0,0;
    reactions[85] = Reaction(in, out, drugPressure*drugKillRates(5));
    // REACTION 86 -> C3 ->
    in  <<= 0,0,0,0,0,0,1,0;
    out <<= 0,0,0,0,0,0,0,0;
    reactions[86] = Reaction(in, out, drugPressure*drugKillRates(6));
    // REACTION 87 -> C3+ ->
    in  <<= 0,0,0,0,0,0,0,1;
    out <<= 0,0,0,0,0,0,0,0;
    reactions[87] = Reaction(in, out, drugPressure*drugKillRates(7));
    //////////////////////////////////////////////////////////////////////////////////////////////
    int Nreactions = 88;
    for (int i=0; i<Nreactions; i++){
        for(int k=0;k<numOfSpecies;k++){
            epsMat(i,k)=0;
        }
    }
    for (int r=0; r<Nreactions; r++) {
        propensities(r) = reactions[r].propensity(x); //calculates all the propensities
        for (int j=0; j<numOfSpecies;j++){
            epsMat(r,j) = reactions[r].delta(j);
        }
    }
}

void dX(popStates &x , popStates &dxdt , double t ){
    double discSum, reactIdx;
    for (unsigned i=0;i<x.size();++i){
        if(x(i)<1){
            x(i) = 0;
        }
    }
    for (unsigned i=0;i<x.size();++i){
        discSum = 0;
        for(int m=0;m<numReactions_c;m++){
            reactIdx = Mc(m);
            discSum = discSum + propensities(reactIdx)*epsMat(reactIdx,i); // SUM OVER CTS REACTIONS
        }
        dxdt(i) =  discSum;
    }
}
void dg0(popStates &g0, popStates &dg0dt , double t ){
    double discSum, reactIdx;
    discSum = 0;
    for(int m=0;m<numReactions_d;m++){
        reactIdx = Md(m);
        discSum = discSum + propensities(reactIdx); // SUM OVER CTS REACTIONS
    }
    dg0dt(0) =  discSum;
}
void initiateMcMd(){
    int flg_d = 0;
    int flg_c = 0;
    int flg_0 = 0;
    int Nreactions = 88;
    for(int i=0; i<Nreactions; i++){Md(i)=0;}
    for(int i=0; i<Nreactions; i++){Mc(i)=0;}
    for (int m=0;m<Nreactions;m++){
        flg_0 = 0;
        for(int i=0;i <numOfSpecies; i++){
            if(epsMat(m,i)>0){
                //            if(x(i)>0){
                if(x(i)<=Cx*abs(epsMat(m,i))){
                    flg_0 = flg_0 +1;
                }
            }
        }
        if((propensities(m)>Ca) && (flg_0==0)){
            //            cout << "R_"<<m<<" can be approximated"<<endl;
            Mc(flg_c) = m;
            flg_c ++;
        }else{
            //            cout << "R_"<<m<<" must be stochastic"<<endl;
            Md(flg_d) = m;
            flg_d ++;
        }
    }
    
    numReactions_c = flg_c; //number of cts reactions -> indexes kept in Mc
    numReactions_d = flg_d; //number of deterministic reactions -> indexes kept in Md
}

void updateMcMd(){
    int flg_d = 0;
    int flg_c = 0;
    int flg_0 = 0;
    int Nreactions = 88;
    for(int i=0; i<Nreactions; i++){Md(i)=0;}
    for(int i=0; i<Nreactions; i++){Mc(i)=0;}
    for (int m=0;m<Nreactions;m++){
        flg_0 = 0;
        for(int i=0;i <numOfSpecies; i++){
            if(epsMat(m,i)>0){
                if(x(i)<=Cx*abs(epsMat(m,i))){
                    flg_0 = flg_0 +1;
                }
            }
        }
        if((propensities(m)>Ca) && (flg_0==0)){
            //            cout << "R_"<<m<<" can be approximated"<<endl;
            Mc(flg_c) = m;
            flg_c ++;
        }else{
            //            cout << "R_"<<m<<" must be stochastic"<<endl;
            Md(flg_d) = m;
            flg_d ++;
        }
    }
    
    numReactions_c = flg_c; //number of cts reactions -> indexes kept in Mc
    numReactions_d = flg_d; //number of deterministic reactions -> indexes kept in Md
}

template <typename T>
void updateDrugPressure(ublas::matrix<T> &tspanTrt){
    drugPressure = 0;
    for(unsigned period = 0; period<tspanTrt.size1(); ++period){
        ublas::matrix_row<ublas::matrix<T>> tempRange (tspanTrt,period);
        const int timeBegin       = tempRange(0);
        const int timeEnd         = tempRange(1);
        //        cout << period<<" "<<timeBegin<<" "<<timeEnd<<" "<<t<<" "<<endl;
        if((t>=(double)timeBegin) && (t<(double)timeEnd)){
            drugPressure = tempRange(2);
            break;
        }
    }
}

template <typename T>
void updateSamples(ublas::matrix<T> &tspanTrt){
    double lastDayOfExposure = tspanTrt(tspanTrt.size1()-1,0);
    for(int i=0; i<sampleTimes.size(); ++i){
        if(abs(t-(lastDayOfExposure+sampleTimes(i)))<dt){
            samplePops(5*i+0)=t;
            samplePops(5*i+1)=x(0);
            samplePops(5*i+2)=x(1);
            samplePops(5*i+3)=x(2);
            samplePops(5*i+4)=x(3);
        }
    }
}

template <typename T>
void updateInfection(ublas::matrix<T> &tspanInf){
    inf = 0;
    for(unsigned period = 0; period<tspanInf.size1(); ++period){
        ublas::matrix_row<ublas::matrix<T>> tempRange (tspanInf,period);
        const int timeBegin       = tempRange(0);
        const int timeEnd         = tempRange(1);
        if((tempRange(1)-tempRange(0))>0){
            if((t>=(double)timeBegin) && (t<(double)timeEnd)){
                inf = tempRange(2);
                break;
            }
        }else{
            if((t-(double)timeBegin)<dt){
                inf = tempRange(2);
                break;
            }
        }
    }
    if(inf==1){
        for(unsigned i=0; i<numOfSpecies; ++i){
            x(i) = x(i)+inf*infAbundVec(i);
        }
    }
}

void updatec0counter(){
    if(x(1)>0){
        c0counter(0)=0;
    }else{
        c0counter(0)=c0counter(0)+1;
    }
    
    if(x(3)>0){
        c0counter(1)=0;
    }else{
        c0counter(1)=c0counter(1)+1;
    }
}

void updateExtinctBoth(){
    if(updateExt==0){
        if((x(1)<1) && (x(3)<1)){
            extinct(0) = 1;
            extinct(1) = t;
            updateExt  = 1;
        }
    }
}


template <typename T>
void runGutStocHybrid(ublas::matrix<T> &tspanTrt, ublas::matrix<T> &tspanInf, string folderChar, int save){
    stringstream stream;
    string fileName;
    
    stream.str(string());
    stream << "schedule_trtInit";
    fileName = stream.str();
    saveUblasVecAppend(trtInit,folderChar,fileName,"%.1f ",save);
    
    stream.str(string());
    stream << "schedule_trtLen";
    fileName = stream.str();
    saveUblasVecAppend(trtDura,folderChar,fileName,"%.1f ",save);
    
    //    stream.str(string());
    //    stream << "schedule_infInit";
    //    fileName = stream.str();
    //    saveUblasVecAppend(infInit,folderChar,fileName,"%.1f ",save);
    
    //////////////////// SET THE INITIAL CONDITIONS ////////////////////
    odeint::runge_kutta4<popStates> rk4;
    t    = 0.0;
    updateDrugPressure(tspanTrt);updateInfection(tspanInf);updatec0counter();updateSamples(tspanTrt);
    x    = mixVecsAbs(xICSum,resAbs); //IC for x
    ////////////////////////////////////////////////////////////////////
    double r1 = ((double) rand() / (RAND_MAX));
    g0(0)= log(r1); //IC for g0
    calculatePropensities(tspanTrt,tspanInf); // t=0
    updateMcMd();
    while(t<obsTime){
        while((g0(0)<0) && (t+dt<obsTime)){
            rk4.adjust_size(x);
            rk4.do_step(dX, x , t , dt);
            rk4.adjust_size(g0);
            rk4.do_step(dg0, g0 , t , dt);
            t+=dt;updateDrugPressure(tspanTrt);updateInfection(tspanInf);updatec0counter();updateSamples(tspanTrt);
            //            xInt(0) = xInt(0)+x(1)*dt;xInt(1) = xInt(1)+x(3)*dt;updateExtinctBoth();updatec0counter();
            calculatePropensities(tspanTrt,tspanInf);
        }
        // g0 > 0 now
        double r2 = ((double) rand() / (RAND_MAX));
        // pick the discrete reaction idx to fire
        double sumAllDiscPropensities = 0;
        double sumDiscPropensities = 0;
        int it;
        for (int i=0; i<numReactions_d; i++){
            sumAllDiscPropensities = sumAllDiscPropensities + r2*propensities(Md(i));
        }
        for (it=0; it<numReactions_d; it++) {
            sumDiscPropensities = sumDiscPropensities+ propensities(Md(it));
            if (sumAllDiscPropensities<sumDiscPropensities){
                break;
            }
        }
        if(t+dt<obsTime){
            x = implementReaction(tspanTrt,Md(it));
            t+=dt;updateDrugPressure(tspanTrt);updateInfection(tspanInf);updatec0counter();updateSamples(tspanTrt);
            //            xInt(0) = xInt(0)+x(1)*dt;xInt(1) = xInt(1)+x(3)*dt;updateExtinctBoth();updatec0counter();
            r1 = ((double) rand() / (RAND_MAX));
            g0(0) = log(r1);
            updateMcMd();
        }else{
            break;
        }
    }
    
    //    stream.str(string());
    //    stream << "xEnds";
    //    fileName = stream.str();
    //    saveUblasVecAppend(x,folderChar,fileName,"%.5f ",save);
    
    //    stream.str(string());
    //    stream << "xInt";
    //    fileName = stream.str();
    //    saveUblasVecAppend(xInt,folderChar,fileName,"%.3f ",save);
    
    //    stream.str(string());
    //    stream << "ext";
    //    fileName = stream.str();
    //    saveUblasVecAppend(extinct,folderChar,fileName,"%.1f ",save);
    
    stream.str(string());
    stream << "samplePops";
    fileName = stream.str();
    saveUblasVecAppend(samplePops,folderChar,fileName,"%.14f ",save);
    
    stream.str(string());
    stream << "extCounter";
    fileName = stream.str();
    saveUblasVecAppend(c0counter,folderChar,fileName,"%.1f ",save);
}


ublas::matrix<double> scheduleRandomTreatment(int trtTime, int obsTime, int totalTrts){
    
    resizeVec(trtInit,totalTrts);
    resizeVec(trtDura,totalTrts);
    ublas::matrix<int> tspanTrtReturn(2*totalTrts+1,3);
    
    if(totalTrts>0){
        trtLengthVector <<= 3,5,7,10;
        ////////// CHOOSE THE NUMBER OF TREATMENTS RANDOMLY /////////
        int dum  = rand() % trtLengthVector.size();
        int idxLastTrtLen = abs(dum);
        int lastTrtLength = trtLengthVector(idxLastTrtLen);
        
        int guardBand = 1;
        int upperLim  = trtTime-(lastTrtLength+guardBand);
        
        ublas::vector<int> diff(totalTrts-1);
        ublas::vector<int> diffKey(totalTrts-1);
        
        for (unsigned i=0; i<totalTrts; ++i){
            trtInit(i)     = abs(1 + rand() % upperLim);
            dum = rand() % trtLengthVector.size();
            int trtDuraIdx = abs(dum);
            trtDura(i)     = trtLengthVector(trtDuraIdx);
        }
        // CORRECT FOR THE LAST TREATMENT LENGHT (OR DURATION)
        trtDura(totalTrts-1) = lastTrtLength;
        sort(trtInit.begin(),trtInit.end());
        for (unsigned i=0; i<totalTrts-1; ++i){
            diff(i)=trtInit(i+1)-trtInit(i);
            if(diff(i)<(trtDura(i)+guardBand)){
                diffKey(i) = 1;
            }else{
                diffKey(i) = 0;
            }
        }
        double sumDiff = sum(diffKey);
        
        while(sumDiff>0){
            for (unsigned i=0; i<totalTrts; ++i){
                trtInit(i)     = abs(1 + rand() % upperLim);
                dum  = (rand() % trtLengthVector.size());
                int trtDuraIdx = abs(dum);
                trtDura(i)     = trtLengthVector(trtDuraIdx);
            }
            // CORRECT FOR THE LAST TREATMENT LENGHT (OR DURATION)
            trtDura(totalTrts-1) = lastTrtLength;
            sort(trtInit.begin(),trtInit.end());
            for (unsigned i=0; i<totalTrts-1; ++i){
                diff(i)=trtInit(i+1)-trtInit(i);
                if(diff(i)<(trtDura(i)+guardBand)){
                    diffKey(i) = 1;
                }else{
                    diffKey(i) = 0;
                }
            }
            sumDiff = sum(diffKey);
        }
        
        tspanTrtReturn(0,0)=0;
        tspanTrtReturn(0,1)=trtInit(0);
        tspanTrtReturn(0,2)=0;
        
        unsigned k=0;
        for (unsigned i=1; i<2*totalTrts-1; i+=2){
            tspanTrtReturn(i,0) = trtInit(k);
            tspanTrtReturn(i,1) = trtInit(k)+trtDura(k);
            tspanTrtReturn(i,2) = 1;
            tspanTrtReturn(i+1,0) = trtInit(k)+trtDura(k);
            tspanTrtReturn(i+1,1) = trtInit(k+1);
            tspanTrtReturn(i+1,2) = 0;
            k = k + 1;
        }
        tspanTrtReturn(2*totalTrts-1,0) = trtInit(k);
        tspanTrtReturn(2*totalTrts-1,1) = trtInit(k)+trtDura(k);
        tspanTrtReturn(2*totalTrts-1,2) = 1;
        tspanTrtReturn(2*totalTrts,0) = trtInit(k)+trtDura(k);
        //        tspanTrtReturn(2*totalTrts,1) = trtTime;
        tspanTrtReturn(2*totalTrts,1) = obsTime;
        tspanTrtReturn(2*totalTrts,2) = 0;
        
    }else{
        tspanTrtReturn(0,0) = 0;
        tspanTrtReturn(0,1) = obsTime;
        tspanTrtReturn(0,2) = 0;
    }
    return tspanTrtReturn;
    
}

ublas::matrix<double> scheduleRandomInfection(int trtTime, int obsTime, int totalInfs){
    
    resizeVec(infInit,totalInfs);
    ublas::matrix<int> tspanTrtReturn(2*totalInfs+1,3);
    
    if(totalInfs>0){
        int guardBand = 1;
        int upperLim  = obsTime;
        ublas::vector<int> diff(totalInfs-1);
        ublas::vector<int> diffKey(totalInfs-1);
        
        for (unsigned i=0; i<totalInfs; ++i){
            infInit(i) = abs(1 + rand() % upperLim);
        }
        // CORRECT FOR THE LAST TREATMENT LENGHT (OR DURATION)
        sort(infInit.begin(),infInit.end());
        for (unsigned i=0; i<totalInfs-1; ++i){
            diff(i)=infInit(i+1)-infInit(i);
            if(diff(i)<guardBand){
                diffKey(i) = 1;
            }else{
                diffKey(i) = 0;
            }
        }
        double sumDiff = sum(diffKey);
        
        while(sumDiff>0){
            for (unsigned i=0; i<totalInfs; ++i){
                infInit(i)     = abs(1 + rand() % upperLim);
            }
            // CORRECT FOR THE LAST TREATMENT LENGHT (OR DURATION)
            sort(infInit.begin(),infInit.end());
            for (unsigned i=0; i<totalInfs-1; ++i){
                diff(i)=infInit(i+1)-infInit(i);
                if(diff(i)<guardBand){
                    diffKey(i) = 1;
                }else{
                    diffKey(i) = 0;
                }
            }
            sumDiff = sum(diffKey);
        }
        
        tspanTrtReturn(0,0)=0;
        tspanTrtReturn(0,1)=infInit(0);
        tspanTrtReturn(0,2)=0;
        
        unsigned k=0;
        for (unsigned i=1; i<2*totalInfs-1; i+=2){
            tspanTrtReturn(i,0) = infInit(k);
            tspanTrtReturn(i,1) = infInit(k);
            tspanTrtReturn(i,2) = 1;
            tspanTrtReturn(i+1,0) = infInit(k);
            tspanTrtReturn(i+1,1) = infInit(k+1);
            tspanTrtReturn(i+1,2) = 0;
            k = k + 1;
        }
        tspanTrtReturn(2*totalInfs-1,0) = infInit(k);
        tspanTrtReturn(2*totalInfs-1,1) = infInit(k);
        tspanTrtReturn(2*totalInfs-1,2) = 1;
        tspanTrtReturn(2*totalInfs,0) = infInit(k);
        tspanTrtReturn(2*totalInfs,1) = obsTime;
        tspanTrtReturn(2*totalInfs,2) = 0;
        
        
    }else{
        tspanTrtReturn(0,0) = 0;
        tspanTrtReturn(0,1) = obsTime;
        tspanTrtReturn(0,2) = 0;
    }
    return tspanTrtReturn;
    
}

double myPow(double a, double b){
    if(b==0){
        return 0;
    }else{
        return pow(a,b);
    }
}

void updateVectors(ublas::vector<double> &parameterVector){
    double c0   = parameterVector(0);
    double c1   = parameterVector(1);
    double h_intra = myPow(10,parameterVector(2));
    double h_inter = myPow(10,parameterVector(3));
    double gam    = parameterVector(4);
    c       <<= c0,c0,c1,c1,0,0,0,0;
    gama    <<= 0,gam,0,gam,0,0,0,0;
    hVec    <<= h_intra,h_inter,0,0,h_inter,h_intra,0,0,0,0,0,0,0,0,0,0;
}

void calcSS(){
    sol01 = (aVec(1)*aVec(6)*aVec(11)*rVec(6)-aVec(1)*aVec(7)*aVec(10)*rVec(6)-aVec(2)*aVec(5)*aVec(11)*rVec(6)+aVec(2)*aVec(7)*aVec(9)*rVec(6)+aVec(3)*aVec(5)*aVec(10)*rVec(6)-aVec(3)*aVec(6)*aVec(9)*rVec(6)-aVec(1)*aVec(6)*aVec(15)*rVec(4)+aVec(1)*aVec(7)*aVec(14)*rVec(4)+aVec(2)*aVec(5)*aVec(15)*rVec(4)-aVec(2)*aVec(7)*aVec(13)*rVec(4)-aVec(3)*aVec(5)*aVec(14)*rVec(4)+aVec(3)*aVec(6)*aVec(13)*rVec(4)+aVec(1)*aVec(10)*aVec(15)*rVec(2)-aVec(1)*aVec(11)*aVec(14)*rVec(2)-aVec(2)*aVec(9)*aVec(15)*rVec(2)+aVec(2)*aVec(11)*aVec(13)*rVec(2)+aVec(3)*aVec(9)*aVec(14)*rVec(2)-aVec(3)*aVec(10)*aVec(13)*rVec(2)-aVec(5)*aVec(10)*aVec(15)*rVec(0)+aVec(5)*aVec(11)*aVec(14)*rVec(0)+aVec(6)*aVec(9)*aVec(15)*rVec(0)-aVec(6)*aVec(11)*aVec(13)*rVec(0)-aVec(7)*aVec(9)*aVec(14)*rVec(0)+aVec(7)*aVec(10)*aVec(13)*rVec(0))/(aVec(0)*aVec(5)*aVec(10)*aVec(15)-aVec(0)*aVec(5)*aVec(11)*aVec(14)-aVec(0)*aVec(6)*aVec(9)*aVec(15)+aVec(0)*aVec(6)*aVec(11)*aVec(13)+aVec(0)*aVec(7)*aVec(9)*aVec(14)-aVec(0)*aVec(7)*aVec(10)*aVec(13)-aVec(1)*aVec(4)*aVec(10)*aVec(15)+aVec(1)*aVec(4)*aVec(11)*aVec(14)+aVec(1)*aVec(6)*aVec(8)*aVec(15)-aVec(1)*aVec(6)*aVec(11)*aVec(12)-aVec(1)*aVec(7)*aVec(8)*aVec(14)+aVec(1)*aVec(7)*aVec(10)*aVec(12)+aVec(2)*aVec(4)*aVec(9)*aVec(15)-aVec(2)*aVec(4)*aVec(11)*aVec(13)-aVec(2)*aVec(5)*aVec(8)*aVec(15)+aVec(2)*aVec(5)*aVec(11)*aVec(12)+aVec(2)*aVec(7)*aVec(8)*aVec(13)-aVec(2)*aVec(7)*aVec(9)*aVec(12)-aVec(3)*aVec(4)*aVec(9)*aVec(14)+aVec(3)*aVec(4)*aVec(10)*aVec(13)+aVec(3)*aVec(5)*aVec(8)*aVec(14)-aVec(3)*aVec(5)*aVec(10)*aVec(12)-aVec(3)*aVec(6)*aVec(8)*aVec(13)+aVec(3)*aVec(6)*aVec(9)*aVec(12));
    sol02 = -1*(aVec(0)*aVec(6)*aVec(11)*rVec(6)-aVec(0)*aVec(7)*aVec(10)*rVec(6)-aVec(2)*aVec(4)*aVec(11)*rVec(6)+aVec(2)*aVec(7)*aVec(8)*rVec(6)+aVec(3)*aVec(4)*aVec(10)*rVec(6)-aVec(3)*aVec(6)*aVec(8)*rVec(6)-aVec(0)*aVec(6)*aVec(15)*rVec(4)+aVec(0)*aVec(7)*aVec(14)*rVec(4)+aVec(2)*aVec(4)*aVec(15)*rVec(4)-aVec(2)*aVec(7)*aVec(12)*rVec(4)-aVec(3)*aVec(4)*aVec(14)*rVec(4)+aVec(3)*aVec(6)*aVec(12)*rVec(4)+aVec(0)*aVec(10)*aVec(15)*rVec(2)-aVec(0)*aVec(11)*aVec(14)*rVec(2)-aVec(2)*aVec(8)*aVec(15)*rVec(2)+aVec(2)*aVec(11)*aVec(12)*rVec(2)+aVec(3)*aVec(8)*aVec(14)*rVec(2)-aVec(3)*aVec(10)*aVec(12)*rVec(2)-aVec(4)*aVec(10)*aVec(15)*rVec(0)+aVec(4)*aVec(11)*aVec(14)*rVec(0)+aVec(6)*aVec(8)*aVec(15)*rVec(0)-aVec(6)*aVec(11)*aVec(12)*rVec(0)-aVec(7)*aVec(8)*aVec(14)*rVec(0)+aVec(7)*aVec(10)*aVec(12)*rVec(0))/(aVec(0)*aVec(5)*aVec(10)*aVec(15)-aVec(0)*aVec(5)*aVec(11)*aVec(14)-aVec(0)*aVec(6)*aVec(9)*aVec(15)+aVec(0)*aVec(6)*aVec(11)*aVec(13)+aVec(0)*aVec(7)*aVec(9)*aVec(14)-aVec(0)*aVec(7)*aVec(10)*aVec(13)-aVec(1)*aVec(4)*aVec(10)*aVec(15)+aVec(1)*aVec(4)*aVec(11)*aVec(14)+aVec(1)*aVec(6)*aVec(8)*aVec(15)-aVec(1)*aVec(6)*aVec(11)*aVec(12)-aVec(1)*aVec(7)*aVec(8)*aVec(14)+aVec(1)*aVec(7)*aVec(10)*aVec(12)+aVec(2)*aVec(4)*aVec(9)*aVec(15)-aVec(2)*aVec(4)*aVec(11)*aVec(13)-aVec(2)*aVec(5)*aVec(8)*aVec(15)+aVec(2)*aVec(5)*aVec(11)*aVec(12)+aVec(2)*aVec(7)*aVec(8)*aVec(13)-aVec(2)*aVec(7)*aVec(9)*aVec(12)-aVec(3)*aVec(4)*aVec(9)*aVec(14)+aVec(3)*aVec(4)*aVec(10)*aVec(13)+aVec(3)*aVec(5)*aVec(8)*aVec(14)-aVec(3)*aVec(5)*aVec(10)*aVec(12)-aVec(3)*aVec(6)*aVec(8)*aVec(13)+aVec(3)*aVec(6)*aVec(9)*aVec(12));
    sol03 = (aVec(0)*aVec(5)*aVec(11)*rVec(6)-aVec(0)*aVec(7)*aVec(9)*rVec(6)-aVec(1)*aVec(4)*aVec(11)*rVec(6)+aVec(1)*aVec(7)*aVec(8)*rVec(6)+aVec(3)*aVec(4)*aVec(9)*rVec(6)-aVec(3)*aVec(5)*aVec(8)*rVec(6)-aVec(0)*aVec(5)*aVec(15)*rVec(4)+aVec(0)*aVec(7)*aVec(13)*rVec(4)+aVec(1)*aVec(4)*aVec(15)*rVec(4)-aVec(1)*aVec(7)*aVec(12)*rVec(4)-aVec(3)*aVec(4)*aVec(13)*rVec(4)+aVec(3)*aVec(5)*aVec(12)*rVec(4)+aVec(0)*aVec(9)*aVec(15)*rVec(2)-aVec(0)*aVec(11)*aVec(13)*rVec(2)-aVec(1)*aVec(8)*aVec(15)*rVec(2)+aVec(1)*aVec(11)*aVec(12)*rVec(2)+aVec(3)*aVec(8)*aVec(13)*rVec(2)-aVec(3)*aVec(9)*aVec(12)*rVec(2)-aVec(4)*aVec(9)*aVec(15)*rVec(0)+aVec(4)*aVec(11)*aVec(13)*rVec(0)+aVec(5)*aVec(8)*aVec(15)*rVec(0)-aVec(5)*aVec(11)*aVec(12)*rVec(0)-aVec(7)*aVec(8)*aVec(13)*rVec(0)+aVec(7)*aVec(9)*aVec(12)*rVec(0))/(aVec(0)*aVec(5)*aVec(10)*aVec(15)-aVec(0)*aVec(5)*aVec(11)*aVec(14)-aVec(0)*aVec(6)*aVec(9)*aVec(15)+aVec(0)*aVec(6)*aVec(11)*aVec(13)+aVec(0)*aVec(7)*aVec(9)*aVec(14)-aVec(0)*aVec(7)*aVec(10)*aVec(13)-aVec(1)*aVec(4)*aVec(10)*aVec(15)+aVec(1)*aVec(4)*aVec(11)*aVec(14)+aVec(1)*aVec(6)*aVec(8)*aVec(15)-aVec(1)*aVec(6)*aVec(11)*aVec(12)-aVec(1)*aVec(7)*aVec(8)*aVec(14)+aVec(1)*aVec(7)*aVec(10)*aVec(12)+aVec(2)*aVec(4)*aVec(9)*aVec(15)-aVec(2)*aVec(4)*aVec(11)*aVec(13)-aVec(2)*aVec(5)*aVec(8)*aVec(15)+aVec(2)*aVec(5)*aVec(11)*aVec(12)+aVec(2)*aVec(7)*aVec(8)*aVec(13)-aVec(2)*aVec(7)*aVec(9)*aVec(12)-aVec(3)*aVec(4)*aVec(9)*aVec(14)+aVec(3)*aVec(4)*aVec(10)*aVec(13)+aVec(3)*aVec(5)*aVec(8)*aVec(14)-aVec(3)*aVec(5)*aVec(10)*aVec(12)-aVec(3)*aVec(6)*aVec(8)*aVec(13)+aVec(3)*aVec(6)*aVec(9)*aVec(12));
    sol04 = -1*(aVec(0)*aVec(5)*aVec(10)*rVec(6)-aVec(0)*aVec(6)*aVec(9)*rVec(6)-aVec(1)*aVec(4)*aVec(10)*rVec(6)+aVec(1)*aVec(6)*aVec(8)*rVec(6)+aVec(2)*aVec(4)*aVec(9)*rVec(6)-aVec(2)*aVec(5)*aVec(8)*rVec(6)-aVec(0)*aVec(5)*aVec(14)*rVec(4)+aVec(0)*aVec(6)*aVec(13)*rVec(4)+aVec(1)*aVec(4)*aVec(14)*rVec(4)-aVec(1)*aVec(6)*aVec(12)*rVec(4)-aVec(2)*aVec(4)*aVec(13)*rVec(4)+aVec(2)*aVec(5)*aVec(12)*rVec(4)+aVec(0)*aVec(9)*aVec(14)*rVec(2)-aVec(0)*aVec(10)*aVec(13)*rVec(2)-aVec(1)*aVec(8)*aVec(14)*rVec(2)+aVec(1)*aVec(10)*aVec(12)*rVec(2)+aVec(2)*aVec(8)*aVec(13)*rVec(2)-aVec(2)*aVec(9)*aVec(12)*rVec(2)-aVec(4)*aVec(9)*aVec(14)*rVec(0)+aVec(4)*aVec(10)*aVec(13)*rVec(0)+aVec(5)*aVec(8)*aVec(14)*rVec(0)-aVec(5)*aVec(10)*aVec(12)*rVec(0)-aVec(6)*aVec(8)*aVec(13)*rVec(0)+aVec(6)*aVec(9)*aVec(12)*rVec(0))/(aVec(0)*aVec(5)*aVec(10)*aVec(15)-aVec(0)*aVec(5)*aVec(11)*aVec(14)-aVec(0)*aVec(6)*aVec(9)*aVec(15)+aVec(0)*aVec(6)*aVec(11)*aVec(13)+aVec(0)*aVec(7)*aVec(9)*aVec(14)-aVec(0)*aVec(7)*aVec(10)*aVec(13)-aVec(1)*aVec(4)*aVec(10)*aVec(15)+aVec(1)*aVec(4)*aVec(11)*aVec(14)+aVec(1)*aVec(6)*aVec(8)*aVec(15)-aVec(1)*aVec(6)*aVec(11)*aVec(12)-aVec(1)*aVec(7)*aVec(8)*aVec(14)+aVec(1)*aVec(7)*aVec(10)*aVec(12)+aVec(2)*aVec(4)*aVec(9)*aVec(15)-aVec(2)*aVec(4)*aVec(11)*aVec(13)-aVec(2)*aVec(5)*aVec(8)*aVec(15)+aVec(2)*aVec(5)*aVec(11)*aVec(12)+aVec(2)*aVec(7)*aVec(8)*aVec(13)-aVec(2)*aVec(7)*aVec(9)*aVec(12)-aVec(3)*aVec(4)*aVec(9)*aVec(14)+aVec(3)*aVec(4)*aVec(10)*aVec(13)+aVec(3)*aVec(5)*aVec(8)*aVec(14)-aVec(3)*aVec(5)*aVec(10)*aVec(12)-aVec(3)*aVec(6)*aVec(8)*aVec(13)+aVec(3)*aVec(6)*aVec(9)*aVec(12));
}

void initiateSampleTimes(int len){
    resizeVec(sampleTimes,len);
    for(unsigned i=0;i<len;++i){
        sampleTimes(i)=15*i;
    }
}

void initiateSamplePops(int len){
    resizeVec(samplePops,5*len);
    for(unsigned i=0;i<5*len;++i){
        samplePops(i)=0;
    }
}

std::vector<std::vector<double> > readCSV(std::string filename)
{
    std::vector<std::vector<double> > contents;
    std::ifstream in(filename);
    std::string line;
    while(getline(in, line))
    {
        contents.push_back( std::vector<double>());
        std::vector<double>& line_contents = contents[contents.size()-1];
        std::istringstream iss(line);
        do
        {
            double value;
            if( iss >> value )
            {
                line_contents.push_back(value);
                
            }
            else
            {
                line_contents.push_back( 0 );
                iss.clear();
            }
            std::string fieldend;
            getline(iss, fieldend, ',');
        }
        while(!iss.eof());
    }
    return contents;
}

int main( int argc , char **argv )
{
    ////// CREATE FOLDER FOR SAVING DATA //////////
    cout << fixed; cout << setprecision(4); // SET THE PRECISION FOR PRINTING OUT STUFF ON THE SCREEN
    save=1; // save = 0 means no files will be saved, save = 1 means files will be saved
    string saveDirectory   = argv[1]; // directory to save files to (this is given in the bash script)
    int numSims            = atoi(argv[2]); // number of simulations
    int resC0              = atoi(argv[3]); // log10(inital resistance frequency of C_0) (set to 0 in the bash sript -> which means C_0(0) = 0 (not 1e0))
    int resC1              = atoi(argv[4]); // log10(inital resistance frequency of C_1) (set to 6 in the bash sript -> which means C_1(0) = 1e6)
    int h1bin              = atoi(argv[5]); // binary indicator of interphyla conjugation (=1, conjugation happens between C_0 and C_1, if =0, only intraphyla conjugation
    ///////////////////////////// SET RNG ///////////////////////////////
    struct timeval t1;
    gettimeofday(&t1, NULL);
    int seed = (t1.tv_usec * t1.tv_sec);
    seed     = abs(seed);
    srand(seed);
    cout << "Today RNG seed to srand() : "<< seed << endl;
    /////////////////////// SET THE PARAMETERS //////////////////////////
    initiateSampleTimes(25); // allocate vector for the sampling times (this is the T_df in the paper, takes 25 different values, which are 0,15,30,...,360.
    initiateSamplePops(25); // allocate vector for the population abundances at sampling times (t) for [t,C_0,C_0^(+),C_1,C_^(+)]
    mu = 0.0;kappa = 0; // These are mu and kappa parameters in the paper
    aVec <<=-9.29999981645995,0,0,-0.00102786454039388,0,-1.21900593182689,0,-1.22599659767443e-05,2.87504542575169e-06,0,-1.81177486681986,0,-0.0853374375569896,0.0118166995595237,0,-93.2701506295908;
    rVec <<= 0.2720,0.2720,0.7286,0.7286,0.5884,0,0.3412,0;
    for(unsigned i=0;i<aVec.size();++i){
        aVec(i) = aVec(i)/popScale;
    }
    calcSS(); // Calculates the steady state population abundances of the system given the growth and interaction parameters
    ds            <<= 0.288,0.398,0.449,0.395; // death rates of sensitive strains C_0,C_1,C_2,C_3
    dr            <<= 0,0,0,0; // death rates of resistant strains C_0+,C_1+,C_2+,C_3+
    drugKillRates <<= mixVecs(ds,dr);
    xICSum        <<= sol01,sol02,sol03,sol04; // Initial condition is the steady state solution
    parameterVector <<= 0.1571,0.0122,-15.8309,-15.9230,0.0131;// rho_0 (cost of resistance for C_0+),rho_1 (cost of resistance for C_1+),log10(h_intra),log10(h_inter),gamma
    updateVectors(parameterVector);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // is and ir are not relevant to the current model.
    // If the user would like to use random recolonization events then needs to adjust the abundances with these vectors
    is            <<= 0,0,0,0; // sensitive strain abundances during recolonization
    ir            <<= 0,myPow(10,resC1),0,0; // resistant strain abundances during recolonization
    infAbundVec   <<= mixVecs(is,ir);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    resAbs        <<= myPow(10,resC0),myPow(10,resC1),0,0; // Initial abundance of the resistant strains
    trtTime         = 1000; // Treatment time (maximum 1000 days)
    obsTime         = trtTime+360+5; // Time the patient is observed -> 365 days more than the treatment time
    trtCountVector<<= 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20; // Number of treatments that will be applied (parameter N in the paper, varies between 1 and 20)
    resizeVec(trtCountVector,20); // reallocates the vector so that it's size is 20
    infCountVector<<= 0; // Number of recolonizations that will be applied (set to 0) and would be used like the variable "totalTrts" in a loop. But since there is no random recolonizations, this variable is taken out of the loop and set to 0.
    resizeVec(infCountVector,1); // reallocates the vector so that it's size is 1
    totalInfs = 0; // No recolonizations are considered for now, but user can change this to have random recolonization events
    Ca = 4; Cx = 4; // These are the tuning parameters of the hybrid stochastic - deterministic ODE solver.
    for(unsigned k=0; k<trtCountVector.size(); ++k){ // loop over the number of treatments (N in the paper)
        totalTrts  = trtCountVector(k); // variable for the total number of treatments (N)
        for (unsigned s=0; s<numSims; ++s){ // loop over the number of simulations per a given number of treatments N
            tspanInf   = scheduleRandomInfection(trtTime,obsTime,totalInfs); // schedules the recolonization times (which is null, no random recolonizations)
            tspanTrt   = scheduleRandomTreatment(trtTime,obsTime,totalTrts); // schedule the treatment times
            // treatment times are schedules as follows,
            // A matrix is generated with 3 columns. First column indicates the starting time, second column indicates the ending time, and third column is a binary indicator of the application of treatment (=0, no treatment, =1, treatment)
            // As an example : say fot the case N=1, the patient had no treatments between day 0 and 10, and a 5 day treatment between day 10 and 15, the matrix would be
            // [0 10 0
            // 10 15 1]
            folderChar = setDirectoriesFullRandomScheduling(saveDirectory,trtTime,totalTrts,totalInfs,save); // Sets the directory name for saving the files
            runGutStocHybrid(tspanTrt,tspanInf,folderChar,save); // runs the model for all the given settings above
        }
    }
    return 0;
}
