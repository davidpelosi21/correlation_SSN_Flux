#ifndef __Vec__
#define __Vec__
#include <cstdio>
#include <string>
using namespace std;

class Vec{
public:
    Vec(int i=1); //default constructor
    Vec(int, double); //set N elements equal to value
    Vec(int,  double*); //set N elements from array
    Vec(const Vec&); // copy constructor
    //////////
    ~Vec(); //destructor
    /////////
    void SetEntries (int, double*);
    int size() const; //Vec size
    double dot(const Vec&) const; //produto interno
    void swap(int, int); //swap Vec elements
    void Print(string opt = "0") ; //class dump
    /////////
    double& operator[](int);
    double operator[](int) const;//Vec is declared as const
    /////////
    Vec& operator=(const Vec&);
    Vec operator+=(const Vec&);
    Vec operator+(const Vec&) const;
    Vec operator-=(Vec&);
    Vec operator-(const Vec&) const;
    Vec operator*(const Vec&) const; //x1x2,y1y2,z1z2
    Vec operator*(double) const;//Vec.operator*(k) = Vec*scalar
    Vec operator-();
    friend Vec operator* (double,  Vec&);
    static Vec ReadFile(string);
    
private:
    int N;
    double* entries=NULL;
    
};

#endif
