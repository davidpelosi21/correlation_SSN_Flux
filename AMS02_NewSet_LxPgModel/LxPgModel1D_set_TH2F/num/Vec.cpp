#include "Vec.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
using namespace std;

Vec::Vec(int i)
{
    if(i==0){
        entries = NULL;
        N = 0;
    }
    else{
        N = i;
        entries = new double[N];
        for(int a=0; a<N; ++a) entries[a] = 1;
    }
}
Vec::Vec(int i, double val): N(i) {
    entries = new double[N];
    for(int a=0; a<N; a++)
        entries[a] = val;
}

Vec::Vec(int i,  double* v): N(i){
    entries = new double[N];
    for(int a=0; a<N; a++){
        entries[a]=v[a];
    }
}

Vec::Vec(const Vec& vec) : N(vec.size())  {
    entries = new double[N];
    for(int a = 0; a<N; a++)
        entries[a]=vec[a];
}

Vec::~Vec(){
    if(entries)
        delete[] entries;
}


int Vec::size() const{
    return N;
}

void Vec::SetEntries(int n, double* v){
    delete[] this->entries;
    N=n;
    entries = new double[N];
    for(int a=0; a<N;a++)
        entries[a]=v[a];
}


double Vec::dot(const Vec& v) const {
    if(this->size()!=v.size()){
        printf("Operação Abortada. Vetores possuem tamanhos diferentes.");
        return 0;
    }
    
    double sum=0;
    for(int i=0; i<N; i++)
        sum+=(entries[i]*v[i]);
    return sum;
}

void Vec::swap(int i, int f){
    double aux;
    
    aux=entries[i];
    entries[i]=entries[f];
    
    entries[f]=aux;
}

void Vec::Print(string opt){
    if(opt=="0"){
        for(int i=0; i<N;i++)
            printf(" %8.3lf ", entries[i]);
        printf("\n");
    }
    
    if(opt=="1"){
        printf("\n");
        for(int i=0; i<N;i++)
            printf("%10.3f\n", entries[i]);
        printf("\n");
    }
}


double Vec::operator[](int i) const{
    return entries[i];
}

double& Vec::operator[](int i){
    return entries[i];
}

Vec& Vec::operator=(const Vec& v){
    if(this!=&v){
        N=v.N;
        delete [] entries;
        entries = new double[N];
        for(int a=0; a<N; a++)
            entries[a]=v.entries[a];
    }
    
    return *this;
}

Vec Vec::operator+=(const Vec& v){
    for(int a=0; a<N; a++)
        entries[a]+=v.entries[a];
    
    return *this;
}

Vec Vec::operator+(const Vec& v) const{
    
    if(this->size()!=v.size()){
        printf("Operação Abortada. Vetores Somados possuem tamanhos diferentes.");
        return (*this);
    }
    
    Vec res;
    res = *this;
    for(int a=0; a<N; a++)
        res.entries[a]+=v.entries[a];
    
    return res;
}

Vec Vec::operator-=( Vec& v){
    for(int a=0; a<N; a++)
        entries[a]-=v.entries[a];
    
    return *this;
}


Vec Vec::operator-(const Vec& v) const
{
    if(this->size()!=v.size()){
        printf("Operação Abortada. Vetores Subtraídos possuem tamanhos diferentes.");
        return (*this);
    }
    
    Vec res;
    res = *this;
    for(int a=0; a<N; a++)
        res.entries[a]-=v.entries[a];
    
    return res;
}

Vec Vec::operator*(const Vec& v) const
{
    if(this->size()!=v.size()){
        printf("Operação Abortada. Vetores Multiplicados possuem tamanhos diferentes.");
        return (*this);
    }
    
    Vec res;
    res = *this;
   	for(int a=0;a <N; a++)
        res.entries[a]*=v.entries[a];
    
   	return res;
}

Vec Vec::operator*(double m) const
{
    Vec res;
    res = *this;
    for(int i=0; i<N; ++i)
        res.entries[i] = entries[i]*m;
    return res;
}

Vec Vec::operator-(){
    for(int i=0; i<N; i++)
        entries[i]=-entries[i];
    return *this;
}

Vec operator* (double m, Vec& v){
    int a = v.size();
    Vec res(v);
    for(int i=0;i<a;i++)
        res[i]*=m;
    
    return res;
}

Vec Vec::ReadFile(string FILE){
    ifstream Data (FILE.c_str());
    double aux;
    vector<double> k;
    
    if (!Data.is_open()) {
        printf("\nERRO: Ficheiro não aberto\n");
        exit(1);
    }
    
    while (Data >> aux)
        k.push_back(aux);
    
    double* b = new double[k.size()];
    
    for(int i=0; i<k.size(); i++)
        b[i]=k[i];
    
    return Vec(k.size(),b);
}


