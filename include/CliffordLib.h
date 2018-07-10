//VERSION ORIGINAL DE JULIO C. ZAMORA 15-11-04
#ifndef CLIFFORDLIB_HPP_
#define CLIFFORDLIB_HPP_

#include <iostream>
#include <math.h>
#include <vector>

#define Error2 0.000001
#define bits(x)(((Bx(x)+(Bx(x)>>4))&0x0F0F0F0F)%255)
#define Bx(x)((x)-(((x)>>1)&0x77777777)-(((x)>>2)&0x33333333)-(((x)>>3)&0x11111111))

//BLADES
class Blade{
	
	public:
	Blade();
	Blade(int tBase);
	Blade(int tBase,float tVal);
   ~Blade();
	
	float mVal;
	int mBase; 
	};
//DEFINITION
typedef std::vector<Blade> Clifford;
typedef std::vector<Blade>::iterator idx; 

//VECTORIAL BASIS 
//extern Clifford e0,e1,e2,e3,e4,e5;
extern Clifford E,Ic;
//extern Clifford e,eo;

//VECTORIAL BASIS 
extern Clifford e0,e1,e2,e3,e4,e5,e6,e7,e8,e9;
extern Clifford e,eo;
extern Clifford ex,ey,ez;
extern Clifford eox,eoy,eoz;
extern Clifford Esd,Eds,Isd,Ids;
extern Clifford Ie;


void operator<=(double &, Clifford);
void operator<=(float &, Clifford);
//  SCALAR PRODUCT
//Clifford operator*(double,Clifford&);
Clifford operator*(float, Clifford);
//  SCALAR PRODUCT
//Clifford operator*(Clifford&,double);
Clifford operator*(Clifford, float);

// *  PRODUCTO CLIFFORD
Clifford operator*(Clifford, Clifford);

//  ^ PRODUCTO WEDGE
Clifford operator^(Clifford, Clifford);

//  | PRODUCTO PUNTO
Clifford operator|(Clifford, Clifford);

//  | PRODUCTO PUNTO same blade
float operator||(Clifford, Clifford);

//  | PRODUCTO ZAMORA
Clifford operator%(Clifford, Clifford);


//  +  SUMA
Clifford operator+(Clifford, Clifford);
//Clifford operator+(double,Clifford&);
Clifford operator+(float, Clifford);
//Clifford operator+(Clifford&,double);
Clifford operator+(Clifford, float);

//  -  RESTA
Clifford operator-(Clifford, Clifford);
//Clifford operator-(double,Clifford&);
Clifford operator-(float, Clifford);
//Clifford operator-(Clifford&,double);
Clifford operator-(Clifford, float);

//divicion
//Clifford operator/(Clifford&,double);
Clifford operator/(Clifford, float);
//Clifford operator/(double,Clifford&);
Clifford operator/(float, Clifford);

//Clifford operator/(Clifford&,double);
Clifford operator/(Clifford, Clifford);

//unarios
Clifford operator~(Clifford A);
Clifford operator!(Clifford A);
Clifford operator*(Clifford A);

//FUNCTIONS
void Geometry(int ip,int iq);
void InitConformal();

void cDump(Clifford Q);
void cDumpPointG41(Clifford Q);
void cDumpG41(Clifford Q);
int sign(int A,int B);
int sign2(int A,int B);
int sign3(int A,int B);
void CloseG(void);

bool isScalar(Clifford A);
Clifford inverse(Clifford A);
float abs(Clifford A);

//CONFORMAL MAPPING
void EtoC(Clifford &);
Clifford  E2toC(float x, float y);
void CtoE(Clifford &);
void Eto63(Clifford &p);
void Eto63_2(Clifford &p);
bool Zero(Clifford &p);
void InitConformal();

//MATH TOOLS
Clifford exp(Clifford&);
Clifford sqrt(Clifford& P);
//void cTabla();
//void cValida();
int signo(int A,int B);
int signo2(int A,int B);

#endif
