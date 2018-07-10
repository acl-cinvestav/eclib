//VERSION ORIGINAL DE JULIO C. ZAMORA 15-11-04
//Contador de bits de orden 1  ZAMORA  3-05-05

#include "CliffordLib.h"

int tp,tq,dim;
float *Re;
int signs[512][512];
//Clifford e0,e1,e2,e3,e4,e5,E,Ie,Ic,e,eo;
//VECTORIAL BASIS 
Clifford e0,e1,e2,e3,e4,e5,e6,e7,e8,e9;
Clifford e,eo;
Clifford ex,ey,ez;
Clifford eox,eoy,eoz;
Clifford Esd,Eds,Isd,Ids;
Clifford Ie;
Clifford E,Ic;

//BLADES______
Blade::Blade() 
{
  mBase=0;
  mVal=0.0;
}; 
//_____________________
Blade::Blade(int tBase) 
{
  mBase=tBase;
  mVal=1.0;
}; 
//_________________________________
Blade::Blade(int tBase,float tVal) 
{
  mBase=tBase;
  mVal=tVal;
};
//_____________
Blade::~Blade()
{
}
//___________SUMA DE CLIFFORDS_______________
Clifford operator+(Clifford A, Clifford B) {
	Clifford C; 
	idx a=A.begin();
	idx b=B.begin();
	while(a!=A.end()&&b!=B.end())
	{
		//switch
		int d=(a->mBase)-(b->mBase);
		if(d<0)
		{C.push_back(Blade(a->mBase,a->mVal));a++;}
		else if(d>0)
		{C.push_back(Blade(b->mBase,b->mVal));b++;}
		else if(0==d)
		{
			float Val=a->mVal+b->mVal;
			if(fabs(Val)>Error2){C.push_back(Blade(a->mBase,Val));}
			a++;b++;
		}
	}
	while(b!=B.end())
	{	
		if(fabs(b->mVal)>Error2){C.push_back(Blade(b->mBase,b->mVal));}
		b++;
	}
	while(a!=A.end())
	{	
		if(fabs(a->mVal)>Error2){C.push_back(Blade(a->mBase,a->mVal));}
		a++;
	}
  return C;
};
//___________RESTA DE CLIFFORDS_______________
Clifford operator-(Clifford A, Clifford B) {
	Clifford C; 
	idx a=A.begin();
	idx b=B.begin();
	while(a!=A.end()&&b!=B.end())
	{
		//switch
		int d=(a->mBase)-(b->mBase);
		if(d<0)
		{C.push_back(Blade(a->mBase,a->mVal));a++;}
		else if(d>0)
		{C.push_back(Blade(b->mBase,-b->mVal));b++;}
		else if(0==d)
		{
			float Val=a->mVal-b->mVal;
			if(fabs(Val)>Error2){C.push_back(Blade(a->mBase,Val));}
			a++;b++;
		}
	}
	while(b!=B.end())
	{	
		if(fabs(b->mVal)>Error2){C.push_back(Blade(b->mBase,-b->mVal));}
		b++;
	}
	while(a!=A.end())
	{	
		if(fabs(a->mVal)>Error2){C.push_back(Blade(a->mBase,a->mVal));}
		a++;
	}
  return C;
};

//______________________________________
//EXTRAE LA PARTE ESCALAR DEL MULTIVECTOR
//Y SE LA ASIGNA A UN DOUBLE
void operator<=(double &tv, Clifford A)
{
	if(A.capacity())
	{
		idx a=A.begin();
		tv=a->mVal;
	}
	else{tv=0.0;}
} 
void operator<=(float &tv, Clifford A)
{
	idx a=A.begin();
	if( A.size() )
	{
		tv=a->mVal;
	}
	else{tv=0.0;}
} 

bool isScalar(Clifford A)
{

	idx a = A.begin();

	bool R = true;

	while(a != A.end() )
	{

		if(fabs(a->mVal)> Error2 && a->mBase != 0){return false;}
		a++;
	}
	return R;
}

Clifford inverse(Clifford A)
{
	float sk; sk<=A;

	if(isScalar(A))
	{
		Clifford result;
		result.push_back(Blade(0, 1.0 / sk));
		return result;
	}

	Clifford xb = !A;
	Clifford temp = A|xb;
	Clifford inv = inverse(temp);
	return (xb*inv);
}

float abs(Clifford A)
{
	float sk; sk<=A;

	if(isScalar(A))
	{
		return sqrt(fabs(sk));
	}

	Clifford xb = !A;
	Clifford temp = A|xb;
	return abs(temp);


}

//_________PRODUCTO POR ESCALAR_________
Clifford operator*(float tv, Clifford A)
{
	Clifford R;
	idx a=A.begin();
	if(fabs(tv)>Error2)
	{
		while(a!=A.end())
		{	
		R.push_back(Blade(a->mBase,tv*a->mVal));
		a++;
		}
	return R;
	}
	else
	{
	R.push_back(Blade());
	}
	return R;
};

//_______PRODUCTO POR ESCALAR______________
Clifford operator*(Clifford A, float tv)
{
 return(tv*A);
};

//_______________________________________
Clifford operator/(Clifford A, float tv)
{
 return((1/tv)*A);
};

Clifford operator/(float tv, Clifford A)
{
	float n;
	n <= (A*(!A));
	return((tv/n)*(!A));
};

//_______PRODUCTO CLIFFORD_________________
Clifford operator*(Clifford A, Clifford B)
{
    Clifford R; 
	for(int i=0;i<dim;i++)Re[i]=0;

    idx a=A.begin();
	while(a!=A.end())
	{	
		int am=a->mBase;
		idx b=B.begin();
		while(b!=B.end())
		{	
			int bm=b->mBase;
			float prd=(a->mVal)*(b->mVal);
			if(fabs(prd)>Error2)
			{
				Re[am^bm]+=signs[am][bm]*prd;
			}
			b++;
		}
		a++;
	}

    
	for(int i=0;i<dim;i++)
	{if(fabs(Re[i])>Error2){R.push_back(Blade(i,Re[i]));}}
	
	//delete Re;
	return R;
}
//______________PRODUCTO WEDGE_____________
Clifford operator^(Clifford A, Clifford B)
{
   Clifford R; 
   //double *Re=new double[dim];
   for(int i=0;i<dim;i++)Re[i]=0;


   idx a=A.begin();
	while(a!=A.end())
	{
		int am=a->mBase;
		idx b=B.begin();
		while(b!=B.end())
		{	
			int bm=b->mBase;
			
			if((am&bm)==0)
			{
				float prd=(a->mVal)*(b->mVal);
				if(fabs(prd)>Error2)
					Re[am^bm]+=signs[am][bm]*prd;

			}
			b++;
		}
		a++;
	}

	for(int i=0;i<dim;i++)
	{if(fabs(Re[i])>Error2){R.push_back(Blade(i,Re[i]));}}
  //delete Re;
  return R;
}
//_____________DOT PRODUCT_________________ 
Clifford operator|(Clifford A, Clifford B)
{
   Clifford R;
   //double *Re=new double[dim];
   for(int i=0;i<dim;i++)Re[i]=0;

    idx a=A.begin();
	while(a!=A.end())
	{
		idx b=B.begin();
		int am=a->mBase;
		while(b!=B.end())
		{
			int bm=b->mBase;
			
			if(((am&bm)==am)||((am&bm)==bm))
			{
				float prd=(a->mVal)*(b->mVal);
				if(fabs(prd)>Error2)
					Re[am^bm]+=signs[am][bm]*prd;
			}
			b++;
		}
		a++;
	}
	
   for(int i=0;i<dim;i++)
   {if(fabs(Re[i])>Error2){R.push_back(Blade(i,Re[i]));}}
  
  //delete Re;
  return R;
}
//_____________JULIO PRODUCT_________________ 
Clifford operator%(Clifford A, Clifford B)
{
   Clifford R;
   //double *Re=new double[dim];
   for(int i=0;i<dim;i++)Re[i]=0;

    idx a=A.begin();
	while(a!=A.end())
	{
		idx b=B.begin();
		int am=a->mBase;
		while(b!=B.end())
		{
			int bm=b->mBase;
			int q=am&bm;
			float c=(a->mVal)*(b->mVal);
			if((c>Error2||c<(-Error2))&&((q!=am)&&(q!=bm)&&(q!=0)))
			{
				//int sig=sign2(am,bm);
				//Re[am^bm]+=sig*c;
				Re[am^bm]+=signs[am][bm]*c;
			}
			b++;
		}
		a++;
	}
	
   for(int i=0;i<dim;i++)
   {if(Re[i]!=0){R.push_back(Blade(i,Re[i]));}}
  
  //delete Re;
  return R;
}
//____________________
void EtoC(Clifford &p)
{
	p=p+0.5*(p|p)*e+eo;
}
//____________________
void CtoE(Clifford &p)
{
  float n;n<=(p|e);
  p=(Ie|p|Ie)/n;
}

Clifford E2toC(float x, float y)
{
	//e1.insert(e1.begin(),Blade(1,x));
	//e2.insert(e2.begin(),Blade(2,y));
	Clifford R;
	R.push_back(Blade(1,x));
	R.push_back(Blade(2,y));
	float d=(x*x+y*y);
	R.push_back(Blade(8,(d-1.0)/2.0));
	R.push_back(Blade(16,(d+1.0)/2.0));
	return R;
}
void Eto63(Clifford &p)
{
   float x;x<=(p|e1);
   float y;y<=(p|e2);
   float z;z<=(p|e3);
   p=p+0.5*(x*x*ex+y*y*ey+z*z*ez)+eo;
}

void Eto63_2(Clifford &p)
{
   float x;x<=(p|e1);
   float y;y<=(p|e2);
   float z;z<=(p|e3);
   p=p+0.5*(x*x*ex+y*y*ey+z*z*ez)+eo;//+x*y*(ex+ey)+x*z*(ex+ez)+z*y*(ex+ez)+eo;
}
//__________________________

void InitConformal()
{
	//dimention of algebra
	tp=4; tq=1;dim=(1<<(tp+tq));
	Re=new float[dim];

	for(int a=0;a<dim;a++)
	{
		for(int b=0;b<dim;b++)
		{
			signs[a][b]=sign2(a,b);	
		}
		
	}
	//init basis
	e0.clear();
	e1.clear();
	e2.clear();
	e3.clear();
	e4.clear();
	Ie.clear();
	 E.clear();
	Ic.clear();
	e0.insert(e0.begin(),Blade(0));
	e1.insert(e1.begin(),Blade(1));
	e2.insert(e2.begin(),Blade(2));
	e3.insert(e3.begin(),Blade(4));
	e4.insert(e4.begin(),Blade(8));
	e5.insert(e5.begin(),Blade(16));
	Ie.insert(Ie.begin(),Blade(7));
	 E.insert( E.begin(),Blade(24));
	Ic.insert(Ic.begin(),Blade(31));
	
	e=e4+e5;
	eo=(e5-e4)/2.0;



}

void Geometry(int ip,int iq)
{
//dimencion of algebra
	tp=ip; tq=iq;
	
	dim=(1<<(tp+tq));
	Re=new float[dim];

	for(int a=0;a<dim;a++)
	{
		for(int b=0;b<dim;b++)
		{
			signs[a][b]=sign2(a,b);	
		}
		
	}

	//init basis
	e0.clear();
	e1.clear();
	e2.clear();
	e3.clear();
	e4.clear();
	e5.clear();
	e6.clear();
	e7.clear();
	e8.clear();
	e9.clear();

	Ie.clear();
	Esd.clear();
	Eds.clear();

	eox.clear();
	eoy.clear();
	eoz.clear();
	ex.clear();
	ey.clear();
	ez.clear();
	e.clear();
	eo.clear();

	e0.insert(e0.begin(),Blade(0));
	e1.insert(e1.begin(),Blade(1));
	e2.insert(e2.begin(),Blade(2));
	e3.insert(e3.begin(),Blade(4));
	e4.insert(e4.begin(),Blade(8));
	e5.insert(e5.begin(),Blade(16));
	e6.insert(e6.begin(),Blade(32));
	e7.insert(e7.begin(),Blade(64));
	e8.insert(e8.begin(),Blade(128));
	e9.insert(e9.begin(),Blade(256));

	
	ex=e7+e4;
    ey=e8+e5;
    ez=e9+e6;
	
    eox=0.5*(e7-e4);
    eoy=0.5*(e8-e5);
    eoz=0.5*(e9-e6);
	
    e=ex+ey+ez;
    eo=eox+eoy+eoz;
  
    Ie=e1^e2^e3;
    Esd=ex^ey^ez^eo;
    Eds=eox^eoy^eoz^e;
    Isd=Ie^Esd;
    Ids=Ie^Eds;

  return ;
	
}
//_______________________
Clifford exp(Clifford& P)
{
  Clifford tP=P,R;
  R.push_back(Blade(0,1.0));
  R=R+P;  
  for(float n=2;n<=6;n++)
  {	
	tP = P * tP/n;
	if(tP.size())R=R+tP;
	else{n=6;}
  }
  return R;
}	

//_______________________
Clifford sqrt(Clifford& P)
{
	//TODO
  return P;
}


//_______________________
bool Zero(Clifford& P)
{
    idx b=P.begin();
	if(fabs(b->mVal)>1004.3)
				return false;

	/*
	while(b!=P.end())
	{
		if(fabs(b->mVal)>4.3)
				return false;
		b++;
	}*/
	return true;
}	

//  | PRODUCTO PUNTO same blade
float operator||(Clifford A, Clifford B)
{
	Clifford virtual_dot = A|B;
	return (*virtual_dot.begin()).mVal;
}



Clifford operator/(Clifford A, Clifford B)
{
	return A*inverse(B);
}


//SALIDA DE DATOS CLIFFORD
void cDump(Clifford Q)
{
 idx i=Q.begin();
 while(i!=Q.end())
	    {
		//SIGNO DEL BLADE	
		if ((*i).mVal > 0)printf("+");

		//VALOR DEL BLADE	
		printf("%f",(i->mVal));
		//NOMBRE DEL BLADE
        if(i->mBase){printf("e");;
		for(register int j=0;j<(tp+tq);j++)
		{
			if(((i->mBase)>>j)&1)
			{
				printf("%d",(j+1));
			}
		}
		}
	    i++;
	    }
		
   printf("\n");
 return ;
};

void cDumpG41(Clifford Q)
{

}

void cDumpPointG41(Clifford Q)
{
	idx i=Q.begin();
	float fe4, fe5;
	 while(i!=Q.end())
	 {
		 if(i->mBase == 8)
		 {
			 fe4 = i->mVal;
			 i++;

		 }else if(i->mBase == 16)
		 {
			 fe5 = i->mVal;
			 i++;
			 continue;
		 }else
		 {
		//SIGNO DEL BLADE
		if ((*i).mVal > 0)printf("+");

		//VALOR DEL BLADE
		printf("%f",(i->mVal));
		//NOMBRE DEL BLADE
	    if(i->mBase)
	    {
	    	printf("e");
		for(register int j=0;j<(tp+tq);j++)
		{
			if(((i->mBase)>>j)&1)
			{
				printf("%d",(j+1));
				}
			}
	    }
	    i++;
	 }
	 }
	 float temp = (fe5-fe4);
	 if (temp > 0)printf("+");
	 printf("%feo",temp);
	 temp = (fe5+fe4)/2.0;
	 if (temp > 0)printf("+");
	 printf("%fe",temp);

	 printf("\n");
	 return ;
}
//_________REVERSE_________
Clifford operator~(Clifford A)
{	float s=1;
	Clifford R;
	idx a=A.begin();
	while(a!=A.end())
	{	
		int r;
		r=bits(a->mBase);
		r=r*(r-1);
		if((r>>1)&1){s=-1;}
		R.push_back(Blade(a->mBase,s*a->mVal));
		a++;
	}
	return R;
}
//_________CONJUGADO_________
Clifford operator!(Clifford A)
{	float s=1;
	Clifford R;
	idx a=A.begin();
	while(a!=A.end())
	{	
		int r;
		r=bits(a->mBase);
		r=r*(r+1);
		if((r>>1)&1){s=-1;}
		R.push_back(Blade(a->mBase,s*a->mVal));
		a++;
	}
	return R;
}

/*
//____________________________
Clifford operator*(Clifford& A) 
{
	Clifford R=A*Ic;
	return R;
}
*/
Clifford operator+(float tv,Clifford A)
{
	Clifford R;R.push_back(Blade(0,tv));
	return (R+A);
}
Clifford operator+(Clifford A,  float tv)
{
	return (tv+A);
}
Clifford operator-(float tv,Clifford A)
{
	Clifford R;R.push_back(Blade(0,tv));
	return (R-A);
}
Clifford operator-(Clifford A,float tv)
{
	Clifford R;R.push_back(Blade(0,-1*tv));
	return (A+R);
}
/*
int sign(int A,int B)
{
 int bm=B,am=A,sig=1,n=0;
 //SIGNO DEL BLADE RESULTANTE
 unsigned int NegBases=((bm>>tp)&(am>>tp));	
 for(register int i=0;i<(tp+tq);i++)
	{
		if(((bm)>>i)&1)
    	{
			unsigned int tA=((am)>>(i+1));
			for(register int j=0;j<(tp+tq-i);j++)
			{
				if((tA>>j)&1){n++;}
			}	  
		}
		if(((NegBases>>i)&1)&&(i<tq)){n++;}	 
	}
	if(n&1){sig=-1;}
	return sig;
}
*/
/*
//SIGNO DEL BLADE RESULTANTE
int sign2(int A,int B)
{
	int sig=1;
	int tA2=((A>>tp)&(B>>tp));
	for(register int i=0;i<(tp+tq);i++)
	{
		if((B>>i)&1)
			tA2=tA2^(A>>(i+1));
	}
	if(bits(tA2)&1){sig=-1;}

return sig;
}
*/
//SIGNO DEL BLADE RESULTANTE
int sign2(int A,int B)
{
	int sig=1;
	int tA2=(A&B)>>tp;
	for(register int i=0;i<(tp+tq);i++)
	{
		A>>=1;
		if((B>>i)&1)
			tA2^=A;
	}
	if(bits(tA2)&1){sig=-1;}

return sig;
}

void CloseG()
{
delete [] Re;
}



//investigacion
int signo(int A,int B)
{
	int sig=1;
	int tA2=((A>>tp)&(B>>tp));
	//int tA2=0; 
	for(register int i=0;i<(tp+tq);i++)
	{
		if((B>>i)&1)
			tA2=tA2^(A>>(i+1));
	}

	if(bits(tA2)&1){sig=-1;}

//return bits(tA2)&1;
return sig;
}
//-----------------------------------------
int signo3(int A,int B)
{
	
	int tA2=0,n=(tp+tq); 
	for(int i=1;i<=n;i++)
		tA2^=(A>>i);

	tA2&=B;
	tA2^=(A&B)>>tp;
	tA2^=(tA2>>8);
	tA2^=(tA2>>4);
	tA2^=(tA2>>2);
	tA2^=(tA2>>1);
		
	if(tA2&1)return -1;
	return 1;
}

int signo2(int A,int B)
{
	//int sig=1;
	//int tA2=((A>>tp)&(B>>tp));
	A>>=1;
	//B&=0xF;
	for(register int i=0;i<3;i++)
		A^=A>>1;
	return bits((A&B))&1;
}

/*
//Imprimir tabla
void cTabla()
{
	for(int a=0;a<32;a+=2)
	{
		for(int b=0;b<32;b++)
		{
			afxDump<<" "<<signo2(a,b);	
		}
		afxDump<<"\n";	
	}
  return ;
};
*/
/*
void cValida()
{
	for(int a=0;a<128;a++)
	{
		for(int b=0;b<128;b++)
		{
			if(signo3(a,b)!=sign2(a,b))afxDump<<"("<<a<<","<<b<<")";	
		}
	}
  return ;
};
*/
