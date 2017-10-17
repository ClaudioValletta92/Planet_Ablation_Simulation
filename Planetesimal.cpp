//in this file we write all th classes and the methods used
#include "Planetesimal.h"		
#include <cmath>
#include <iostream>
#include <math.h>
using namespace std;

Planetesimal::Planetesimal() {
	m_Xmu=0.; 
	m_Pden=0.;
	m_Break=0.;
	m_E0=0.; 
	m_Ef=0.;
	m_v=0;
	m_vx=0;
	m_vy=0;

}

Planetesimal::Planetesimal(double rirat,char model,double size,double vx, double vy, double v,double x,double y,double r ) {
	m_size=size; 
	m_v=v;
	m_vx=vx;
	m_vy=vy;
	m_x=x;
	m_y=y;
	m_r=r;
	
	if (model=='W') 
		{
			m_Xmu=18.; 
			m_Pden=1.;
			m_Break=pow(10,6);
			m_E0=2.8*pow(10,8); 
			m_Ef=3.30*pow(10,9);
			m_rirat=0.;
		}
	else if (model=='I')
		{
			m_Xmu=56.; 
			m_Pden=7.8;
			m_Break=pow(10,8);
			m_E0=8.26*pow(10,10); 
			m_Ef=2.69*pow(10,9);
			m_rirat=0.;
		}
		else if (model=='R')
		{
			m_Xmu=50.; 
			m_Pden=3.4;
			m_Break=pow(10,8);
			m_E0=8.08*pow(10,10); 
			m_Ef=1.43*pow(10,10);
			m_rirat=0.;
		}
		else if (model=='M')
		{
			m_Xmu=18.; 
			m_Pden=2.;
			m_Break=pow(10,6);
			m_E0=2.8*pow(10,10); 
			m_Ef=3.3*pow(10,9)/(1+rirat);
			m_rirat=rirat;
		}
		
		else cout<<"Lettera sbagliata, stai attento"<<endl;
			
}
double Planetesimal::tmelt(){

		if (m_model=='W') 
		{
			return 273.;
		}
		else if (m_model=='I')
		{
			return 1.8*pow(10,3);
		}
		else if (m_model=='R')
		{
			return 1.8*pow(10,3);
		}
		else if (m_model=='M')
		{
			return 273.;
		}
}

double Planetesimal::cp(){

		if (m_model=='W') 
		{
			return 4.19*pow(10,7);
		}
		else if (m_model=='I')
		{
			return 6.91*pow(10,6);
		}
		else if (m_model=='R')
		{
			return 8.95*pow(10,6);
		}
		else if (m_model=='M')
		{
			return (4.19*pow(10,7)+m_rirat*8.95*pow(10,6))/(1.+m_rirat);
		}
}
double Planetesimal::tcrit(){
		
		if (m_model=='W') 
		{
			return 648.;
		}
		else if (m_model=='I')
		{
			return 4*pow(10,3);
		}
		else if (m_model=='R')
		{
			return 4*pow(10,3);
		}
		else if (m_model=='M')
		{
			return 648.;
		}
}

void Planetesimal::setsize(double size){
	m_size=size;
}
void Planetesimal::setrirat(double rirat){
	m_rirat=rirat;
}
void Planetesimal::setvx(double vx){
	m_vx=vx;
}
void Planetesimal::setvy(double vy){
	m_vy=vy;
}
void Planetesimal::setx(double x){
	m_x=x;
}
void Planetesimal::sety(double y){
	m_y=y;
}
char Planetesimal::getmodel(){
	return m_model;
}
void Planetesimal::setr(){
	m_r=sqrt(pow(m_x,2)+pow(m_y,2));
}
double Planetesimal::getx(){
	return m_x;
}
double Planetesimal::gety(){
	return m_y;
}
double Planetesimal::getr(){
	return m_r;
}
void Planetesimal::setv(){
	m_v=sqrt(pow(m_vx,2)+pow(m_vy,2));
}
double Planetesimal::getv(){
	return m_v;
}
double Planetesimal::getvx(){
	return m_vx;
}
double Planetesimal::getvy(){
	return m_vy;
}
double Planetesimal::getsize(){
	return m_size;
}
double Planetesimal::getxmu() {
	return m_Xmu;
}
double Planetesimal::getpden() {
	return m_Pden;
}
double Planetesimal::getstrength() {
	return m_Break;
}
double Planetesimal::getE0() {
	return m_E0;
}
double Planetesimal::getEf() {
	return m_Ef;
}
double Planetesimal::getrirat() {
	return m_rirat;
}
double Planetesimal::getmass() {
	
	return 4./3.*M_PI*m_Pden*pow(m_size,3);
}
double Starting_Point::con1(){
	return sqrt(6.673*pow(10,-8)*m_atem);
}
Layer::Layer(){
	m_r=0;
	m_rho=0;
	m_p=0;
	m_t=0;
	m_xx=0;
}
Layer::Layer(double r, double rho,double p, double t, double xx){
	m_r=r;
	m_rho=rho;
	m_p=p;
	m_t=t;
	m_xx=xx;
}
double Layer::getr(){
	return m_r;
}
double Layer::getrho(){
	return m_rho;
}
void Layer::print(){
	cout<<m_r<<" "<<m_rho<<" "<<m_p<<" "<<m_t<<" "<<m_xx<<endl;
}
double Layer::mg(){
	return 6.673*pow(10,-8)*m_xx;
}
double Layer::getp(){
	return m_p;
}
double Layer::gett(){
	return m_t;
}
double Layer::getxx(){
	return m_xx;
}
void Layer::setr(double r){
	m_r=r;
}
void Layer::setrho(double rho){
	m_rho=rho;
}
void Layer::setp(double p){
	m_p=p;
}
void Layer::sett(double t){
	m_t=t;
}
void Layer::setxx(double xx){
	m_xx=xx;
}
double Layer::Visc(){
	return 4.3*pow(10,-6)*sqrt(m_t);
}

Starting_Point::Starting_Point(){
	m_corem=0;
	m_v0=0;
	m_imp=0;
	m_fp=0;
	m_rmax=0;
	m_atem=0;
}

Starting_Point::Starting_Point(double corem, double v0,double imp,double fp, double rmax, double atem){
	m_corem=corem;
	m_v0=v0;
	m_imp=imp;
	m_fp=fp;
	m_rmax=rmax;
	m_atem=atem;
}
void Starting_Point::setcorem(double corem){
	m_corem=corem;
}
void Starting_Point::setv0(double v0){
	m_v0=v0;
}
double Planetesimal::deltam(double deltasize){
	return -12.57*m_size*m_size*deltasize*m_Pden;
}
void Starting_Point::setimp(double imp){
	m_imp=imp;
}
void Starting_Point::setfp(double fp){
	m_fp=fp;
}
void Starting_Point::setrmax(double rmax){
	m_rmax=rmax;
}
void Starting_Point::setatem(double atem){
	m_atem=atem;
}
double Starting_Point::getcorem(){
	return m_corem;
}
double Starting_Point::getv0(){
	return m_v0;
}
double Starting_Point::getimp(){
	return m_imp;
}
double Starting_Point::getfp(){
	return m_fp;
}
double Starting_Point::getrmax(){
	return m_rmax;
}
double Starting_Point::getatem(){
	return m_atem;
}
double Starting_Point::getexlam(){
	return (m_rmax+getsemia())/getsemia();
}
double Starting_Point::getecca(){
	return sqrt(1.+pow(m_v0,4)*m_imp*m_imp/getmg0()/getmg0());
}
double Starting_Point::getsemia(){
	return getmg0()/m_v0/m_v0;
}
double Starting_Point::getvx(){
	
	return sqrt(m_v0*m_v0+2.*getmg0()/m_rmax);
	
}
double Starting_Point::getfgrav0(){
	return -6.673*pow(10,-8)*m_atem;
}
double Starting_Point::getmg0(){
	return 6.673*pow(10,-8)*m_atem;
}
	
double Mach_Number(Planetesimal Halley,Layer boh){
	return Halley.getv()/sqrt(1.4*boh.getp()/boh.getrho());
}
double Reynold_Number(Planetesimal Halley,Layer boh){

	return Halley.getv()*boh.getrho()*Halley.getsize()/(boh.Visc());
}
double Energy(Planetesimal Halley, Layer boh){
	return 0.5*(pow(Halley.getvx(),2)+pow(Halley.getvy(),2))-boh.getxx()*6.673*pow(10,-8)/Halley.getr();
}
double dk(Planetesimal Halley,Layer boh){
double dkt=0;
if(Reynold_Number(Halley,boh)<=1)
	{ 
		return 1;
	}
        if(1<Reynold_Number(Halley,boh)<pow(10,3))
	{
		dkt=6/sqrt(Reynold_Number(Halley,boh));
		if(Mach_Number(Halley,boh)>1)
		{
			dkt=1.1-log(Reynold_Number(Halley,boh))/6.;
		}
	}
 	if(pow(10,3)<Reynold_Number(Halley,boh)<pow(10,5))
	{
		dkt=0.2;
		if(Mach_Number(Halley,boh)>1)
		{
			dkt=0.5;
		}
	}
	if(Reynold_Number(Halley,boh)>pow(10,5))
	{
		dkt=0.15;
		if(Mach_Number(Halley,boh)>1)
		{
			dkt=0.5;
		}
	
       }
       
       return dkt;
       
};
double Drag(Planetesimal Halley, Layer boh,double radmax){
	
		
	double  fpi,dkn,psi,f;
	if (Halley.getr()>radmax)
	{
		cout<<"No forza d'attrito"<<endl;
		return 0;
	}
	if(Reynold_Number(Halley,boh)<=1)
	{ 
		
        	fpi=1.019*pow(10,9)*Halley.getpden(); 
        	dkn=1./(Halley.getsize()*fpi); 
        	psi=1.+dkn*(1.25+.42*exp(-.87/dkn)); 
        	f=18.85*boh.Visc()*Halley.getv()*Halley.getsize()/psi; 
		cout<<"forza drag "<<f<<endl;

		
		return f;
	}


       f=dk(Halley,boh)*Halley.getsize()*Halley.getsize()*3.14159*boh.getrho()*Halley.getv()*Halley.getv();
	
	
       return f;

       
};

Runge_Kutta::Runge_Kutta(){
	m_t=0;
	m_dx=0;
	m_dy=0;
	m_dvx=0;
	m_dvy=0;
}
Runge_Kutta::Runge_Kutta(double t, double x,double y,double vx, double vy){
	m_t=t;
	m_dx=0;
	m_dy=0;
	m_dvx=0;
	m_dvy=0;
}	
double Runge_Kutta::Planetesimal_Gravitationalplusdrag_Equation_xaxe(Planetesimal Halley, Layer boh, Starting_Point Initial,double x,double 									     y, double vx, double vy, double vx1, double vy1)
{
	cout<<" la prima parte della forza è "<<Halley.getmass()*Initial.getfgrav0()*x/pow(x*x+y*y,1.5)<<" la seconda "<<-Halley.getmass()*Drag(Halley,boh,Initial.getrmax())*abs(x)/sqrt(x*x+y*y)*vx/(abs(vx)*Halley.getmass())<<endl;
	
	
	return Initial.getfgrav0()*x/pow(x*x+y*y,1.5)-Drag(Halley,boh,Initial.getrmax())/(vx1*vx1+vy1*vy1)*abs(x)*(vx*vx+vy*vy)/sqrt(x*x+y*y)*vx/(abs(vx)*Halley.getmass());
}
void Runge_Kutta::sett(double t)
{
	m_t=t;
}
double Runge_Kutta::Planetesimal_Gravitationalplusdrag_Equation_yaxe(Planetesimal Halley, Layer boh, Starting_Point Initial,double x,double 									     y, double vx, double vy,double vx1, double vy1)
{
	
	return Initial.getfgrav0()*y/pow(x*x+y*y,1.5)-Drag(Halley,boh,Initial.getrmax())/(vx1*vx1+vy1*vy1)*(vx*vx+vy*vy)*abs(y)/sqrt(x*x+y*y)*vy/(abs(vy)*Halley.getmass());
}
int Runge_Kutta::Solve(Planetesimal Halley, Layer boh, Starting_Point Initial,double precision,int counter){
	double x0,x1,x2,x3,x4=0;
	double vx0,vx1,vx2,vx3,vx4=0;
	double y0,y1,y2,y3,y4=0;
	double vy0,vy1,vy2,vy3,vy4=0;
	double dxdouble,dydouble, dvxdouble, dvydouble;
	double dx2,dy2,dvx2,dvy2;
	double dx1,dy1,dvx1,dvy1;
	double dxnormal,dynormal,dvxnormal,dvynormal;
	counter=0;
	do
	{
	x0=Halley.getx();
	vx0=Halley.getvx();
	vy0=Halley.getvy();
	y0=Halley.gety();
	cout<<"xo "<<x0<<" yo "<<y0<<" vxo "<<vx0<<"vyo "<<vy0<<endl; 
	vx1=m_t*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,boh,Initial,x0,y0,vx0,vy0,vx0,vy0);
        vy1=m_t*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,boh,Initial,x0,y0,vx0,vy0,vx0,vy0);
	x1=m_t*vx0;
	y1=m_t*vy0;   
	cout<<"x1 "<<x1<<" y1 "<<y1<<" vx1 "<<vx1<<" vy1 "<<vy1<<endl; 
	vx2=m_t*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,boh,Initial,x0+x1/2.,y0+y1/2.,vx0+vx1/2.,vy0+vy1/2.,vx0,vy0);
        vy2=m_t*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,boh,Initial,x0+x1/2.,y0+y1/2.,vx0+vx1/2.,vy0+vy1/2.,vx0,vy0);
     	x2=m_t*(vx0+vx1/2.);
        y2=m_t*(vy0+vy1/2.);
	cout<<"x2 "<<x2<<" y2 "<<y2<<" vx2 "<<vx2<<" vy2 "<<vy2<<endl;
        vx3=m_t*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,boh,Initial,x0+x2/2.,y0+y2/2.,vx0+vx2/2.,vy0+vy2/2.,vx0,vy0);
        vy3=m_t*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,boh,Initial,x0+x2/2.,y0+y2/2.,vx0+vx2/2.,vy0+vy2/2.,vx0,vy0);
	x3=m_t*(vx0+vx2/2.);
        y3=m_t*(vy0+vy2/2.);
  	vx4=m_t*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,boh,Initial,x0+x3,y0+y3,vx0+vx3,vy0+vy3,vx0,vy0);
        vy4=m_t*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,boh,Initial,x0+x3,y0+y3,vx0+vx3,vy0+vy3,vx0,vy0);
 	x4=m_t*(vx0+vx3);
	y4=m_t*(vy0+vy3);
	dxdouble=(x1+2.*x2+2.*x3+x4)/6; 
       	dydouble=(y1+2.*y2+2.*y3+y4)/6;
  	dvxdouble=(vx1+2.*vx2+2.*vx3+vx4)/6; 
       	dvydouble=(vy1+2.*vy2+2.*vy3+vy4)/6;

	//ora ho rifaccio con la metà del passo 
	
	vx1=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,boh,Initial,x0,y0,vx0,vy0,vx0,vy0);
        vy1=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,boh,Initial,x0,y0,vx0,vy0,vx0,vy0);
	x1=m_t/2.*vx0;
	y1=m_t/2.*vy0;   
	cout<<"x1 "<<x1<<" y1 "<<y1<<" vx1 "<<vx1<<" vy1 "<<vy1<<endl; 
	vx2=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,boh,Initial,x0+x1/2.,y0+y1/2.,vx0+vx1/2.,vy0+vy1/2.,vx0,vy0);
        vy2=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,boh,Initial,x0+x1/2.,y0+y1/2.,vx0+vx1/2.,vy0+vy1/2.,vx0,vy0);
     	x2=m_t/2.*(vx0+vx1/2.);
        y2=m_t/2.*(vy0+vy1/2.);
	cout<<"x2 "<<x2<<" y2 "<<y2<<" vx2 "<<vx2<<" vy2 "<<vy2<<endl;
        vx3=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,boh,Initial,x0+x2/2.,y0+y2/2.,vx0+vx2/2.,vy0+vy2/2.,vx0,vy0);
        vy3=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,boh,Initial,x0+x2/2.,y0+y2/2.,vx0+vx2/2.,vy0+vy2/2.,vx0,vy0);
	x3=m_t/2.*(vx0+vx2/2.);
        y3=m_t/2.*(vy0+vy2/2.);
  	vx4=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,boh,Initial,x0+x3,y0+y3,vx0+vx3,vy0+vy3,vx0,vy0);
        vy4=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,boh,Initial,x0+x3,y0+y3,vx0+vx3,vy0+vy3,vx0,vy0);
 	x4=m_t/2.*(vx0+vx3);
	y4=m_t/2.*(vy0+vy3);
	dx1=(x1+2.*x2+2.*x3+x4)/6; 
       	dy1=(y1+2.*y2+2.*y3+y4)/6;
  	dvx1=(vx1+2.*vx2+2.*vx3+vx4)/6; 
       	dvy1=(vy1+2.*vy2+2.*vy3+vy4)/6;
	x0=x0+dx1;
	y0=y0+dy1;
	vx0=vx0+dvx1;
	vy0=vy0+dvy1;
	
	//ecco fatto primo
	//mo faccio il secondo
	vx1=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,boh,Initial,x0,y0,vx0,vy0,vx0,vy0);
        vy1=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,boh,Initial,x0,y0,vx0,vy0,vx0,vy0);
	x1=m_t/2.*vx0;
	y1=m_t/2.*vy0;   
	cout<<"x1 "<<x1<<" y1 "<<y1<<" vx1 "<<vx1<<" vy1 "<<vy1<<endl; 
	vx2=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,boh,Initial,x0+x1/2.,y0+y1/2.,vx0+vx1/2.,vy0+vy1/2.,vx0,vy0);
        vy2=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,boh,Initial,x0+x1/2.,y0+y1/2.,vx0+vx1/2.,vy0+vy1/2.,vx0,vy0);
     	x2=m_t/2.*(vx0+vx1/2.);
        y2=m_t/2.*(vy0+vy1/2.);
	cout<<"x2 "<<x2<<" y2 "<<y2<<" vx2 "<<vx2<<" vy2 "<<vy2<<endl;
        vx3=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,boh,Initial,x0+x2/2.,y0+y2/2.,vx0+vx2/2.,vy0+vy2/2.,vx0,vy0);
        vy3=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,boh,Initial,x0+x2/2.,y0+y2/2.,vx0+vx2/2.,vy0+vy2/2.,vx0,vy0);
	x3=m_t/2.*(vx0+vx2/2.);
        y3=m_t/2.*(vy0+vy2/2.);
  	vx4=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_xaxe(Halley,boh,Initial,x0+x3,y0+y3,vx0+vx3,vy0+vy3,vx0,vy0);
        vy4=m_t/2.*Planetesimal_Gravitationalplusdrag_Equation_yaxe(Halley,boh,Initial,x0+x3,y0+y3,vx0+vx3,vy0+vy3,vx0,vy0);
 	x4=m_t/2.*(vx0+vx3);
	y4=m_t/2.*(vy0+vy3);
	dx2=(x1+2.*x2+2.*x3+x4)/6; 
       	dy2=(y1+2.*y2+2.*y3+y4)/6;
  	dvx2=(vx1+2.*vx2+2.*vx3+vx4)/6; 
       	dvy2=(vy1+2.*vy2+2.*vy3+vy4)/6;
	cout<<"sono quiiiii"<<endl;
	counter++;
	cout<<counter<<endl;
	dxnormal=dx1+dx2;
	dynormal=dy1+dy2;
	if(abs(sqrt(dxdouble*dxdouble+dydouble*dydouble)-sqrt(dxnormal*dxnormal+dynormal*dynormal))/15.>precision)
	{
		m_t=m_t/2.;
	}
	
	} while(abs(sqrt(dxdouble*dxdouble+dydouble*dydouble)-sqrt(dxnormal*dxnormal+dynormal*dynormal))/15.>precision);
	
	
	
	
	//m_dx=(x1+2.*x2+2.*x3+x4)/6; 
       	//m_dy=(y1+2.*y2+2.*y3+y4)/6;
  	//m_dvx=(vx1+2.*vx2+2.*vx3+vx4)/6; 
       	//m_dvy=(vy1+2.*vy2+2.*vy3+vy4)/6;
	m_dx=dxdouble; 
       	m_dy=dydouble;
  	m_dvx=dvxdouble; 
       	m_dvy=dvydouble;

return counter;
}

double Runge_Kutta::getdx(){
	return m_dx;
}
double Runge_Kutta::getdy(){
	return m_dy;
}
double Runge_Kutta::getdvx(){
	return m_dvx;
}
double Runge_Kutta::getdvy(){
	return m_dvy;
}
double Runge_Kutta::getds(){
	return sqrt(pow(m_dx,2)+pow(m_dy,2));
}
double Runge_Kutta::getdv(){
	return sqrt(pow(m_dvx,2)+pow(m_dvy,2));
}
double Ablate(Planetesimal Halley, Layer boh, double t,char model,double rmax){
	double Esput,Eabl,con,energy,tlow,thigh,T,pvap,evap,R;	
	
	double Q=0.5*dk(Halley,boh);
	if(Halley.getr()>rmax)
	{
		Q=0;
		
	}
	double Eint=5.67*pow(10,-5)*pow(boh.gett(),4);
	if(Halley.getr()>rmax)
	{
		Eint=Eint*pow(rmax/Halley.getr(),2);
	}
	if (model=='M')
	{
		
		Esput=Q*boh.getrho()*pow(Halley.getv(),3);
		Eabl=Eint+0.25*Esput-5.67*pow(10,-5)*pow(Halley.tmelt(),4);
		if (Eabl<0)
			{
				Eabl=0;
			}
		return -Eabl*t/Halley.getpden()/(Halley.cp()*Halley.tmelt()+Halley.getEf()); 
	}
	
	//ora faccio la vaporizzazione
	con=Q*boh.getrho()*pow(Halley.getv(),3);
	energy=Eint+con/4.;
	tlow=10.;
	thigh=7.*pow(10,4);
	do
	{
		T=sqrt(tlow*thigh);
		if (Halley.getmodel()=='W') 
			{
				cout<<"water"<<endl;
				pvap=pow(10.,(-2104.2/T+5.5901))*pow(10,6);
			}
		else if (Halley.getmodel()=='I')
			{
				cout<<"iron"<<endl;
				pvap=pow(10.,(12.509-20014./T));
				if(T>3840.)
				{
					pvap=61.9*exp(3.301D-3*T); 
				}
			
			}
		else if (Halley.getmodel()=='R')
			{
				cout<<"rock"<<endl;
				pvap=pow(10.,(13.176-24605./T));
			}
		else if (Halley.getmodel()=='M')
			{
				cout<<"mixture"<<endl;
				pvap=pow(10.,(-2104.2/T+5.5901))*pow(10,6);
			}
		evap=pvap*sqrt(Halley.getxmu()/(5.22*pow(10,8)*T));
		R=5.67*pow(10,-5)*pow(T,4)+evap*Halley.getE0();
		if(R>=energy)
		{
			
			thigh=T;
		} 
        	if(R<energy)
		{
			
			tlow=T;
		}
	}
	while(abs(1.-R/energy)>pow(10,-2));
	

	if(T>Halley.tcrit())
	T=Halley.tcrit();
	evap=(Eint+.25*Q*boh.getrho()*pow(Halley.getv(),3)-5.67*pow(10,-5)*pow(T,4))/Halley.getE0();
	if(evap<0)
	{
		evap=0;
	}
	if(Halley.getmodel()=='M'){
		return -evap*t/Halley.getpden()*3.;
	}
	return -evap*t/Halley.getpden()*1.;


		
}
double Time_Step(Starting_Point initial, Planetesimal Halley){
	
	return initial.getfp()*6.28*pow(Halley.getr(),1.5)/initial.con1();
} 

double Planetesimal::angolar_momentum(){
	return m_x*m_vy-m_y*m_vx;
}

double Starting_Point::impacty(){
	return m_imp*(getexlam()-1.)/sqrt(getexlam()*getexlam()-1.);
}
double pdyn(Planetesimal Halley,Layer boh)
{
	return Halley.getv()*Halley.getv()*boh.getrho()/2.;
}
double rdyn(Planetesimal Halley,Layer boh)
{
	return sqrt(10.*pdyn(Halley,boh)/(25.13*6.673*pow(10,-8)))/Halley.getpden();
}
