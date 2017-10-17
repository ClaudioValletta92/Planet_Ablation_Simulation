
#ifndef __planetesimal_h__ 
#define __planetesimal_h__ 
//Definizions of all the classes and the functions
;class Planetesimal {
public:
	//costruttore di default
	Planetesimal(); 
/* mean molecular weight xmu,planetesimal density pden,material strength strength, energy of vaporization E0, heat of fusion particle Ef, size radius of planetesimal  */
	
	Planetesimal(double rirat,char model,double size,double v_x,double v_y,double v, double x,double y, double r); //costruttore completo
	//singoli costruttori e metodi per ottenere il valore delle variabili
	void setsize(double size);
	void setrirat(double rirat);
	void setvx(double vx);
	void setvy(double vy);
	void setv();
	void setx(double x);
	void sety(double y);
	void setr();
	void print();
	double getxmu();
	char getmodel();
	double getpden();
	double getstrength();
	double getE0();
	double getEf();
	double getsize();
	double getrirat();
	double getmass();
	double getvx();
	double getvy();
	double getv();
	double getx();
	double gety();
	double getr();
	double tmelt();
	double cp();
	double deltam(double deltasize);

	double tcrit();
	double angolar_momentum();
	
private:
//le variabili della classe, definite private
	double m_Xmu, m_Pden, m_Break, m_E0, m_Ef, m_size, m_rirat, m_vx, m_vy, m_v, m_x, m_y, m_r;
	char m_model;
}; 

class Layer {
public:
	Layer();
	Layer(double r, double rho,double p, double t, double xx);
	void setr(double r);
	void setrho(double rho);
	void setp(double p);
	void sett(double t);
	void setxx(double xx);
	double getr();
	double getrho();
	double getp();
	double gett();
	double getxx();
	double Visc();
	double mg();
	void print();
	
private:
	double m_r, m_rho, m_p, m_t, m_xx;

};

class Starting_Point {
public:
	Starting_Point();
	Starting_Point(double corem,double v0,double imp, double fp, double rmax, double atem);
	void setcorem(double corem);
	void setv0(double v0);
	void setimp(double imp);
	void setfp(double fp);
	void setrmax(double rmax);
	void setatem(double atem);
	double getcorem();
	double getv0();
	double getimp();
	double getfp();
	double getrmax();
	double getatem();
	double getfgrav0();
	double getmg0();
	double getvx();
	double getsemia();
	double getecca();
	double getexlam();
	double con1();	
	double impacty();

private:
	double m_corem, m_v0, m_imp, m_fp, m_rmax, m_atem;
};

class Runge_Kutta{
	public:
	Runge_Kutta();
	Runge_Kutta(double t,double x,double y,double vx,double vy);
	double Planetesimal_Gravitationalplusdrag_Equation_xaxe(Planetesimal Halley, Layer boh, Starting_Point Initial,double x, double y, 									double vx, double vy,double vx1, double vy1);
	double Planetesimal_Gravitationalplusdrag_Equation_yaxe(Planetesimal Halley, Layer boh, Starting_Point Initial, double x, double y, 									double vx, double vy,double vx1, double vy1);
	int Solve(Planetesimal Halley, Layer boh, Starting_Point Initial, double precision,int counter);
	double getdx();
	double getdy();
	double getdvx();
	double getdvy();
	double getds();
	double getdv();
	void sett(double t);
	private:
	double m_t, m_dx, m_dy, m_dvx, m_dvy;
	
};
	
double Mach_Number(Planetesimal Halley,Layer boh);
double Reynold_Number(Planetesimal Halley,Layer boh);
double dk(Planetesimal Halley,Layer boh);
double Energy(Planetesimal Halley, Layer boh);
double Drag(Planetesimal Halley, Layer boh,double radmax);
double Ablate(Planetesimal Halley, Layer boh, double t,char model, double rmax);
double Time_Step(Starting_Point initial, Planetesimal Halley);
double pdyn(Planetesimal Halley,Layer boh);
double rdyn(Planetesimal Halley,Layer boh);

#endif 
