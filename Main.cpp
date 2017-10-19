using namespace std							//usual declaration in c++
#include "Planetesimal.h"						//file with the definition of the classes and functions
#include <iostream>
#include <fstream>
#include <cmath>
/*in this program I compute the ablation of a planetesimal due to gas drag of the planet atmosphere. the atmosphere is divided in layers, each layer with its own thermodynamical quantities. In the first section we take the value of the layers from a file and put it into an array. In the second section we compute where the planetesimal is. In the third section we compute the motion of the planetesimal and the ablation by the atmosphere with the functions and classes defined in the included file. Finally we see if we can exit from the code that happens if
-the planetesimal is completely ablated
-the planetesimal is broken
-the planetesimal reached the core
-the planetesimal is in a loop completely contained in the innermost layer of the atmosphere
What do I mean by loop? the planetesimal can be stuck in a loop with the radius going up and down and if the ablation is too slow the computational time is too long, so if the loop is completely contained in the innermost layer we exit from the program and say that all the planetesimal is ablated in the innermost layer
*/
int main ()
{
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------FIRST PART OF THE CODE--------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
	int Layer_Number;						//the number of the layers used
	char model,ablation;						//il modello per il planetesimal e per l'ablazione/vaporizzazione
	double dr,dden,dtemp,dpress,dxm;
	double rirat, size, v_x, v_y, v, x, y, rad;
	double r, rho, p, t,xx;
	double fdrag;
	double edep;
	double time=0;
	double dedr; 
	double rdep;
	double looppoint;
	int counter;
	double impact;

	double depint=0;
	double fgravx, fdragx;
	double fgravy, fdragy;
	double newx,newy,newvx,newvy;
	double oldr,newr;
	double difference1,difference2;
	double ratio;
	Layer local;
	Runge_Kutta SolveEqMotion;
	double corem, v0, imp, fp, rmax, atem;
	double rockm=0;
	ifstream in("Build.dat");					//file for the initial parameter
	ofstream out;							//files for the output
	out.open("Trajectory.dat", ios::out);
	
	

	in>>Layer_Number;
	Layer *atmosphere= new Layer[Layer_Number]();			//array of object of the class atmosphere
	in>>size>>ablation>>model>>rirat;
	in>>corem>>v0>>imp>>fp>>rmax>>atem;
	
	for(int a=Layer_Number-1;a>=0;a--)				//i put the right values in the atmosphere 
	{
		in>>r>>rho>>p>>t>>xx;
		atmosphere[a].setr(r);
		atmosphere[a].setrho(rho);
		atmosphere[a].setp(p);
		atmosphere[a].sett(t);
		atmosphere[a].setxx(xx);
		cout<<"layer number "<<104-a<<endl;
		atmosphere[a].print();
	}
do{		
//	Starting_Point initial(corem,v0,pow(10,11),fp,rmax,atem);	//object of the class initial, it contains all the information of the initial status of the calculation
	Starting_Point initial(corem,v0,impact,fp,rmax,atem);
	Planetesimal Halley(rirat,model,size,initial.getvx(),-0.1,sqrt(pow(initial.getvx(),2)+0.1*0.1),-rmax,initial.getimp(),sqrt(impact*impact+rmax*rmax));						//object of the class planetesimal
	looppoint=rmax;					
	
	in.close();
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------SECOND PART OF THE CODE--------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------



	for(int c=1;c<1000000000000000000;c++)					//the principal for cycle, where we calculate the ablation and the equation of motion
	{
	int i=Layer_Number-1;
	while(atmosphere[i].getr()<Halley.getr())			//here i calculate where is the planetesimal, in which layer
	{	
		if(i==0){
				break;
			}
		i--;
	}
	
	dr=atmosphere[i].getr()-Halley.getr();				//here i calculate what is the environment araound the planetesimal, so i comput the derivatives in all the thermodynamics quantities and i calculate the temperature, pressure and so on that are present around the planetesimal
	
	
	if(dr<0)
	{
		dr=0;
	}

	cout<<" the i is "<<i<<endl;
	
	dden=(-atmosphere[i].getrho()+atmosphere[i+1].getrho())/(-atmosphere[i].getr()+atmosphere[i+1].getr());
	dtemp=(-atmosphere[i].gett()+atmosphere[i+1].gett())/(-atmosphere[i].getr()+atmosphere[i+1].getr());
	dpress=(-atmosphere[i].getp()+atmosphere[i+1].getp())/(-atmosphere[i].getr()+atmosphere[i+1].getr());
	dxm=(-atmosphere[i].getxx()+atmosphere[i+1].getxx())/(-atmosphere[i].getr()+atmosphere[i+1].getr());	
	
	local.setr(Halley.getr());
	local.setrho(atmosphere[i].getrho()-dden*dr);
	local.sett(atmosphere[i].gett()-dtemp*dr);
	local.setp(atmosphere[i].getp()-dpress*dr);
	local.setxx(atmosphere[i].getxx()-dxm*dr);
	oldr=Halley.getr();
	

	cout<<"local layer"<<endl;
	local.print();		
	cout<<"that must be between"<<endl;
	atmosphere[i+1].print();
	cout<<"and"<<endl;				
	atmosphere[i].print();
/*
	fgravx=initial.getfgrav0()*Halley.getx()/pow(Halley.getx()*Halley.getx()+Halley.gety()*Halley.gety(),1.5);
	fdragx=-Drag(Halley,local,initial.getrmax())*abs(Halley.getx())/sqrt(Halley.getx()*Halley.getx()+Halley.gety()*Halley.gety())*Halley.getvx()/(abs(Halley.getvx())*Halley.getmass());

	fgravy=initial.getfgrav0()*Halley.gety()/pow(Halley.getx()*Halley.getx()+Halley.gety()*Halley.gety(),1.5);
	fdragy=-Drag(Halley,local,initial.getrmax())*abs(Halley.gety())/sqrt(Halley.getx()*Halley.getx()+Halley.gety()*Halley.gety	())*Halley.getvy()/(abs(Halley.getvy())*Halley.getmass());
*/


//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------THIRD PART OF THE CODE--------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
	
	//now i solve the equation of motion, so I calculate the drag force and than i solve the equation of motion
	fdrag=Drag(Halley,local,initial.getrmax());
	SolveEqMotion.sett(Time_Step(initial, Halley));	
	SolveEqMotion.Solve(Halley,local,initial,0.1*abs(atmosphere[i+1].getr()-atmosphere[i].getr()),counter);
	cout<<"il contatore qua e "<<counter<<endl;
	//out<<counter<<endl;
	
	
	cout<<"before solving the equation of motion x is "<<Halley.getx()<<endl;
	cout<<" and y is "<<Halley.gety()<<endl;
	cout<<"before solving the equation of motion vx is "<<Halley.getvx()<<" and vy is "<<Halley.getvy()<<endl;
	
	newx=Halley.getx()+SolveEqMotion.getdx();
	newvx=Halley.getvx()+SolveEqMotion.getdvx();
	newy=Halley.gety()+SolveEqMotion.getdy();
	newvy=Halley.getvy()+SolveEqMotion.getdvy();
	Halley.setx(newx);
	Halley.sety(newy);
	Halley.setvx(newvx);
	Halley.setvy(newvy);
	Halley.setr();
	Halley.setv();
	
	newr=Halley.getr();
	//now i check if we are stuck in a loop
	if(c%2!=0)
	{
		difference1=newr-oldr;
	}
	if(c%2==0)
	{
		difference2=newr-oldr;
	}
	ratio=difference2/difference1;
	//done
	

	
	edep=fdrag*SolveEqMotion.getds();
	//now that i solved the equation of motion i calculate the ablation of the planetesimal
	Halley.setsize(Halley.getsize()+Ablate(Halley,local,SolveEqMotion.gett(),ablation,rmax));
	rockm=rockm+Halley.deltam(SolveEqMotion.getds())*Halley.getrirat()/(1+Halley.getrirat());
	
	if(local.gett()>2.3*pow(10,3))
	{
		edep=edep-8.08*pow(10,10)*rockm;
          	rockm=0.;
	}
	edep=edep-Halley.deltam(SolveEqMotion.getds())*Halley.getE0()/(1.+Halley.getrirat())+rockm*(local.mg()/Halley.getr()-local.mg()/oldr);

//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------FOURTH PART OF THE CODE--------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------

	if(Halley.getsize()<=0)
	{
		cout<<"The planetesimal is ablated, congrats"<<endl;
		return 0;
	}
	time=time+SolveEqMotion.gett();
	if(ratio<0)
	{
		
		if(looppoint<atmosphere[Layer_Number-2].getr()&&Halley.getr()<atmosphere[Layer_Number-2].getr())
		{
			cout<<"all the loop is in the innermost layer!"<<endl;
			return 0;
		}
		looppoint=Halley.getr();
		//out<<"looppoint is "<<looppoint<<endl;
	}
	if(Halley.getr()<atmosphere[Layer_Number-1].getr())
	{
		cout<<"The planetesimal reached the solid part of the planet"<<endl;
		edep=edep+(Halley.getmass()+rockm)*Halley.getv()*Halley.getv()/2.;
        	depint=depint+edep;
        	dedr=edep/SolveEqMotion.getds();
		return 0;
	}
	
	rdep=(Halley.getr()+oldr)/2.;
	depint=depint+edep;

	if(pdyn(Halley,local)>Halley.getstrength()&&rdyn(Halley,local)>Halley.getsize())
	{
		cout<<"Everything broke up"<<endl;
		edep=Halley.getmass()*(Halley.getv()*Halley.getv()/2.-(Halley.getE0()+8.08*pow(10,10)*Halley.getrirat())/(1.+Halley.getrirat()))-rockm*8.08*pow(10,10);
		depint=depint+edep;
        	dedr=edep/SolveEqMotion.getds();
		return 0;
	}

	//out<<Halley.getr()<<" "<<Halley.getsize()<<endl;
	//out<<fgravx<<" "<<fdragx<<" "<<fgravy<<" "<<fdragy<<endl;
	//out<<Energy(Halley,local)<<endl;
	//out<<fgravx/fdragx<<" "<<fgravy/fdragy<<endl;
	out<<Halley.getx()<<" "<<Halley.gety()<<endl;

	}
	impact=impact+1*pow(10,11);	
	}while(impact<pow(10,14));
		

	
	
			
		


return 0;
}
