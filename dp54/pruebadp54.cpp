#include <iostream>
#include "dp54-sim.h"
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>

void dp54(std::vector<double> & pos, std::vector<double> & vel, const double tini, const double tend)
{
  std::vector<double> k1, k2, k3, k4, k5, k6, k7, velaux, posaux, vaux, paux;
  k1.resize(6);
  k2.resize(6);
  k3.resize(6);
  k4.resize(6);
  k5.resize(6);
  k6.resize(6);
  k7.resize(6);
  vaux.resize(vel.size());
  paux.resize(vel.size());
  velaux.resize(vel.size());
  posaux.resize(pos.size());
  const double eps=0.001;
  double dt=T/6000.0;
  double dt2=0;
  int nt=0;
  
  //calculo
  for (int tt=0;tt<13400;tt++)
    {
      double t = tini + dt*nt;

       std::copy(vel.begin(), vel.end(), velaux.begin());
       std::copy(pos.begin(), pos.end(), posaux.begin());

       // k1
       for(int ii = 0; ii < vel.size(); ++ii) {
	 k1[ii+3] = dt*compute(pos, vel,t, ii+1);
	 k1[ii] = dt*compute(pos, vel,t, ii+4);
       }
       // k2 aux
      for(int ii = 0; ii < vel.size(); ++ii) {
	vaux[ii] = vel[ii] + k1[ii+3]/5.0;
	paux[ii] = pos[ii] + k1[ii]/5.0;
      }
      //k2
      for(int ii = 0; ii < vel.size(); ++ii) {
	k2[ii+3] = dt*compute(paux,vaux,t+dt/5, ii+1);
	k2[ii] = dt*compute(paux,vaux,t+dt/5, ii+4);
      }
      // k3 aux
      for(int ii = 0; ii < vel.size(); ++ii) {
	vaux[ii] = vel[ii] + 3.0*k1[ii+3]/40.0 + 9.0*k2[ii+3]/40.0;
	paux[ii] = pos[ii] + 3.0*k1[ii]/40.0 + 9*k2[ii]/40.0;
      }
      //k3
      for(int ii = 0; ii < vel.size(); ++ii) {
	k3[ii+3] = dt*compute(paux,vaux,3*(t+dt)/10, ii+1);
	k3[ii] = dt*compute(paux,vaux,3*(t+dt)/10, ii+4);
      }
      // k4 aux
      for(int ii = 0; ii < vel.size(); ++ii) {
	vaux[ii] = vel[ii] + 44*k1[ii+3]/45 - 56*k2[ii+3]/15 + 32*k3[ii+3]/9;
	paux[ii] = pos[ii] + 44*k1[ii]/45 - 56*k2[ii]/15 + 32*k3[ii]/9;
      }
      //k4
      for(int ii = 0; ii < vel.size(); ++ii) {
	k4[ii+3] = dt*compute(paux,vaux,4*(t+dt)/5, ii+1);
	k4[ii] = dt*compute(paux,vaux,4*(t+dt)/5, ii+4);
      }
      // k5 aux
      for(int ii = 0; ii < vel.size(); ++ii) {
	vaux[ii] = vel[ii] + 38744.0*k1[ii+3]/13122 - 25360.0*k2[ii+3]/2187 + 64448.0*k3[ii+3]/6561 - 212.0*k4[ii+3]/729;
	paux[ii] = pos[ii] + 38744*k1[ii]/13122 - 25360*k2[ii]/2187 + 64448*k3[ii]/6561 - 212*k4[ii]/729;
      }
      //k5
      for(int ii = 0; ii < vel.size(); ++ii) {
	k5[ii+3] = dt*compute(paux,vaux,8*(t+dt)/9, ii+1);
	k5[ii] = dt*compute(paux,vaux,8*(t+dt)/9, ii+4);
      }      
      // k6 aux
      for(int ii = 0; ii < vel.size(); ++ii) {
	vaux[ii] = vel[ii] + 9017*k1[ii+3]/3168 - 355*k2[ii+1]/33 + 46732*k3[ii+3]/5247 + 49*k4[ii+3]/176 - 5103*k5[ii+3]/18656;
	paux[ii] = pos[ii] + 9017*k1[ii]/3168 - 355*k2[ii]/33 + 46732*k3[ii]/5247 + 49*k4[ii]/176 - 5103*k5[ii]/18656;
      }
      //k6
      for(int ii = 0; ii < vel.size(); ++ii) {
	k6[ii+3] = dt*compute(paux, vaux, t + dt, ii+1);
	k6[ii] = dt*compute(paux, vaux, t + dt, ii+4);
      }
      // k7 aux
      for(int ii = 0; ii < vel.size(); ++ii) {
	vaux[ii] = vel[ii] + 35*k1[ii+3]/384 + 500*k3[ii+3]/1113 + 125*k4[ii+3]/192 - 2187*k5[ii+3]/6784 + 11*k6[ii+3]/84;
	paux[ii] = pos[ii] + 35*k1[ii]/384 + 500*k3[ii]/1113 + 125*k4[ii]/192 - 2187*k5[ii]/6784 + 11*k6[ii]/84;
      }
      //k7
      for(int ii = 0; ii < vel.size(); ++ii) {
	k7[ii+3] = dt*compute(paux, vaux, t + dt, ii+1);
	k7[ii] = dt*compute(paux, vaux, t + dt, ii+4);
      }
      
      // new vel, pos
      for(int ii = 0; ii < vel.size(); ++ii) {
	vel[ii] = vel[ii] + (35*k1[ii+3]/384 + 0 + 500*k3[ii+3]/1113 + 125*k4[ii+3]/192 - 2187*k5[ii+3]/6884) + 11*k6[ii+3]/84;
	pos[ii] = pos[ii] + (35*k1[ii]/384 + 0 + 500*k3[ii]/1113 + 125*k4[ii]/192 - 2187*k5[ii]/6884) + 11*k6[ii]/84;	
      }
      // vel,pos aux
      for(int ii = 0; ii < vel.size(); ++ii) {
	velaux[ii] = velaux[ii] + (5179*k1[ii+3]/57600 + 0 + 7571*k3[ii+3]/16695 + 393*k4[ii+3]/640 - 92097*k5[ii+3]/339200) + 187*k6[ii+3]/2100 + k7[ii+3]/40;
      }
      for(int ii = 0; ii < vel.size(); ++ii) {
	posaux[ii] = posaux[ii] + (5179*k1[ii]/57600 + 0 + 7571*k3[ii]/16695 + 393*k4[ii]/640 - 92097*k5[ii]/339200) + 187*k6[ii]/2100 + k7[ii]/40;
      }
      //new dt
      dt2=dtnew(E(vel,velaux,pos,posaux),dt);
      dt=dt2;
      
      std::cout << t  << " ";
      print(pos);
      print(vel);
      std::cout << "\n";
      
      
      nt++;
    }
  std::cout << nt << "\n";
}  
  




  
double compute(const std::vector<double> & pos, const std::vector<double> & vel, const double t, const int id)
{
  double r = rvar(pos,vel);
  double s = svar(pos,vel);
  
  if(id==1){ //retorna funcion para calcular la velocidad en x
    return (pos[0] + 2*vel[1] - ((1-u)*(pos[0]+u)/(r*r*r)) - (u*(pos[0]-1+u)/(s*s*s)));
  }
  if(id==2){ //retorna funcion para calcular la velocidad en y
    return (pos[1] - 2*vel[0] - ((1-u)*pos[1]/(r*r*r)) - (u*pos[1]/(s*s*s)));
  }
  if(id==3){ //retorna funcion para calcular la velocidad en z
    return (((u-1)*pos[2]/(r*r*r)) - (u*pos[2]/(s*s*s)));
  }
  if(id==4){ //retorna la velocidad en x
    return vel[0];
  }
  if(id==5){ //retorna la velocidad en y
    return vel[1];
  }
  if(id==6){ //retorna la velocidad en z
    return vel[2];
  }
  else{ //advierte de un error
    std::cerr << "Error!!!!" << id << std::endl;
    exit(1);
  }
}

void initial_condition(std::vector<double> & pos, std::vector<double> & vel)
{
  //posiciones
  pos = {0.994, 0.0, 0.0};
  
  //velocidades
  vel = {0.0, -2.0015851063790825224053786224, 0.0};
}

void print(const std::vector<double> & v)
{
  for (int i = 0; i < v.size(); ++i)
    {
      std::cout << v[i] << " ";
    }
}

double rvar(const std::vector<double> & pos, const std::vector<double> & vel)
{
  return std::sqrt((pos[0]+u)*(pos[0]+u)+(pos[1]*pos[1])+(pos[2]*pos[2]));
}

double svar(const std::vector<double> & pos, const std::vector<double> & vel)
{
  return std::sqrt((pos[0]+u-1)*(pos[0]+u-1)+(pos[1]*pos[1])+(pos[2]*pos[2]));
}

 double dtnew(double e,double dt)
{
  return dt* std::min(5.0,std::max(0.2,(std::pow(0.38,1/5))*(std::pow(1/e,1/5))));
}
 
double E(const std::vector<double> & vel, const std::vector<double> & velaux,const std::vector<double> & pos, const std::vector<double> & posaux)//medida conjunta del error
{
  double sum=0;
  for(int ii=0;ii<vel.size();ii++)
    {
      sum+=((vel[ii]-velaux[ii])/((vel[ii]-velaux[ii])+(std::max(std::abs(vel[ii]),std::abs(velaux[ii])))))*((vel[ii]-velaux[ii])/((vel[ii]-velaux[ii])+(std::max(std::abs(vel[ii]),std::abs(velaux[ii]))))) + ((pos[ii]-posaux[ii])/((pos[ii]-posaux[ii])+(std::max(std::abs(pos[ii]),std::abs(posaux[ii])))))*((pos[ii]-posaux[ii])/((pos[ii]-posaux[ii])+(std::max(std::abs(pos[ii]),std::abs(posaux[ii]))))); //Ai con o sin valor absoluto?
      
    }
  return sqrt(sum/3);
}

