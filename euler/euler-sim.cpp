#include <iostream>
#include "euler-sim.h"
#include <cmath>
#include <cstdlib>
#include <vector>

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

double compute(const std::vector<double> & pos, const std::vector<double> & vel, const double t, const int id)
{
  double r = rvar(pos,vel);
  double s = svar(pos,vel);
  
  if(id==1){ //returna funcion para calcular la velocidad en x
    return (pos[0] + 2*vel[1] - ((1-u)*(pos[0]+u)/(r*r*r)) - (u*(pos[0]-1+u)/(s*s*s)));
  }
  if(id==2){ //returna funcion para calcular la velocidad en y
    return (pos[1] - 2*vel[0] - ((1-u)*pos[1]/(r*r*r)) - (u*pos[1]/(s*s*s)));
  }
  if(id==3){ //returna funcion para calcular la velocidad en z
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

void euler(std::vector<double> & pos, std::vector<double> & vel, const double tini, const double tend)
{
  std::vector<double> posaux(pos.size());
  std::vector<double> velaux(vel.size());
  double dt = T/(24000.0);
  int N = (int)((tend-tini)/dt);
  for (int i = 0; i < N; ++i)
    {
      double t = tini +  dt*i;
      std::copy(vel.begin(), vel.end(), velaux.begin());
      std::copy(pos.begin(), pos.end(), posaux.begin()); //extra
      for (int j = 0; j < (vel.size()) ; ++j)
	{
	  vel[j] = vel[j] + dt*compute(posaux, velaux, t, j+1); 
	}
      
      for (int k = 0; k < (pos.size()) ; ++k)
	{
	  pos[k] = pos[k] + dt*compute(posaux, velaux, t, k+4);
	}
      
      std::cout << t << " ";
      print(pos);
      print(vel);
      std::cout << "\n";
    }
}



