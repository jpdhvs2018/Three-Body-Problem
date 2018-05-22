#include "dp54-sim.h"
#include <vector>
#include <iostream>
int main()
{
  const double ti=0.0;
  const double te=T;//17.0652165601579625588917206249;
  std::vector<double> position;
  std::vector<double> velocity;
  initial_condition(position, velocity);
  //dormand-prince 5(4)(position, velocity, ti, te);
  dp54(position, velocity, ti, te);
  return 0;
}
