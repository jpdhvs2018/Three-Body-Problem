#include "rk4-sim.h"
#include <vector>
#include <iostream>
int main()
{
  const double ti=0.0;
  const double te=17.0652165601579625588917206249;
  std::vector<double> position;
  std::vector<double> velocity;
  initial_condition(position, velocity);
  rk4(position, velocity, ti, te);
  return 0;
}
