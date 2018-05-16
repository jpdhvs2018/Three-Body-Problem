#include "euler-sim.h"
#include <vector>
#include <iostream>

int main()
{
  const double ti=0.0;
  const double te=1000.0;
  std::vector<double> position;
  std::vector<double> velocity;
  initial_condition(position, velocity);
  euler(position, velocity, ti, te);
  return 0;
}
