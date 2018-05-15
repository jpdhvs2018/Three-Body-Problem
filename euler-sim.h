include  "particle.h"

//function declaration
void initial_condition(particle ball);

void euler(const double step, const particle ball, const double tini, const double tend); // imprime los valores en un intervalo

void compute(const particle ball, const double t, const double id); // returna valores para la posicion, velocidad, aceleracion

double rvar(const particle ball);

double svar(const particle ball);
