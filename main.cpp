// file main.cpp

#include "ODDtool.h"

using namespace std;

// user-supplied history function 
void myHistory(double t, vector<double>& result) { 
   result[0] = abs(cos(t)) + 1; 
   result[1] = 1 + sin(t); 
}

// Two coupled Mackey-Glass systems
void MG2Rhs(double t, const vector<double>& x, ODD_delayed_values& xd, const vector<double>& p, vector<double>& result) {
   result[0] = p[0] * xd[0][1] / (1 + pow(xd[0][1], p[1])) - p[2] * x[0];
   result[1] = p[0] * xd[0][0] / (1 + pow(xd[0][0], p[1])) - p[2] * x[1];
}

int main() { 
   // create an integrator 
   ODD_integrator<CKStepper> integrator(MG2Rhs, myHistory);

   // start integration over time interval specified in parameter file (ODD_parameters.txt) 
   integrator.integrate(); 
}
