#include <iostream>
#include <fstream>
#include <cmath>
#include "JI.h"

using namespace std;

int main() {

    int twist = 3;

    JI ji_Instance = JI();
    //Kinematics is QQ, then Bjorken x, then momentum transfer t, then E beam energy
    double kinematics[4] = {1.82 /*QQ*/, 0.34 /*x_b*/, -0.17 /*t*/, 5.75 /*E_b*/};
    ji_Instance.SetKinematics(kinematics);
    ji_Instance.setAlphaBeta(0.5, 0.5);

    ofstream output;
    output.open("Fig4_hU_tilde_Testing_Twist3.txt");

    for(double phi = 0; phi < 361; phi++){
        ji_Instance.setValuesOf_h(phi, twist);
        output << phi << " " << ji_Instance.get_h_U_tilde() << endl;
    }
    cout << "test" << endl;

}