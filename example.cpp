#include <qpOASES.hpp>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>
#include <algorithm>
#include <yaml.h>


using namespace std;
using namespace Eigen;

/** Example for qpOASES main function using the SQProblem class. */
int main( )
{
	USING_NAMESPACE_QPOASES

	/* Setup data of first QP. */
	// real_t H[3*3] = { 1.0, 0.0, 0.0, 0.5, 0, 0, 1, 0, 0 };
	// real_t A[1*3] = { 1.0, 1.0 ,1};
	// real_t g[2] = { 1.5, 1.0 };
	// real_t lb[3] = { 0.5, -2.0 ,-5};
	// real_t ub[3] = { 5.0, 2.0 ,5};
	// real_t lbA[1] = { -1.0 };
	// real_t ubA[1] = { 2.0 };

    MatrixXf H(7,7);
    H.fill(0);
    MatrixXf H_topleft(3,3);
    H_topleft <<    1.62e+08, 1.35e+07,   720000,
                    1.35e+07,  1.2e+06,    72000, 
                    720000,    72000,     5760;
    H.block<3,3>(0,0) = H_topleft;
    cout<<"H:"<<endl<<H<<endl;
    real_t H_r[H.size()];
    float* hp = H.data();
        for(int i=0;i<H.size();i++) {
        H_r[i] = *hp++;
    }
    
    MatrixXf A(7,6);
    A <<        0,     0,     0, 15625, 18750, 18750,
                0,     0,     0,  3125,  3125,  2500,
                0,     0,     0,   625,   500,   300,
                0,     0,     0,   125,    75,    30,
                0,     0,     2 ,   25,    10,     2,
                0,     1 ,    0,     5,     1,     0,
                1,     0,     0,     1,     0,     0;
    real_t A_r[A.size()];
    float* ap = A.data();
    for(int i=0;i<A.size();i++) {
        A_r[i] = *ap++;
    }

    real_t lba[] = {5,0,0,10,0,0};
    real_t uba[] = {5,0,0,10,0,0};
    real_t g[] = {0,0,0,0,0,0,0};

	/* Setup data of second QP. */
	// real_t H_new[2*2] = { 1.0, 0.5, 0.5, 0.5 };
	// real_t A_new[1*2] = { 1.0, 5.0 };
	// real_t g_new[2] = { 1.0, 1.5 };
	// real_t lb_new[2] = { 0.0, -1.0 };
	// real_t ub_new[2] = { 5.0, -0.5 };
	// real_t lbA_new[1] = { -2.0 };
	// real_t ubA_new[1] = { 1.0 };

	int_t nWSR = 100;

	/* Setting up SQProblem object. */
	QProblem example(7,5);
	returnValue rv = example.init( H_r,g,A_r,nullptr,nullptr,lba,uba, nWSR,0 );
    cout<<"rv: "<<rv<<endl;
	real_t xOpt[7];
	Options options;
    options.setToMPC();
    example.setOptions(options);
	example.getPrimalSolution(xOpt);
    for(int i=1;i<=7;i++) {
        cout<<xOpt[i-1]<<" ";
        if(i%7 == 0)
            cout<<endl;
    }

//repeat

	/* Solve second QP. */
	QProblem example1(7, 5);
	Options options1;
	options1.setToMPC();
    example1.setOptions(options1);	
	// nWSR = 10;
	returnValue rv1 = example1.init( H_r,g,A_r,nullptr,nullptr,lba,uba, nWSR,0 );
    cout<<"rv: "<<rv1<<endl;

	/* Get and print solution of second QP. */
	real_t xOpt1[7];
	example1.getPrimalSolution( xOpt1 );
    for(int i=1;i<=7;i++) {
        cout<<xOpt1[i-1]<<" ";
        if(i%7 == 0)
            cout<<endl;
    }

	return 0;
}


/*
 *	end of file
 */
