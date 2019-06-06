#include <qpOASES.hpp>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>
#include <algorithm> 

using namespace Eigen;
using namespace std;

MatrixXf getPosTimeVec(double t) {
    MatrixXf posTimes(1,7);
    posTimes << pow(t, 6), pow(t, 5), pow(t, 4), pow(t, 3), pow(t, 2), pow(t, 1), 1;
    return posTimes;
}

MatrixXf getVelTimeVec(double t) {
    MatrixXf velTimes(1,7);
    velTimes << 6*pow(t, 5), 5*pow(t, 4), 4*pow(t, 3), 3*pow(t, 2), 2*t, 1, 0;
    return velTimes; 
}

MatrixXf getAccTimeVec(double t) {
    MatrixXf accTimes(1,7);
    accTimes << 30*pow(t, 4), 20*pow(t, 3), 12*pow(t, 2), 6*t, 2, 0, 0;
    return accTimes;
}

MatrixXf getHblock(double t0, double t1) {
    int n = 7;
    MatrixXf hblock(n,n);
    hblock.fill(0);
    hblock(0,0) = 2*72*360*(pow(t1,5)-pow(t0,5));
    hblock(1,1) = 2*40*120*(pow(t1,3)-pow(t0,3));
    hblock(2,2) = 2*576*(t1-t0);
    hblock(1,0) = 30*720*(pow(t1,4)-pow(t0,4));
    hblock(0,1) = 30*720*(pow(t1,4)-pow(t0,4));
    hblock(2,1) = 24*120*(pow(t1,2)-pow(t0,2));
    hblock(1,2) = 24*120*(pow(t1,2)-pow(t0,2));
    hblock(0,2) = 16*360*(pow(t1,3)-pow(t0,3));
    hblock(2,0) = 16*360*(pow(t1,3)-pow(t0,3));
    return hblock;
}

int main() {
USING_NAMESPACE_QPOASES

	vector<Vector3d> posList;
	Vector3d p1, p2, p3;
	p1 << 0,0,2;
	p2 << 5,5,2;
	p3 << 0,10,2;
	posList.push_back(p1);
	posList.push_back(p2);
	posList.push_back(p3);

	double t1, t2, t3;
	t1 = 0;
	t2 = 5;
	t3 = 10;
	vector<double> tList;
	tList.push_back(t1);
	tList.push_back(t2);
	tList.push_back(t3);


	const int n = 7; //number of coefficients (degree + 1) 
    int K = 1; //number of drones
    const int M = 1; //number of splines
    const int D = 1; //dimensions
    int nx = n*M*K*D; //number of decision variables
    const int wpts = posList.size();

// construct Hessian
    MatrixXf H(K*M*D*n, K*M*D*n);
    for(int k=0;k<K;k++) {
        MatrixXf mstacked(M*D*n, M*D*n);
        for(int m=0;m<1;m++) {
            double t0 = tList[0];
            double t1 = tList[2];
            MatrixXf dstacked(D*n,D*n);
            dstacked.fill(0);
            for(int d=0;d<D;d++) {
                MatrixXf hblock = getHblock(t0, t1);
                dstacked.block<7,7>(d*7,d*7) = hblock;
            }
            mstacked.block<7*D,7*D>(7*D*m,7*D*m) = dstacked;
        }
        H.block<M*D*7,M*D*7>(n*D*M*k,n*D*M*k) = mstacked;
    }

    real_t H_r[H.size()];
    MatrixXf H_t = H.transpose();
    cout<<"H: "<<H.size()<<" "<<H.rows()<<" "<<H.cols()<<endl<<H<<endl;
    float* hp = H_t.data();
    for(int i=0;i<H.size();i++) {
        H_r[i] = *hp++;
    }

    MatrixXf A(7,7);
    MatrixXf tPos = getPosTimeVec(0);
    A.block<1,7>(0,0) = tPos;
    MatrixXf tVel = getVelTimeVec(0);
    A.block<1,7>(1,0) = tVel;
    MatrixXf tAcc = getAccTimeVec(0);
    A.block<1,7>(2,0) = tAcc;
    tPos = getPosTimeVec(5);
    A.block<1,7>(3,0) = tPos;
    tPos = getPosTimeVec(10);
    A.block<1,7>(4,0) = tPos;
    tVel = getVelTimeVec(10);
    A.block<1,7>(5,0) = tVel;
    tAcc = getAccTimeVec(10);
    A.block<1,7>(6,0) = tAcc;

    real_t A_r[A.size()];
    MatrixXf A_t = A.transpose();
    cout<<"A: "<<A.size()<<" "<<A.rows()<<" "<<A.cols()<<endl;
    float* ap = A_t.data();
    for(int i=0;i<A.size();i++) {
        A_r[i] = *ap++;
    }

    real_t g[7];
    fill(g, g+7, 0);

    real_t lba[] = {0,0,0,5,10,0,0};
    real_t uba[] = {0,0,0,5,10,0,0};
    real_t ub[7], lb[7];
    fill(lb,lb+7,-50);
    fill(ub,ub+7,50);

	/* Setting up QProblem object. */
	QProblem example(7,7);

	Options options;
	example.setOptions(options);

	int_t nWSR = 100;
	example.init( H_r,g,A_r,lb,ub,lba,uba, nWSR );

	real_t xOpt[7];
	example.getPrimalSolution(xOpt);
    for(int i=0;i<7;i++) {
        cout<<xOpt[i]<<endl;
    }
	return 0;


}