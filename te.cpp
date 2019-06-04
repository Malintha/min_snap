#include <qpOASES.hpp>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>
#include <algorithm> 

using namespace Eigen;
using namespace std;
USING_NAMESPACE_QPOASES

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
    hblock(3,0) = 2*30*720*(pow(t1,4)-pow(t0,4));
    hblock(4,1) = 2*24*120*(pow(t1,2)-pow(t0,2));
    hblock(5,2) = 2*16*360*(pow(t1,3)-pow(t0,3));
    return hblock;
}

int main() {
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

	double max_vel = 4;
	double max_acc = 5;
	const int n = 7; //number of coefficients (degree + 1) 
    int K = 1; //number of drones
    const int M = 1; //number of splines
    const int D = 3; //dimensions
    int nx = n*M*K*D; //number of decision variables
    const int wpts = posList.size();
// construct Hessian
    MatrixXf H(K*M*D*n, K*M*D*n);
    for(int k=0;k<K;k++) {
        MatrixXf mstacked(M*D*n, M*D*n);
        for(int m=0;m<M;m++) {
            double t0 = tList[m];
            double t1 = tList[m+1];
            MatrixXf dstacked(D*n,D*n);
            dstacked.fill(0);
            for(int d=0;d<D;d++) {
                MatrixXf hblock = getHblock(t0, t1);
                dstacked.block<7,7>(d*7,d*7) = hblock;
            }
            mstacked.block<7*3,7*3>(7*3*m,7*3*m) = dstacked;
        }
        H.block<M*D*7,M*D*7>(n*D*M*k,n*D*M*k) = mstacked;
    }

//construct A and constraint matrices
// TODO: change the first param of rows allocator of block() to const
    const int constraints = (posList.size() + 4)*D*K;
    MatrixXf A(constraints, nx);
    vector<double> lb, ub;
    for(int k=0;k<K;k++) {
        MatrixXf dstacked(constraints/K, n*D);
        dstacked.fill(0);
        for(int d=0;d<D;d++) {
            int dConstraints = constraints/(K*D);
            MatrixXf mstacked(dConstraints, n);
            mstacked.fill(0);
            for(int m=0;m<wpts;m++) {
                double t = tList[m];
                int mdConstraints;
                m==0 || m==wpts-1 ? mdConstraints = 3: mdConstraints = 1;
                MatrixXf mdblock(mdConstraints, n);
                mdblock.fill(0);
                //position row
                MatrixXf tPos = getPosTimeVec(t);
                mdblock.block<1,7>(0,0) = tPos;
                if(m==0 || m==posList.size()-1) {
                    //add velocity row
                    MatrixXf tVel = getVelTimeVec(t);
                    mdblock.block<1,7>(1,0) = tVel;
                    //acceleration row
                    MatrixXf tAcc = getAccTimeVec(t);
                    mdblock.block<1,7>(2,0) = tAcc;
                }
                int rowIdx;
                m == 0? rowIdx = 0 : rowIdx = m+2;
                mstacked.block<3, 7>(rowIdx, 0) = mdblock;
            }
            int dColIdx = n*d;
            int dRowIdx = dConstraints*d;
            dstacked.block<7,n>(dRowIdx, dColIdx) = mstacked;
        }
        int kRowIdx = dstacked.rows();
        int kColIdx = dstacked.cols();
        A.block<21,n*D>(k*kRowIdx, k*kColIdx) = dstacked;
    }
    cout<<"A: "<<A.cols()<<" "<<A.rows()<<endl<<A<<endl;


//   cout << "Block in the middle" << endl;
//   cout << m1.block<2,2>(1,1) << endl << endl;
//   for (int i = 1; i <= 3; ++i)
//   {
//     cout << "Block of size " << i << "x" << i << endl;
//     cout << m1.block(0,0,i,i) << endl << endl;
//   }

    // MatrixXf dstacked(14,14);
    // dstacked.fill(0);
    // MatrixXf hblock1 = getHblock(0, 1);
    // MatrixXf hblock2 = getHblock(1,2);
    // dstacked.block<7,7>(0,0) = hblock1;
    // dstacked.block<7,7>(7,7) = hblock2;
    // cout<<dstacked<<endl;
}