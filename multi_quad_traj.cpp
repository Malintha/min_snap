#include <qpOASES.hpp>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>
#include <algorithm>
#include <yaml.h>

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
	Vector3d p1, p2, p3, p4;
	p1 << 5,0,2.5;
	p2 << 10,5,2.5;
	p3 << 5,10,2.5;
	posList.push_back(p1);
	posList.push_back(p2);
	// posList.push_back(p3);

	double t1, t2, t3, t4;
	t1 = 0;
	t2 = 5;
	t3 = 10;

	vector<double> tList;
	tList.push_back(t1);
	tList.push_back(t2);
	// tList.push_back(t3);

	double max_vel = 4; //maximum velocity m/s
	double max_acc = 5; //maximum acceleration m/s/s
	const int n = 7; //number of coefficients (degree + 1) 
    const int K = 1; //number of drones
    const int M = 1; //number of splines
    const int D = 1; //dimensions
    int nx = n*K*D; //number of decision variables
    const int nwpts = posList.size();
    int nChecks = 0; //number of interim checks for velocity and acceleration
    const int nc = (posList.size() + 4 + 2*nChecks*(posList.size()-1))*D*K;

    MatrixXf A(nc, nx);
    vector<double> lba;
    vector<double> uba;
    cout<<"nx: "<<nx<<" nc: "<<nc<<endl;

// construct Hessian
    MatrixXf H(K*M*D*n, K*M*D*n);
    for(int k=0;k<K;k++) {
        MatrixXf mstacked(M*D*n, M*D*n);
        for(int m=0;m<M;m++) {
            double t0 = tList[0];
            double t1 = tList[tList.size()-1];
            MatrixXf dstacked(D*n,D*n);
            dstacked.fill(0);
            for(int d=0;d<D;d++) {
                MatrixXf hblock = getHblock(t0, t1);
                dstacked.block<n,n>(d*7,d*7) = hblock;
            }
            mstacked.block<7*D,7*D>(7*D*m,7*D*m) = dstacked;
        }
        H.block<M*D*7,M*D*7>(n*D*M*k,n*D*M*k) = mstacked;
    }

//construct A and constraint matrices
    for(int k=0;k<K;k++) {
        MatrixXf dstacked(nc/(K), n*D);
        dstacked.fill(0);
        for(int d=0;d<D;d++) {
            int dConstraints = nc/(K*D);
            MatrixXf mstacked(dConstraints, n);
            mstacked.fill(0);
            vector<double> lba_d;
            vector<double> uba_d;
            for(int m=0;m<nwpts;m++) {
                vector<double> lba_m;
                vector<double> uba_m;
                int ccount = 0;
                double t = tList[m];
                int mdConstraints;
                m==0 || m== nwpts-1 ? mdConstraints = 3: mdConstraints = 1;
                if(nwpts > 0) {
                    m < nwpts - 1 ? mdConstraints += 2*nChecks : mdConstraints = mdConstraints;
                }
                MatrixXf mdblock(mdConstraints, n);
                mdblock.fill(0);

                //position equality
                MatrixXf tPos = getPosTimeVec(t);
                mdblock.block<1,n>(ccount++,0) = tPos;

                lba.push_back(posList[m][d]);
                uba.push_back(posList[m][d]);
                if(m==0 || m==posList.size()-1) {
                    //add velocity equality
                    MatrixXf tVel = getVelTimeVec(t);
                    mdblock.block<1,n>(ccount++,0) = tVel;
                    lba.push_back(0);
                    uba.push_back(0);
                    //acceleration equality
                    MatrixXf tAcc = getAccTimeVec(t);
                    mdblock.block<1,n>(ccount++,0) = tAcc;
                    lba.push_back(0);
                    uba.push_back(0);
                }
                //add interim velocity and acceleration limits
                if(m < nwpts-1) {
                    double t1 = tList[m+1];
                    for(int i=1;i <= nChecks;i++) {
                        double tCheck = t + ((double)i/(nChecks+1))*(t1 - t);
                        MatrixXf tVel = getVelTimeVec(tCheck);
                        mdblock.block<1,n>(ccount++, 0) = tVel;                        
                        lba.push_back(-max_vel);
                        uba.push_back(max_vel);

                        MatrixXf tAcc = getAccTimeVec(tCheck);
                        mdblock.block<1,n>(ccount++,0) = tAcc;

                        lba.push_back(-max_acc);
                        uba.push_back(max_acc);
                    }
                }          

                int rowIdx;
                int mcrows = mdblock.rows();
                m == 0? rowIdx = 0 : rowIdx = 2+m*(2*nChecks + 1);
                // cout<<"mdblock"<<endl<<mdblock<<endl;
                mstacked.block(rowIdx, 0, mcrows, n) = mdblock;
                // cout<<"mstacked"<<endl<<mstacked<<endl;
            }
            int dColIdx = n*d;
            int dRowIdx = dConstraints*d;
            int dcrows = mstacked.rows();
            int dccols = mstacked.cols();
            dstacked.block(dcrows*d, dccols*d, dcrows, dccols) = mstacked;
        }

        cout<<"dStacked: "<<dstacked.rows()<<" "<<dstacked.cols()<<endl;
        int kRowIdx = dstacked.rows()*k;
        int kColIdx = dstacked.cols()*k;
        A.block(kRowIdx, kColIdx,nc/K,n*D) = dstacked;
        // A = dstacked;
    }

    real_t lb_r[lba.size()];
    real_t ub_r[uba.size()];
    cout<<"lba"<<endl;
    for(int i=0;i<lba.size();i++) {
        cout<<lba[i]<<endl;
    }
    cout<<"uba: "<<endl;
    for(int i=0;i<uba.size();i++) {
        cout<<lba[i]<<endl;
    }

    copy(lba.begin(), lba.begin()+nc,lb_r);
    copy(uba.begin(), uba.begin()+nc,ub_r);

    real_t A_r[A.size()];
    MatrixXf A_t = A.transpose();
    float* ap = A_t.data();
    for(int i=0;i<A.size();i++) {
        A_r[i] = *ap++;
    }

    real_t H_r[H.size()];
    MatrixXf H_t = H.transpose();
    float* hp = H_t.data();
    for(int i=0;i<H.size();i++) {
        H_r[i] = *hp++;
    }
    
    cout<<"lba: "<<lba.size()<<endl<<"uba: "<<uba.size()<<endl<<endl;

    real_t g[] = {0,0,0,0,0,0,0};
	QProblem qp(7,6);
	int_t nWSR = 100;
	qp.init( H_r,g,A_r,nullptr,nullptr,lb_r,ub_r, nWSR,0 );
	real_t xOpt[nx];
    Options options;
    options.setToMPC();
	qp.setOptions( options );
	qp.getPrimalSolution(xOpt);
    for(int i=1;i<=nx;i++) {
        cout<<xOpt[i-1]<<" ";
        if(i%7 == 0)
            cout<<endl;
    }
    
    // cout<<"H: "<<H.rows()<<" "<<H.cols()<<endl<<H_t<<endl;;
    // cout<<"A: "<<A.rows()<<" "<<A.cols()<<endl<<A_t<<endl;

//repeat
    // real_t lba_1[] = {7,1,0,12,2,0};
    // real_t uba_1[] = {7,1,0,12,2,0};
	QProblem qp1(7,6);
	int_t nWSR1 = 100;

	qp1.init( H_r,g,A_r,nullptr,nullptr,lb_r,ub_r, nWSR1,0);
    Options options1;
    options1.setToMPC();
	qp1.setOptions( options1 );
	real_t xOpt1[nx];
	qp1.getPrimalSolution( xOpt1 );

    for(int i=1;i<=nx;i++) {
        cout<<xOpt1[i-1]<<" ";
        if(i%7 == 0)
            cout<<endl;
    }

	return 0;
}