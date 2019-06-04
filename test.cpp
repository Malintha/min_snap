#include <qpOASES.hpp>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>
#include <algorithm> 

using namespace Eigen;
using namespace std;
USING_NAMESPACE_QPOASES

vector<double> getPosTimeVec(double t) {
    vector<double> d_vec(7);
    double data[7] = {pow(t, 6), pow(t, 5), pow(t, 4), pow(t, 3), pow(t, 2), pow(t, 1), 1};
    for(int i=0;i<7;i++) {
        d_vec[i] = data[i];
    }
    return d_vec;
}

vector<double> getVelTimeVec(double t) {
    vector<double> d_vec(7);
    double data[7] = {6*pow(t, 5), 5*pow(t, 4), 4*pow(t, 3), 3*pow(t, 2), 2*t, 1, 0};
    for(int i=0;i<7;i++) {
        d_vec[i] = data[i];
    }
    return d_vec; 
}

vector<double> getAccTimeVec(double t) {
    double data[7] = {30*pow(t, 4), 20*pow(t, 3), 12*pow(t, 2), 6*t, 2, 0, 0};
    vector<double> d_vec(7);
    for(int i=0;i<7;i++) {
        d_vec[i] = data[i];
    }
    return d_vec;
}

void printmat(vector<double> *H) {
    int n = sqrt(H->size());
    cout<<"Dim: "<<n<<"x"<<n<<endl;
    for(int i=0;i<H->size();i++) {
        int c = i%n;
        int r = i/n;
        cout<<H->at(i)<<" ";
        if(c==n-1) {
            cout<<endl;
        }
    }
    cout<<endl;
}

void printmat(vector<double> *H, int cols) {
    int rows = H->size()/cols;
    cout<<"Dim: "<<rows<<"x"<<cols<<endl;

    for(int i=0;i<H->size();i++) {
        int c = i%cols;
        int r = i/rows;
        cout<<H->at(i)<<" ";
        if(c==cols-1) {
            cout<<endl;
        }
    }
    cout<<endl;
}

vector<double> getHblock(double t0, double t1) {
    int n = 7;
    vector<double> hblock_(n*n);
    fill(hblock_.begin(), hblock_.begin()+n*n-1, 0);
    for(int i=0;i<n*n;i++) {
        int r = i/n;
        int c = i%n;
        if(r==0 && c==0) {
            hblock_[i] = 2*72*360*(pow(t1,5)-pow(t0,5));
        }
        else if(r==1 && c==1) {
            hblock_[i] = 2*40*120*(pow(t1,3)-pow(t0,3));
        }
        else if(r == 2 && c==2) {
            hblock_[i] = 2*576*(t1-t0);
        }
        else if(r==3 && c==0) {
            hblock_[i] = 2*30*720*(pow(t1,4)-pow(t0,4));
        }
        else if(r==4 && c==1) {
            hblock_[i] = 2*24*120*(pow(t1,2)-pow(t0,2));
        }
        else if(r==5 && c==2) {
            hblock_[i] = 2*16*360*(pow(t1,3)-pow(t0,3));
        }
    }
    return hblock_;
}

void reSizeMat(vector<double> *A, int prevDim, int newDim) {
    vector<double> B(newDim*newDim);
    fill(B.begin(),B.begin() + newDim*newDim, 0);
    for(int i=0;i<newDim*newDim;i++) {
        int r = i/newDim;  
        int c = i%newDim;
        if(r < prevDim && c < prevDim) {
            B[i] = A->at(r*prevDim + c);
        }
    }
    *A = B;
}

void reSizeMat(vector<double> *A, int prevCols, int newRows, int newCols) {
    int prevRows;
    prevCols == 0 ? prevRows = 0 : prevRows = A->size()/prevCols;
    vector<double> B(newRows*newCols);
    fill(B.begin(),B.begin() + newCols*newRows, 0);
    for(int i=0;i<newCols*newRows;i++) {
        int r = i/newCols;
        int c = i%newCols;    
        if(r < prevRows && c < prevCols) {
            B[i] = A->at(r*prevCols + c);
        }
    }
    *A = B;
}

void blockDiag(vector<double> *H, int HCols, vector<double> Hn, int HnCols) {
    int HnRows = Hn.size()/HnCols;
    int HRows;
    HCols == 0 ? HRows = 0 : HRows = H->size()/HCols;
    // int currDim = sqrt(H->size());
    int nblocks = HRows/HnRows;
    int newRows = HnRows*(nblocks+1);
    int newCols = HnCols*(nblocks+1);
    reSizeMat(H, HCols, newRows, newCols);
    // printmat(H, newCols);
    double* fp = H->data();
    for(int i=0;i<(newRows*newCols);i++) {
        int c = i%newCols;
        int r = i/newCols;
        if(r > HRows -1 || c > HCols -1) {
            if(r<HRows && c > HCols-1) {
                fp[i] = 0;
            }
            else if(r>HRows-1 && c < HCols) {
                fp[i] = 0;
            }
            else {
                int r_n = r%HnRows;
                int c_n = c%HnRows;
                fp[i] = (Hn[r_n*HnRows + c_n]);
            }
        }
    }
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
	int n = 7; //number of coefficients (degree + 1) 
    int K = 1; //number of drones
    int M = 1; //number of splines
    int D = 3; //dimensions
    int nx = n*M*K*D; //number of decision variables
    int wpts = posList.size();
//construct Hessian
    vector<double> kstacked;
    for(int k=0;k<K;k++) {
        vector<double> mstacked;
        for(int m=0;m<M;m++) {
            double t0 = tList[m];
            double t1 = tList[m+1];
            vector<double> dstacked;
            for(int d=0;d<D;d++) {
                vector<double> hblock = getHblock(t0, t1);
                blockDiag(&dstacked, d*n, hblock, n);
            }
        blockDiag(&mstacked, m*D*n, dstacked, D*n);
        }
    blockDiag(&kstacked,k*M*D*n, mstacked, M*D*n);
    }
    // cout<<"H: "<<endl;
    // printmat(&kstacked);

//construct A and constraint matrices
    vector<double> A;
    vector<double> lb, ub;
    for(int k=0;k<K;k++) {
        vector<double> mstacked;
        for(int m=0;m<wpts;m++) {
            double t = tList[m];
            vector<double> dstacked;
            for(int d=0;d<D;d++) {
                vector<double> dblock;
                //position row
                vector<double> tPos = getPosTimeVec(t);
                for(int i=0;i<tPos.size();i++) {
                    dblock.push_back(tPos[i]);
                }
                if(m==0 || m==posList.size()-1) {
                    //add velocity row
                    vector<double> tVel = getVelTimeVec(t);
                    for(int i=0;i<tVel.size();i++) {
                        dblock.push_back(tVel[i]);
                    }
                    //acceleration row
                    vector<double> tAcc = getAccTimeVec(t);
                    for(int i=0;i<tAcc.size();i++) {
                        dblock.push_back(tAcc[i]);
                    }
                }
                blockDiag(&dstacked,n*D,dblock,n);
            }
            blockDiag(&mstacked, m*D*n, dstacked, D*n);
        }
        // blockDiag(&A, k*wpts*D*n,mstacked, wpts*D*n);
    }
    // printmat(&A, K*wpts*D*n);

    // vector<double> test(10);
    // fill(test.begin(),test.begin()+10,1);
    // vector<double> test1(10);
    // fill(test1.begin(),test1.begin()+10,2);
    // printmat(&test, 5);
    // printmat(&test1, 5);
    // blockDiag(&test, 5, test1, 5);
    // printmat(&test, 10);

}
