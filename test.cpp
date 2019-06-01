#include <qpOASES.hpp>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>
#include <algorithm> 

using namespace Eigen;
using namespace std;
USING_NAMESPACE_QPOASES

void reSizeMat(vector<int> *A, int prevDim, int newDim) {
    vector<int> B(newDim*newDim);
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

void printmat(vector<int> *H) {
    int n = sqrt(H->size());
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

void blockDiag(vector<int> *H, real_t *Hn, int HnRows) {
    int HnLen = HnRows*HnRows;
    int currDim = sqrt(H->size());
    int nblocks = currDim/HnRows;
    int newDim = HnRows*(nblocks+1);
    reSizeMat(H, currDim, newDim);

    // if(nblocks == 0) {
    //     copy(Hn, Hn+HnLen, H->begin());
    // }
    // else {
    for(int i=0;i<(newDim*newDim);i++) {
        int c = i%newDim;
        int r = i/newDim;
        if(r > currDim -1 || c > currDim -1) {
            if(r<currDim && c > currDim-1) {
                H->at(i) = 0;
            }
            else if(r>currDim-1 && c < currDim) {
                H->at(i) = 0;
            }
            else {
                int r_n = r%HnRows;
                int c_n = c%HnRows;
                H->at(i) = (Hn[r_n*HnRows + c_n]);
            // }
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

	float t1, t2, t3;
	t1 = 0;
	t2 = 5;
	t3 = 10;
	vector<float> tList;
	tList.push_back(t1);
	tList.push_back(t2);
	tList.push_back(t3);

	float max_vel = 4;
	float max_acc = 5;

	int n = 7; //number of coefficients (degree + 1) 
    int K = 1; //number of drones
    int M = posList.size() - 1; //number of splines
    int DIM = 3; //dimensions

    int nx = n*M*K*DIM; //number of decision variables
    int nConstraints = 4 + DIM*posList.size(); //init and final velocity + acceleration = 0

    // real_t H[nx*nx], A[1*nx], eqConstraints[nConstraints];

    real_t Hn1[n*n],Hn2[n*n],Hn3[n*n], Hd[n*n*DIM*DIM];
    fill(Hn1,Hn1+n*n,1);

    vector<int> Hv;  
    blockDiag(&Hv, Hn1, n);
    printmat(&Hv);
    blockDiag(&Hv, Hn1, n);
    printmat(&Hv);
    blockDiag(&Hv, Hn1, n);
    printmat(&Hv);

}
