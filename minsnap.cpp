#include <qpOASES.hpp>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;
USING_NAMESPACE_QPOASES

void blockDiag(real_t *A, real_t *B) {
	//implement the block diagonal matrix concatenation of 2 arrays 
}

void set_H(real_t *H, int n, int K, int M, int DIM, vector<float> tList) {
	real_t Hk[(n*n)*DIM*M];
	for(int k=0;k<K;k++) {
		real_t Hm[n*DIM*M];
		for(int m=0;m<M;m++) {
			float t0 = tList[m];
			float t1 = tList[m+1];
			real_t Hd[n*n*DIM];
			for(int d=0;d<DIM;d++) {
				real_t Hn[n*n];
				for(int i=0;i<n;i++) {
					int r = i/n;
					int c = i%n;
					if(r==0 && c==0) {
						Hn[i] = 72*360*(pow(t1,5) - pow(t0,5));
					}
					else if(r==1 && c ==1) {
						Hn[i] = 40*120*(pow(t1,3) - pow(t0,3));
					}
					else if(r==2 && c==2) {
						Hn[i] = 576*(t1-t0);
					}
					else if(r=3 && c ==0) {
						Hn[i] = 30*720*(pow(t1,4) - pow(t0,4));
					}
					else if(r==4 && c==1) {
						Hn[i] = 24*120*(pow(t1,2) - pow(t0,2));
					}
					else if(r==5 && c==2) {
						Hn[i] = 16*360*(pow(t1,3) - pow(t0,3));
					}
					else {
						Hn[i] = 0;
					}
				}
			blockDiag(Hd, Hn);
			}
		}


	}
}

void set_g() {

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

    real_t H[nx*nx], A[1*nx], eqConstraints[nConstraints];

	set_H(H, n, K, M, DIM);


}
