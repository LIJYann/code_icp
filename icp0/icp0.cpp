// icp0.cpp : définit le point d'entrée pour l'application console.
//

#include "stdafx.h"
#include<algorithm>
#include<vector>
#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Eigenvalues>

using namespace Eigen;
using namespace std;

typedef struct Point3D
{
	float x,y,z;
};

typedef struct Rotation
{
	float x1,x2,x3,y1,y2,y3,z1,z2,z3;
};

float sq(float x) {return x*x;}

void CalculateMeanPoint3D(vector<Point3D> &P, Point3D &mean)
{
	vector<Point3D>::iterator it;
	mean.x = 0;
	mean.y = 0;
	mean.z = 0;
	for(it=P.begin(); it!=P.end(); it++)
	{
		mean.x += it->x;
		mean.y += it->y;
		mean.z += it->z;
	}
	mean.x = mean.x/P.size();
	mean.y = mean.y/P.size();
	mean.z = mean.z/P.size();
}

void MovingPointSet(vector<Point3D> &P, Point3D &T)
{
	vector<Point3D>::iterator it;
	for(it=P.begin(); it!=P.end(); it++)
	{
		it->x -= T.x;
		it->y -= T.y;
		it->z -= T.z;
	}
}

void FindCorrespondingPoint(vector<Point3D> &P, vector<Point3D> &Q, vector<Point3D> &X, float &E){
	if (P.size()>Q.size()) {FindCorrespondingPoint(Q, P, X, E);return;}
	vector<Point3D>::iterator itp, itq;
	E=0;
	for (itp=P.begin(); itp!=P.end(); itp++){
		float mindist = -1;
		int num,count=0;
		for (itq=Q.begin(); itq!=Q.end(); itq++){
			count++;
			float dist=sq(itq->x-itp->x)+sq(itq->y-itp->y)+sq(itq->z-itp->z);
			if (mindist<0 || mindist>dist) {mindist=dist;num=count;}
		}
		E += mindist/P.size();
		X.push_back(*itp);
	}
}

void CalculateRotation(vector<Point3D> &P, vector<Point3D> &X, Rotation &R){
	//build covariance matrix between P and X
	float cov[9]={0,0,0,0,0,0,0,0,0};
	vector<Point3D>::iterator itx,itp;	for (itp=P.begin(),itx=X.begin();itp!=P.end();itp++,itx++){
		cov[0] += (itp->x)*(itx->x)/P.size();//(1,1)
		cov[1] += (itp->x)*(itx->y)/P.size();//(1,2)
		cov[2] += (itp->x)*(itx->z)/P.size();//(1,3) 
		cov[3] += (itp->y)*(itx->x)/P.size();//(2,1)
		cov[4] += (itp->y)*(itx->y)/P.size();//(2,2)
		cov[5] += (itp->y)*(itx->z)/P.size();//(2,3)
		cov[6] += (itp->z)*(itx->x)/P.size();//(3,1)
		cov[7] += (itp->z)*(itx->y)/P.size();//(3,2)
		cov[8] += (itp->z)*(itx->z)/P.size();//(3,3)
	}

	//build 4*4 symetric matrix
	float Q[16];
	//first line
	Q[0] = cov[0]+cov[4]+cov[8];
	Q[1] = cov[5]-cov[7];
	Q[2] = cov[6]-cov[2];
	Q[3] = cov[1]-cov[3];
	//second line
	Q[4] = cov[5]-cov[7];
	Q[5] = cov[0]-cov[4]-cov[8];
	Q[6] = cov[1]+cov[3];
	Q[7] = cov[6]+cov[2];
	//third line
	Q[8] = cov[6]-cov[2];
	Q[9] = cov[1]+cov[3];
	Q[10] = cov[4]-cov[0]-cov[8];
	Q[11] = cov[5]+cov[7];
	//last line
	Q[12] = cov[1]-cov[3];
	Q[13] = cov[2]+cov[6];
	Q[14] = cov[5]+cov[7];
	Q[15] = cov[8]-cov[0]-cov[4];
	
	Matrix4f A;
	A << Q[0], Q[1], Q[2], Q[3], Q[4], Q[5], Q[6], Q[7], Q[8], Q[9], Q[10], Q[11], Q[12], Q[13], Q[14], Q[15];
	EigenSolver<Matrix4f> es(A);

	Matrix4f D = es.pseudoEigenvalueMatrix();
	Matrix4f V = es.pseudoEigenvectors();

	float biggestValue=D(0,0);
	int pos=0;
	for (int i=1; i<4; i++){
        if (biggestValue<D(i,i)){biggestValue=D(i,i);pos=i;}
	}
	float q[]={V(0,pos), V(1,pos), V(2,pos), V(3,pos)};
	
	//calculate rotation matrix with unit quaternion
	R.x1 = q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];//(1,1)
	R.y1 = 2*(q[1]*q[2]-q[0]*q[3]);//(1,2)
	R.z1 = 2*(q[1]*q[3]+q[0]*q[2]);//(1,3)
	R.x2 = 2*(q[1]*q[2]+q[0]*q[3]);//(2,1)
	R.y2 = q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];//(2,2)
	R.z2 = 2*(q[2]*q[3]-q[0]*q[1]);//(2,3)
	R.x3 = 2*(q[1]*q[3]-q[0]*q[2]);//(3,1)
	R.y3 = 2*(q[2]*q[3]+q[0]*q[1]);//(3,2)
	R.z3 = q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];//(3,3)
}


int _tmain(int argc, _TCHAR* argv[])
{
	
	FILE *fp=NULL;
	float f;
	Point3D xyz,T,mean;
	Rotation R;
	vector<Point3D> P,X,Q;

	char* F_PATH="C:\Users\li\Desktop\Localization\Data\bunny.asc"; 
	//int iLength ;  
	//iLength = WideCharToMultiByte(CP_ACP, 0, argv[1], -1, NULL, 0, NULL, NULL);      
	//WideCharToMultiByte(CP_ACP, 0, argv[1], -1, F_PATH, iLength, NULL, NULL);
	fp=fopen(F_PATH,"r");
	while(fscanf(fp,"%f",&f)){
		xyz.x=f;printf("%f/n",f);
		fscanf(fp,"%f",&f);
		xyz.y=f;
		fscanf(fp,"%f",&f);
		xyz.z=f;
		P.push_back(xyz);
	}
	if (fclose(fp)) {printf("error, file fails to close./n");}
	fp=NULL;

	/*F_PATH="C:\Users\li\Desktop\Localization\Data\bunny_perturbed.asc";
	//iLength = WideCharToMultiByte(CP_ACP, 0, argv[1], -1, NULL, 0, NULL, NULL);      
	//WideCharToMultiByte(CP_ACP, 0, argv[2], -1, F_PATH, iLength, NULL, NULL);
	fp=fopen(F_PATH,"r");
	while(fscanf(fp,"%f",&f)){
		xyz.x=f;
		fscanf(fp,"%f",&f);
		xyz.y=f;
		fscanf(fp,"%f",&f);
		xyz.z=f;
		Q.push_back(xyz);
	}
	if (fclose(fp)) {printf("error, file fails to close./n");}
	fp=NULL;*/


	return 0;
}

