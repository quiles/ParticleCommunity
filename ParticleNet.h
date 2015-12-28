/*
 * PhaseSeparation.h
 *
 *  Created on: Set/10/2015
 *      Author: quiles
 */

#ifndef PHASESEPARATION_H_
#define PHASESEPARATION_H_


#include <iostream>
#include <vector>
#include "Snap.h"


using namespace std;

class TParticle{
public:
    int node_id;
    TUNGraph::TNodeI node;
    float x, y, z;
    int degree;
	float dxA, dyA, dzA;
	float dxR, dyR, dzR;
//	float erro;
	int index; // community id obtained by the algorithm
	int indexReal; // real community id
};

class TCentroid{
public:
    float x, y, z;
    float error;
    int nparticles; // store the number of associated particles
    int comm_id;
};


class TParticleNet {
private:
    PUNGraph Network;
    vector <TCentroid> Centroids;
    vector <TParticle> Particles;

    float alpha; // attraction strength
    float beta; // repulsion strength
    float gamma; // repulsion decai, constant at 1.0
    float eta; // numerical integration step (delta T), constant at 1.0
    float addCentroidThreshold;
    float mergeCentroidThreshold;
    
    int nextComId;

    bool toAddCentroid; // flag to add a new centroid
    bool toRemoveCentroid; // flag to remove a centroid
    int centroid2remove;
    int idCentroidMaxError; // store the id of the centroid with the largest error

    void assignCentroids();
    void computeCentroids();
    void removeCentroid();
    void addCentroid();
    
public:

    TParticleNet(const char *filename);
    ~TParticleNet();

    // format: 1 -> LFR files
    //         2 -> SNAP files
    void LoadCommunities(const char *filename, int format);

    PUNGraph GetNetwork(){
        return Network;
    };
    void SetModelParameters(float a, float b, float g){
        alpha=a; beta=b; gamma=g; eta=1.0;
    };
    void SetDetectionParameters(float a, float m){
        addCentroidThreshold=a; mergeCentroidThreshold=m;
    };
    
    void RunByStep();
    int RunModel();
    
    void SaveParticlePosition(const char *filename);
    void CommunityDetection();
    void SaveCommunities(const char *filename);
    float NMI();
    
//	int iteracao;
//	int stepAddCentroid; // define o momento de inser��o de novos centroids
//	int centroidToRemove;
//	int countDown;
//	float fATotal,fRTotal;
//	int maxError; // armazena o indice do v�rtice com maior erro em rela��o ao centroid
//	float accErrorCentroid; // armazena o erro m�dio acumulado de todos os centroids
//
//	int startCentroid; // inicio do c�lculo dos centroids
//	int startCountDown;
//	int startStepAdd;
//	int startMinComSize;
//
//	bool calcCentroids();
//	float detecta();
//	void calcDegree();
//	void resetCentroids();
//	void addCentroids();
//	void removeCentroids(int i);
//	void mergeCentroids();

//	void Reset();
//
//	int Run5(float eta, float at, float rep, float decai);
//	int Run6(float eta, float at, float rep, float decai);
//	int Run6L(float eta, float at, float rep, float decai, int it);
//	void SavePosition(string fname);
//	void SaveCentroids(string fname);
//	void GetPrecisionMI(float &pMI, float &pP, const char *atributo);
//	float GetFTotal();
//	void GetNumC(int &numCom, int &numCen);
//
//	void SetParametersDetection(int ce, int cd, int sa, int cs);
//	void SetThresholdToAdd(float threshold);
};

#endif /* PHASESEPARATION_H_ */


