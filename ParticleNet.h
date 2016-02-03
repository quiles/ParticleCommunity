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

#undef min
#undef max

using namespace std;

class TCentroid{
public:
    float x, y, z;
    float error;
    int nparticles; // store the number of associated particles
    int comm_id;
    float totalR; // sum of the repulsion inside the cluster
    float totalA; // sum of the attraction inside the cluster
};

class TParticle{
public:
    int node_id;
//    TUNGraph::TNodeI node;
    float x, y, z;
    int degree;
	float dxA, dyA, dzA;
	float dxR, dyR, dzR;
//	float erro;
    int cluster_id;
	TCentroid *index; // community id obtained by the algorithm / pointer to the centroid
	int indexReal; // real community id
    vector <int> indexRealH; // real community id
};


class TParticleNet {
private:
    PUNGraph Network;
    vector <TCentroid> Centroids;
    vector <TParticle> Particles;
    vector <int> RefNodes;
    
    float alpha; // attraction strength
    float beta; // repulsion strength
    float gamma; // repulsion decai, constant at 1.0
    float eta; // numerical integration step (delta T), constant at 1.0
    float addCentroidThreshold;
    float mergeCentroidThreshold;
//    int centroidTransient;
    
    int nextComId;
    int numCommunities;
    vector<int> numCommunitiesH;
    int numClusters;
    

    bool toAddCentroid; // flag to add a new centroid
    bool toRemoveCentroid; // flag to remove a centroid
    int centroid2remove;
    int idCentroidMaxError; // store the id of the centroid with the largest error
    float accError;
    float RR, oldRR, oldRR2;
    
    
    void assignCentroids();
    void computeCentroids();
    void removeCentroid();
    void addCentroid();
    void mergeCentroids();
    void fineTuning(int t);
    bool ExpandCluster(TParticle *p, int clusterID);
    vector<TParticle*> GetRegion(TParticle *p);
    void ShrinkCentroids();

    
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
    void RunForNewNodes(int steps);

    int RunModel(int maxIT, float minDR, bool verbose);
    
    void SaveParticlePosition(const char *filename);
    void SaveCommunities(const char *filename);
    void SaveCentroids(const char *filename);
    float NMI();
    float NMI2();
    float NMIH(int nivel);
    void ResetParticles();
    int getNumCommunities();
    int getNumParticles();
    void printCentroids();
    float printCentroidsError();
    
    void ReloadNetwork(const char *filename);
    void AddNode(int node_id);
    void AddLink(int i, int j);
    void DeleteNode(int node_id);
    void DeleteLink(int i,int j);
    
    void CommunityDetection2(); // DBScan Algorithm - under development....
    int CommunityDetection3(); // current clustering approach
    
};

#endif /* PHASESEPARATION_H_ */


