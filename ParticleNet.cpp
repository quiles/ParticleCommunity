

#include "ParticleNet.h"


//#define RAND rand()

//#include <math.h>
#include <time.h>
//#include <stdlib.h>
//#include <stdio.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
//#include <stdio.h>
#include <iomanip>


TParticleNet::TParticleNet(const char *filename){
	TCentroid centroid;
    TParticle particle;
    
    Network = TSnap::LoadEdgeList<PUNGraph>(filename,0,1);
    
    for (TUNGraph::TNodeI NI = Network->BegNI(); NI < Network->EndNI(); NI++) {
        particle.x = (float)(rand()%2000 - 1000) / 10000.0;
        particle.y = (float)(rand()%2000 - 1000) / 10000.0;
        particle.z = (float)(rand()%2000 - 1000) / 10000.0;
        particle.index = 0;
        particle.indexReal = 0;
        particle.node_id = NI.GetId();
        particle.node = NI;
        particle.degree = NI.GetDeg();
        Particles.push_back(particle);
    }

	centroid.x = 0.0;
	centroid.y = 0.0;
    centroid.z = 0.0;
	centroid.comm_id = 1;
    nextComId = 2; // id of the next detected community (used to identify the centroids/communities);
	Centroids.push_back(centroid);
}

TParticleNet::~TParticleNet() {
}

void TParticleNet::LoadCommunities(const char *filename, int format){
    ifstream file;
    string line;
    stringstream ss;
    int node_id;
    int com;
    vector<TParticle>::iterator it;
    
    file.open(filename, ifstream::in);
    
    switch (format) {
        case 1: { // LFR
            int teste=0;
            while (file >> node_id && file >> com){
                for (it = Particles.begin() ; it != Particles.end() && it->node_id != node_id ; ++it);
//                cout << node_id << " - " << com << endl;
                it->indexReal = com;
                teste++;
            }
        }; break;
        case 2: { // SNAP
            com = 1;
            while (getline(file, line)){
                ss.clear();
                ss << line;
                while (ss >> node_id){
                    for (it = Particles.begin() ; it != Particles.end() && it->node_id != node_id ; ++it) {
                    }
//                    cout << node_id << " - " << it->node_id << endl;
                    it->indexReal = com;
                }
                com++;
            }
        }; break;
    }
    
    file.close();
}

void TParticleNet::RunByStep(){
    float sumX, sumY, sumZ, r, dist;
    vector<TParticle>::iterator i,j;
    TUNGraph::TNodeI nodeIt;
    
    for (i = Particles.begin() ; i != Particles.end(); ++i){
        i->dxA = 0.0;
        i->dyA = 0.0;
        i->dzA = 0.0;
        i->dxR = 0.0;
        i->dyR = 0.0;
        i->dzR = 0.0;
    }

    for (i = Particles.begin() ; i != Particles.end()-1; ++i){
        nodeIt = Network->GetNI(i->node_id);
        for (j = i+1 ; j != Particles.end(); ++j){
//            if (i->node.IsNbrNId(j->node_id)) { // i and j are connected
            if (nodeIt.IsNbrNId(j->node_id)) { // i and j are connected
                sumX = alpha*(j->x - i->x);
                sumY = alpha*(j->y - i->y);
                sumZ = alpha*(j->z - i->z);
                i->dxA += sumX;
                i->dyA += sumY;
                i->dzA += sumZ;
                j->dxA -= sumX;
                j->dyA -= sumY;
                j->dzA -= sumZ;
            }
            else { // otherwise
                r = pow(i->x - j->x,2) + pow(i->y - j->y,2) + pow(i->z - j->z,2);
                r = sqrt(r);
                dist = exp(-gamma*r)*beta;
                sumX = dist*(j->x - i->x)/r;
                sumY = dist*(j->y - i->y)/r;
                sumZ = dist*(j->z - i->z)/r;
                i->dxR += sumX;
                i->dyR += sumY;
                i->dzR += sumZ;
                j->dxR -= sumX;
                j->dyR -= sumY;
                j->dzR -= sumZ;
            }
        }
    }
    
    for (i = Particles.begin() ; i != Particles.end(); ++i){
        if (i->degree != 0){
            i->x += (eta*(i->dxA - i->dxR)) / i->degree;
            i->y += (eta*(i->dyA - i->dyR)) / i->degree;
            i->z += (eta*(i->dzA - i->dzR)) / i->degree;
//            i->x += (eta*(i->dxA - i->dxR));
//            i->y += (eta*(i->dyA - i->dyR));
//            i->z += (eta*(i->dzA - i->dzR));
        }
    }
    
}

int TParticleNet::RunModel(){
    float sumX, sumY, sumZ, r, dist;
    vector<TParticle>::iterator i,j;
    TUNGraph::TNodeI nodeIt;
    float R, oldR;
    int steps=0;
    
    do {
//        assignCentroids();
//        if (toAddCentroid || toRemoveCentroid){
//            if (toRemoveCentroid) {
//                removeCentroid();
//                				cout << "R" << endl;
//            }
//            if (toAddCentroid) {
//                addCentroid();
//                				cout << "A" << endl;
//            }
//            assignCentroids();
//            //			mergeCentroids();
//        }
//        computeCentroids();

        for (i = Particles.begin() ; i != Particles.end(); ++i){
            i->dxA = 0.0;
            i->dyA = 0.0;
            i->dzA = 0.0;
            i->dxR = 0.0;
            i->dyR = 0.0;
            i->dzR = 0.0;
        }
        
        for (i = Particles.begin() ; i != Particles.end()-1; ++i){
            nodeIt = Network->GetNI(i->node_id);
            for (j = i+1 ; j != Particles.end(); ++j){
                if (nodeIt.IsNbrNId(j->node_id)) { // i and j are connected
                    sumX = alpha*(j->x - i->x);
                    sumY = alpha*(j->y - i->y);
                    sumZ = alpha*(j->z - i->z);
                    i->dxA += sumX;
                    i->dyA += sumY;
                    i->dzA += sumZ;
                    j->dxA -= sumX;
                    j->dyA -= sumY;
                    j->dzA -= sumZ;
                }
                else { // otherwise
                    r = pow(i->x - j->x,2) + pow(i->y - j->y,2) + pow(i->z - j->z,2);
                    r = sqrt(r);
                    dist = exp(-gamma*r)*beta;
                    sumX = dist*(j->x - i->x)/r;
//cout << fabs(i->x) << endl;
                    sumY = dist*(j->y - i->y)/r;
                    sumZ = dist*(j->z - i->z)/r;
                    i->dxR += sumX;
                    i->dyR += sumY;
                    i->dzR += sumZ;
                    j->dxR -= sumX;
                    j->dyR -= sumY;
                    j->dzR -= sumZ;
                }
            }
        }
        
        oldR = R;
        R=0.0;
        for (i = Particles.begin() ; i != Particles.end(); ++i){
            if (i->degree != 0){
                i->x += (eta*(i->dxA - i->dxR)) / i->degree;
                i->y += (eta*(i->dyA - i->dyR)) / i->degree;
                i->z += (eta*(i->dzA - i->dzR)) / i->degree;
//                i->x += (eta*(i->dxA - i->dxR));
//                i->y += (eta*(i->dyA - i->dyR));
//                i->z += (eta*(i->dzA - i->dzR));

                R += fabs(i->dxR) + fabs(i->dyR) + fabs(i->dzR);
//                cout << fabs(i->dxR) << endl;
            }
        }
        ++steps;
        cout << "Repulsao: " << R << " Diferenca: " << fabs(R-oldR) << endl;
    } while (fabs(R - oldR) > 1.0);

    return steps;
}


void TParticleNet::assignCentroids(){
    float dist, dist2, maxError=0.0;
    vector<TCentroid>::iterator c, c_assigned;
    vector<TParticle>::iterator p;
    
    for (c=Centroids.begin() ; c!=Centroids.end(); ++c){
        c->error = 0.0;
        c->nparticles = 0;
    }
    

    
    for (p=Particles.begin() ; p!=Particles.end(); ++p){
        c = Centroids.begin();
        dist = pow(p->x - c->x,2) + pow(p->y - c->y,2) + pow(p->z - c->z,2);
        c_assigned = c;
        ++c;
        while (c!=Centroids.end()){
            dist2 = pow(p->x - c->x,2) + pow(p->y - c->y,2) + pow(p->z - c->z,2);
            if (dist2 < dist){
                c_assigned = c;
                dist = dist2;
            }
            ++c;
        }
        p->index = c_assigned->comm_id;
        c_assigned->error += dist;
        c_assigned->nparticles++;
        if (dist > maxError) {
            idCentroidMaxError = c_assigned - Centroids.begin(); // particle's index in the vector
            maxError = dist;
        }
    }

    addCentroidThreshold = 0.5;
    
    toAddCentroid = false;
    toRemoveCentroid = false;
    for (c=Centroids.begin() ; c!=Centroids.end() ; ++c){
        if (c->nparticles) c->error /= c->nparticles;
        if (c->error == 0.0) {
            toRemoveCentroid = true;
            centroid2remove = c - Centroids.begin(); // centroid's index in the vector
        }
        else if (c->error > addCentroidThreshold){
            toAddCentroid = true;
//            if (c->error < 2.0*thresholdToAdd) idCentroidMaxError = c - Centroids.begin();
//            // if the centroid error is small, it means that a group of particles is dividing
//            // otherwise, it represents a group without a previsous assigned centroid.
//            // Whether the group is dividing, the new centroid must be add inside this group.
        }
    }
}

void TParticleNet::computeCentroids(){

    vector<TCentroid>::iterator c;
    vector<TParticle>::iterator p;

    
    for (c=Centroids.begin() ; c!=Centroids.end() ; ++c){
        c->x = c->y = c->z = 0.0;
        c->nparticles = 0;
        
        for (p=Particles.begin() ; p!=Particles.end() ; ++p){
            if (p->index == c->comm_id){
                c->x += p->x;
                c->y += p->y;
                c->z += p->z;
                c->nparticles++;
            }
        }
        if (c->nparticles){
            c->x /= (float) c->nparticles;
            c->y /= (float) c->nparticles;
            c->z /= (float) c->nparticles;
        }
        else {
            c->x = 0.0;
            c->y = 0.0;
            c->z = 0.0;
        }
    }
}

               
void TParticleNet::addCentroid(){
    TCentroid centroid;
    
    centroid.x = Centroids[idCentroidMaxError].x + (float)(rand()%100) / 100.0;
    centroid.y = Centroids[idCentroidMaxError].y + (float)(rand()%100) / 100.0;
    centroid.z = Centroids[idCentroidMaxError].z + (float)(rand()%100) / 100.0;
    centroid.comm_id = nextComId++;
    
    Centroids.push_back(centroid);
    toAddCentroid = false;
}

void TParticleNet::removeCentroid(){
    Centroids.erase(Centroids.begin()+centroid2remove);
    toRemoveCentroid = false;
}


void TParticleNet::CommunityDetection(){

    vector<TCentroid>::iterator c;

    toAddCentroid = false;
    toRemoveCentroid = false;
    
    //    mergeCentroids();
    do {
        assignCentroids();
        if (toAddCentroid || toRemoveCentroid){
            if (toRemoveCentroid) {
                removeCentroid();
            }
            if (toAddCentroid) {
                addCentroid();
            }
            assignCentroids();
        }
        computeCentroids();
        cout << "Number of Seeds: " << Centroids.size() << endl;
        
        
        for (c=Centroids.begin() ; c!=Centroids.end() ; ++c){
            cout << "Centroid: " << c - Centroids.begin() << " " << c->x << " " << c->y << " " << c->z << endl;
        }
        
    } while (toAddCentroid || toRemoveCentroid);
}

bool myfunction (TParticle i, TParticle j) { return (i.node_id<j.node_id); }


void TParticleNet::SaveCommunities(const char *filename){
    ofstream file;

    sort (Particles.begin(), Particles.end(), myfunction);
    
    vector<TParticle>::iterator i;
    
    file.open(filename, ofstream::out);
    for (i = Particles.begin() ; i != Particles.end(); ++i){
        file << i->index << "\t" << i->node_id << endl;
    }
    file.close();
    
}

float TParticleNet::NMI(){
    
}

void TParticleNet::SaveParticlePosition(const char *filename){
    ofstream file;
    vector<TParticle>::iterator i;
    
    file.open(filename, ofstream::out);
    for (i = Particles.begin() ; i != Particles.end(); ++i){
        file << i->x << "\t" << i->y << "\t" << i->z << "\t" << i->index << "\t" << i->indexReal << endl;
    }
    file.close();
}




