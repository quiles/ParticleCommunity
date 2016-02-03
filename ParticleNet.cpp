

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
#include <algorithm>



TParticleNet::TParticleNet(const char *filename){
	TCentroid centroid;
    TParticle particle;
    
    Network = TSnap::LoadEdgeList<PUNGraph>(filename,0,1);
    
    centroid.x = 0.0;
    centroid.y = 0.0;
    centroid.z = 0.0;
    addCentroidThreshold = 0.5;
    centroid.comm_id = 1;
    numCommunities = 1;
    nextComId = 2; // id of the next detected community (used to identify the centroids/communities);
    Centroids.push_back(centroid);
    
    oldRR2 = 0.0;
    oldRR = 0.0;
    RR = 0.0;
    
    numClusters = 0;

    for (TUNGraph::TNodeI NI = Network->BegNI() ; NI < Network->EndNI(); NI++) {
        particle.x = (float)(rand()%2000 - 1000) / 10000.0;
        particle.y = (float)(rand()%2000 - 1000) / 10000.0;
        particle.z = (float)(rand()%2000 - 1000) / 10000.0;
        particle.index = NULL;
        particle.indexReal = 0;
        particle.node_id = NI.GetId();
//        particle.node = NI;
        particle.degree = NI.GetDeg();
        particle.cluster_id = 0;
        Particles.push_back(particle);
        RefNodes.push_back(NI.GetId());
    }

}

TParticleNet::~TParticleNet() {
}

void TParticleNet::ReloadNetwork(const char *filename){
    TParticle particle;
    vector<TParticle>::iterator it;
    ifstream file;
    
    cout << " nodes: " << Network->GetNodes() << " # of particles: " << Particles.size() <<  endl;

    //Still need to find a better and optimized way to update the network....
    
    Network = TSnap::LoadEdgeList<PUNGraph>(filename,0,1);

    // Case 1: find new nodes to add the new particles
    for (TUNGraph::TNodeI NI = Network->BegNI(); NI < Network->EndNI(); NI++) {
        if (find(RefNodes.begin(), RefNodes.end(), NI.GetId()) == RefNodes.end()){
            particle.x = (float)(rand()%2000 - 1000) / 10000.0;
            particle.y = (float)(rand()%2000 - 1000) / 10000.0;
            particle.z = (float)(rand()%2000 - 1000) / 10000.0;
            particle.index = NULL;
            particle.indexReal = 0;
            particle.node_id = NI.GetId();
            particle.degree = NI.GetDeg();
            particle.cluster_id = 0;
            Particles.push_back(particle);
//            cout << "A new particle associated to node: " << particle.node_id << " was inserted\n";
            RefNodes.push_back(NI.GetId());
        }
    }
    
    // Case 2: delete particles related to deleted nodes
    for (it = Particles.begin() ; it != Particles.end(); ++it){
        if (!Network->IsNode(it->node_id)){
//            cout << "Particle associated to node: " << it->node_id << " was deleted\n";
            RefNodes.erase(find(RefNodes.begin(), RefNodes.end(), it->node_id));
            Particles.erase(it);
        }
    }
}


void TParticleNet::AddNode(int node_id){
    TParticle particle;
    Network->AddNode(node_id);
    particle.x = (float)(rand()%2000 - 1000) / 10000.0;
    particle.y = (float)(rand()%2000 - 1000) / 10000.0;
    particle.z = (float)(rand()%2000 - 1000) / 10000.0;
    particle.index = NULL;
    particle.indexReal = 0;
    particle.node_id = node_id;
    particle.degree = 0;
    particle.cluster_id = 0;
    Particles.push_back(particle);
}

void TParticleNet::AddLink(int i, int j){
    vector<TParticle>::iterator it;
    if (Network->IsNode(i) && Network->IsNode(j)){
        Network->AddEdge(i,j);
        for (it = Particles.begin() ; it != Particles.end(); ++it){
            if (it->node_id == i) it->degree++;
            if (it->node_id == j) it->degree++;
        }
    }
}

void TParticleNet::DeleteNode(int node_id){
    vector<TParticle>::iterator it;
    if (Network->IsNode(node_id)) {
        Network->DelNode(node_id);
        for (it = Particles.begin() ; it != Particles.end(); ++it){
            if (it->node_id == node_id){
                Particles.erase(it);
                return;
            }
        }
    }
}

void TParticleNet::DeleteLink(int i, int j){
    vector<TParticle>::iterator it;
    if (Network->IsNode(i) && Network->IsNode(j)){
        Network->DelEdge(i,j);
        for (it = Particles.begin() ; it != Particles.end(); ++it){
            if (it->node_id == i) it->degree--;
            if (it->node_id == j) it->degree--;
        }
    }
}


void TParticleNet::ResetParticles() {
    TCentroid centroid;
    TParticle particle;

    Particles.clear();
    Centroids.clear();
    
    
    for (TUNGraph::TNodeI NI = Network->BegNI(); NI < Network->EndNI(); NI++) {
        particle.x = (float)(rand()%2000 - 1000) / 10000.0;
        particle.y = (float)(rand()%2000 - 1000) / 10000.0;
        particle.z = (float)(rand()%2000 - 1000) / 10000.0;
        particle.index = NULL;
        particle.indexReal = 0;
        particle.node_id = NI.GetId();
//        particle.node = NI;
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


void TParticleNet::LoadCommunities(const char *filename, int format){
    ifstream file;
    string line;
    stringstream ss;
    int node_id;
    int com;
    vector<TParticle>::iterator it;
    
    numCommunities = 0;
    
    file.open(filename, ifstream::in);

    
    switch (format) {
        case 1: { // LFR
            int teste=0;
            cout << "Loading Loading: " << file.is_open()<< endl;
            while (file >> node_id && file >> com){
                for (it = Particles.begin() ; it != Particles.end() && it->node_id != node_id ; ++it);
                it->indexReal = com;
                teste++;
                if (com > numCommunities) numCommunities = com;
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
                    it->indexReal = com;
                    if (com > numCommunities) numCommunities = com;
                }
                com++;
            }
        }; break;
    }
    
    file.close();
}

// still under development -> it aims to reduce the number of steps when new nodes are added.
// Some iterations performed only on the new nodes can make the convergence quick...

void TParticleNet::RunForNewNodes(int steps){
    float sumX, sumY, sumZ, r, dist;
    vector<TParticle>::iterator i,j;
    TUNGraph::TNodeI nodeIt;
    int s;
    
    
    for (s=0 ; s<steps ; s++){
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
                r = pow(i->x - j->x,2) + pow(i->y - j->y,2) + pow(i->z - j->z,2);
                r = sqrt(r);
                if (r==0) r = 0.0001;
                
                if (nodeIt.IsNbrNId(j->node_id)) { // i and j are connected
                    sumX = (j->x - i->x)/r;
                    sumY = (j->y - i->y)/r;
                    sumZ = (j->z - i->z)/r;
                    i->dxA += sumX;
                    i->dyA += sumY;
                    i->dzA += sumZ;
                    j->dxA -= sumX;
                    j->dyA -= sumY;
                    j->dzA -= sumZ;
                }
                else { // otherwise
                    dist = exp(-r);
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
                i->x += ((alpha*i->dxA - beta*i->dxR)) / (float)i->degree;
                i->y += ((alpha*i->dyA - beta*i->dyR)) / (float)i->degree;
                i->z += ((alpha*i->dzA - beta*i->dzR)) / (float)i->degree;
            }
        }
    }
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
            r = pow(i->x - j->x,2) + pow(i->y - j->y,2) + pow(i->z - j->z,2);
            r = sqrt(r);
            if (r==0) r = 0.0001;

            if (nodeIt.IsNbrNId(j->node_id)) { // i and j are connected
                sumX = (j->x - i->x)/r;
                sumY = (j->y - i->y)/r;
                sumZ = (j->z - i->z)/r;
                i->dxA += sumX;
                i->dyA += sumY;
                i->dzA += sumZ;
                j->dxA -= sumX;
                j->dyA -= sumY;
                j->dzA -= sumZ;
            }
            else { // otherwise
                dist = exp(-r);
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
            i->x += ((alpha*i->dxA - beta*i->dxR)) / (float)i->degree;
            i->y += ((alpha*i->dyA - beta*i->dyR)) / (float)i->degree;
            i->z += ((alpha*i->dzA - beta*i->dzR)) / (float)i->degree;
        }
    }
}

int TParticleNet::RunModel(int maxIT, float minDR, bool verbose){
    float sumX, sumY, sumZ, r, dist;
    vector<TParticle>::iterator i,j;
    TUNGraph::TNodeI nodeIt;
//    float oldRR, oldRR2, RR=0, oldACC, deltaOld=0.0;
    int steps=0;
    float value;
    
    do {
        
        for (i = Particles.begin() ; i != Particles.end(); ++i){
            i->dxA = 0.0;
            i->dyA = 0.0;
            i->dzA = 0.0;
            i->dxR = 0.0;
            i->dyR = 0.0;
            i->dzR = 0.0;
        }
        // Algorithm CORE (O(n^2))
        for (i = Particles.begin() ; i != Particles.end()-1; ++i){
            nodeIt = Network->GetNI(i->node_id);
            for (j = i+1 ; j != Particles.end(); ++j){
                r = pow(i->x - j->x,2) + pow(i->y - j->y,2) + pow(i->z - j->z,2);
                r = sqrt(r);
                if (r==0) r = 0.0001;
                if (nodeIt.IsNbrNId(j->node_id)) { // i and j are connected
                    sumX = (j->x - i->x)/r;
                    sumY = (j->y - i->y)/r;
                    sumZ = (j->z - i->z)/r;
                    i->dxA += sumX;
                    i->dyA += sumY;
                    i->dzA += sumZ;
                    j->dxA -= sumX;
                    j->dyA -= sumY;
                    j->dzA -= sumZ;
                }
                else { // otherwise -> Most of the computation time is wasted in this part (repulsion)
                    dist = exp(-r);
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
        
//        oldACC = accError;
        oldRR2 = oldRR;
        oldRR = RR;
        RR = 0.0;
        for (i = Particles.begin() ; i != Particles.end(); ++i){
            if (i->degree != 0){
                sumX = ((alpha*i->dxA - beta*i->dxR)) / (float)i->degree;
                sumY = ((alpha*i->dyA - beta*i->dyR)) / (float)i->degree;
                sumZ = ((alpha*i->dzA - beta*i->dzR)) / (float)i->degree;
                i->x += sumX;
                i->y += sumY;
                i->z += sumZ;
                RR += (fabs(i->dxR) + fabs(i->dyR) + fabs(i->dzR)) / (float)i->degree;
            }
        }
        
        RR = RR*0.1 + oldRR*0.9; // used to make the evolution of RR more stable (due to the numerical integration step DT=1.0.
        ++steps;
    
        value = fabs(oldRR2-RR); // / (float)Particles.size();
        if (verbose) cout << steps << " " << RR << " " << value << endl;
    } while (steps<maxIT && value > minDR);
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
    accError = 0.0;

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

        p->index = &(*c_assigned); //->comm_id;
        c_assigned->error += dist;
//        accError2 += pow(dist,2);
//        accError += dist;
        c_assigned->nparticles++;
        if (dist > maxError) {
            idCentroidMaxError = c_assigned - Centroids.begin(); // Centroid's index associated to the farthest particle. When needed, the new centroid is added nearby this one.
            maxError = dist;
        }
    }
    
    for (c=Centroids.begin() ; c!=Centroids.end(); ++c){
        accError += c->error;
    }
    

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
            if (p->index == &(*c)){ //->comm_id){
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
//            cout << "D";
            c->x = 0.0;
            c->y = 0.0;
            c->z = 0.0;
        }
    }
}

               
void TParticleNet::addCentroid(){
    TCentroid centroid;
    vector<int>::iterator i;

    
    centroid.x = Centroids[idCentroidMaxError].x + (float)(rand()%10) / 100.0;
    centroid.y = Centroids[idCentroidMaxError].y + (float)(rand()%10) / 100.0;
    centroid.z = Centroids[idCentroidMaxError].z + (float)(rand()%10) / 100.0;
    
    centroid.comm_id = nextComId++;
    numCommunities++;
    
    Centroids.push_back(centroid);
    toAddCentroid = false;
//    cout << "A";
}

void TParticleNet::removeCentroid(){
    Centroids.erase(Centroids.begin()+centroid2remove);
    toRemoveCentroid = false;
    numCommunities--;
//    cout << "R";
}

void TParticleNet::mergeCentroids(){
    vector<TCentroid>::iterator c1,c2,c_remove;
    float min=99999, dist;

    for (c1=Centroids.begin() ; c1!=(Centroids.end()-1) ; ++c1){
        for (c2=c1+1 ; c2!=(Centroids.end()) ; ++c2){
            dist = pow(c1->x - c2->x, 2) + pow(c1->y - c2->y,2) + pow(c1->z - c2->z,2);
            if (dist < min){
                min = dist;
                c_remove = c1;
            }
//            dist = pow(dist, 0.5);
//            if (dist < 0.9) {
//                
//            }
        }
    }
    if (sqrt(min)<0.9) {
//        cout << "M";
        Centroids.erase(c_remove);
        numCommunities--;
    }
}

int TParticleNet::CommunityDetection3(){
    
    vector<TCentroid>::iterator c;
    vector<TParticle>::iterator p;
    float oldAcc;
    float diffError;
    int steps=0;
    
    toAddCentroid = false;
    toRemoveCentroid = false;
    
    
    do {
        oldAcc = accError;
        mergeCentroids();
        assignCentroids();
        if (toAddCentroid || toRemoveCentroid){
            if (toRemoveCentroid) removeCentroid();
            if (toAddCentroid) addCentroid();
            assignCentroids();
        }
        computeCentroids();
        accError = accError*0.1 + oldAcc*0.9; // used to make the evolution of the centroid error more stable...
        diffError = fabs(oldAcc - accError); // / (float)Centroids.size();
        ++steps;
    } while (diffError>0.01);
    return steps;
}


bool myfunction (TParticle i, TParticle j) { return (i.node_id<j.node_id); }

void TParticleNet::SaveCentroids(const char *filename){
    ofstream file;
        
    vector<TCentroid>::iterator i;

    int aux;
    for (aux=1, i=Centroids.begin() ; i!=Centroids.end(); ++i, aux++){
        i->comm_id = aux;
    }
    
    file.open(filename, ofstream::out);
    file << "% id_centroids; position; # of associated particles; total repulsion; total attraction; sum of dists " << endl;
    for (i = Centroids.begin() ; i != Centroids.end(); ++i){
        file << i->x << "\t"
             << i->y << "\t"
             << i->z << "\t"
             << (int) i->comm_id << "\t"
             << i->nparticles << "\t"
             << i->totalR << "\t"
             << i->totalA << "\t"
             << i->error << endl;
    }
    file.close();
}

    
void TParticleNet::SaveCommunities(const char *filename){
    ofstream file;
    
    int aux;
    vector<TCentroid>::iterator c;
    for (aux=1, c=Centroids.begin() ; c!=Centroids.end(); ++c, aux++){
        c->comm_id = aux;
    }
    sort (Particles.begin(), Particles.end(), myfunction);
    
    vector<TParticle>::iterator i;
    
    file.open(filename, ofstream::out);
    for (i = Particles.begin() ; i != Particles.end(); ++i){
//        file << i->index << "\t" << i->node_id << endl;
        file << i->node_id <<  "\t" << i->index->comm_id << endl;
    }
    file.close();
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////


float TParticleNet::NMI(){
    
    int i,j;
    int **labelCount;
    int *numCount, *sumI, *sumJ;
    float precision, Si=0,Sj=0;
    int p, valor, comUser;
    int N, numCentroids;
    
    
    if (Centroids.size()<=1) return 0;
    
    vector<TCentroid>::iterator c;
    vector<TParticle>::iterator part;
    for (i=0, c=Centroids.begin() ; c!=Centroids.end(); ++c, i++){
        c->comm_id = i;
    }
    
//    	cout << "Com: " << numCommunities << endl;
//    	cout << "Cen: " << numCentroids << endl;
    
    numCentroids = Centroids.size();
    N = Particles.size();
    
    if (numCentroids > numCommunities) comUser = numCentroids;
    else comUser = numCommunities;
    
    labelCount = new int*[comUser];
    sumI = new int[comUser];
    sumJ = new int[comUser];
    for (i=0 ; i<comUser ; i++){
        sumI[i] = sumJ[i] = 0.0;
        labelCount[i] = new int[comUser];
        for (j=0 ; j<comUser ; j++){
            labelCount[i][j] = 0;
        }
    }
    
    for (part = Particles.begin() ; part != Particles.end() ; ++part){
        labelCount[part->index->comm_id][part->indexReal-1]++;
    }

    
    // sorting by rows
    for (j=0 ; j<comUser-1 ; j++){
        p = j;
        valor = labelCount[j][j];
        for (i=j+1 ; i<comUser ; i++){
            if (labelCount[i][j] > valor){
                p = i;
                valor = labelCount[i][j];
            }
        }
        if (p!=j) {
            numCount = labelCount[j];
            labelCount[j] = labelCount[p];
            labelCount[p] = numCount;
        }
    }
    
    // sorting by columns
    for (j=0 ; j<comUser-1 ; j++){
        p = j;
        valor = labelCount[j][j];
        for (i=j+1 ; i<comUser ; i++){
            if (labelCount[j][i] > valor){
                p = i;
                valor = labelCount[j][i];
            }
        }
        if (p!=j) {
            for (i=0 ; i<comUser ; i++){
                valor = labelCount[i][j];
                labelCount[i][j] = labelCount[i][p];
                labelCount[i][p] = valor;
            }
        }
    }
    

//    cout << "Matrix:\n";
    for (j=0 ; j<comUser ; j++){
        for (i=0 ; i<comUser ; i++){
            sumJ[j] += labelCount[j][i];
            sumI[i] += labelCount[j][i];
//            cout << labelCount[j][i] << "\t";
        }
//        cout << ";" << endl;
    }

    comUser = numCommunities;
    
    precision = 0.0;
    Si = 0.0;
    Sj = 0.0;
    int pNormal=0;
    for (i=0 ; i<comUser ; i++){
        pNormal += labelCount[i][i];
        for (j=0 ; j<comUser ; j++){
            if (labelCount[j][i]){
                precision += (float)labelCount[j][i] * log( (float)(labelCount[j][i]*N)  /  (float)(sumI[i]*sumJ[j]) );
            }
        }
        if (sumI[i]) Si += (float)sumI[i]*log((float)sumI[i]/(float)N);
        if (sumJ[i]) Sj += (float)sumJ[i]*log((float)sumJ[i]/(float)N);
    }
    
    precision *= -2;
    precision /= (Si+Sj);
    
//    for (i=0 ; i<comUser ; i++)
//        delete[] labelCount[i];
//    delete[] labelCount;
//    delete[] sumI;
//    delete[] sumJ;
    
    return precision;
}

float TParticleNet::NMI2(){
    vector<TParticle>::iterator part;
    int i,j;
    int **labelCount;
    int *numCount, *sumI, *sumJ;
    float precision, Si=0,Sj=0;
    int p, valor, comUser;
    int N;
    
    if (!numClusters) return 0;
    
//    numClusters--;
    
    N = Particles.size();

    if (numClusters > numCommunities) comUser = numClusters;
    else comUser = numCommunities;
    
//    cout << "#Clusters: " << numClusters << endl;

    labelCount = new int*[comUser];
    sumI = new int[comUser];
    sumJ = new int[comUser];
    for (i=0 ; i<comUser ; i++){
        sumI[i] = sumJ[i] = 0.0;
        labelCount[i] = new int[comUser];
        for (j=0 ; j<comUser ; j++){
            labelCount[i][j] = 0;
        }
    }

//    cout << "#Clusters: " << numClusters << endl;
    
    
    for (part = Particles.begin() ; part != Particles.end() ; ++part){
//        cout << part->node_id << " -> " << part->cluster_id << endl;
//        labelCount[part->cluster_id-1][part->indexReal-1]++;
        labelCount[part->cluster_id][part->indexReal-1]++;
    }
    
    
    // sorting by rows
    for (j=0 ; j<comUser-1 ; j++){
        p = j;
        valor = labelCount[j][j];
        for (i=j+1 ; i<comUser ; i++){
            if (labelCount[i][j] > valor){
                p = i;
                valor = labelCount[i][j];
            }
        }
        if (p!=j) {
            numCount = labelCount[j];
            labelCount[j] = labelCount[p];
            labelCount[p] = numCount;
        }
    }
    
    // sorting by columns
    for (j=0 ; j<comUser-1 ; j++){
        p = j;
        valor = labelCount[j][j];
        for (i=j+1 ; i<comUser ; i++){
            if (labelCount[j][i] > valor){
                p = i;
                valor = labelCount[j][i];
            }
        }
        if (p!=j) {
            for (i=0 ; i<comUser ; i++){
                valor = labelCount[i][j];
                labelCount[i][j] = labelCount[i][p];
                labelCount[i][p] = valor;
            }
        }
    }
    
    
//    cout << "Matrix:\n";
    for (j=0 ; j<comUser ; j++){
        for (i=0 ; i<comUser ; i++){
            sumJ[j] += labelCount[j][i];
            sumI[i] += labelCount[j][i];
//            cout << labelCount[j][i] << "\t";
        }
//        cout << ";" << endl;
    }
    
    comUser = numCommunities;
    
    precision = 0.0;
    Si = 0.0;
    Sj = 0.0;
    int pNormal=0;
    for (i=0 ; i<comUser ; i++){
        pNormal += labelCount[i][i];
        for (j=0 ; j<comUser ; j++){
            if (labelCount[j][i]){
                precision += (float)labelCount[j][i] * log( (float)(labelCount[j][i]*N)  /  (float)(sumI[i]*sumJ[j]) );
            }
        }
        if (sumI[i]) Si += (float)sumI[i]*log((float)sumI[i]/(float)N);
        if (sumJ[i]) Sj += (float)sumJ[i]*log((float)sumJ[i]/(float)N);
    }
    
    precision *= -2;
    precision /= (Si+Sj);
    
    //    for (i=0 ; i<comUser ; i++)
    //        delete[] labelCount[i];
    //    delete[] labelCount;
    //    delete[] sumI;
    //    delete[] sumJ;
    
    return precision;
}

float TParticleNet::NMIH(int nivel){
    
    int i,j;
    int **labelCount;
    int *numCount, *sumI, *sumJ;
    float precision, Si=0,Sj=0;
    int p, valor, comUser;
    int N, numCentroids;
    
    
    ShrinkCentroids();
    
    if (Centroids.size()<=1) return 0;
    
    vector<TCentroid>::iterator c;
    vector<TParticle>::iterator part;
    for (i=0, c=Centroids.begin() ; c!=Centroids.end(); ++c, i++){
        c->comm_id = i;
    }
    
    //    	cout << "Com: " << numCommunities << endl;
    //    	cout << "Cen: " << numCentroids << endl;
    
    numCentroids = Centroids.size();
    N = Particles.size();
    
    numCommunities = numCommunitiesH[nivel];
    
//    cout << "\n\n# Com: " << numCommunities << endl;
    
    if (numCentroids > numCommunities) comUser = numCentroids;
    else comUser = numCommunities;
    
    labelCount = new int*[comUser];
    sumI = new int[comUser];
    sumJ = new int[comUser];
    for (i=0 ; i<comUser ; i++){
        sumI[i] = sumJ[i] = 0.0;
        labelCount[i] = new int[comUser];
        for (j=0 ; j<comUser ; j++){
            labelCount[i][j] = 0;
        }
    }
    
    for (part = Particles.begin() ; part != Particles.end() ; ++part){
        labelCount[part->index->comm_id][part->indexRealH[nivel]-1]++;
//        cout << "Par: " << part->node_id << " com ac: " << part->index->comm_id << " com re: " << part->indexRealH[nivel] << endl;
    }
    
    
    // sorting by rows
    for (j=0 ; j<comUser-1 ; j++){
        p = j;
        valor = labelCount[j][j];
        for (i=j+1 ; i<comUser ; i++){
            if (labelCount[i][j] > valor){
                p = i;
                valor = labelCount[i][j];
            }
        }
        if (p!=j) {
            numCount = labelCount[j];
            labelCount[j] = labelCount[p];
            labelCount[p] = numCount;
        }
    }
    
    // sorting by columns
    for (j=0 ; j<comUser-1 ; j++){
        p = j;
        valor = labelCount[j][j];
        for (i=j+1 ; i<comUser ; i++){
            if (labelCount[j][i] > valor){
                p = i;
                valor = labelCount[j][i];
            }
        }
        if (p!=j) {
            for (i=0 ; i<comUser ; i++){
                valor = labelCount[i][j];
                labelCount[i][j] = labelCount[i][p];
                labelCount[i][p] = valor;
            }
        }
    }
    
    
//        cout << "Matrix:\n";
    for (j=0 ; j<comUser ; j++){
        for (i=0 ; i<comUser ; i++){
            sumJ[j] += labelCount[j][i];
            sumI[i] += labelCount[j][i];
//                        cout << labelCount[j][i] << "\t";
        }
//                cout << ";" << endl;
    }
//    cout << endl;
    
    comUser = numCommunities;
    
    precision = 0.0;
    Si = 0.0;
    Sj = 0.0;
    int pNormal=0;
    for (i=0 ; i<comUser ; i++){
        pNormal += labelCount[i][i];
        for (j=0 ; j<comUser ; j++){
            if (labelCount[j][i]){
                precision += (float)labelCount[j][i] * log( (float)(labelCount[j][i]*N)  /  (float)(sumI[i]*sumJ[j]) );
            }
        }
        if (sumI[i]) Si += (float)sumI[i]*log((float)sumI[i]/(float)N);
        if (sumJ[i]) Sj += (float)sumJ[i]*log((float)sumJ[i]/(float)N);
    }
    
    precision *= -2;
    precision /= (Si+Sj);
    
    //    for (i=0 ; i<comUser ; i++)
    //        delete[] labelCount[i];
    //    delete[] labelCount;
    //    delete[] sumI;
    //    delete[] sumJ;
    
    return precision;
}



//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////


// fixed << setprecision(2)

void TParticleNet::SaveParticlePosition(const char *filename){
    ofstream file;
    vector<TParticle>::iterator i;
    
    file.open(filename, ofstream::out);
    file << "% Particle's file -> a snapshot of the particle space\n";
    file << "% Simulation Parameters -> alpha: " << alpha << " beta: " << beta << endl;
    file << "% x, y, z, ground truth community id, assigned community id, node_id\n";

    for (i = Particles.begin() ; i != Particles.end(); ++i){
        if (i->index)
            file << fixed << setprecision(2) <<
                    i->x << "\t" << i->y << "\t" << i->z << "\t" <<
                    i->indexReal << "\t" <<
                    i->index->comm_id << "\t" <<
//                    i->cluster_id << "\t" <<
                    i->node_id << endl;
        else
            file << fixed << setprecision(2) <<
                    i->x << "\t" << i->y << "\t" << i->z << "\t" <<
                    i->indexReal << "\t" <<
                    "-1" << "\t" <<
//                    i->cluster_id << "\t" <<
                    i->node_id << endl;
    }
    file.close();
}

int TParticleNet::getNumCommunities(){
    return Centroids.size();
}

int TParticleNet::getNumParticles(){
    return Particles.size();
}

void TParticleNet::fineTuning(int t){
    int i;
    for (i=0 ; i<t ; i++){
        assignCentroids();
        computeCentroids();
    }
}

void TParticleNet::ShrinkCentroids(){
    vector<TCentroid>::iterator i;
    int aux;
    for (aux=1, i=Centroids.begin() ; i!=Centroids.end(); ++i, aux++){
        i->comm_id = aux;
    }
}

void TParticleNet::printCentroids(){
    vector<TCentroid>::iterator i;
    float E=0.0, TR=0.0, TA=0.0;
    
    int aux;
    for (aux=1, i=Centroids.begin() ; i!=Centroids.end(); ++i, aux++){
        i->comm_id = aux;
    }
    
    for (i = Centroids.begin() ; i != Centroids.end(); ++i){
//        cout << "C: " << i->comm_id << " E: " << i->error << " A: " << i->totalA << " R: " << i->totalR << endl;
        E += i->error;
        TA += i->totalA;
        TR += i->totalR;
//        file << i->x << "\t"
//        << i->y << "\t"
//        << i->z << "\t"
//        << (int) i->comm_id << "\t"
//        << i->nparticles << "\t"
//        << i->totalR << "\t"
//        << i->totalA << "\t"
//        << i->error << endl;
    }
//    file.close();
    E /= Centroids.size();
    TA /= Centroids.size();
    TR /= Centroids.size();
    cout << " " << E << " " << TA << " " << TR << endl;
//    cout << " " << E << " " << Centroids.size() << endl;

}

float TParticleNet::printCentroidsError(){
    return accError;
}


////////////////
////////////////
////////////////
////////////////
////////////////
////////////////
////////////////
////////////////
////////////////
////////////////
////////////////


#define EPS 1.0
#define MINPTS 4

vector<TParticle*> TParticleNet::GetRegion(TParticle *p){
    float dist;
    vector<TParticle>::iterator n;
    vector<TParticle*> region;
    
    for (n=Particles.begin() ; n!=Particles.end(); ++n){
        if (&(*n)==p) continue;
        dist = pow(n->x - p->x,2) + pow(n->y - p->y,2) + pow(n->z - p->z,2);
        dist = sqrt(dist);
        if (dist <= EPS) region.push_back(&(*n));
//        cout << dist << "\t";
    }
    return region;
}

bool TParticleNet::ExpandCluster(TParticle *p, int clusterID){
    
    vector<TParticle*> seeds = GetRegion(p);
    
    if (seeds.size() < MINPTS)
    {
        p->cluster_id = -1;
        return false;
    }
    else // all points in seeds are density reachable from point 'p'
    {
        p->cluster_id = clusterID;
        for (int s=0; s<seeds.size(); s++) {
            seeds[s]->cluster_id = clusterID;
//            cout << " " << seeds[s]->node_id;
        }
        while (seeds.size() > 0)
        {
            TParticle* currentP = seeds.back();
            seeds.pop_back();
            vector<TParticle*> seeds2 = GetRegion(currentP);
            if (seeds2.size() >= MINPTS)
            {
                for (int i=0; i<seeds2.size(); i++)
                {
                    TParticle* currentP2 = seeds2.back();
                    if (currentP2->cluster_id == 0  || currentP2->cluster_id == -1)
                    {
                        if (currentP2->cluster_id == 0) seeds.push_back(currentP2);
                        currentP2->cluster_id = clusterID;
//                        cout << " " << currentP2->node_id;
                    }
                }
            }
        }
        return true;
    }
}

void TParticleNet::CommunityDetection2(){
    
    vector<TParticle>::iterator i;
    int clusterID=1;

    for (i=Particles.begin() ; i!=Particles.end(); ++i){
        i->cluster_id = 0;
    }
    for (i = Particles.begin() ; i != Particles.end(); ++i){
//        cout << "Seed: " << i->node_id;
        if (i->cluster_id == 0){
            if (ExpandCluster(&(*i),clusterID)) clusterID++;
        }
//        cout << endl;
    }
    numClusters = clusterID;

//    for (i=Particles.begin() ; i!=Particles.end(); ++i){
//        cout << "Particle: " << i->cluster_id << endl;
//    }
    for (i=Particles.begin() ; i!=Particles.end(); ++i){
        if (i->cluster_id == -1) i->cluster_id = 0;
    }


    
}
