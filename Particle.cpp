
#include "Snap.h"

#undef min
#undef max

#include "ParticleNet.h"
#include <iostream>
#include <time.h>
#include <string>
#include <stdlib.h>

using namespace std;

void PrintGStats(const char s[], PUNGraph Graph) {
    cout << "Graph: " << s << " nodes: " << Graph->GetNodes()
         << " edges: " << Graph->GetEdges() << endl;
}

long int numCom=0;
float alpha=1.0;
float beta=0.50, betaM=1.0, betaS=0.1;
bool saveStates=false;
bool timeVarying=false;
bool snapFormat=false;
bool searchBeta=false;
bool verbose=false;
int steps=1000;
float minDR=1.0;
string fName, fNameCom;
int maxSteps=1000;


void Message(int i){
    if (i==0){
        cout << "Community Detection - Particle's Model\n\n";
        cout << "\tUsage: $ ./Particle filename.dat [options]\n\n";
        cout << "\t[options]: \n";
        cout << "\t\t-a value -> attractive parameter (alpha) [default: 1.0]\n";
        cout << "\t\t-b value -> repulsive parameter (beta) [default: 0.1]\n";
        cout << "\t\t-tr value ->  define the stop condition - theta_R  [default: 1.0]\n";
        cout << "\t\t-max value ->  define the maximum number of steps (iterations)  [default: 1000]\n";
//        cout << "\t\t-ss number_of_steps -> save states (files.par & files.cen)\n";
//        cout << "\t\t-sf -> load the input file using the snapFormat\n";
        cout << "\t\t-sb initial_beta final_beta step -> searching beta\n";
        cout << "\t\t-v -> verbose\n";
        cout << "\n\n";
        cout << "\tExample: $ ./Particle rede001.dat -b 0.3 -ss 100\n\n";
    }
    else if (i==1) cout << "error: no valid input files\n\n";
}

void Model0_debug(){
    TParticleNet *Model=NULL;
    string saveName;
    int it, st;
    clock_t ini,end;
    Model = new TParticleNet(fName.c_str());

    Model->LoadCommunities(fNameCom.c_str(),1); // 1 for LFR community file format / 2 for SNAP community file format

    Model->SetModelParameters(alpha, beta, 1.0);
    ini = clock();
    it = Model->RunModel(maxSteps,minDR,false);
    st = Model->CommunityDetection3();
    end = clock();
    cout <<
            Model->getNumParticles() << " " <<
            Model->getNumCommunities() << " " <<
            it << " " <<
            st << " " <<
            Model->NMI() << " " <<
            ((float)(end-ini))/CLOCKS_PER_SEC <<
            endl;
}

void ModelSearchBeta_debug(){
    TParticleNet *Model=NULL;
    string saveName;
    int it, st;
    char out[256];
    Model = new TParticleNet(fName.c_str());
    Model->LoadCommunities(fNameCom.c_str(),1);

    Model->SetModelParameters(1.0, beta, 1.0);

    sprintf(out,"time_0.par");
    saveName = fName;
    saveName.replace(fName.size()-3,3,out);
    Model->SaveParticlePosition(saveName.c_str());
    
    for (float b=beta ; b<betaM ; b+=betaS){
//        cout << "Beta: " << b << endl;
//        cout << b << " 1 1 1 1 1" << endl;
        Model->SetModelParameters(1.0, b, 1.0);
        it = Model->RunModel(maxSteps,minDR,false);
        st = Model->CommunityDetection3();
        cout << "\n\n" << b << " " <<  // beta
                it << " " << // model iterations
                st << " " << // clustering iterations
                Model->getNumCommunities() << " " << // # of communities detected
                Model->printCentroidsError() << " " << // sum of the error of the centroids
//                Model->NMI() << // NMI
                endl << endl << endl;
//        cout << "Steps: " << it << " #Com: " << Model->getNumCommunities() << " NMI: " << Model->NMI() << endl;
//        cout << "ErroR Centroids: " << Model->printCentroidsError() << endl << endl;
        sprintf(out,"beta_%.3f.com",beta);
        saveName = fName;
        saveName.replace(fName.size()-3,3,out);
        Model->SaveCommunities(saveName.c_str());
    }
    
    //    cout << itt << " NMI: " << Model->NMI() << " #Com: " << Model->getNumCommunities() << endl;
}

void ModelSearchBeta(){
    TParticleNet *Model=NULL;
    string saveName;
    int it, st;
    char out[256];
    Model = new TParticleNet(fName.c_str());
    Model->SetModelParameters(1.0, beta, 1.0);
    
    Model->LoadCommunities(fNameCom.c_str(),1);
    
    sprintf(out,"time_0.par");
    saveName = fName;
    saveName.replace(fName.size()-3,3,out);
    cout << "Initial state: " << saveName << endl;
    Model->SaveParticlePosition(saveName.c_str());
    
    for (float b=beta ; b<betaM ; b+=betaS){
        Model->SetModelParameters(1.0, b, 1.0);
        cout << "Model running - [beta: " << b << "]\n";
        it = Model->RunModel(maxSteps,minDR,verbose);
        st = Model->CommunityDetection3();
        cout << "Total steps: " << it << endl;
        cout << "# of communities detected: " << Model->getNumCommunities() << " " << "NMI: " << Model->NMI() << endl;
        cout << "Accumulated centroid error: " << Model->printCentroidsError() << endl;

        sprintf(out,"beta_%.3f.par",b);
        saveName = fName;
        saveName.replace(fName.size()-3,3,out);
        Model->SaveParticlePosition(saveName.c_str());
        cout << "Current state file: " << saveName << endl;

        sprintf(out,"beta_%.3f.com",b);
        saveName = fName;
        saveName.replace(fName.size()-3,3,out);
        Model->SaveCommunities(saveName.c_str());
        cout << "Current community structure: " << saveName << endl << endl;
    }
}


void Model0(){
    TParticleNet *Model=NULL;
    string saveName;
    char out[256];
    int it, st;
    clock_t ini,end;
    
    ini = clock();
    Model = new TParticleNet(fName.c_str());

    sprintf(out,"time_0.par");
    saveName = fName;
    saveName.replace(fName.size()-3,3,out);
    Model->SaveParticlePosition(saveName.c_str());
    cout << "Initial state: " << saveName << endl;
    
    Model->SetModelParameters(alpha, beta, 1.0);
    cout << "Model running...\n";
    it = Model->RunModel(maxSteps,minDR,verbose);
    cout << "Detecting clusters...\n";
    st = Model->CommunityDetection3();
    end = clock();

    cout << "Total steps: " << it << endl;
    cout << "# of communities detected: " << Model->getNumCommunities() << endl;
    cout << "Elapsed time (s): " << ((float)(end-ini))/CLOCKS_PER_SEC << endl;
    
    sprintf(out,"time_final.par");
    saveName = fName;
    saveName.replace(fName.size()-3,3,out);
    Model->SaveParticlePosition(saveName.c_str());
    cout << "Final state: " << saveName << endl;

    sprintf(out,"com");
    saveName = fName;
    saveName.replace(fName.size()-3,3,out);
    Model->SaveCommunities(saveName.c_str());
    cout << "Community structure: " << saveName << endl;
    
}


int main(int argc,char *argv[]){
    int i=0;
    srand (time(NULL));
    
    if (argc <= 1) Message(0);
    else {
        FILE *stream;
        stream = fopen(argv[1], "r");
        if (!stream){ Message(1); return 0;}
        fclose(stream);
        fName = argv[1];
        for (i=2 ; i<argc ; i++){
            if (strcmp(argv[i],"-c") == 0){
                if (++i>=argc) break;
                numCom = strtol(argv[i],NULL, 10);
            }
            else if (strcmp(argv[i],"-max") == 0){
                if (++i>=argc) break;
                maxSteps = strtol(argv[i],NULL, 10);
            }
            else if (strcmp(argv[i],"-tr") == 0){
                if (++i>=argc) break;
                minDR = (float)strtod(argv[i],NULL);
            }
            else if (strcmp(argv[i],"-a") == 0){
                if (++i>=argc) break;
                alpha = (float)strtod(argv[i],NULL);
            }
            else if (strcmp(argv[i],"-cf") == 0){
                if (++i>=argc) break;
                fNameCom = argv[i];
            }
            else if (strcmp(argv[i],"-b") == 0){
                if (++i>=argc) break;
                beta = (float) strtod(argv[i],NULL);
            }
            else if (strcmp(argv[i],"-ss") == 0){
                saveStates = true;
                if (++i>=argc) break;
                steps = strtol(argv[i],NULL,10);
                if (steps == 0) --i;
            }
            else if (strcmp(argv[i],"-tv") == 0){
                timeVarying = true;
            }
            else if (strcmp(argv[i],"-sf") == 0){
                snapFormat = true;
            }
            else if (strcmp(argv[i],"-v") == 0){
                verbose = true;
            }
            else if (strcmp(argv[i],"-sb") == 0){
                searchBeta = true;
                if (++i>=argc) break;
                beta = (float) strtod(argv[i],NULL);
                if (++i>=argc) break;
                betaM = (float) strtod(argv[i],NULL);
                if (++i>=argc) break;
                betaS = (float) strtod(argv[i],NULL);
            }
        }
//        cout << numCom << " beta " << beta << " alpha " << alpha << endl;
        if (searchBeta) ModelSearchBeta();
        else Model0();

//        if (saveStates) ModelStep(fName,fNameCom,alpha,beta, steps);
//        else Model0(fName,fNameCom,alpha,beta);
    }

    return 0;
}
