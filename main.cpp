
#include "Snap.h"
#include "ParticleNet.h"
#include <iostream>
#include <time.h>

using namespace std;

void PrintGStats(const char s[], PUNGraph Graph) {
    cout << "Graph: " << s << " nodes: " << Graph->GetNodes()
         << " edges: " << Graph->GetEdges() << endl;
}

int main(){
//    PUNGraph Graph2;
    TParticleNet *Model;
    int i;
    char file[256];
    float t_start;
  
    t_start = clock();
//    Model = new TParticleNet("./datasets/com-dblp.ungraph.txt");
//    Model = new TParticleNet("./datasets/com-youtube.ungraph.txt");
//    Model = new TParticleNet("./datasets/rede1.txt");
//    Model = new TParticleNet("./datasets/net_128_m0.05_r0.dat");
    Model = new TParticleNet("./datasets/net_128_m0.40_r0.dat");
//    Model = new TParticleNet("./datasets/network_10.dat");
    cout << "Tempo de carga [rede]: " << (double)((clock()-t_start)/CLOCKS_PER_SEC) << endl;

//    t_start = clock();
//    Model->LoadCommunities("./datasets/com-dblp.all.cmty.txt",2);
//    Model->LoadCommunities("./datasets/com-youtube.all.cmty.txt",2);
//    Model->LoadCommunities("./datasets/comu1.txt",1);
//    Model->LoadCommunities("./datasets/com_128_m0.05_r0.dat",1);
//    Model->LoadCommunities("./datasets/community_10.dat",1);
//    cout << "Tempo de carga [comunidades]: " << (double)((clock()-t_start)/CLOCKS_PER_SEC) << endl;

    t_start = clock();
    
    PrintGStats("Dados do Grafo",Model->GetNetwork());
    
    Model->SetModelParameters(1.0, 0.5, 1.0);

    t_start = (float)clock();
//    for (i=0 ; i<50 ; i++){
//        t_start = (float)clock();
//        Model->RunByStep();
////        sprintf(file,"position%d.txt",i);
////        Model->SaveParticlePosition(file);
//        cout << "Processing time: " << (clock() - t_start) / CLOCKS_PER_SEC << endl;
//    }
////    cout << "Processing time: " << (clock() - t_start) / CLOCKS_PER_SEC << endl;
    
    int itt = Model->RunModel();
    cout << "\n\nSteps: " << itt << endl;
    Model->CommunityDetection();
    Model->SaveCommunities("output.txt");

    cout << "Tempo do modelo: " << (double)((clock()-t_start)/CLOCKS_PER_SEC) << endl;

    
    
//    Model->SaveParticlePosition("position0.txt");

    
//    net = TSnap::LoadEdgeList<PNEANet>("Rede001.dat",0,1);
//    PrintGStats("ManipulateNodesEdges:Graph2",Net);
//    
//    Graph2 = TSnap::LoadEdgeList<PUNGraph>("network_100.dat",0,1);
//    Graph2 = TSnap::LoadEdgeList<PNEANet>("Rede001.dat",0,1);
    return 0;
}
