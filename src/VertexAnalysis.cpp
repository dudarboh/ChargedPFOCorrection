#include "VertexAnalysis.hpp"


VertexAnalysis aVertexAnalysis;

VertexAnalysis::VertexAnalysis():marlin::Processor("VertexAnalysis"){
    _description = "Write interesting parameters into root file";
}

void VertexAnalysis::init(){
    _file.reset( new TFile("rename.root", "RECREATE") );
    _tree.reset( new TTree("VertexAnalysis", "VertexAnalysis") );

    _tree->Branch("primary_vtx_old", &_primVtxOld);
    _tree->Branch("primary_vtx_new", &_primVtxNew);
    _tree->Branch("primary_vtx_mc", &_primVtxMc);

    _h1 = new TH1F("h1", "N prim vertices,N vtx; Events", 10, 0, 10);
    _h2 = new TH1F("h2", "N Build Up vertices,N vtx; Events", 10, 0, 10);
    _h3 = new TH1F("h3", "N Build Up V0 vertices,N vtx; Events", 10, 0, 10);
    _h4 = new TH1F("h4", "N PFOs in prim vtx,N PFOs; N prim vtx", 60, 0, 60);
    _h5 = new TH1F("h5", "N PFOs in build up vtx, N PFOs; N build up vtx", 10, 0, 10);
    _h6 = new TH1F("h6", "N PFOs in build up v0 vtx, N PFOs; N build up v0 vtx", 10, 0, 10);
    _h7 = new TH1F("h7", "nPFO == nTracks for prim vtx, nPFO == nTracks; N prim vtx", 10, 0, 10);
    _h8 = new TH1F("h8", "nPFO == nTracks for build up vtx, nPFO == nTracks; N build up vtx", 10, 0, 10);
    _h9 = new TH1F("h9", "nPFO == nTracks for build up v0 vtx, nPFO == nTracks; N build up v0 vtx", 10, 0, 10);

    // _tree->Branch("nChargedPFOs", &_nChargedPFOs);
    // _tree->Branch("nPrimaryVerticesOld", &_nPrimaryVerticesOld);
    // _tree->Branch("nBuildUpVerticesOld", &_nBuildUpVerticesOld);
    // _tree->Branch("nBuildUpV0VerticesOld", &_nBuildUpV0VerticesOld);
    // 
    // _tree->Branch("nPrimaryVerticesNew", &_nPrimaryVerticesNew);
    // _tree->Branch("nBuildUpVerticesNew", &_nBuildUpVerticesNew);
    // _tree->Branch("nBuildUpV0VerticesNew", &_nBuildUpV0VerticesNew);

}

void VertexAnalysis::processEvent( EVENT::LCEvent* event ){
    std::cout<<"Event "<<++_nEvt<<std::endl;

    LCCollection* vtxCol = event->getCollection("BuildUpVertex");
    LCRelationNavigator navRecoToMc( event->getCollection("RecoMCTruthLink") );
    LCRelationNavigator navTrackToMc( event->getCollection("MarlinTrkTracksMCTruthLink") );




    for(int i=0; i<vtxCol->getNumberOfElements(); ++i){
        Vertex* vertex = static_cast<Vertex*> (vtxCol->getElementAt(i));
        std::vector <ReconstructedParticle*> rps = vertex->getAssociatedParticle()->getParticles();
        for(int j=0; j<rps.size(); ++j){
            ReconstructedParticle* rp = rps[j];
            if (rp->getTracks().size() != 1) break;
            Track* track = rp->getTracks()[0];
            // std::cout<<*track<<std::endl;

            const std::vector<LCObject*>& mcs = navRecoToMc.getRelatedToObjects( rp );
            for(int k=0; k<rps.size(); ++k){
                MCParticle* mc = static_cast<MCParticle*> (mcs[k]);
                std::cout<<*mc<<std::endl;
            }

            // Track* track = rp->getTracks()[0];
            
        }
    }




    // LCCollection* primaryVerticesNew = event->getCollection("PrimaryVertex_2");



    // LCCollection* pfos = event->getCollection("PandoraPFOs");
    // LCRelationNavigator pfoToMc( event->getCollection("RecoMCTruthLink") );
    // 
    // 
    // LCCollection* mcs = event->getCollection("MCParticle");
    // // UTIL::LCTOOLS::printMCParticles(mcs);
    // MCParticle* mc = static_cast<MCParticle*> (mcs->getElementAt(0));
    // _primVtxMc.SetCoordinates( mc->getVertex() );
    // 
    // 
    // Vertex* primVtxOld = static_cast<Vertex*> (primaryVerticesOld->getElementAt(0));
    // _primVtxOld.SetCoordinates( primVtxOld->getPosition() );
    // 
    // Vertex* primVtxNew = static_cast<Vertex*> (primaryVerticesNew->getElementAt(0));
    // _primVtxNew.SetCoordinates( primVtxNew->getPosition() );
    // 
    // _tree->Fill();
}

void VertexAnalysis::end(){
    _file->Write();
}

EVENT::MCParticle* VertexAnalysis::getMcMaxWeight(UTIL::LCRelationNavigator pfoToMc, EVENT::ReconstructedParticle* pfo){
    MCParticle* mc = nullptr;
    const std::vector <LCObject*>& mcs = pfoToMc.getRelatedToObjects(pfo);
    const std::vector <float>& mcWeights = pfoToMc.getRelatedToWeights(pfo);
    // weight //10000 = energy fraction of this MC particle in the cluster
    // weight % 10000 = number of tracker hits in per-mile of this MC particle in the track
    if (mcs.size() == 0) return mc;
    //I simply check the larger weight, meaning highest energy fraction. If two particles have equal energy fraction, choose with more track hits..
    int maxW = std::max_element(mcWeights.begin(), mcWeights.end()) - mcWeights.begin();
    mc = static_cast <MCParticle*> ( mcs[maxW] );
    return mc;
}
