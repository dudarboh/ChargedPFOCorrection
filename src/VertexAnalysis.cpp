#include "VertexAnalysis.hpp"

using namespace std;

VertexAnalysis aVertexAnalysis;

VertexAnalysis::VertexAnalysis():marlin::Processor("VertexAnalysis"){
    _description = "Write interesting parameters into root file";

    registerProcessorParameter("refit_option",
                              "If true works with verticies after track refitting",
                              _refitOpt,
                              bool(false));

}

void VertexAnalysis::init(){

    std::string filename;
    if (_refitOpt) filename = "after_refit.root";
    else filename = "before_refit.root";
    _file.reset( new TFile(filename.c_str(), "RECREATE") );
    _tree2.reset( new TTree("VertexAnalysis2", "VertexAnalysis") );
    _tree2->Branch("n_build_up_vertices", &_nBuildUpVertices);
    _tree.reset( new TTree("VertexAnalysis", "VertexAnalysis") );
    _tree->Branch("reco_vtx_pos", &_vtxRecoPos);
    _tree->Branch("mc_vtx_pos", &_vtxMcPos);
    _tree->Branch("decay_parent_type", &_decayParentType);
    _tree->Branch("n_confused_tracks", &_nConfusedTracks);
    _tree->Branch("n_missed_tracks", &_nMissedTracks);
    _tree->Branch("n_tracks", &_nTracks);
}

void VertexAnalysis::processEvent( EVENT::LCEvent* event ){
    std::cout<<"Event "<<++_nEvent<<std::endl;

    std::string vtxColName;
    if (_refitOpt) vtxColName = "BuildUpVertex_2";
    else vtxColName = "BuildUpVertex";
    LCCollection* vtxCol = event->getCollection(vtxColName);
    _nBuildUpVertices = vtxCol->getNumberOfElements();
    _tree2->Fill();

    LCRelationNavigator navRecoToMc( event->getCollection("RecoMCTruthLink") );
    XYZVector ip;
    ip.SetCoordinates( static_cast<MCParticle*> (event->getCollection("MCParticle")->getElementAt(0))->getVertex() );

    for(int i=0; i<_nBuildUpVertices; ++i){
        Vertex* recoVertex = static_cast<Vertex*> (vtxCol->getElementAt(i));
        _vtxRecoPos.SetCoordinates(recoVertex->getPosition());

        //searching for a true vertex
        std::vector<XYZVector> candidates; // true vertex candidates
        std::vector<int> candidatesWeights; // n prongs (tracks) from this vertex
        //create a vector of all possible true vertex positions from reconstructed prongs
        vector<ReconstructedParticle*> pfos = recoVertex->getAssociatedParticle()->getParticles();
        for( size_t j=0; j<pfos.size(); ++j ){
            ReconstructedParticle* pfo = pfos[j];
            ReconstructedParticle* defaultPfo;
            if (_refitOpt) defaultPfo = getDefaultPfo(event, pfo);
            else defaultPfo = pfo;
            //track weights in relations are based on old tracks/before refit! Does it have an impact?
            MCParticle* mc = getMcMaxTrackWeight(defaultPfo, navRecoToMc);
            XYZVector candidate;
            candidate.SetCoordinates( mc->getVertex() );
            bool foundInTheList = std::find(candidates.begin(), candidates.end(), candidate) != candidates.end();
            if (foundInTheList){
                int idx = std::find(candidates.begin(), candidates.end(), candidate) - candidates.begin();
                candidatesWeights[idx]++;
            }
            else{
                candidates.push_back(candidate);
                candidatesWeights.push_back(1);
            }
        }
        //take vertex which has the most attached prongs (weights) to it.
        int maxWeightIdx = std::max_element(candidatesWeights.begin(), candidatesWeights.end()) - candidatesWeights.begin();
        int maxWeight = candidatesWeights[maxWeightIdx];
        // check if there are multiple verticies with highest number of prongs
        if (std::count(candidatesWeights.begin(), candidatesWeights.end(), maxWeight) == 1){
            //only single max value -- easy case just take vertex at this index
            _vtxMcPos = candidates[maxWeightIdx];
        }
        else{
            //multiple highest weight cases. Loop and take the furthest from the IP
            double dToIp = 0.;
            for(size_t j=0; j <candidates.size(); ++j){
                if( candidatesWeights[j] == maxWeight && ( (candidates[j] - ip).r() > dToIp ) ){
                    dToIp = (candidates[j] - ip).r();
                    _vtxMcPos = candidates[j];
                }
            }
        }

        //get number of confused tracks
        //get list of reco prongs from the true vertex
        std::vector<MCParticle*> recoProngsFromTrueVertex;
        for( size_t j=0; j<pfos.size(); ++j ){
            ReconstructedParticle* pfo = pfos[j];
            ReconstructedParticle* defaultPfo;
            if (_refitOpt) defaultPfo = getDefaultPfo(event, pfo);
            else defaultPfo = pfo;
            MCParticle* mc = getMcMaxTrackWeight(defaultPfo, navRecoToMc);
            //check if this prong from the true vertex
            XYZVector pos;
            pos.SetCoordinates( mc->getVertex() );
            if (pos == _vtxMcPos) recoProngsFromTrueVertex.push_back(mc);
        }
        int nTotalRecoProngs = recoVertex->getAssociatedParticle()->getParticles().size();
        int nTrueRecoProngs = recoProngsFromTrueVertex.size();
        _nConfusedTracks = nTotalRecoProngs - nTrueRecoProngs;

        //get meson type K/D/B or just other
        MCParticle* trueProng = recoProngsFromTrueVertex[0];
        XYZVector prongVertex;
        prongVertex.SetCoordinates( trueProng->getVertex() );
        MCParticle* parent = nullptr;
        int pdg;
        if ( trueProng->getParents().size() > 0 ){                
            parent = trueProng->getParents()[0];
            XYZVector parentVertex;
            parentVertex.SetCoordinates( parent->getVertex() );
            while ( parentVertex == prongVertex && parent->getPDG() != 92){
                parent = parent->getParents()[0];
                parentVertex.SetCoordinates( parent->getVertex() );
            }
            pdg = parent->getPDG();
        }
        else{
            _decayParentType = 0;
            _nMissedTracks = 0;
            // cout<<"****** Analyzing vertex # "<<i<<" ***********"<<endl;
            // cout<<"Reconstructed position: "<<_vtxRecoPos.r()<<" mm"<<endl;
            // cout<<"True position: "<<_vtxMcPos.r()<<" mm"<<endl;
            // cout<<"D to IP: "<<(_vtxMcPos - ip).r()<<" mm"<<endl;
            // cout<<"total reconstructed prongs: "<<nTotalRecoProngs<<endl;
            // cout<<"total TRUE reconstructed prongs: "<<nTrueRecoProngs<<endl;
            // cout<<"Confusing tracks: "<<_nConfusedTracks<<endl;
            // cout<<"Decay parent: "<<_decayParentType<<endl;
            // cout<<"total TRUE prongs: "<<nTrueRecoProngs<<endl;
            // cout<<"Missed tracks: "<<_nMissedTracks<<endl;
            _tree->Fill();

            continue;
        }

        vector<int> strangeMesons = {130, 310, 311, 321, 9000311, 9000321, 10311, 10321, 100311, 100321, 9010311, 9010321, 9020311, 9020321, 313, 323, 10313, 10323, 20313, 20323, 100313, 100323, 9000313, 9000323, 30313, 30323, 315, 325, 9000315, 9000325, 10315, 10325, 20315, 20325, 9010315, 9010325, 9020315, 9020325, 317, 327, 9010317, 9010327, 319, 329, 9000319, 9000329};
        vector<int> charmMesons = {411, 421, 10411, 10421, 413, 423, 10413, 10423, 20413, 20423, 415, 425, 431, 10431, 433, 10433, 20433, 435};
        vector<int> bottomMesons = {511, 521, 10511, 10521, 513, 523, 10513, 10523, 20513, 20523, 515, 525, 531, 10531, 533, 10533, 20533, 535, 541, 10541, 543, 10543, 20543, 545};
        if( std::find(strangeMesons.begin(), strangeMesons.end(), std::abs(pdg)) != strangeMesons.end() ){
            _decayParentType = 3;
        }
        else if ( std::find(charmMesons.begin(), charmMesons.end(), std::abs(pdg)) != charmMesons.end() ){
            _decayParentType = 4;
        }
        else if ( std::find(bottomMesons.begin(), bottomMesons.end(), std::abs(pdg)) != bottomMesons.end() ){
            _decayParentType = 5;
        }
        else{
            _decayParentType = 0;
        }

        //get number of missed tracks
        int nTotalProngs = 0;
        XYZVector decayVertex;
        decayVertex.SetCoordinates( parent->getEndpoint() );
        vector<MCParticle*> decayChain;
        fillDecayChainDown(parent, decayChain);
        for(size_t j=0; j <decayChain.size(); ++j){
            MCParticle* daughter = decayChain[j];
            XYZVector daughterVertex;
            daughterVertex.SetCoordinates( daughter->getVertex() );
            if (daughterVertex  == decayVertex && daughter->getCharge() != 0 && daughter->getGeneratorStatus() <= 1) nTotalProngs++;
        }

        _nMissedTracks = nTotalProngs - nTrueRecoProngs;

        // cout<<"****** Analyzing vertex # "<<i<<" ***********"<<endl;
        // cout<<"Reconstructed position: "<<_vtxRecoPos.r()<<" mm"<<endl;
        // cout<<"True position: "<<_vtxMcPos.r()<<" mm"<<endl;
        // cout<<"D to IP: "<<(_vtxMcPos - ip).r()<<" mm"<<endl;
        // cout<<"total reconstructed prongs: "<<nTotalRecoProngs<<endl;
        // cout<<"total TRUE reconstructed prongs: "<<nTrueRecoProngs<<endl;
        // cout<<"Confusing tracks: "<<_nConfusedTracks<<endl;
        // cout<<"Decay parent: "<<_decayParentType<<endl;
        // cout<<"total TRUE prongs: "<<nTotalProngs<<endl;
        // cout<<"Missed tracks: "<<_nMissedTracks<<endl;


        _tree->Fill();
    }
    
  
}

void VertexAnalysis::end(){
    _file->Write();
}

EVENT::MCParticle* VertexAnalysis::getMcMaxTrackWeight(EVENT::ReconstructedParticle* pfo, UTIL::LCRelationNavigator nav){
    const vector<LCObject*>& mcs = nav.getRelatedToObjects(pfo);
    const vector<float>& weights = nav.getRelatedToWeights(pfo);
    //get index of highest TRACK weight MC particle
    int i = std::max_element(weights.begin(), weights.end(), [](float a, float b){return (int(a)%10000)/1000. < (int(b)%10000)/1000.;}) - weights.begin();
    return static_cast<MCParticle*> ( mcs[i] );
}

void VertexAnalysis::fillDecayChainDown(EVENT::MCParticle* mc, std::vector<EVENT::MCParticle*>& decayChain){
    // stop iterating down at particles not created by generator
    if(mc->getGeneratorStatus() == 0) return;
    decayChain.push_back(mc);
    const vector<MCParticle*> daughters = mc->getDaughters();
    for(auto daughter : daughters){
        bool foundInTheList = std::find(decayChain.begin(), decayChain.end(), daughter) != decayChain.end();
        if ( !foundInTheList ) fillDecayChainDown(daughter, decayChain);
    }
}

ReconstructedParticle* VertexAnalysis::getDefaultPfo(EVENT::LCEvent* event, EVENT::ReconstructedParticle* pfo){
    LCCollection* pfos1 = event->getCollection("PandoraPFOs");
    LCCollection* pfos2 = event->getCollection("updatedPandoraPFOs");
    for(int i=0; i < pfos2->getNumberOfElements(); ++i){
        ReconstructedParticle* pfoTmp = static_cast<ReconstructedParticle*> ( pfos2->getElementAt(i) );
        if (pfo == pfoTmp) return static_cast<ReconstructedParticle*> ( pfos1->getElementAt(i) );
    }
    return nullptr;
}
