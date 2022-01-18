#ifndef VertexAnalysis_h
#define VertexAnalysis_h 1

#include <string>
#include <vector>
#include "marlin/Processor.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/Vector3D.h"
#include "EVENT/LCObject.h"
#include "EVENT/MCParticle.h"
#include "UTIL/LCRelationNavigator.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/LCCollection.h"
#include "EVENT/Vertex.h"
#include "UTIL/Operators.h"
#include "UTIL/LCTOOLS.h"
#include "TH1F.h"


using namespace lcio;
using ROOT::Math::XYZVector;
using ROOT::Math::XYZVectorF;


class VertexAnalysis : public marlin::Processor {
    public:

        VertexAnalysis(const VertexAnalysis&) = delete;
        VertexAnalysis& operator=(const VertexAnalysis&) = delete;
        marlin::Processor* newProcessor() { return new VertexAnalysis; }

        VertexAnalysis();
        void init();
        void processEvent(EVENT::LCEvent* evt);
        void end();

        EVENT::MCParticle* getMcMaxWeight(UTIL::LCRelationNavigator pfoToMc, EVENT::ReconstructedParticle* pfo);

    private:
        int _nEvent;

        std::unique_ptr<TFile> _file;
        std::unique_ptr<TTree> _tree;

        // Track parameters at IP for single track PFOs
        // int _pdg;
        // double _omegaOld;
        // double _tanLOld;
        // double _phiOld;
        // double _d0Old;
        // double _z0Old;
        // double _omegaNew;
        // double _tanLNew;
        // double _phiNew;
        // double _d0New;
        // double _z0New;
        
        
        XYZVectorF _primVtxOld;
        XYZVectorF _primVtxNew;
        XYZVector _primVtxMc;
        int _nChargedPFOs;
        int _nPrimaryVerticesOld;
        int _nBuildUpVerticesOld;
        int _nBuildUpV0VerticesOld;
        int _nPrimaryVerticesNew;
        int _nBuildUpVerticesNew;
        int _nBuildUpV0VerticesNew;
        int _nEvt;
        TH1F* _h1;
        TH1F* _h2;
        TH1F* _h3;
        TH1F* _h4;
        TH1F* _h5;
        TH1F* _h6;
        TH1F* _h7;
        TH1F* _h8;
        TH1F* _h9;

};

#endif
