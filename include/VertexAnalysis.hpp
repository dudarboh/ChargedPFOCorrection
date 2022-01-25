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

        EVENT::MCParticle* getMcMaxTrackWeight(EVENT::ReconstructedParticle* pfo, UTIL::LCRelationNavigator nav);
        void fillDecayChainDown(EVENT::MCParticle* mc, std::vector<EVENT::MCParticle*>& decayChain);
        ReconstructedParticle* getDefaultPfo(EVENT::LCEvent* event, EVENT::ReconstructedParticle* pfo);


    private:
        int _nEvent;
        bool _refitOpt;

        std::unique_ptr<TFile> _file;
        std::unique_ptr<TTree> _tree;
        std::unique_ptr<TTree> _tree2;
        int _nBuildUpVertices;
        XYZVectorF _vtxRecoPos;
        XYZVector _vtxMcPos;
        int _decayParentType;
        int _nConfusedTracks;
        int _nMissedTracks;
        int _nTracks;

};

#endif
