#ifndef ChargedPFOCorrection_h
#define ChargedPFOCorrection_h 1
#include <marlin/Processor.h>
#include <marlin/Global.h>
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/LCIterator.h"
#include "UTIL/Operators.h"
#include <UTIL/BitField64.h>
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>
#include "EVENT/LCStrVec.h"
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include "IMPL/LCCollectionVec.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "lcio.h"
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TStyle.h"
#include "TCanvas.h"

// using namespace lcio ;

class ChargedPFOCorrection : public marlin::Processor{
	public:
		marlin::Processor*  newProcessor(){ return new ChargedPFOCorrection; }
		ChargedPFOCorrection();
		~ChargedPFOCorrection() = default;

		void init();
		void processEvent( EVENT::LCEvent *event );
        void end();

		int getTrackPDG(EVENT::LCEvent *event , EVENT::Track* track);
		int getTrackIndex(EVENT::LCCollection* trackCollection, EVENT::Track* track);
		TLorentzVector getTrackFourMomentum(EVENT::Track* track , double mass);
		std::vector<float> updateChargedPFOCovMat(EVENT::Track* track , double mass);
    private:
        int _nEvt;
        double _bField;
};

#endif
