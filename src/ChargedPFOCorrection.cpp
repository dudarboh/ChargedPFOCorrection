#include "ChargedPFOCorrection.hpp"
#include <iostream>
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/LCObject.h"
#include "IMPL/LCCollectionVec.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "IMPL/ParticleIDImpl.h"
#include "UTIL/PIDHandler.h"
#include "marlin/VerbosityLevels.h"
#include <GeometryUtil.h>
#include "TLorentzVector.h"

using namespace lcio;

ChargedPFOCorrection aChargedPFOCorrection;

ChargedPFOCorrection::ChargedPFOCorrection():marlin::Processor("ChargedPFOCorrection"){
    _description = "ChargedPFOCorrection creates new RECONSTRUCTEDPARTICLE collection that PFOs are updated using tracks refitted with true mass for protons and kaons";
}

void ChargedPFOCorrection::init(){
	_bField = MarlinUtil::getBzAtOrigin();
}

void ChargedPFOCorrection::processEvent( EVENT::LCEvent* event ){
	LCCollection* pfos = event->getCollection("PandoraPFOs");
	LCCollection* tracksCol = event->getCollection("MarlinTrkTracks");
	LCCollection* tracksKaonRefitCol = event->getCollection("MarlinTrkTracksKaon");
	LCCollection* tracksProtonRefitCol = event->getCollection("MarlinTrkTracksProton");
	LCCollection* mcs = event->getCollection("MCParticle");
	LCCollectionVec* outputPfosCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

	for(int i=0; i<pfos->getNumberOfElements(); ++i){
		ReconstructedParticle* pfo = static_cast<ReconstructedParticle*>( pfos->getElementAt(i) );
        int charge = pfo->getCharge();
        int pfoPDG = pfo->getType();
        const std::vector<Track*>& tracks = pfo->getTracks();
        int nTracks = tracks.size();

        //Substitute ALL Kaon/Proton tracks
        std::vector<Track*> outputTracks;
        for(int j=0; j<nTracks; ++j){
            Track* track = tracks[j];
            int pdg = getTrackPDG(event, track);
            int idx = getTrackIndex(tracksCol, track);
            if ( std::abs(pdg) == 321 ) track = static_cast<Track*> ( tracksKaonRefitCol->getElementAt( idx ) );
            else if( std::abs(pdg) == 2212 ) track = static_cast<Track*> ( tracksProtonRefitCol->getElementAt( idx ) );
            outputTracks.push_back(track);
        }

        // create new PFO substituting tracks and 4mom and cov matrix in simple cases. Everything else just copy from the old pfo
        ReconstructedParticleImpl* outputPFO = new ReconstructedParticleImpl();
        for(unsigned int j=0; j<pfo->getTracks().size(); ++j) outputPFO->addTrack(outputTracks[j]);

        if (nTracks == 1){
            Track* track = tracks[0];
            int pdg = getTrackPDG(event, track);
            if ( std::abs(pdg) == 321 || std::abs(pdg) == 2212 ){
                double mass;
                if ( std::abs(pdg) == 321 ) mass = 0.493677;
                else if ( std::abs(pdg) == 2212 ) mass = 0.938272088;
                TLorentzVector fourMomentum = getTrackFourMomentum(track, mass);
                const double momentum[3] = {fourMomentum.Px(), fourMomentum.Py(), fourMomentum.Pz()};
                std::vector<float> covMatrix = updateChargedPFOCovMat(track, mass);
                outputPFO->setType(pdg);
                outputPFO->setMomentum(momentum);
                outputPFO->setEnergy(fourMomentum.E());
                outputPFO->setCovMatrix(covMatrix);
                outputPFO->setMass(mass);
            }
            else{
                outputPFO->setType(pfo->getType());
                outputPFO->setMomentum(pfo->getMomentum());
                outputPFO->setEnergy(pfo->getEnergy());
                outputPFO->setCovMatrix(pfo->getCovMatrix());
                outputPFO->setMass(pfo->getMass());
            }
        }
        else{
            outputPFO->setType(pfo->getType());
            outputPFO->setMomentum(pfo->getMomentum());
            outputPFO->setEnergy(pfo->getEnergy());
            outputPFO->setCovMatrix(pfo->getCovMatrix());
            outputPFO->setMass(pfo->getMass());            
        }
        outputPFO->setCharge(pfo->getCharge());
        outputPFO->setReferencePoint(pfo->getReferencePoint());
        outputPFO->setParticleIDUsed(pfo->getParticleIDUsed());
        outputPFO->setGoodnessOfPID(pfo->getGoodnessOfPID());
        outputPFO->setStartVertex(pfo->getStartVertex());
        for(unsigned int j=0; j<pfo->getParticles().size(); ++j) outputPFO->addParticle(pfo->getParticles()[j]);
        for(unsigned int j=0; j<pfo->getClusters().size(); ++j) outputPFO->addCluster(pfo->getClusters()[j]);
        for(unsigned int j=0; j<pfo->getParticleIDs().size(); ++j){
            ParticleIDImpl* inPID = static_cast<ParticleIDImpl*>( pfo->getParticleIDs()[j] );
            ParticleIDImpl* outPID = new ParticleIDImpl;
            outPID->setType(inPID->getType());
            outPID->setPDG(inPID->getPDG());
            outPID->setLikelihood(inPID->getLikelihood());
            outPID->setAlgorithmType(inPID->getAlgorithmType()) ;
            for(unsigned int k=0; k<inPID->getParameters().size(); ++k) outPID->addParameter( inPID->getParameters()[k] );
            outputPFO->addParticleID( outPID );
        }

        outputPfosCol->addElement( outputPFO );
    }
    event->addCollection( outputPfosCol , "updatedPandoraPFOs" );
}

void ChargedPFOCorrection::end(){}


int ChargedPFOCorrection::getTrackIndex(EVENT::LCCollection* tracks, EVENT::Track* track){
	for (int i=0; i<tracks->getNumberOfElements(); ++i){
		Track* pionTrack = static_cast<EVENT::Track*>( tracks->getElementAt(i) );
		if ( pionTrack == track ) return i;
    }
    return -1;
}


int ChargedPFOCorrection::getTrackPDG(EVENT::LCEvent* event, EVENT::Track* track){
    LCRelationNavigator trackToMcNav( event->getCollection("MarlinTrkTracksMCTruthLink") );
    const std::vector<LCObject*>& mcParticles = trackToMcNav.getRelatedToObjects(track);
    const std::vector<float>&  mcWeights = trackToMcNav.getRelatedToWeights(track);

    int nMCParticles = mcParticles.size();
    for(int i=0; i<nMCParticles; ++i){
        MCParticle* mc = static_cast<MCParticle*> (mcParticles[i]);
        
    }
    //can't assign PDG if no MC particles linked
    if( nMCParticles == 0 ) return 0;

    int mcIdx = std::max_element( mcWeights.begin(), mcWeights.end() ) - mcWeights.begin();
	MCParticle* mc = static_cast<MCParticle*> (mcParticles[mcIdx]);

    return mc->getPDG();
}


TLorentzVector ChargedPFOCorrection::getTrackFourMomentum(EVENT::Track* track, double mass){
	double pt = 0.299792458e-3 * _bField / std::abs( track->getOmega() );
	double px = pt*std::cos( track->getPhi() );
	double py = pt*std::sin( track->getPhi() );
	double pz = pt*track->getTanLambda();
	double energy = std::sqrt(mass*mass + px*px + py*py + pz*pz);
	return TLorentzVector(px, py, pz, energy);
}


std::vector<float> ChargedPFOCorrection::updateChargedPFOCovMat(EVENT::Track* track, double mass){
    //	Obtain covariance matrix on (px,py,pz,E) from the covariance matrix on helix parameters.
    //			Dpx/DTanLambda		Dpy/DTanLambda		Dpz/DTanLambda		DE/DTanLambda
    //	 J =	Dpx/DPhi		      Dpy/DPhi		      Dpz/DPhi		      DE/DPhi
    //			Dpx/DOmega		     Dpy/DOmega		     Dpz/DOmega		     DE/DOmega
    //
    //			0			         0			     Pt	          Pz.Pt/E
    //	 J =	-Py			        Px			     0			    0
    //			-Px/Omega		-Py/Omega		-Pz/Omega		-P2/E.Omega
    //
    //	Order in the covariance matrix on helix parameters:
    //			TanLambda.TanLambda	       TanLambda.phi		   TanLambda.Omega
    //	Cov =	phi.TanLambda		         phi.phi			      phi.Omega
    //			Omega.TanLambda		        Omega.phi		         Omega.Omega
	int rows = 3; // n rows jacobian
	int columns	= 4; // n columns jacobian

	double omega = track->getOmega();
    double pt = 0.299792458e-3 * _bField / std::abs( track->getOmega() );
	double px = pt*std::cos( track->getPhi() );
	double py = pt*std::sin( track->getPhi() );
    double pz = pt*track->getTanLambda();
    double p = std::sqrt(pt*pt + pz*pz);

	double energy = std::sqrt(mass*mass + p*p);

    double jacobianByRows[rows*columns] = {
        0, 0, pt, pz*pt/energy,
		-py, px, 0,	0,
		-px/omega, -py/omega, -pz/omega, -p*p/(energy*omega)
    };
    // Create Jacobian ("F" fill by columns, "C" fill by rows)
    TMatrixD jacobian(rows, columns, jacobianByRows, "C");

    // we ignore d0, z0 parameters and using only phi, omega, tanL
    std::vector<float> flatTrackCovMat = track->getCovMatrix();
	double trackCovMatByRows[rows*rows] = {
        flatTrackCovMat[14], flatTrackCovMat[11], flatTrackCovMat[12],
		flatTrackCovMat[11], flatTrackCovMat[2], flatTrackCovMat[4],
		flatTrackCovMat[12], flatTrackCovMat[4], flatTrackCovMat[5]
	};
	TMatrixD trackCovMat(rows, rows, trackCovMatByRows, "C");

    //4momentum cov matrix and its diagonal 10 elements flattened into array
	TMatrixD covMat4Mom(4, 4);
    covMat4Mom.Mult( TMatrixD(jacobian, TMatrixD::kTransposeMult, trackCovMat), jacobian );

	std::vector<float> flatCovMat4Mom;
	flatCovMat4Mom.push_back( covMat4Mom(0,0) ); // x-x
	flatCovMat4Mom.push_back( covMat4Mom(1,0) ); // y-x
	flatCovMat4Mom.push_back( covMat4Mom(1,1) ); // y-y
	flatCovMat4Mom.push_back( covMat4Mom(2,0) ); // z-x
	flatCovMat4Mom.push_back( covMat4Mom(2,1) ); // z-y
	flatCovMat4Mom.push_back( covMat4Mom(2,2) ); // z-z
	flatCovMat4Mom.push_back( covMat4Mom(3,0) ); // e-x
	flatCovMat4Mom.push_back( covMat4Mom(3,1) ); // e-y
	flatCovMat4Mom.push_back( covMat4Mom(3,2) ); // e-z
	flatCovMat4Mom.push_back( covMat4Mom(3,3) ); // e-e

	return flatCovMat4Mom;
}
