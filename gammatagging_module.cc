////////////////////////////////////////////////////////////////////////
// Class:       gammatagging
// Module Type: producer 
// File:        gammatagging_module.cc
// Author: Erin Conley (eec29@duke.edu)
// Description: Creates an art::Assns between the MC truth process name
//				and the recob::Hit information in an event. Used in 
//				event display truth studies.
////////////////////////////////////////////////////////////////////////

#ifndef gammatagging_Module
#define gammatagging_Module

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "art/Framework/Core/EventSelector.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Assns.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardata/Utilities/PtrMaker.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larsim/MCCheater/BackTracker.h"

// C++ includes
#include <string>
#include <cstdlib>
#include <cmath>

//ROOT includes
#include "TVector3.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TH2.h"

#include "gammatagging.h"

namespace gammatagging {

class gammatagging: public art::EDProducer {

public:
	explicit gammatagging(fhicl::ParameterSet const &p);
	virtual ~gammatagging();
	
	void beginJob();
	void endJob();
    void produce (art::Event& evt) override;
    void reconfigure(fhicl::ParameterSet const& pset) override;	
	template <class T> std::vector<T> uniqueValues(std::vector<T> v);
    void computeBetaAndFill(std::vector<TVector3> pos, TH1D *hist);
    void sumAndFillDouble(std::vector<double> v, TH1D *hist);
    void sumAndFillTVector(std::vector<TVector3> v, TH1D *hist);
	TVector3 normalize(TVector3 v);
	int checkHit(geo::WireID hitID, std::vector<geo::WireID> hitsOnTrack);
	TVector3 makeVector(TVector3 pos, TVector3 electron);
	int findLongestTrack(std::vector<art::Ptr<recob::Track>> tracks);

private:
	std::string prim;
	std::string brem;
	
    //backtracker object
	art::ServiceHandle<cheat::BackTracker> bt;
    //labels
    std::string HitClusterLabel;
    std::string MCParticleLabel;
    std::string fDetSimProducerLabel;
    std::string HitSpacePointLabel;
    //All Truth
    TH1D* fVectorLengthBrem;
    TH1D* fVectorLengthDeex;
    TH1D* fDotProductBrem;
    TH1D* fDotProductDeex;
    TH1D* fBetaBrem;
    TH1D* fBetaDeex;
    TH1D* fVectorLengtheIoni;
    TH1D* fDotProducteIoni;
    TH1D* fBetaeIoni;
    TH1D* fVectorLengthComb;
    TH1D* fDotProductComb;
    TH1D* fBetaComb;
    //All Vectors; compare to recob to determine back/forward
    TH1D* fVectorLengthAll;
    TH1D* fVectorLengthAllB;
    TH1D* fBetaAll;
    TH1D* fBetaAllB;
    
    //Reconstruction information
    TH1D* fVectorLengthReco;
    TH1D* fDotProductReco;
    TH1D* fBetaReco;

    //Hits with tracks/no tracks study
    TH1D* fVectorLengthT;
    TH1D* fVectorLengthNT;
    TH1D* fDotProductT;
    TH1D* fDotProductNT;
    TH1D* fBetaT;
    TH1D* fBetaNT;

    
    
	
};

    
gammatagging::gammatagging(fhicl::ParameterSet const& pset)
{
    // Read in the parameters from the .fcl file.
    this->reconfigure(pset);
    //objects to add to the events
    produces< std::vector<std::string> >();
        //must include name in () to include product instance name
    produces< art::Assns<std::string, recob::Hit> >();
    
}

gammatagging::~gammatagging(){}

void gammatagging::reconfigure(fhicl::ParameterSet const& pset){
    //grab the labels from the FCL file
    HitClusterLabel = pset.get<std::string>("HitClusterLabel");
    MCParticleLabel= pset.get<std::string>("MCParticleLabel");
    fDetSimProducerLabel = pset.get< std::string >("DetSimLabel");
    HitSpacePointLabel = pset.get< std::string >("HitSpacePointLabel");
    
}

void gammatagging::beginJob(){
    //declare service handle
    art::ServiceHandle<art::TFileService> tfs;
    //strings
    prim = "primary";
    brem = "eBrem";
    //histograms
    //truth summed vector lengths
    fVectorLengthBrem = tfs->make< TH1D >("brem_truth_lengths","Normalized truth vector lengths for brem energy depositions;Vector Magnitude (cm);Number of events",60,0,1);
    fVectorLengthDeex = tfs->make< TH1D >("deex_truth_lengths","Normalized truth vector lengths for de-excitation energy depositions;Vector Magnitude (cm);Number of events",60,0,1);
    //truth all summed vectors; forward and backward
    fVectorLengthAll = tfs->make< TH1D >("all_truth_lengths", "All normalized truth vector lengths;Vector Magnitude (cm);Number of vectors", 60, 0, 1);
    fVectorLengthAllB = tfs->make< TH1D >("all_truth_lengths_B", "All normalized truth vector lengths (from end of electron track);Vector Magnitude (cm);Number of events", 60, 0, 1);
    fBetaAll = tfs->make< TH1D >("beta_all", "#beta_{14} for all vectors", 50, -1, 5);
    fBetaAllB = tfs->make< TH1D >("beta_allB", "#beta_{14} for all vectors (backwards)", 50, -1, 5);
    //truth summed dot products
    fDotProductBrem = tfs->make< TH1D >("brem_truth_dot","Sum of dot products between truth brem energy deposition vectors and electron track vector;Dot product (square cm);Number of dot products",20,-1,1);
    fDotProductDeex = tfs->make< TH1D >("deex_truth_dot","Sum of dot products between truth de-ex energy deposition vectors and electron track vector;Dot product (square cm);Number of dot products",20,-1,1);
	//isotropy condition
	fBetaBrem = tfs->make< TH1D >("brem_truth_beta", "#beta_{14} for all truth bremsstrahlung energy depositions;#beta_{14};Number of events", 100, -1, 5);
	fBetaDeex = tfs->make< TH1D >("deex_truth_beta", "#beta_{14} for all truth de-excitation energy depositions;#beta_{14};Number of events", 100, -1, 5);

	//test
	fVectorLengtheIoni = tfs->make< TH1D >("eioni_truth_lengths","Normalized truth vector lengths for ionization electrons;Vector magnitude (cm);Number of events",60,0,1);
	fDotProducteIoni = tfs->make< TH1D >("eioni_truth_dot","Sum of dot products between truth ionization electron vectors and electron track vector;Dot product (cm^{2});Number of dot products",20,-1,1);
	fBetaeIoni = tfs->make< TH1D >("eioni_truth_beta","#beta_{14} for all truth ionization electrons;#beta_{14};Number of events",100,-1,5);
	//comb
	fVectorLengthComb = tfs->make< TH1D >("comb_truth_lengths","Normalized truth vector lengths for ionization electrons and brem energy depositions;Vector magnitude (cm);Number of events",60,0,1);
	fDotProductComb = tfs->make< TH1D >("comb_truth_dot","Sum of dot products between truth ionization electron/brem depositions vectors and electron track vector;Dot product (cm^{2});Number of dot products",20,-1,1);
	fBetaComb = tfs->make< TH1D >("comb_truth_beta","#beta_{14} for all truth ionization electrons and brem energy depositions;#beta_{14};Number of events",100,-1,5);
	
	//more test
	fVectorLengthReco = tfs->make< TH1D >("length_reco", "Normalized vector Lengths for reconstructed information (HitToXYZ);Vector magnitude (cm);Number of events", 60, 0, 1);
	fDotProductReco = tfs->make< TH1D >("dot_reco", "Dot products for normalized reconstructed information (HitToXYZ);Dot products (cm^{2});Number of events", 20,-1,1);
	fBetaReco = tfs->make< TH1D >("beta_reco", "#beta_{14} for reconstructed information (HitToXYZ);#beta_{14};Number of events", 50, -1, 5);

	//track/no track study
	fVectorLengthT = tfs->make< TH1D >("reco_length_track", "Normalized vector lengths for hits associated with tracks;Vector Magnitude (cm);Number of events", 60,0,1);
	fVectorLengthNT = tfs->make< TH1D >("reco_length_noTrack", "Normalized vector lengths for hits not associated with tracks;Vector Magnitude (cm);Number of events", 60,0,1);
	fDotProductT = tfs->make< TH1D >("reco_dot_track", "Dot products for normalized vectors for hits associated with tracks;Dot product;Number of dot products",20,-1,1);
	fDotProductNT = tfs->make< TH1D >("reco_dot_noTrack", "Dot products for normalized vectors for hits not associated with tracks;Dot product;Number of dot products",20,-1,1);
	fBetaT = tfs->make< TH1D >("reco_beta_track", "#beta_{14} for hits associated with tracks;#beta_{14};Number of hits",50,-1,5);
	fBetaNT = tfs->make< TH1D >("reco_beta_noTrack", "#beta_{14} for hits not associated with tracks;#beta_{14};Number of hits",50,-1,5);	
}

void gammatagging::endJob(){}

void gammatagging::produce(art::Event& evt)
{

/////////////////////////////////////////////////////////////////////////////////
////////////////////////// CREATE DATA PRODUCT //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//create handle for the hits
auto hitListHandle = evt.getValidHandle<std::vector<recob::Hit>>(HitClusterLabel);

//declare unique_ptrs to be used in the association
std::unique_ptr<std::vector<std::string>> stringcol(new std::vector<std::string>);
std::unique_ptr<art::Assns<recob::Hit, std::string>> assn(new art::Assns<recob::Hit, std::string>);

// helper to create art pointer of the future string collection data product
lar::PtrMaker<std::string> stringPtr(evt, *this); //formerly had ", label" in the parantheses to get product instance name

auto const& backtracker = *bt; // speed up by accessing service directly

    int numBrem = 0; int numDeex = 0; int numOther = 0;
for (std::size_t iHit = 0; iHit < hitListHandle->size(); ++iHit) {
    // create a art pointer to the element of existing hit collection data product
    art::Ptr<recob::Hit> hitPtr(hitListHandle, iHit);

    std::vector<sim::TrackIDE> trackides = bt->HitToTrackID(hitPtr);
    //had to change back to bt as the other variable crashed the program
    
    //go through the tracks
    for(std::size_t i = 0; i < trackides.size(); ++i){
        auto const tid = trackides[i].trackID;
        //get truth information from the backtracker
        const simb::MCParticle* part = backtracker.TrackIDToParticle(tid);
        //fill vector of processes
        if(part->PdgCode() == 22){ //specifically for gammas
            if(prim.compare(part->Process()) != 0) stringcol->emplace_back(part->Process());
            else stringcol->emplace_back("primaryG");
 
            //keep a count of how many gammas are in the event
            if(prim.compare(part->Process()) == 0) numDeex++;
            else if(part->Process() == "eBrem") numBrem++;
            else if(prim.compare(part->Process()) != 0) numOther++; //not the same, so add to other count
		}
        else stringcol->emplace_back("NULL"); //rest of them don't matter
        
        // add an association between the hit and the string we just added
        assn->addSingle(hitPtr, stringPtr(stringcol->size() - 1));
        
    } // for track IDEs
} // for hits

    //put the data product in the event
    evt.put(std::move(stringcol)); //must include label at end to include product instance name
    evt.put(std::move(assn));
    
    //std::cout << "EVENT " << evt.id().event() << ", NUM BREMS " << numBrem << ", NUM DEEX " << numDeex << ", NUM OTHER " << numOther << std::endl;

////////////////////////////////////////////////////////////////////////
/////////////////////////// VECTOR SUM STUDY: TRUTH ////////////////////
////////////////////////////////////////////////////////////////////////

    auto particleHandle
    = evt.getValidHandle<std::vector<simb::MCParticle>>(MCParticleLabel);
    
    //delcare objects to hold positions, dot products
    std::vector<TVector3> deexPos; 
    std::vector<TVector3> bremPos; 
    std::vector<double> bremDot; 
    std::vector<double> deexDot; 
    std::vector<TVector3> ePos;
    std::vector<double> eDot;
    //hold all positions and stuff
    std::vector<TVector3> allPos;
    std::vector<TVector3> allPosB;
    //vectors to hold beginning/end of truth electron track
    TVector3 electron; TVector3 electronEnd;  
    //vectors to hold daughter IDs, IDs of brems/de-ex gammas
    std::vector<int> bremDaughtTruth; std::vector<int> deexDaughtTruth;
    std::vector<int> bremIDTruth; std::vector<int> deexIDTruth;
    std::vector<int> eDaughtTruth; std::vector<int> eIDTruth;
    std::vector<int> combDaughtTruth;
    
    //loop over particles to get the electron position
	for (simb::MCParticle const& part: *particleHandle) {   
		//std::cout << "Particle " << part.TrackId() << ", " << part.PdgCode() << ", " << part.Process() << ", " << part.NumberDaughters() << std::endl;		
				
		if(prim.compare(part.Process()) == 0){ //primary particles
            if(part.PdgCode() == 22){ //DE-EXCITATION GAMMAS	
				unsigned int numDaught = part.NumberDaughters();
				for(unsigned int i = 0; i < numDaught; ++i){
					deexDaughtTruth.push_back(part.Daughter(i));
					deexIDTruth.push_back(part.TrackId());
					//now we know the track ID of gamma mothers for every blip
				}
			}
			else if(part.PdgCode() == 11){ 
				//primary electron
				electron.SetXYZ( part.Vx(), part.Vy(), part.Vz() );
				electronEnd.SetXYZ(part.EndX(),part.EndY(),part.EndZ());
				unsigned int numDaught = part.NumberDaughters();
				for(unsigned int i = 0; i < numDaught; ++i){
					eDaughtTruth.push_back(part.Daughter(i));
					combDaughtTruth.push_back(part.Daughter(i));
					eIDTruth.push_back(part.TrackId());
				}
			}
		}
		else{ //non-primary particles
            if(brem.compare(part.Process())==0){ //BREMSSTRAHLUNG GAMMAS
				unsigned int numDaught2 = part.NumberDaughters();
				for(unsigned int i = 0; i < numDaught2; ++i){
					bremDaughtTruth.push_back(part.Daughter(i));
					combDaughtTruth.push_back(part.Daughter(i));
					bremIDTruth.push_back(part.TrackId());
				}
			}
		}                
    }//end loop over MCParticles
    
    TVector3 electronnorm(normalize(electron));
    TVector3 electronV( makeVector(electronEnd, electron) );
    TVector3 electronVnorm(normalize(electronV));
    
    //std::cout << "Event " << evt.id().event() << " electron vertex begin: (" << electron.X() << ", " << electron.Y() << ", " << electron.Z() << ")" << std::endl;
    //std::cout << "Event " << evt.id().event() << " electron vertex end: (" << electronEnd.X() << ", " << electronEnd.Y() << ", " << electronEnd.Z() << ")" << std::endl;
    
    std::vector<TVector3> combPos;
    std::vector<double> combDot;
    
    //second loop over MCParticles
    for (simb::MCParticle const& part: *particleHandle) {
        
        if(prim.compare(part.Process()) != 0){
			//all non-primary particles, which is what I want for truth
			TVector3 v(part.Vx(),part.Vy(),part.Vz());
			TVector3 tmpAll(makeVector(v,electron));
			TVector3 tmpAllB(makeVector(v,electronEnd));
			TVector3 tmpAllnorm(normalize(tmpAll));
			TVector3 tmpAllBnorm(normalize(tmpAllB));
			allPos.push_back(tmpAllnorm);			
			allPosB.push_back(tmpAllBnorm);
		}
		//electron ionization daughters
		for(unsigned int i = 0; i < eDaughtTruth.size(); ++i){
			if(part.TrackId() == eDaughtTruth[i]){
				if(part.PdgCode()==11 && prim.compare(part.Process())!=0){
					TVector3 v(part.Vx(),part.Vy(),part.Vz());
					TVector3 tmp(makeVector(v,electron));
					TVector3 tmpnorm(normalize(tmp));	
					ePos.push_back(tmpnorm);
					eDot.push_back(tmpnorm.Dot(electronVnorm));				
				}
			}
		}
		
		for(unsigned int i = 0; i < combDaughtTruth.size(); ++i){
			if(part.TrackId() == combDaughtTruth[i]){
				if(part.PdgCode() == 11 && prim.compare(part.Process()) != 0){
					TVector3 v(part.Vx(),part.Vy(),part.Vz());
					TVector3 tmp(makeVector(v,electron));
					TVector3 tmpnorm(normalize(tmp));	
					combPos.push_back(tmpnorm);
					combDot.push_back(tmpnorm.Dot(electronVnorm));						
				}
			}
		}
        
        //de-ex daughters
        for(unsigned int i = 0; i < deexDaughtTruth.size(); ++i){
			if(part.TrackId() == deexDaughtTruth[i]){
				if(part.PdgCode() == 11 && prim.compare(part.Process()) != 0){
					TVector3 v(part.Vx(),part.Vy(),part.Vz());
					TVector3 tmp(makeVector(v,electron));
					TVector3 tmpnorm(normalize(tmp));
					deexPos.push_back(tmpnorm);
					deexDot.push_back(tmpnorm.Dot(electronVnorm)); 
				}//end check for non-primary electron
			}//end check for particle being de-ex daughter
		}//end loop over all de-ex daughters

        for(unsigned int i = 0; i < bremDaughtTruth.size(); ++i){
			if(part.TrackId() == bremDaughtTruth[i]){
				if(part.PdgCode() == 11 && prim.compare(part.Process()) != 0){
					//Brem gamma blips
					TVector3 v(part.Vx(),part.Vy(),part.Vz());
					TVector3 tmp2(makeVector(v,electron));
					TVector3 tmp2norm(normalize(tmp2));     
					bremPos.push_back(tmp2norm);
					double dotprod( tmp2norm.Dot(electronVnorm) );
					bremDot.push_back(dotprod);
					
					//ePos.push_back(tmp2norm);
					//eDot.push_back(dotprod);
					
				}//end check for non-primary electron
			}//end check for particle being brem daughter
		}//end loop over all bren daughters
		   
    }//end loop over all particles
    
    if(deexPos.size() > 1) computeBetaAndFill(deexPos, fBetaDeex);
    if(bremPos.size() > 1) computeBetaAndFill(bremPos, fBetaBrem);
	if(allPos.size() > 1) computeBetaAndFill(allPos, fBetaAll);
	if(allPosB.size() > 1) computeBetaAndFill(allPosB, fBetaAllB);
	if(ePos.size() > 1) computeBetaAndFill(ePos, fBetaeIoni);
	if(combPos.size() > 1) computeBetaAndFill(combPos, fBetaComb);

	
	//if(combPos.size() > 1) computeBetaAndFill(combPos, fBetaComb);
	
	deexIDTruth = uniqueValues(deexIDTruth);
	
    //sum vectors, fill hists
    sumAndFillTVector(bremPos, fVectorLengthBrem);
    if(deexIDTruth.size() > 1) sumAndFillTVector(deexPos, fVectorLengthDeex);
    sumAndFillDouble(bremDot, fDotProductBrem);
    sumAndFillDouble(deexDot, fDotProductDeex);
    sumAndFillTVector(allPos, fVectorLengthAll);
    sumAndFillTVector(allPosB, fVectorLengthAllB);
    sumAndFillTVector(ePos, fVectorLengtheIoni);
    sumAndFillDouble(eDot, fDotProducteIoni);
    sumAndFillDouble(combDot, fDotProductComb);
    sumAndFillTVector(combPos, fVectorLengthComb);

    ////////////////////////////////////////////////////////////////////
    /////////////////// RECONSTRUCTION INFORMATION /////////////////////
    ////////////////////////////////////////////////////////////////////
	
	//grab the hit IDs, use BackTracker to get the XYZ positions
	std::vector<geo::WireID> hitWireIDs;
	std::vector<TVector3> hitXYZs;
	for(std::size_t iHit = 0; iHit < hitListHandle->size(); ++iHit){
		art::Ptr<recob::Hit> hitPtr(hitListHandle, iHit);
		hitWireIDs.push_back(hitPtr->WireID());
		std::vector<sim::TrackIDE> trackides = bt->HitToTrackID(hitPtr);
		std::vector<double> xyz; TVector3 tmp;
		for(std::size_t i = 0; i < trackides.size(); ++i) xyz = bt->HitToXYZ(hitPtr);
		if(xyz.size()>0){ 
			tmp.SetXYZ(xyz[0],xyz[1],xyz[2]);
			hitXYZs.push_back(tmp);
		}
	}
	
	//define track stuff
	art::Handle< std::vector<recob::Track> > trackcol;
	evt.getByLabel(HitSpacePointLabel, trackcol);
	art::FindMany<recob::Hit> fmh(trackcol, evt, HitSpacePointLabel);
	
	std::vector<art::Ptr<recob::Track>> tracks;
	art::fill_ptr_vector(tracks, trackcol);
	
	std::vector<geo::WireID> hitWireIDsTrack;
	std::vector<TVector3> trStart;
	std::vector<TVector3> trEnd;
	std::vector<TVector3> hitXYZtrack;
	
	for(size_t i = 0; i < tracks.size(); ++i){
		const recob::Track *tr = tracks[i].get();
		trStart.push_back(tr->Vertex());
		trEnd.push_back(tr->End());
		
		std::vector<const recob::Hit*> hits = fmh.at(i);
		
		for(size_t j = 0; j < hits.size(); ++j){
			const recob::Hit *h = hits[j];
			hitWireIDsTrack.push_back(h->WireID());
			
			for(size_t k = 0; k < hitWireIDs.size(); ++k){
				if(h->WireID()==hitWireIDs[k]){
					//associated with track
					if(hitXYZs.size()>0) hitXYZtrack.push_back(hitXYZs[k]);
					break;
				}//end check for wire IDs on track
			}//end loop over wire IDs
		}//end loop over hits
	}//end loop over tracks
	
	std::vector<TVector3> recoPos;
	std::vector<double> recoDot;
	std::vector<TVector3> recoPosB;
	std::vector<double> recoDotB;
	
	if(tracks.size() > 0){
		int longest(findLongestTrack(tracks));
		
		TVector3 eR(trStart[longest].X(), trStart[longest].Y(), trStart[longest].Z());
		TVector3 eRV(makeVector(trEnd[longest], eR));
		TVector3 eRVnorm(normalize(eRV));
		
		TVector3 eRB(trEnd[longest].X(), trEnd[longest].Y(), trEnd[longest].Z());
		TVector3 eRVB(makeVector(trStart[longest], eRB));
		TVector3 eRVBnorm(normalize(eRVB)); 
		//so this only works if there is a track AND if BackTracker worked correctly
		//only use the first track for now?
		for(size_t i = 0; i < hitXYZs.size(); ++i){
			TVector3 tmp2(makeVector(hitXYZs[i],eR));
            TVector3 tmp2norm(normalize(tmp2));
			recoPos.push_back(tmp2norm);
			//normalize
			double dotprod(tmp2norm.Dot(eRVnorm));
			recoDot.push_back(dotprod);		
			//BACKWARDS
			TVector3 tmp2b(makeVector(hitXYZs[i],eRB));	
            TVector3 tmp2bnorm(normalize(tmp2b));
			recoPosB.push_back(tmp2bnorm);
			//normalize
			double dotprodb(tmp2bnorm.Dot(eRVBnorm));
			recoDotB.push_back(dotprodb);				
		}
		
	}
	
	if(recoPos.size() > 1) computeBetaAndFill(recoPos, fBetaReco);
	sumAndFillTVector(recoPos, fVectorLengthReco);
	sumAndFillDouble(recoDot, fDotProductReco);

	std::vector<geo::WireID> hitWireIDsNotTrack;
	std::vector<TVector3> hitXYZnotTrack;
	for(size_t allHits =0; allHits < hitWireIDs.size(); ++allHits){
		geo::WireID hit1 = hitWireIDs[allHits]; 
		int check(checkHit(hit1, hitWireIDsTrack));
		if(check == 0){ 
			hitWireIDsNotTrack.push_back(hit1);
			if(hitXYZs.size()>0) hitXYZnotTrack.push_back(hitXYZs[allHits]);
		}
	}   
	    
	hitWireIDs = uniqueValues(hitWireIDs);
	hitWireIDsTrack = uniqueValues(hitWireIDsTrack);
	hitWireIDsNotTrack = uniqueValues(hitWireIDsNotTrack);
	
	std::vector<TVector3> tPos;
	std::vector<TVector3> ntPos;
	std::vector<double> tDot;
	std::vector<double> ntDot;
	
	if(tracks.size() > 0){
		int longest(findLongestTrack(tracks));
		
		TVector3 eR(trStart[longest].X(), trStart[longest].Y(), trStart[longest].Z());
		TVector3 eRV(makeVector(trEnd[longest], eR));
		TVector3 eRVnorm(normalize(eRV));
		
		//so this only works if there is a track AND if BackTracker worked correctly
		//only use the first track for now?
		for(size_t i = 0; i < hitXYZnotTrack.size(); ++i){
			TVector3 tmp2(makeVector(hitXYZnotTrack[i],eR));
            TVector3 tmp2norm(normalize(tmp2));
			ntPos.push_back(tmp2norm);
			ntDot.push_back(tmp2norm.Dot(eRVnorm));					
					
		}
		for(size_t i = 0; i < hitXYZtrack.size(); ++i){
			TVector3 tmp2(makeVector(hitXYZtrack[i],eR));
            TVector3 tmp2norm(normalize(tmp2));
			tPos.push_back(tmp2norm);
			tDot.push_back(tmp2norm.Dot(eRVnorm));				
		}
		
	}
	
	if(tPos.size() > 1) computeBetaAndFill(tPos, fBetaT);
	if(ntPos.size() > 1) computeBetaAndFill(ntPos, fBetaNT);	
	sumAndFillTVector(tPos, fVectorLengthT);
	sumAndFillTVector(ntPos, fVectorLengthNT);
	sumAndFillDouble(tDot, fDotProductT);
	sumAndFillDouble(ntDot, fDotProductNT);	
	   
    
} //END OF MAIN CODE

////////////////////////////////////////////////////////////////////////
//////////////////////BEGIN OTHER FUNCTIONS/////////////////////////////
////////////////////////////////////////////////////////////////////////

TVector3 gammatagging::makeVector(TVector3 pos, TVector3 electron){
	TVector3 v(pos.X()-electron.X(), pos.Y()-electron.Y(), pos.Z()-electron.Z());
	return v;
}

template <class T> std::vector<T> gammatagging::uniqueValues(std::vector<T> v){
	std::sort(v.begin(), v.end());
    v.erase(unique(v.begin(), v.end()), v.end());        
    return v;
}

void gammatagging::computeBetaAndFill(std::vector<TVector3> pos, TH1D *hist){
	double beta1(0.0); double beta4(0.0);
	
	double normConst(2.0/(pos.size()*(pos.size()-1)));
    	for(unsigned int i = 0; i < pos.size()-1; ++i){
			for(unsigned int j = i+1; j < pos.size(); ++j){
				TVector3 part1(pos[i].X(), pos[i].Y(), pos[i].Z());
				TVector3 part2(pos[j].X(), pos[j].Y(), pos[j].Z());
				double angle(part1.Angle(part2)); //in radians		
				double c(cos(angle));
				double leg1(c);
				double leg4(0.125*(35.0*pow(c,4)-30.0*pow(c,2)+3.0));
				beta1+=normConst*leg1;
				beta4+=normConst*leg4;					
		}
	}

	double beta14(beta1+4.0*beta4);
	//fill histogram
	hist->Fill(beta14);
}

void gammatagging::sumAndFillDouble(std::vector<double> v, TH1D *hist){
	double toFill(0.0);
    for(size_t i = 0; i < v.size(); ++i) toFill += v[i];
    //fill histogram
    if(v.size() > 0) hist->Fill(toFill/v.size());
}

void gammatagging::sumAndFillTVector(std::vector<TVector3> v, TH1D *hist){
	TVector3 toFill;
	for(size_t i = 0; i < v.size(); ++i) toFill += v[i];
	//fill histogram	
	if(v.size() > 0) hist->Fill(toFill.Mag()/v.size());
}

TVector3 gammatagging::normalize(TVector3 v){
	TVector3 norm(v.X()/v.Mag(), v.Y()/v.Mag(), v.Z()/v.Mag());
	return norm;
}

int gammatagging::checkHit(geo::WireID hitID, std::vector<geo::WireID> hitsOnTrack){
	int onTrack(0);
	for(size_t i = 0; i < hitsOnTrack.size(); ++i){
		geo::WireID h = hitsOnTrack[i];
		if(h == hitID){
			onTrack = 1;
			break;
		}
	}
	return onTrack;
}

int gammatagging::findLongestTrack(std::vector<art::Ptr<recob::Track>> tracks){
	std::vector<double> trackLengths;
	//loop through tracks, grab all the track lengths
	for(size_t i = 0; i < tracks.size(); ++i){
		const recob::Track *tr = tracks[i].get();
		trackLengths.push_back(tr->Length());
	}
	//find max track length
	auto longest = std::max_element(trackLengths.begin(), trackLengths.end());
	//put the index into an integer
	int index = std::distance(trackLengths.begin(), longest);
	
	//return the index of the longest track
	return index;
}

}  // namespace gammatagging

namespace gammatagging{
	DEFINE_ART_MODULE(gammatagging)
}



#endif // gammatagging_module

