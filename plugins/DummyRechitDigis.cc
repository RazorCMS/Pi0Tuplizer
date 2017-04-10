// -*- C++ -*-
//
// Package:    ECALlite/DummyRechitDigis
// Class:      DummyRechitDigis
// 
/**\class DummyRechitDigis DummyRechitDigis.cc ECALlite/DummyRechitDigis/plugins/DummyRechitDigis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joshua Robert Hardenbrook
//         Created:  Mon, 01 Jun 2015 06:58:24 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

//
// class declaration
//

class DummyRechitDigis : public edm::EDProducer {
public:
  explicit DummyRechitDigis(const edm::ParameterSet&);
  ~DummyRechitDigis();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
      
  // ----------member data ---------------------------
  edm::InputTag tag_barrelHitProducer_;
  edm::InputTag tag_endcapHitProducer_;
  edm::InputTag tag_preshowerHitProducer_;
  const std::string barrelRecHitCollection_;
  const std::string endcapRecHitCollection_;
  const std::string preshowerRecHitCollection_;
  edm::InputTag tag_barrelDigiProducer_;
  edm::InputTag tag_endcapDigiProducer_;
  const std::string barrelDigiCollection_;
  const std::string endcapDigiCollection_;
  const bool doDigi_;
};

DummyRechitDigis::DummyRechitDigis(const edm::ParameterSet& iConfig):
  tag_barrelHitProducer_       (iConfig.getParameter< edm::InputTag > ("barrelHitProducer")),
  tag_endcapHitProducer_       (iConfig.getParameter< edm::InputTag > ("endcapHitProducer")),
  tag_preshowerHitProducer_       (iConfig.getParameter< edm::InputTag > ("preshowerHitProducer")),
  barrelRecHitCollection_  (iConfig.getUntrackedParameter<std::string>("barrelRecHitCollection")),
  endcapRecHitCollection_  (iConfig.getUntrackedParameter<std::string>("endcapRecHitCollection")),
  preshowerRecHitCollection_  (iConfig.getUntrackedParameter<std::string>("preshowerRecHitCollection")),
  tag_barrelDigiProducer_      (iConfig.getParameter< edm::InputTag > ("barrelDigis")),
  tag_endcapDigiProducer_      (iConfig.getParameter< edm::InputTag > ("endcapDigis")),
  barrelDigiCollection_ (iConfig.getUntrackedParameter<std::string>("barrelDigiCollection")),
  endcapDigiCollection_ (iConfig.getUntrackedParameter<std::string>("endcapDigiCollection")),
  doDigi_ (iConfig.getUntrackedParameter<bool>("doDigi"))
{
  if(doDigi_) { 
    produces< EBDigiCollection >(barrelDigiCollection_);
    produces< EEDigiCollection >(endcapDigiCollection_);
  }
  else {
    if(tag_barrelHitProducer_.label() !="") produces< EcalRecHitCollection >(barrelRecHitCollection_);
    if(tag_endcapHitProducer_.label() !="") produces< EcalRecHitCollection >(endcapRecHitCollection_);
    if(tag_preshowerHitProducer_.label() !="") produces< EcalRecHitCollection >(preshowerRecHitCollection_);
  }
}

DummyRechitDigis::~DummyRechitDigis(){ }

void DummyRechitDigis::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  //std::cout << "\n-----------New Event ----------------\n " << std::endl;
   using namespace edm;   
   // build an empty collection

   // fake rechits
   // handle to try to fill
   edm::Handle<EcalRecHitCollection> barrelRecHitsHandle;
   edm::Handle<EcalRecHitCollection> endcapRecHitsHandle;
   edm::Handle<EcalRecHitCollection> preshowerRecHitsHandle;
   // dummy collection to Put()
   std::unique_ptr< EcalRecHitCollection > rechits_temp( new EcalRecHitCollection);
   std::unique_ptr< EcalRecHitCollection > rechits_temp2( new EcalRecHitCollection);
   std::unique_ptr< EcalRecHitCollection > rechits_temp3( new EcalRecHitCollection);

   //   EcalRecHit::EcalRecHit(const DetId& id, float energy, float time, uint32_t flags, uint32_t flagBits)
   // add one rechit
   EcalRecHit zero_rechit(0,0,0,0,0);
   EcalRecHitCollection zero_collection;
   zero_collection.push_back(zero_rechit);
   *rechits_temp  = zero_collection;
   *rechits_temp2 = zero_collection;
   *rechits_temp3 = zero_collection;

   std::unique_ptr< EcalRecHitCollection > rechits_eb( new EcalRecHitCollection);
   std::unique_ptr< EcalRecHitCollection > rechits_ee( new EcalRecHitCollection);
   std::unique_ptr< EcalRecHitCollection > rechits_es( new EcalRecHitCollection);

   // fake digis
   // handle to try to fill
   Handle<EBDigiCollection> digisEBHandle;
   Handle<EEDigiCollection> digisEEHandle;
   // dummy collection to Put()
   std::unique_ptr<EBDigiCollection> outputEBDigiCollection( new EBDigiCollection );
   std::unique_ptr<EEDigiCollection> outputEEDigiCollection( new EEDigiCollection );

   //Digi zero_digi;
   EBDigiCollection ebfakecol;
   EEDigiCollection eefakecol;
   //ebfakecol.push_back(zerodigi);
   //eefakecol.push_back(zerodigi);

   // fake empty collections
   std::unique_ptr<EBDigiCollection> fakeEBDigiCollection( new EBDigiCollection );
   *fakeEBDigiCollection = ebfakecol;
   std::unique_ptr<EEDigiCollection> fakeEEDigiCollection( new EEDigiCollection) ;
   *fakeEEDigiCollection = eefakecol;


   if(!doDigi_) { 
     bool foundEBRechit = true;
     bool foundEERechit = true;
     bool foundESRechit = true;

     // if you dont find the barrel rechits youre looking for, put in a fake one
     if(tag_barrelHitProducer_.label() !="")
     {
     try { 
       iEvent.getByLabel(tag_barrelHitProducer_, barrelRecHitsHandle);
       *rechits_eb = *(barrelRecHitsHandle.product());       
     }
     catch(cms::Exception& ex) { foundEBRechit = rechits_eb->size() > 0;}     
     // if you found the collection put it back into the event
     iEvent.put( foundEBRechit ? std::move(rechits_eb) : std::move(rechits_temp), barrelRecHitCollection_);
     }

     // if you dont find the endcap rechits youre looking for, put in a fake one
     if(tag_endcapHitProducer_.label() !="")
     {
     try {
       iEvent.getByLabel(tag_endcapHitProducer_, endcapRecHitsHandle);
       *rechits_ee = *(endcapRecHitsHandle.product());
     } 
     catch (cms::Exception& ex){ foundEERechit = rechits_ee->size() > 0;
     } 

     iEvent.put( foundEERechit ? std::move(rechits_ee) : std::move(rechits_temp2), endcapRecHitCollection_);
     }
     // if you dont find the preshower rechits youre looking for, put in a fake one
     if(tag_preshowerHitProducer_.label() !="")
     {
     try {
       iEvent.getByLabel(tag_preshowerHitProducer_, preshowerRecHitsHandle);
       *rechits_es = *(preshowerRecHitsHandle.product());
     } 
     catch (cms::Exception& ex){ foundESRechit = rechits_es->size() > 0;
     } 

     iEvent.put( foundESRechit ? std::move(rechits_es) : std::move(rechits_temp3), preshowerRecHitCollection_);
     }

   } // end dummy rechits
   
   // Build fake digi collections
   else {    
     bool foundEBDigi = true;
     bool foundEEDigi = true;

     if(tag_barrelDigiProducer_.label() !="")
     {
     try { // barrel digis
       iEvent.getByLabel(tag_barrelDigiProducer_, digisEBHandle);
       *outputEBDigiCollection = *(digisEBHandle.product());
     }
     catch (cms::Exception& ex) {
       foundEBDigi = outputEBDigiCollection->size() > 0;
     }     

     // insert the EB collection
     iEvent.put(foundEBDigi ? std::move(outputEBDigiCollection) : std::move(fakeEBDigiCollection), barrelDigiCollection_);
     }

     if(tag_endcapDigiProducer_.label() !="")
     {
     try { // endcap digis
       iEvent.getByLabel(tag_endcapDigiProducer_, digisEEHandle);
       *outputEEDigiCollection = *(digisEEHandle.product());
     }
     catch (cms::Exception& ex) {
       foundEEDigi = outputEEDigiCollection->size() > 0;
     }
     //     std::cout << "Putting Real EE collection?  " << foundEEDigi <<  std::endl;
     iEvent.put(foundEEDigi ? std::move(outputEEDigiCollection) : std::move(fakeEEDigiCollection), endcapDigiCollection_);
     } 
   }
   //std::cout << "-----------End Dummy Rechits ---------- " << std::endl;
}

void DummyRechitDigis::beginJob(){ }
void DummyRechitDigis::endJob(){ }
void DummyRechitDigis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DummyRechitDigis);
