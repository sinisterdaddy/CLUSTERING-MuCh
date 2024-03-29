//main changes in module by module
// and ExecClusteringPeaks
//clusterformation and neighbour finding
 
#include "CbmMuchFindHitsGem.h"

#include "CbmDigiManager.h"
#include "CbmMuchCluster.h"
#include "CbmMuchGeoScheme.h"
#include "CbmMuchModule.h"
#include "CbmMuchModuleGem.h"
#include "CbmMuchPad.h"
#include "CbmMuchPixelHit.h"

#include "FairRootManager.h"

#include "TClonesArray.h"
#include "TFile.h"
#include "TMath.h"
#include "TStopwatch.h"
//#include "CbmTimeSlice.h"
#include "CbmMuchAddress.h"
#include "CbmMuchBeamTimeDigi.h"
#include "CbmMuchDigi.h"

#include <algorithm>
#include <iomanip>
#include <iostream>

using std::cout;
using std::endl;
using std::fixed;
using std::left;
using std::multimap;
using std::right;
using std::setw;
using std::vector;

// -------------------------------------------------------------------------
CbmMuchFindHitsGem::CbmMuchFindHitsGem(const char* digiFileName, Int_t flag)
  : FairTask("MuchFindHitsGem", 1)
  , fDigiFile(digiFileName)
  , fFlag(flag)
  , fAlgorithm(3)
  , fClusterSeparationTime(26.)
  , fMaxClusterSeprationTime(500.)
  , fThresholdRatio(0.1)
  , fEvent(0)
  , fNofTimeslices(0)
  ,
  //fDigis(NULL),
  fEvents(NULL)
  , fClusterCharges()
  , fLocalMax()
  , fClusterPads()
  , fNeighbours()
  , fClusters(new TClonesArray("CbmMuchCluster", 1000))
  , fHits(new TClonesArray("CbmMuchPixelHit", 1000))
  , fGeoScheme(CbmMuchGeoScheme::Instance())
  , fDigiIndices()
  , fFiredPads()
  ,
  // fDaq(),
  // fTimeSlice(NULL),
  // fDigiData(),
  fuClusters(0)
{
}

// -----   Private method Init   -------------------------------------------
InitStatus CbmMuchFindHitsGem::Init()
{
  FairRootManager* ioman = FairRootManager::Instance();
  //if (fDaq) fTimeSlice = (CbmTimeSlice*) ioman->GetObject("TimeSlice.");
  //else      fDigis     = (TClonesArray*) ioman->GetObject("MuchDigi");

  // --- Digi Manager for reading digis which were stored in vector
  fDigiManager = CbmDigiManager::Instance();
  if (bBeamTimeDigi) fDigiManager->UseMuchBeamTimeDigi();
  fDigiManager->Init();

  // fDigis will not be used now. Just for checking. Need to remove
  /*fDigis     = (TClonesArray*) ioman->GetObject("MuchDigi");
  if (! fDigis)
    fDigis     = (TClonesArray*) ioman->GetObject("MuchBeamTimeDigi");
  if (! fDigis)
    fDigis     = (TClonesArray*) ioman->GetObject("CbmMuchBeamTimeDigi");
  if (! fDigis)
    fDigis     = (TClonesArray*) ioman->GetObject("CbmMuchDigi");
  if (! fDigis)
    LOG(info) << "MuchFindHitsGem: No MuchDigi or MuchBeamTimeDigi or CbmMuchDigi or CbmMuchBeamTimeDigi exist";
    */

  // Implementation of Event by event execution after TimeSlice Mode and Event Building should have one Event branch.

  fEvents = dynamic_cast<TClonesArray*>(FairRootManager::Instance()->GetObject("Event"));
  if (!fEvents) fEvents = dynamic_cast<TClonesArray*>(FairRootManager::Instance()->GetObject("CbmEvent"));

  if (!fEvents) {
    LOG(info) << GetName() << ": No event branch present.";
    // return kfatal;
  }
  else {
    fEventMode = kTRUE;
    //fMode=kCbmEvent;
    LOG(info) << GetName() << "TimeSlice: Event-by-event mode after Event Building selected. ";
  }

  ioman->Register("MuchCluster", "Cluster in MUCH", fClusters, IsOutputBranchPersistent("MuchCluster"));
  ioman->Register("MuchPixelHit", "Hit in MUCH", fHits, IsOutputBranchPersistent("MuchPixelHit"));

  // Initialize GeoScheme
  /// Save old global file and folder pointer to avoid messing with FairRoot
  TFile* oldFile     = gFile;
  TDirectory* oldDir = gDirectory;

  TFile* file         = new TFile(fDigiFile);
  LOG_IF(fatal, !file) << "Could not open file " << fDigiFile;
  TObjArray* stations = file->Get<TObjArray>("stations");
  LOG_IF(fatal, !stations) << "TObjArray stations not found in file " << fDigiFile;
  file->Close();
  file->Delete();
  /// Restore old global file and folder pointer to avoid messing with FairRoot
  gFile      = oldFile;
  gDirectory = oldDir;

  fGeoScheme->Init(stations, fFlag);
  return kSUCCESS;
}
// -------------------------------------------------------------------------

// -----   Task execution   ------------------------------------------------
void CbmMuchFindHitsGem::Exec(Option_t*)
{
  TStopwatch timer;
  timer.Start();
  // fDigiData.clear();
  // Removing SetDaq functionality as Cluster and Hit Finder algorithm is same for both the Time Based and Event Based mode.
  //if (fDaq) ;
  //fDigiData = fTimeSlice->GetMuchData();
  // else {
  //LOG(debug)<<"Start Reading digi from a module ";

  /*for (Int_t iDigi = 0; iDigi < fDigiManager->GetNofDigis(ECbmModuleId::kMuch); iDigi++) {
    //Reading digi from CbmDigiManager which stors digis in vector
    const auto * digi;
    if(!bBeamTimeDigi) digi = (CbmMuchDigi*) fDigiManager->Get<CbmMuchDigi>(iDigi);
    else digi = (CbmMuchBeamTimeDigi*) fDigiManager->Get<CbmMuchBeamTimeDigi>(iDigi);
    //const CbmMuchDigi* digi = (CbmMuchDigi*) fDigis->At(iDigi);
    CbmMuchModule* module = fGeoScheme->GetModuleByDetId(digi->GetAddress()); //AZ
    //std::cout << module->GetDetectorType() << std::endl; //AZ
    if (module->GetDetectorType() == 2) continue; //AZ - skip 2-D straws
    fDigiData.push_back(*digi);
  }*/
  //}

  // Clear output array
  if (fHits) fHits->Clear();
  if (fClusters) fClusters->Delete();  //Clear(); // Delete because of memory escape

  fuClusters = 0;
  // --- Time-slice mode: process entire array

  if (!fEventMode) {
    //if ( fMode == kCbmTimeslice ){
    ProcessData(nullptr);
    LOG(info) << setw(20) << left << GetName() << ": processing time is " << timer.RealTime()
              << " Time Slice Number is " << fNofTimeslices << " digis "
              << fDigiManager->GetNofDigis(ECbmModuleId::kMuch)
              //<< "s digis " <<  fDigis->GetEntriesFast()
              << " clusters " << fClusters->GetEntriesFast() << " total hits " << fHits->GetEntriesFast();
  }
  // --- Event mode: loop over events
  else {
    assert(fEvents);
    Int_t nEvents = fEvents->GetEntriesFast();
    for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
      CbmEvent* event = dynamic_cast<CbmEvent*>(fEvents->At(iEvent));
      assert(event);
      Int_t nDigis =
        (event ? event->GetNofData(ECbmDataType::kMuchDigi) : fDigiManager->GetNofDigis(ECbmModuleId::kMuch));
      //: fDigis->GetEntriesFast() );
      //if (event) LOG(debug)<<" Timeslice "<< fNofTimeslices <<" event : " << event->GetNumber() <<" nDigi : " << nDigis;
      ProcessData(event);
      LOG(debug) << setw(20) << left
                 << GetName()
                 //<< ": Processing Time for an event is " << timer.RealTime()
                 << ": Time slice " << fNofTimeslices << " with " << nEvents << (nEvents == 1 ? " event" : " events")
                 << " and processing event nubmer " << iEvent << " digis "
                 << nDigis
                 //<< "s digis " <<  fDigis->GetEntriesFast()
                 << " and created cluster " << event->GetNofData(ECbmDataType::kMuchCluster) << " and created hit "
                 << event->GetNofData(ECbmDataType::kMuchPixelHit);
    }  //# events
    LOG(info) << setw(20) << left << GetName() << ": Processing Time is " << timer.RealTime() << ": Time slice "
              << fNofTimeslices << " with " << nEvents << (nEvents == 1 ? " event" : " events") << "s digis "
              << fDigiManager->GetNofDigis(ECbmModuleId::kMuch)
              //<< "s digis " <<  fDigis->GetEntriesFast()
              << " and event wise total "
              << " clusters " << fClusters->GetEntriesFast() << " total hits " << fHits->GetEntriesFast();

  }  //? event mode
  fNofTimeslices++;
}
// -------------------------------------------------------------------------

// -----   Public method Exec   --------------------------------------------
void CbmMuchFindHitsGem::ProcessData(CbmEvent* event)
{
  TStopwatch EventTimer;
  EventTimer.Start();

  fEvent++;
  //LOG(debug3)<<" Start creating cluster ";
  // Find clusters
  FindClusters(event);
  Int_t NuOfClusterInEvent = (event ? event->GetNofData(ECbmDataType::kMuchCluster) : fClusters->GetEntriesFast());

  for (Int_t clusterIndex = 0; clusterIndex < NuOfClusterInEvent; ++clusterIndex) {
    UInt_t iCluster         = (event ? event->GetIndex(ECbmDataType::kMuchCluster, clusterIndex) : clusterIndex);
    CbmMuchCluster* cluster = (CbmMuchCluster*) fClusters->At(iCluster);
    switch (fAlgorithm) {
        // Hit
      case 0:
        // One hit per pad
      case 1: {
        // One hit per cluster
        CreateHits(cluster, iCluster, event);
        break;
      }
      case 2: {
        // Simple cluster deconvolution
        ExecClusteringSimple(cluster, iCluster, event);
        break;
      }
      case 3: {
        ExecClusteringPeaks(cluster, iCluster, event);
        break;
      }
      default: {
        Fatal("CbmMuchFindHitsGem::Exec:", "The algorithm index does not exist.");
        break;
      }
    }
  }
  fDigiIndices.clear();
  fFiredPads.clear();
}
// -------------------------------------------------------------------------
// -----   Private method FindClusters  ------------------------------------
void CbmMuchFindHitsGem::FindClusters(CbmEvent* event)
{

  Int_t nDigis = (event ? event->GetNofData(ECbmDataType::kMuchDigi) : fDigiManager->GetNofDigis(ECbmModuleId::kMuch));
  //: fDigis->GetEntriesFast() );
  if (event)
    LOG(debug2) << " Timeslice " << fNofTimeslices << " event : " << event->GetNumber() << " nDigi : " << nDigis;
  if (nDigis < 0) return;
  if (fAlgorithm == 0) {
    for (Int_t iDigi = 0; iDigi < nDigis; iDigi++) {
      UInt_t digiIndex = (event ? event->GetIndex(ECbmDataType::kMuchDigi, iDigi) : iDigi);
      fDigiIndices.clear();
      fDigiIndices.push_back(digiIndex);
      //const CbmMuchDigi* digi = (CbmMuchDigi*) fDigiManager->Get<CbmMuchDigi>(digiIndex);
      const CbmMuchDigi* digi;
      if (!bBeamTimeDigi) digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchDigi>(digiIndex));
      else
        digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchBeamTimeDigi>(digiIndex));
      //const CbmMuchDigi* digi = static_cast<const CbmMuchDigi*>(fDigis->At(digiIndex));
      CbmMuchCluster* cluster = new ((*fClusters)[fuClusters++]) CbmMuchCluster();
     //cluster address will be full address of the first digi of the cluster
      cluster->SetAddress(digi->GetAddress());
      cluster->AddDigis(fDigiIndices);
      // --- In event-by-event mode after event building: register clusters to event using ECbmDataType::kMuchCluster
      //Uncomment below code
      if (event) { event->AddData(ECbmDataType::kMuchCluster, fuClusters - 1); }  //? Event object
    }
    return;
  }

  vector<CbmMuchModuleGem*> modules = fGeoScheme->GetGemModules();

  // Clear array of digis in the modules
  for (UInt_t m = 0; m < modules.size(); m++)
    modules[m]->ClearDigis();

  // Fill array of digis in the modules. Digis are automatically sorted in time

  for (Int_t iDigi = 0; iDigi < nDigis; iDigi++) {
    UInt_t digiIndex = (event ? event->GetIndex(ECbmDataType::kMuchDigi, iDigi) : iDigi);
    //const CbmMuchDigi* digi = (CbmMuchDigi*) fDigiManager->Get<CbmMuchDigi>(digiIndex);
    //const auto * digi;
    const CbmMuchDigi* digi;
    if (!bBeamTimeDigi) digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchDigi>(digiIndex));
    else
      digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchBeamTimeDigi>(digiIndex));
    //const CbmMuchDigi* digi =static_cast<const CbmMuchDigi*>(fDigis->At(digiIndex));
    Double_t time = digi->GetTime();
    //    Double_t chanid = digi->GetChannelId();
    UInt_t address = digi->GetAddress();
    //    UInt_t adc = digi->GetAdc();
    fGeoScheme->GetModuleByDetId(address)->AddDigi(time, digiIndex);
  }
   

  for (UInt_t m = 0; m < modules.size(); m++) {
    CbmMuchModuleGem* module = modules[m];
    assert(module);
    multimap<Double_t, Int_t> digis = modules[m]->GetDigis();
    multimap<Double_t, Int_t>::iterator it, itmin, itmax;
    // Split module digis into time slices according to fClusterSeparationTime
    vector<multimap<Double_t, Int_t>::iterator> slices;
    Double_t tlast = -10000;


    
    
    // Declare a set to store addresses within the current time slice
    std::set<UInt_t> addressesInCurrentSlice;

    for (it = digis.begin(); it != digis.end(); ++it) {
        Double_t t = it->first;
	// Condition 1: Create a new time slice if the time gap exceeds fClusterSeparationTime
	if (t > tlast + fClusterSeparationTime) {
		slices.push_back(it);
		tlast = t;
		continue;
	}
 	tlast = t;
      //Condition 2
      // Extract the address from the current digi
      const CbmMuchDigi* currentDigi;
      if (!bBeamTimeDigi) currentDigi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchDigi>(it->second));
      else currentDigi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchBeamTimeDigi>(it->second));
      Int_t address = currentDigi->GetAddress();
	 LOG(debug2) << "Address: " << address;
        // Check if the address is already in the set for the current time slice
        if (addressesInCurrentSlice.count(address) > 0) {
            LOG(debug2) << "Duplicate address detected in present cluster slice. Creating a new cluster slice.";
	    // Create a new cluster slice since a digi with the same address is encountered
            slices.push_back(it);
            tlast = t;
            // Clear the set for the new time slice
            addressesInCurrentSlice.clear();
            continue;
        }
	tlast = t;
	LOG(debug2) << "After if block - addressesInCurrentSlice size: " << addressesInCurrentSlice.size();
        // Add the current address to the set for the current time slice
        addressesInCurrentSlice.insert(address);
    
	//condition3
	// Check iterator validity before dereferencing
        if (it == digis.begin()) {
            // First digi, no previous digi to check
            slices.push_back(it);
            addressesInCurrentSlice.clear();
            
        } else {        

        // Calculate the time difference between the first digi and the current digi
        Double_t timeDiffFirstToCurrent = t - slices.back()->first;

        // Debug: Print current digi's time and time difference
        LOG(debug2) << "Current Digi Time: " << t ;
        LOG(debug2) << "Time Diff (First to Current): " << timeDiffFirstToCurrent ;

        if (timeDiffFirstToCurrent > fMaxClusterSeprationTime) {
            if(std::next(it) == digis.end()){ 
           	slices.push_back(it);
            	continue;
            }
            Double_t timeDiffCurrentToNext = std::next(it)->first - t;
            Double_t timeDiffCurrentToPrev = t - std::prev(it)->first;

            // Debug: Print time differences with next and previous digis
            LOG(debug2) << "Time Diff (Current to Next): " << timeDiffCurrentToNext;
            LOG(debug2) << "Time Diff (Current to Previous): " << timeDiffCurrentToPrev;

            // Condition 3.1: Create a new time slice if the time difference between the current digi
            // and the digi next to it is less than the time difference between the current digi
            // and the digi previous to it
            if (timeDiffCurrentToNext < timeDiffCurrentToPrev) {
		slices.push_back(it);
            	tlast = t;
	           // Clear the set for the new time slice
	            addressesInCurrentSlice.clear();
                LOG(debug2) << "Condition 1 Met: Creating New Slice" << std::endl;
	            continue;
                // Debug: Print a message when Condition 1 is met
            }
            // Condition 3.2: Create a new time slice from the first digi to the current digi
            else {
                // Remove the previous iterator from slices before creating the new slice
                //slices.pop_back();
                it++;
                slices.push_back(it);
                tlast = it->first;
		addressesInCurrentSlice.clear();
		addressesInCurrentSlice.insert(it->second);
                // Debug: Print a message when Condition 2 is met
                LOG(debug2) << "Condition 3.2 Met: Creating New Slice" << std::endl;
                continue;
            }
        }
	tlast = t;
	}  // ... continue processing digis ...
	
	tlast = t;
   }
  slices.push_back(it);
    for (UInt_t s = 1; s < slices.size(); s++) {
      fFiredPads.clear();
      for (it = slices[s - 1]; it != slices[s]; it++) {
        Int_t iDigi = it->second;
        //const CbmMuchDigi* digi = (CbmMuchDigi*) fDigiManager->Get<CbmMuchDigi>(iDigi);
        const CbmMuchDigi* digi;
        //const auto * digi;
        if (!bBeamTimeDigi) digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchDigi>(iDigi));
        else
          digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchBeamTimeDigi>(iDigi));
        //const CbmMuchDigi* digi = static_cast<const CbmMuchDigi*>(fDigis->At(iDigi));
        CbmMuchPad* pad = module->GetPad(digi->GetAddress());
        pad->SetDigiIndex(iDigi);
        fFiredPads.push_back(pad);
      }
      // Loop over fired pads in a time slice 
      for (UInt_t p = 0; p < fFiredPads.size(); p++) {
        fDigiIndices.clear();
        CreateCluster(fFiredPads[p]);
        if (fDigiIndices.size() == 0) continue;
        //const CbmMuchDigi* digi = static_cast<const CbmMuchDigi*>(fDigis->At(fDigiIndices.front()));
        //const CbmMuchDigi* digi = (CbmMuchDigi*) fDigiManager->Get<CbmMuchDigi>(fDigiIndices.front());
        //const auto * digi;
        const CbmMuchDigi* digi;
        if (!bBeamTimeDigi)
          digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchDigi>(fDigiIndices.front()));
        else
          digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchBeamTimeDigi>(fDigiIndices.front()));

        CbmMuchCluster* cluster = new ((*fClusters)[fuClusters++]) CbmMuchCluster();
        cluster->SetAddress(digi->GetAddress());
        cluster->AddDigis(fDigiIndices);
        // --- In event-by-event mode after event building: register clusters to event using ECbmDataType::kMuchCluster
        if (event) { event->AddData(ECbmDataType::kMuchCluster, fuClusters - 1); }  //? Event object
      }
    }
  }
}


// -------------------------------------------------------------------------


// -----   Private method CreateCluster  -----------------------------------
void CbmMuchFindHitsGem::CreateCluster(CbmMuchPad* pad)
{
  Int_t digiIndex = pad->GetDigiIndex();
  if (digiIndex < 0) return;
  fDigiIndices.push_back(digiIndex);
  pad->SetDigiIndex(-1);
  vector<CbmMuchPad*> neighbours = pad->GetNeighbours();
  for (UInt_t i = 0; i < neighbours.size(); i++)
    CreateCluster(neighbours[i]);
}
// -------------------------------------------------------------------------


// -----   Private method ExecClusteringSimple  ----------------------------
void CbmMuchFindHitsGem::ExecClusteringSimple(CbmMuchCluster* cluster, Int_t iCluster, CbmEvent* event)
{
  //const CbmMuchDigi* digi = (CbmMuchDigi*) fDigiManager->Get<CbmMuchDigi>(cluster->GetDigi(0));
  //const auto * digi;
  const CbmMuchDigi* digi;
  if (!bBeamTimeDigi) digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchDigi>(cluster->GetDigi(0)));
  else
    digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchBeamTimeDigi>(cluster->GetDigi(0)));
  //CbmMuchDigi* digi = static_cast<CbmMuchDigi*>(fDigis->At(cluster->GetDigi(0)));
  CbmMuchModule* m         = fGeoScheme->GetModuleByDetId(digi->GetAddress());
  CbmMuchModuleGem* module = (CbmMuchModuleGem*) m;
  //  Int_t iStation = CbmMuchAddress::GetStationIndex(digi->GetAddress());

  Int_t maxCharge = 0;
  for (Int_t iDigi = 0; iDigi < cluster->GetNofDigis(); iDigi++) {
    Int_t digiIndex = cluster->GetDigi(iDigi);
    //digi = (CbmMuchDigi*) fDigiManager->Get<CbmMuchDigi>(digiIndex);
    if (!bBeamTimeDigi) digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchDigi>(digiIndex));
    else
      digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchBeamTimeDigi>(digiIndex));
    //digi = static_cast<CbmMuchDigi*> (fDigis->At(digiIndex));
    Int_t charge = digi->GetAdc();
    if (charge > maxCharge) maxCharge = charge;
  }

  UInt_t threshold = UInt_t(fThresholdRatio * maxCharge);

  // Fire pads which are higher than threshold in a cluster
  fFiredPads.clear();
  for (Int_t iDigi = 0; iDigi < cluster->GetNofDigis(); iDigi++) {
    Int_t digiIndex = cluster->GetDigi(iDigi);
    //digi = (CbmMuchDigi*) fDigiManager->Get<CbmMuchDigi>(digiIndex);
    if (!bBeamTimeDigi) digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchDigi>(digiIndex));
    else
      digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchBeamTimeDigi>(digiIndex));
    //digi = static_cast<CbmMuchDigi*> (fDigis->At(digiIndex));
    if (digi->GetAdc() <= threshold) continue;
    CbmMuchPad* pad = module->GetPad(digi->GetAddress());
    pad->SetDigiIndex(digiIndex);
    fFiredPads.push_back(pad);
  }
  for (UInt_t p = 0; p < fFiredPads.size(); p++) {
    fDigiIndices.clear();
    CreateCluster(fFiredPads[p]);
    if (fDigiIndices.size() == 0) continue;
    CbmMuchCluster cl;
    cl.AddDigis(fDigiIndices);
    CreateHits(&cl, iCluster, event);
  }
}
// -------------------------------------------------------------------------


// -------------------------------------------------------------------------
void CbmMuchFindHitsGem::ExecClusteringPeaks(CbmMuchCluster* cluster, Int_t iCluster, CbmEvent* event)
{
  Int_t nDigis = cluster->GetNofDigis();
  if (nDigis <= 2) {
    CreateHits(cluster, iCluster, event);
    return;
  }
  fClusterCharges.clear();
  fClusterPads.clear();
  fLocalMax.clear();
  //  for (Int_t i=0;i<fNeighbours.size();i++) fNeighbours[i].clear();
  fNeighbours.clear();

  // Fill cluster map
  // Populating Cluster Data:
  //Iterates over the digis in the cluster and retrieves information (address, pad, adc) from each digi.
  //Fills vectors (fClusterPads, fClusterCharges, fLocalMax) with relevant data.
  for (Int_t i = 0; i < nDigis; i++) {
    Int_t iDigi = cluster->GetDigi(i);
    //const CbmMuchDigi* digi = (CbmMuchDigi*) fDigiManager->Get<CbmMuchDigi>(iDigi);
    //const auto * digi;
    const CbmMuchDigi* digi;
    if (!bBeamTimeDigi) digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchDigi>(iDigi));
    else
      digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchBeamTimeDigi>(iDigi));
    //const CbmMuchDigi* digi = static_cast<const CbmMuchDigi*>(fDigis->At(iDigi));
    UInt_t address           = digi->GetAddress();
    CbmMuchModuleGem* module = (CbmMuchModuleGem*) fGeoScheme->GetModuleByDetId(address);
    CbmMuchPad* pad          = module->GetPad(address);
    Int_t adc                = digi->GetAdc();
    fClusterPads.push_back(pad);
    fClusterCharges.push_back(adc);
    fLocalMax.push_back(1);
  }

  // Fill neighbours
  //Finding Neighbours:
  //Iterates over the pads in the cluster and identifies their neighbors.
  //Populates fNeighbours vector with indices of neighboring pads.
  for (Int_t i = 0; i < nDigis; i++) {
    CbmMuchPad* pad                = fClusterPads[i];
    vector<CbmMuchPad*> neighbours = pad->GetNeighbours();
    vector<Int_t> selected_neighbours;
    for (UInt_t ip = 0; ip < neighbours.size(); ip++) {
      CbmMuchPad* p                    = neighbours[ip];
      vector<CbmMuchPad*>::iterator it = find(fClusterPads.begin(), fClusterPads.end(), p);
      if (it == fClusterPads.end()) continue;
      selected_neighbours.push_back(it - fClusterPads.begin());
    }
    fNeighbours.push_back(selected_neighbours);
  }

  // Flag local maxima
  //Local Maxima Identification:
  //Flags local maxima in the cluster based on charge comparisons with neighboring pads.
  for (Int_t i = 0; i < nDigis; i++) {
    Int_t c = fClusterCharges[i];
    for (UInt_t n = 0; n < fNeighbours[i].size(); n++) {
      Int_t in = fNeighbours[i][n];
      Int_t cn = fClusterCharges[in];
      if (cn < c) fLocalMax[in] = 0;
    }
  }

  // Fire pads corresponding to local maxima
  fFiredPads.clear();
  for (Int_t i = 0; i < nDigis; i++) {
    CbmMuchPad* pad = fClusterPads[i];
    pad->SetDigiIndex(cluster->GetDigi(i));
    if (fLocalMax[i] == 0) continue;
    fFiredPads.push_back(pad); //fFiredPads is vector for local maximas.
  }

  // Create clusters
  //Depending on the number of local maxima, it either creates hits directly or splits the cluster and adds new clusters to fClusters.
  if(fFiredPads.size()<=1){ //for only one local maxima.
     fDigiIndices.clear();
    for (Int_t i = 0; i < nDigis; i++) {
    CbmMuchPad* pad = fClusterPads[i];
    pad->SetDigiIndex(-1);}//setting all pads as not fired.
    CreateHits(cluster, iCluster, event);
 }
 //if more than one  Local Maxima than splitting cluster accordingly and incorporating in fClusters array.
  else{	  
  Int_t NoOfSplittedClusters = 0;
  for (UInt_t p = 0; p < fFiredPads.size(); p++) {
    fDigiIndices.clear();

    CbmMuchPad* pad                = fFiredPads[p];
    fDigiIndices.push_back(pad->GetDigiIndex());
    pad->SetDigiIndex(-1);
    vector<CbmMuchPad*> neighbours = pad->GetNeighbours();
    vector<Int_t> selected_neighbours;
    for (UInt_t ip = 0; ip < neighbours.size(); ip++) {
      CbmMuchPad* np                    = neighbours[ip];
      if(np->GetDigiIndex()>0) {fDigiIndices.push_back(np->GetDigiIndex());
      pad->SetDigiIndex(-1);}//setting all pads as not fired.
    }
  
    if (fDigiIndices.size() == 0) continue;
    CbmMuchCluster cl;
    cl.AddDigis(fDigiIndices); 
    cl.SetAddress(fFiredPads[p]->GetAddress());
    //we have to add above cluster into the list of cluster
    NoOfSplittedClusters++;
  if (fClusters) {
    Int_t newIndex = iCluster + NoOfSplittedClusters;

    if (newIndex >= 0 && newIndex <= fClusters->GetEntriesFast()) {
        LOG(debug3) << "Debug: Adding cluster at index " << newIndex;
        if(NoOfSplittedClusters==1){
        fClusters->AddAt(&cl, iCluster);}
        else{
        fClusters->AddAtAndExpand(&cl,(iCluster + NoOfSplittedClusters-1));
        }
    } else {   
        LOG(debug3) << "Warning: Attempting to add at an out-of-bounds index.";
    }
}else {
    LOG(debug3) << "Error: fClusters is not a valid pointer." ;
}
 
  }

  //adding all the clusters of this event due to splitting clusters will be  added in the event
   
	  if (event) {
	  	vector <int> IndexOfClusters;
		  Int_t ClusterInEvent =  event->GetNofData(ECbmDataType::kMuchCluster);
		  for(int i =0;i<ClusterInEvent;i++){
		  IndexOfClusters.push_back(event->GetIndex(ECbmDataType::kMuchCluster,i));
	 	 }
	  
	  	for(int i=NoOfSplittedClusters-1;i>0;i--){
	  	IndexOfClusters.push_back(event->GetIndex(ECbmDataType::kMuchCluster,fClusters->GetEntriesFast()-i));
		  } 
          	event->ClearData(ECbmDataType::kMuchCluster);
	  	for(auto i =0;i<static_cast<int>(IndexOfClusters.size());i++){
	  	event->AddData(ECbmDataType::kMuchCluster,IndexOfClusters.at(i));
	  	}
	  } 

	if (fClusters) {
	    Int_t newIndex = iCluster + NoOfSplittedClusters;
	    if (newIndex >= 0 && newIndex < fClusters->GetEntriesFast()) {
	    
		CbmMuchCluster* NewFirstSplittedCluster = (CbmMuchCluster*) fClusters->At(iCluster);
		
		if (NewFirstSplittedCluster) {
		    LOG(debug3)<< "Debug: Accessing cluster at index " << iCluster;
		    //creating hits
		    CreateHits(NewFirstSplittedCluster, iCluster, event);
		    
		    LOG(debug3) << "Debug: CreateHits completed for cluster at index " << iCluster;
		} else {
		    LOG(debug3)<< "Error: Accessed cluster at index " << iCluster << " is not a valid pointer.";
		}
	    } else {
		LOG(debug3) << "Warning: Attempting to access an out-of-bounds cluster at index " << iCluster;
			}
	} else {
	    LOG(debug3) << "Error: fClusters is not a valid pointer.";

}
}
}
// -------------------------------------------------------------------------


// -----   Private method CreateHits  --------------------------------------
void CbmMuchFindHitsGem::CreateHits(CbmMuchCluster* cluster, Int_t iCluster, CbmEvent* event)
{
  Int_t nDigis  = cluster->GetNofDigis();
  if (nDigis>2) cout << "found more than 2 digid in a cluister for creatinf hits" << endl;
  Double_t sumq = 0, sumx = 0, sumy = 0, sumdx2 = 0, sumdy2 = 0, sumdxy2 = 0,
           sumdt2 = 0;  // , sumt =0 // not used FU 22.03.23
  Double_t q = 0, x = 0, y = 0, t = 0, z = 0, dx = 0, dy = 0, dxy = 0, dt = 0;
  Double_t nX = 0, nY = 0, nZ = 0;
  Int_t address            = 0;
  Int_t planeId            = 0;
  CbmMuchModuleGem* module = NULL;

  Double_t tmin = -1;

  for (Int_t i = 0; i < nDigis; i++) {
    Int_t iDigi = cluster->GetDigi(i);
    //const CbmMuchDigi* digi = (CbmMuchDigi*) fDigiManager->Get<CbmMuchDigi>(iDigi);
    //const auto * digi;
    const CbmMuchDigi* digi;
    if (!bBeamTimeDigi) digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchDigi>(iDigi));
    else
      digi = static_cast<const CbmMuchDigi*>(fDigiManager->Get<CbmMuchBeamTimeDigi>(iDigi));
    //const CbmMuchDigi* digi = static_cast<const CbmMuchDigi*>(fDigis->At(iDigi));
    if (i == 0) {
      address = CbmMuchAddress::GetElementAddress(digi->GetAddress(), kMuchModule);
      planeId = fGeoScheme->GetLayerSideNr(address);
      module  = (CbmMuchModuleGem*) fGeoScheme->GetModuleByDetId(address);
      z       = module->GetPosition()[2];
    }

    CbmMuchModule* moduleDet = fGeoScheme->GetModuleByDetId(digi->GetAddress());  /// added

    CbmMuchPad* pad = module->GetPad(digi->GetAddress());
    x               = pad->GetX();
    y               = pad->GetY();

    /// -----------------------------Drift time correction for GEM and RPC hits------------------------
    Double_t gemDriftTimeCorrc = 15.0;  // Drift time mean for GEM is 15 ns. Drift vel =100um/ns and drift width 3mm.
    Double_t rpcDriftTimeCorrc = 8.33;  // Drift time mean for RPC is 8.33 ns. Drift vel =120um/ns and drift width 2mm.
    Double_t gemHitTimeError = 4.0;  // Hit time error for GEM = 4.0 as residual dist width is 4, to make pull width 1.
    Double_t rpcHitTimeError =
      2.3;  // Hit time error for RPC = 2.3 ns, as residual dist width is 2.3. This is made to make pull dist width ~1
    Double_t timeShift = 0.5;  // this is added because residual time dist is shifted by -0.5 from 0.

    if (moduleDet->GetDetectorType() == 3)  ///GEM
    {
      if (fFlag == 0) t = digi->GetTime() - gemDriftTimeCorrc + timeShift;
      else
        t = digi->GetTime();  // Not correcting Drift Time for mCBM data
      dt = gemHitTimeError;
    }
    if (moduleDet->GetDetectorType() == 4)  ////RPC
    {
      t  = digi->GetTime() - rpcDriftTimeCorrc + timeShift;
      dt = rpcHitTimeError;
    }

    if (tmin < 0) tmin = t;
    if (tmin < t) tmin = t;
    q   = digi->GetAdc();
    dx  = pad->GetDx();
    dy  = pad->GetDy();
    dxy = pad->GetDxy();
    //dt  = 4.; // digi->GetDTime(); //TODO introduce time uncertainty determination
    sumq += q;
    sumx += q * x;
    sumy += q * y;
    // sumt += q * t; // not used FU 22.03.23
    sumdx2 += q * q * dx * dx;
    sumdy2 += q * q * dy * dy;
    sumdxy2 += q * q * dxy * dxy;
    sumdt2 += q * q * dt * dt;
    //std::cout<<" i "<<i<<" q "<<q<<" sumq "<<sumq<<std::endl;
  }

  x = sumx / sumq;
  y = sumy / sumq;
  //  t   = sumt/sumq;
  t          = tmin;
  dx         = sqrt(sumdx2 / 12) / sumq;
  dy         = sqrt(sumdy2 / 12) / sumq;
  dxy        = sqrt(sumdxy2 / 12) / sumq;
  dt         = sqrt(sumdt2) / sumq;
  Int_t iHit = fHits->GetEntriesFast();

  //------------------------------Added by O. Singh 11.12.2017 for mCbm ---------------------------
  Double_t tX = 18.5, tY = 80.0;
  nX = x + tX;  // Ajit + OS + Apar -> For miniMUCH setup in March 2019
  nY = y + tY;  // Ajit + OS + Apar -> For miniMUCH setup in March 2019
  nZ = z;

  if (fFlag == 1) {
    new ((*fHits)[iHit]) CbmMuchPixelHit(address, nX, nY, nZ, dx, dy, 0, dxy, iCluster, planeId, t, dt);  //mCbm
  }
  else {
    new ((*fHits)[iHit]) CbmMuchPixelHit(address, x, y, z, dx, dy, 0, dxy, iCluster, planeId, t, dt);  //Cbm
  }
  //Adding CbmMuchPixelHit entries in the CbmEvent
  if (event) { event->AddData(ECbmDataType::kMuchPixelHit, iHit); }  //? Event object
}
// ---------------------------------------------------------------------------------

ClassImp(CbmMuchFindHitsGem)
