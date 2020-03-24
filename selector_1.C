#define selector_1_cxx
// The class definition in selector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("selector.C")
// root> T->Process("selector.C","some options")
// root> T->Process("selector.C+")
//


#include "selector_1.h"
#include <TH2.h>
#include <TStyle.h>
#include <math.h>



void selector_1::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   file = new TFile("output_files/output.root","RECREATE");
//   file = new TFile("output_files/lxplus_output.root","RECREATE");
/*
   hist_pt_1 = new TH1F("pT", "n1==1", 100, 0., 1000.);
   hist_eta_1 = new TH1F("eta", "n1==1", 100, -5., 5.);
   hist_phi_1 = new TH1F("phi", "n1==1", 100, -4.,4.);
   hist_E_1 = new TH1F("E", "n1==1", 100, 0., 10e5);
*/
   hist_pt_2 = new TH1F("pT_n1==2", "n1==2", 100, 0., 1000.);
   hist_eta_2 = new TH1F("eta_n1==2", "n1==2", 100, -5., 5.);
   hist_phi_2 = new TH1F("phi_n1==2", "n1==2", 100, -4.,4.);
   hist_E_2 = new TH1F("E_n1==2", "n1==2", 100, 0., 10e5);

   hist_ip2d_pb = new TH1F("ip2d_pb", "pb", 100., -0.1, 1.);
   hist_ip2d_pc = new TH1F("ip2d_pc", "pc", 100., -0.1, 1.);
   hist_ip2d_pu = new TH1F("ip2d_pu", "pu", 100., -0.1, 1.);
   hist_ip2d_llr = new TH1F("ip2d_llr", "ip2d_b_jets", 100., -20., 100.);
   hist_ip3d_pb = new TH1F("ip3d_pb", "pb", 100., -0.1, 1.);
   hist_ip3d_pc = new TH1F("ip3d_pc", "pc", 100., -0.1, 1.);
   hist_ip3d_pu = new TH1F("ip3d_pu", "pu", 100., -0.1, 1.);
   hist_ip3d_llr = new TH1F("ip3d_llr", "ip3d_b_jets", 100., -20., 100.);
   hist_dl1_pb = new TH1F("dl1_pb","pb",100.,-0.1,1.);
   hist_dl1_pc = new TH1F("dl1_pc","pc",100.,-0.1,1.);
   hist_dl1_pu = new TH1F("dl1_pu","pu",100.,-0.1,1.);
   hist_dl1rnn_pb = new TH1F("dl1rnn_pb","pb",100.,-0.1,1.);
   hist_dl1rnn_pc = new TH1F("dl1rnn_pc","pc",100.,-0.1,1.);
   hist_dl1rnn_pu = new TH1F("dl1rnn_pu","pu",100.,-0.1,1.);

   hist_ip2d_llr_inB = new TH1F("ip2d_llr_inB", "ip2d_b_jets", 100., -20., 100.);
   hist_ip3d_llr_inB = new TH1F("ip3d_llr_inB", "ip3d_b_jets", 100., -20., 100.);
   hist_ip2d_llr_inC = new TH1F("ip2d_llr_inC", "ip2d_c_jets", 100., -20., 100.);
   hist_ip3d_llr_inC = new TH1F("ip3d_llr_inC", "ip3d_c_jets", 100., -20., 100.);
   hist_ip2d_llr_exB = new TH1F("ip2d_llr_exB", "ip2d_b_jets", 100., -20., 100.);
   hist_ip3d_llr_exB = new TH1F("ip3d_llr_exB", "ip3d_b_jets", 100., -20., 100.);
   hist_ip2d_llr_exC = new TH1F("ip2d_llr_exC", "ip2d_c_jets", 100., -20., 100.);
   hist_ip3d_llr_exC = new TH1F("ip3d_llr_exC", "ip3d_c_jets", 100., -20., 100.);

   hist_dl1_inC = new TH1F("DL1_inC","DL1_inC",100,-8.,12.);
   hist_dl1_inB = new TH1F("DL1_inB","DL1_inB",100,-8.,12.);
   hist_dl1_exC = new TH1F("DL1_exC","DL1_exC",100,-8.,12.);
   hist_dl1_exB = new TH1F("DL1_exB","DL1_exB",100,-8.,12.);
   /*
   hist_pt_2b = new TH1F("pT", "n2==1", 100, 0., 1000.);
   hist_eta_2b = new TH1F("eta", "n2==1", 100, -5., 5.);
   hist_phi_2b = new TH1F("phi", "n2==1", 100, -4.,4.);
   hist_E_2b = new TH1F("E", "n2==1", 100, 0., 10e5);

   hist_pt_3a = new TH1F("pT", "n1==3", 100, 0., 1000.);
   hist_eta_3a = new TH1F("eta", "n1==3", 100, -5., 5.);
   hist_phi_3a = new TH1F("phi", "n1==3", 100, -4.,4.);
   hist_E_3a = new TH1F("E", "n1==3", 100, 0., 10e5);

   hist_pt_3b = new TH1F("pT", "n3==1", 100, 0., 1000.);
   hist_eta_3b = new TH1F("eta", "n3==1", 100, -5., 5.);
   hist_phi_3b = new TH1F("phi", "n3==1", 100, -4.,4.);
   hist_E_3b = new TH1F("E", "n3==1", 100, 0., 10e5);

   hist_pt_4 = new TH1F("pT", "n1==4", 100, 0., 1000.);
   hist_eta_4 = new TH1F("eta", "n1==4", 100, -5., 5.);
   hist_phi_4 = new TH1F("phi", "n1==4", 100,  -4.,4.);
   hist_E_4 = new TH1F("E", "n1==4", 100, 0., 10e5);
*/


   TString option = GetOption();
}

void selector_1::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

Bool_t selector_1::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetLocalEntry(entry);

//   all TTreeReaderArray<int> and TTreeReaderArray<float> have the same dim, namely dim=*njets, the nuber of jets in each event.

   int nB=0,nBjets=0;
   int n1=0,n2=0,n3=0; //<- the file has max=3 b-tagged particles in a single jet
   std::vector<int> isB,is1B, is2B, is3B;

   m_Ntot++;

   for(int i=0;i<jet_nBHadr.GetSize();i++) {

     nB+=jet_nBHadr[i];
     if(jet_nBHadr[i]>0){
       isB.push_back(i);
       nBjets++;  m_nbjets++;
       if(jet_nCHadr[i]>0)  m_bc_overlap++;
     }
     if(jet_nBHadr[i]==1)  {is1B.push_back(i); n1++;}
//     if(jet_nBHadr[i]==2)  {is2B.push_back(i); n2++;}
//     if(jet_nBHadr[i]==3)  {is3B.push_back(i); n3++;}

   }

   if(nB==0)  {m_noB++;}

/*   if (n1+2*n2+3*n3!=nB) {
       std::cout << "Warning: n1+n2+n3!=nB\t" << n1<<"\t"<<n2<<"\t"<<n3<<"\t"<<nB<<"\n";
   }
*/

   int nC=0,nCjets=0;
   int nC1=0,nC2=0,nC3=0; //<- the file has max=3 b-tagged particles in a single jet
   std::vector<int> isC,is1C,is2C,is3C;

   for(int i=0;i<jet_nCHadr.GetSize();i++) {

     //     nC+=jet_nCHadr[i];
     if(jet_nCHadr[i]>0){
       isC.push_back(i);
       nCjets++;
     }
     if(jet_nCHadr[i]==1)  {is1C.push_back(i); nC1++;}
     //     if(jet_nCHadr[i]==2)  {is2C.push_back(i); nC2++;}
     //     if(jet_nCHadr[i]==3)  {is3C.push_back(i); nC3++;}
   }

double DL1=0;
//inclusive plots
  if(nCjets>0){
    for(std::vector<int>::iterator it = isC.begin(); it != isC.end(); ++it){
        if(jet_ip2d_pb[*it]!=-99){
        hist_ip2d_llr_inC->Fill(jet_ip2d_llr[*it]); //llr is computed as log(pb/pu)
//        if(jet_ip2d_llr[*it]>m_cut){ m_c2d++; }
        }
        if(jet_ip3d_pb[*it]!=-99){
          hist_ip3d_llr_inC->Fill(jet_ip3d_llr[*it]); //llr is computed as log(pb/pu)
//               if(jet_ip3d_llr[*it]>m_cut){ m_c3d++; }
        }
        if(jet_dl1_pb[*it]!=-99){
          DL1=log(jet_dl1_pb[*it]/(m_fc*jet_dl1_pc[*it]+(1-m_fc)*jet_dl1_pu[*it]));
          hist_dl1_inC->Fill(DL1);
        }
      }
  }

   if(nBjets>0){
     for(std::vector<int>::iterator it = isB.begin(); it != isB.end(); ++it){
       if(jet_ip2d_pb[*it]!=-99){
         hist_ip2d_llr_inB->Fill(jet_ip2d_llr[*it]); //llr is computed as log(pb/pu)
//      if(jet_ip2d_llr[*it]>m_cut){ m_c2d++; }
       }
       if(jet_ip3d_pb[*it]!=-99){
         hist_ip3d_llr_inB->Fill(jet_ip3d_llr[*it]); //llr is computed as log(pb/pu)
//      if(jet_ip3d_llr[*it]>m_cut){ m_c3d++; }
       }
       if(jet_dl1_pb[*it]!=-99){
         DL1=log(jet_dl1_pb[*it]/(m_fc*jet_dl1_pc[*it]+(1-m_fc)*jet_dl1_pu[*it]));
         hist_dl1_inB->Fill(DL1);
       }
     }
   }


//exclusive plots
   if(nCjets>0 && nBjets==0){
     for(std::vector<int>::iterator it = isC.begin(); it != isC.end(); ++it){
       if(jet_ip2d_pb[*it]!=-99){
         hist_ip2d_llr_exC->Fill(jet_ip2d_llr[*it]); //llr is computed as log(pb/pu)
//       if(jet_ip2d_llr[*it]>m_cut){ m_c2d++; }
         }
         if(jet_ip3d_pb[*it]!=-99){
           hist_ip3d_llr_exC->Fill(jet_ip3d_llr[*it]); //llr is computed as log(pb/pu)
//       if(jet_ip3d_llr[*it]>m_cut){ m_c3d++; }
         }
         if(jet_dl1_pb[*it]!=-99){
           DL1=log(jet_dl1_pb[*it]/(m_fc*jet_dl1_pc[*it]+(1-m_fc)*jet_dl1_pu[*it]));
           hist_dl1_exC->Fill(DL1);
         }
       }
    }

  if(nBjets>0 && nC1==0){
    for(std::vector<int>::iterator it = isB.begin(); it != isB.end(); ++it){
      if(jet_ip2d_pb[*it]!=-99){
        hist_ip2d_llr_exB->Fill(jet_ip2d_llr[*it]); //llr is computed as log(pb/pu)
//              if(jet_ip2d_llr[*it]>m_cut){ m_c2d++; }
        }
      if(jet_ip3d_pb[*it]!=-99){
        hist_ip3d_llr_exB->Fill(jet_ip3d_llr[*it]); //llr is computed as log(pb/pu)
//          if(jet_ip3d_llr[*it]>m_cut){ m_c3d++; }
      }
      if(jet_dl1_pb[*it]!=-99){
        DL1=log(jet_dl1_pb[*it]/(m_fc*jet_dl1_pc[*it]+(1-m_fc)*jet_dl1_pu[*it]));
        hist_dl1_exB->Fill(DL1);
      }
    }
   }


//select by: n1==1 && n2==0 && n3==0
   if(n1==1 && n2==0 && n3==0){
      m_b++;
/*
       for(std::vector<int>::iterator it = is1B.begin(); it != is1B.end(); ++it){
           hist_pt_1->Fill(jet_pt[*it]*0.001);
           hist_eta_1->Fill(jet_eta[*it]);
           hist_phi_1->Fill(jet_phi[*it]);
           hist_E_1->Fill(jet_E[*it]);
      }
      */
  }




//select by: n1==2 && n2==0 && n3==0
   if(n1==2 && n2==0 && n3==0){
     m_bb++;
       for(std::vector<int>::iterator it = is1B.begin(); it != is1B.end(); ++it){
           m_N++;
           hist_pt_2->Fill(jet_pt[*it]*0.001);
           hist_eta_2->Fill(jet_eta[*it]);
           hist_phi_2->Fill(jet_phi[*it]);
           hist_E_2->Fill(jet_E[*it]);
           if(jet_ip2d_pb[*it]!=-99) {
             //std::cout << "DEBUG:\tip2d_pb=" << jet_ip2d_pb[*it] << " at " << *it << " with nB=" << jet_nBHadr[*it] <<"\n";
             hist_ip2d_pb->Fill(jet_ip2d_pb[*it]);
             hist_ip2d_pc->Fill(jet_ip2d_pc[*it]);
             hist_ip2d_pu->Fill(jet_ip2d_pu[*it]);
             hist_ip2d_llr->Fill(jet_ip2d_llr[*it]); //llr is computed as log(pb/pu)

             if(jet_ip2d_llr[*it]>m_cut){
               m_b2d++;
             }
           }
           if(jet_ip3d_pb[*it]!=-99) {
             hist_ip3d_pb->Fill(jet_ip3d_pb[*it]);
             hist_ip3d_pc->Fill(jet_ip3d_pc[*it]);
             hist_ip3d_pu->Fill(jet_ip3d_pu[*it]);
             hist_ip3d_llr->Fill(jet_ip3d_llr[*it]); //llr is computed as log(pb/pu)

             if(jet_ip3d_llr[*it]>m_cut){
               m_b3d++;
             }
           }
           hist_dl1_pb->Fill(jet_dl1_pb[*it]);
           hist_dl1_pc->Fill(jet_dl1_pc[*it]);
           hist_dl1_pu->Fill(jet_dl1_pu[*it]);
           hist_dl1rnn_pb->Fill(jet_dl1rnn_pb[*it]);
           hist_dl1rnn_pc->Fill(jet_dl1rnn_pc[*it]);
           hist_dl1rnn_pu->Fill(jet_dl1rnn_pu[*it]);
//           std::cout << jet_sv1_m[*it] << "\n"; //-1000 ?
//           std::cout << jet_sv1_efc[*it] << "\n";  //-1 ?
//           std::cout << jet_sv1_n2t[*it] << "\n\n";   //-1 ?
        }
    }
/*
//select by: n1==0 && n2==1 && n3==0
   if(n1==0 && n2==1 && n3==0){
       for(std::vector<int>::iterator it = is2B.begin(); it != is2B.end(); ++it){
           hist_pt_2b->Fill(jet_pt[*it]*0.001);
           hist_eta_2b->Fill(jet_eta[*it]);
           hist_phi_2b->Fill(jet_phi[*it]);
           hist_E_2b->Fill(jet_E[*it]);
        }
    }

//select by: n1==3 && n2==0 && n3==0
   if(n1==3 && n2==0 && n3==0){
       for(std::vector<int>::iterator it = is1B.begin(); it != is1B.end(); ++it){
           hist_pt_3a->Fill(jet_pt[*it]*0.001);
           hist_eta_3a->Fill(jet_eta[*it]);
           hist_phi_3a->Fill(jet_phi[*it]);
           hist_E_3a->Fill(jet_E[*it]);
        }
   }

//select by: n1==0 && n2==0 && n3==1
   if(n1==0 && n2==0 && n3==1){
       for(std::vector<int>::iterator it = is3B.begin(); it != is3B.end(); ++it){
           hist_pt_3b->Fill(jet_pt[*it]*0.001);
           hist_eta_3b->Fill(jet_eta[*it]);
           hist_phi_3b->Fill(jet_phi[*it]);
           hist_E_3b->Fill(jet_E[*it]);
        }
   }

//select by: n1==4 && n2==0 && n3==0
   if(n1==4 && n2==0 && n3==0){
       for(std::vector<int>::iterator it = is1B.begin(); it != is1B.end(); ++it){
           hist_pt_4->Fill(jet_pt[*it]*0.001);
           hist_eta_4->Fill(jet_eta[*it]);
           hist_phi_4->Fill(jet_phi[*it]);
           hist_E_4->Fill(jet_E[*it]);
        }
   }
*/
//inspect tracks
/*
   int cnt_1=0,cnt_2=0;
   for(TTreeReaderArray<std::vector<float, std::allocator<float>>>::iterator it=jet_trk_d0.begin(); it!=jet_trk_d0.end(); ++it){
      std::vector v=*it;
      double m;
      v.size()==0 ? m=-999. : m=std::accumulate(v.begin(), v.end(), 0.0)/v.size();
//       std::cout << m << "\t" << v.size() << "\n";
      cnt_1++;
  }
*/
/*
   std::vector<int> count(*njets);
   int n=0;
   for(TTreeReaderArray<std::vector<int, std::allocator<int>>>::iterator it=jet_trk_pdg_id.begin(); it!=jet_trk_pdg_id.end(); ++it){
      std::vector v=*it;
      int k=0;
      count.at(n)=0;
      for(std::vector<int>::iterator itv=v.begin(); itv!=v.end(); ++itv){
          if(*itv==11){
              k++;
              count[n]=k;
          }
//          std::cout<< *itv << "\n";
      }
      n++;
//    std::cout<< "\n";
   }


//    std::cout << count.size() << "\t" << n << "\n";
   int i=0;
   for(std::vector<int>::iterator it=count.begin(); it!=count.end(); ++it) {
     std::cout << *it << "\t" << jet_nBHadr[i] << "\n";
     if(jet_nBHadr[i]-*it){
       std::cout << "WARNING" << "\n";
     }
     i++;
   }

  std::cout << "\n";


//   if (n1>3 || n2>1 || n3>1) std::cout << "Warning" << "\t"<< n1<<"\t"<<n2<<"\t"<<n3<<"\t"<<nB<<"\n";

//  for(int i=0;i<*njets;i++){
//    std::cout << jet_ip2d_pb[i] << "\n";




//   std::cout << "\n";

   for(TTreeReaderArray<std::vector<float, std::allocator<float>>>::iterator it=jet_trk_vtx_Z.begin(); it!=jet_trk_vtx_Z.end(); ++it){
      std::vector v=*it;
      double m=0.,sum=0.;
      int cnt=0;
      if(v.size()==0)   m=-999.;
      if(v.size()!=0) {
          for(std::vector<float>::iterator itv=v.begin(); itv!=v.end(); ++itv){
              if(*itv!=-999){
                  sum+=*itv;
                  cnt++;
              }
           }
      m=sum/cnt;
      cnt=0;
      }
      std::cout << m << "\t" << v.size() << "\n";
//        std::cout << v.size() << "\n";
      cnt_2++;
   }


//       std::cout << "\n";
//       std::cout << cnt_1 << "\t" << "\n";
//       std::cout << "\n";

*/
   //secondary vtx inspection: jet_sv1_Nvtx is a boolean vector (0,1), with jet_sv1_Nvtx.GetSize()=*njet
   //idea: get the pdgId && jet_sv1_Nvtx==1 distribution


   return kTRUE;
}

void selector_1::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void selector_1::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

//    hist_1->Draw();
//    hist_2->Draw();
/*
    hist_pt_1->Write();
    hist_eta_1->Write();
    hist_phi_1->Write();
    hist_E_1->Write();
*/

    hist_pt_2->Write();
    hist_eta_2->Write();
    hist_phi_2->Write();
    hist_E_2->Write();

    hist_ip2d_pb->Write();
    hist_ip2d_pc->Write();
    hist_ip2d_pu->Write();
    hist_ip2d_llr->Write();
    hist_ip3d_pb->Write();
    hist_ip3d_pc->Write();
    hist_ip3d_pu->Write();
    hist_ip3d_llr->Write();
    hist_dl1_pb->Write();
    hist_dl1_pc->Write();
    hist_dl1_pu->Write();
    hist_dl1rnn_pb->Write();
    hist_dl1rnn_pc->Write();
    hist_dl1rnn_pu->Write();

    hist_dl1_inC->Write();
    hist_dl1_inB->Write();
    hist_dl1_exC->Write();
    hist_dl1_exB->Write();
/*
    Double_t norm_1 = hist_ip2d_llr_inB->GetEntries();
    hist_ip2d_llr_inB->Scale(1/norm_1);
    Double_t norm_2 = hist_ip3d_llr_inB->GetEntries();
    hist_ip3d_llr_inB->Scale(1/norm_2);
    Double_t norm_3 = hist_ip2d_llr_inC->GetEntries();
    hist_ip2d_llr_inC->Scale(1/norm_3);
    Double_t norm_4 = hist_ip3d_llr_inC->GetEntries();
    hist_ip3d_llr_inC->Scale(1/norm_4);
    Double_t norm_5 = hist_ip2d_llr_exB->GetEntries();
    hist_ip2d_llr_exB->Scale(1/norm_5);
    Double_t norm_6 = hist_ip3d_llr_exB->GetEntries();
    hist_ip3d_llr_exB->Scale(1/norm_6);
    Double_t norm_7 = hist_ip2d_llr_exC->GetEntries();
    hist_ip2d_llr_exC->Scale(1/norm_7);
    Double_t norm_8 = hist_ip3d_llr_exC->GetEntries();
    hist_ip3d_llr_exC->Scale(1/norm_8);
*/
    hist_ip2d_llr_inB->Write();
    hist_ip3d_llr_inB->Write();
    hist_ip2d_llr_inC->Write();
    hist_ip3d_llr_inC->Write();
    hist_ip2d_llr_exB->Write();
    hist_ip3d_llr_exB->Write();
    hist_ip2d_llr_exC->Write();
    hist_ip3d_llr_exC->Write();

/*
    hist_pt_2b->Write();
    hist_eta_2b->Write();
    hist_phi_2b->Write();
    hist_E_2b->Write();

    hist_pt_3a->Write();
    hist_eta_3a->Write();
    hist_phi_3a->Write();
    hist_E_3a->Write();

    hist_pt_3b->Write();
    hist_eta_3b->Write();
    hist_phi_3b->Write();
    hist_E_3b->Write();

    hist_pt_4->Write();
    hist_eta_4->Write();
    hist_phi_4->Write();
    hist_E_4->Write();
*/
    file->Close();

    std::cout<< "fraction of events without b:\t" << (double) m_noB/m_Ntot << "\n";
    std::cout<< "fraction of events with one single b:\t" << (double) m_b/m_Ntot << "\n";
    std::cout<< "fraction of events with two single b's:\t" << (double) m_bb/m_Ntot << "\n";
    std::cout<< "total b-c overlap:\t" << (double) m_bc_overlap/m_nbjets << "\n";
    std::cout<< "Score ip2d with cut=" << m_cut << "\t" <<(double) m_b2d/m_N << "\n";
    std::cout<< "Score ip3d with cut=" << m_cut << "\t" <<(double) m_b3d/m_N << "\n";
}
