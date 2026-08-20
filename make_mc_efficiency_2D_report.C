#include <fstream>
#include <iomanip>
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
#include <Math/Vector4D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using namespace ROOT::VecOps;

const double MUON_MASS = 0.1056583745;

void make_mc_efficiency_2D_report() {
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("0.6f");
    gStyle->SetPalette(kRainbow);

     ROOT::EnableImplicitMT(8);

    TString fileName = "/eos/user/y/yingyinw/ntuple_practice/CMSSW_13_0_13/src/UserCode/MultiLepPAT/test/Run_Crab_Script/crab_Official_MC_JpsiPT5_MultiLepPAT_v4/results/mymultilep_*.root";

    ROOT::RDataFrame df("mkcands/X_data", fileName);

    auto df_with_gen = df.Filter("Offical_Pythia_J_MC_xPx.size()>0");

    std::vector<double> yBins = {0.0, 0.3, 0.6, 0.9, 1.2, 1.8, 2.4};
    std::vector<double> ptBins = {
        10.0, 10.5, 11.0, 11.5, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0,
        18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0,
        28.0, 29.0, 30.0, 32.0, 34.0, 36.0, 38.0, 42.0, 46.0, 50.0,
        60.0, 75.0, 95.0, 120.0
    };
    int nYBins = yBins.size() - 1;
    int nPtBins = ptBins.size() - 1;

    // ---- 分母：所有候选，无 cuts ----
    auto df_den = df_with_gen
        .Define("denPt", [](const RVec<float>& px, const RVec<float>& py) {
            RVec<float> pts;
            for (size_t i = 0; i < px.size(); ++i) {
                pts.push_back(static_cast<float>(sqrt(px[i]*px[i] + py[i]*py[i])));
            }
            return pts;
        }, {"mumuonlyPx", "mumuonlyPy"})
        .Define("denAbsY", [](const RVec<float>& px, const RVec<float>& py,
                              const RVec<float>& pz, const RVec<float>& mass) {
            RVec<float> ys;
            for (size_t i = 0; i < px.size(); ++i) {
                ROOT::Math::PxPyPzMVector p4(px[i], py[i], pz[i], mass[i]);
                ys.push_back(static_cast<float>(fabs(p4.Rapidity())));
            }
            return ys;
        }, {"mumuonlyPx", "mumuonlyPy", "mumuonlyPz", "mumuonlyMass"});

    auto hDen = df_den.Histo2D(
        {"hDen", "Reconstructed Candidates (with gen truth);|y|;p_{T} [GeV]",
         nYBins, yBins.data(), nPtBins, ptBins.data()},
        "denAbsY", "denPt"
    );

    // ---- 分子：从分母派生，展开候选，逐步 cuts ----
    auto df_idx = df_den
        .Define("cand_idx", [](unsigned int n) {
            std::vector<unsigned int> idx(n);
            for (unsigned int i=0; i<n; ++i) idx[i]=i;
            return idx;
        }, {"nmumuonly"})
        .Define("mu1Idx", [](const RVec<float>& mu1Idx) {
            return mu1Idx.empty() ? -1 : (int)mu1Idx[0];
        }, {"mumuonlymu1Idx"})
        .Define("mu2Idx", [](const RVec<float>& mu2Idx) {
            return mu2Idx.empty() ? -1 : (int)mu2Idx[0];
        }, {"mumuonlymu2Idx"})
        .Define("vtxCL", [](const RVec<float>& vtxCL) {
            return vtxCL.empty() ? 0.0f : vtxCL[0];
        }, {"mumuonlyVtxCL"})
        .Define("ctau", [](const RVec<float>& ctau) {
            return ctau.empty() ? 0.0f : ctau[0];
        }, {"mumuonlyctau"})
        .Filter([](int idx) { return idx >= 0; }, {"mu1Idx"})
        // 提取 muon 属性
        .Define("mu1Px", [](int i, const RVec<float>& px) { return px[i]; }, {"mu1Idx", "muPx"})
        .Define("mu1Py", [](int i, const RVec<float>& py) { return py[i]; }, {"mu1Idx", "muPy"})
        .Define("mu1Pz", [](int i, const RVec<float>& pz) { return pz[i]; }, {"mu1Idx", "muPz"})
        .Define("mu1Charge", [](int i, const RVec<float>& charge) { return charge[i]; }, {"mu1Idx", "muCharge"})
        .Define("mu1Soft", [](int i, const RVec<int>& soft) { return soft[i]; }, {"mu1Idx", "muIsPatSoftMuon"})
        .Define("mu1TrigMatch", [](int i, const RVec<int>& match) { return match[i]; }, {"mu1Idx", "muIsJpsiTrigMatch"})
        .Define("mu2Px", [](int i, const RVec<float>& px) { return px[i]; }, {"mu2Idx", "muPx"})
        .Define("mu2Py", [](int i, const RVec<float>& py) { return py[i]; }, {"mu2Idx", "muPy"})
        .Define("mu2Pz", [](int i, const RVec<float>& pz) { return pz[i]; }, {"mu2Idx", "muPz"})
        .Define("mu2Charge", [](int i, const RVec<float>& charge) { return charge[i]; }, {"mu2Idx", "muCharge"})
        .Define("mu2Soft", [](int i, const RVec<int>& soft) { return soft[i]; }, {"mu2Idx", "muIsPatSoftMuon"})
        .Define("mu2TrigMatch", [](int i, const RVec<int>& match) { return match[i]; }, {"mu2Idx", "muIsJpsiTrigMatch"})
        .Define("mu1Pt", [](float px, float py) -> float { return static_cast<float>(sqrt(px*px + py*py)); }, {"mu1Px","mu1Py"})
        .Define("mu2Pt", [](float px, float py) -> float { return static_cast<float>(sqrt(px*px + py*py)); }, {"mu2Px","mu2Py"})
        .Define("mu1Eta", [](float px, float py, float pz) -> float {
            ROOT::Math::PxPyPzMVector p(px, py, pz, MUON_MASS);
            return static_cast<float>(p.Eta());
        }, {"mu1Px","mu1Py","mu1Pz"})
        .Define("mu2Eta", [](float px, float py, float pz) -> float {
            ROOT::Math::PxPyPzMVector p(px, py, pz, MUON_MASS);
            return static_cast<float>(p.Eta());
        }, {"mu2Px","mu2Py","mu2Pz"})
        .Define("mu1Phi", [](float px, float py) -> float { return static_cast<float>(atan2(py, px)); }, {"mu1Px","mu1Py"})
        .Define("mu2Phi", [](float px, float py) -> float { return static_cast<float>(atan2(py, px)); }, {"mu2Px","mu2Py"});

    // 构建 dimuon
    auto df_dimu = df_idx
        .Define("dimuP4", [](float px1, float py1, float pz1, float px2, float py2, float pz2) {
            ROOT::Math::PxPyPzMVector m1(px1, py1, pz1, MUON_MASS);
            ROOT::Math::PxPyPzMVector m2(px2, py2, pz2, MUON_MASS);
            return m1 + m2;
        }, {"mu1Px","mu1Py","mu1Pz","mu2Px","mu2Py","mu2Pz"})
        .Define("dimuPt", [](const ROOT::Math::PxPyPzMVector& p) -> float { return static_cast<float>(p.pt()); }, {"dimuP4"})
        .Define("dimuAbsY", [](const ROOT::Math::PxPyPzMVector& p) -> float { return static_cast<float>(fabs(p.Rapidity())); }, {"dimuP4"})
        .Define("dimuMass", [](const ROOT::Math::PxPyPzMVector& p) -> float { return static_cast<float>(p.M()); }, {"dimuP4"})
        .Define("dimuDphi", [](float phi1, float phi2) -> float {
            float dphi = phi1 - phi2;
            if (dphi > M_PI) dphi -= 2*M_PI;
            if (dphi < -M_PI) dphi += 2*M_PI;
            return dphi;
        }, {"mu1Phi","mu2Phi"})
        .Define("passTrigger", [](const RVec<std::string>& names, const RVec<unsigned int>& res) {
            for (size_t t=0; t<names.size(); ++t) {
                if (names[t].find("HLT_DoubleMu4_3_LowMass_v") != std::string::npos && res[t]==1) return true;
            }
            return false;
        }, {"TrigNames","TrigRes"})
        .Define("cut_kin", [](float pt1, float pt2, float eta1, float eta2) {
            return (pt1 > 5.0f && fabs(eta1) < 2.4f && pt2 > 5.0f && fabs(eta2) < 2.4f);
        }, {"mu1Pt","mu2Pt","mu1Eta","mu2Eta"})
        .Define("cut_charge", [](float q1, float q2) { return q1 * q2 < 0; }, {"mu1Charge","mu2Charge"})
        .Define("cut_mass", [](float mass) { return mass >= 2.5f && mass <= 3.5f; }, {"dimuMass"})
        .Define("cut_soft", [](int s1, int s2) { return s1 == 1 && s2 == 1; }, {"mu1Soft","mu2Soft"})
        .Define("cut_vtx", [](float vtx) { return vtx > 0.01f; }, {"vtxCL"})
        .Define("cut_dphi", [](float dphi) { return dphi > 0; }, {"dimuDphi"})
        .Define("cut_match", [](int m1, int m2) { return m1 == 1 && m2 == 1; }, {"mu1TrigMatch","mu2TrigMatch"})
        .Define("cut_ctau", [](float ctau) { return ctau >= -0.03f && ctau <= 0.12f; }, {"ctau"});

    auto df_f = df_dimu
        .Filter("passTrigger", "Trigger")
        .Filter("cut_kin", "Kinematics")
        .Filter("cut_charge", "Charge")
        .Filter("cut_mass", "MassWindow")
        .Filter("cut_soft", "SoftMuon")
        .Filter("cut_vtx", "VertexCL")
        .Filter("cut_dphi", "DPhi")
        .Filter("cut_match", "TrigMatch")
        .Filter("cut_ctau", "CTau");

    auto df_ana_vars = df_f
        .Define("anaPt", [](float pt) { return pt; }, {"dimuPt"})
        .Define("anaAbsY", [](float y) { return y; }, {"dimuAbsY"});

    auto hNum = df_ana_vars.Histo2D(
        {"hNum", "Analysis;|y|;p_{T} [GeV]",
         nYBins, yBins.data(), nPtBins, ptBins.data()},
        "anaAbsY", "anaPt"
    );

    // ---- 打印 cut flow report ----
    auto report = df_f.Report();
    report->Print();

    // ---- 效率计算 ----
    TH2D* hDen_ptr = (TH2D*)hDen.GetPtr();
    TH2D* hNum_ptr = (TH2D*)hNum.GetPtr();

    std::cout << "\nDenominator entries: " << hDen_ptr->GetEntries() << std::endl;
    std::cout << "Numerator entries:   " << hNum_ptr->GetEntries() << std::endl;

    if (hDen_ptr->GetEntries() == 0) {
        std::cerr << "ERROR: Denominator histogram is empty!" << std::endl;
        return;
    }

    TH2D* hEff = (TH2D*)hNum_ptr->Clone("hEff");
    hEff->SetTitle("J/#psi efficiency;|y|;p_{T} [GeV]");

    for (int ix = 1; ix <= hEff->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= hEff->GetNbinsY(); ++iy) {
            double num = hNum_ptr->GetBinContent(ix, iy);
            double den = hDen_ptr->GetBinContent(ix, iy);
            if (den > 0)
                hEff->SetBinContent(ix, iy, num / den);
            else
                hEff->SetBinContent(ix, iy, 0.0);
        }
    }

    double maxEff = hEff->GetMaximum();
    std::cout << "Maximum efficiency: " << maxEff << std::endl;

    std::ofstream fout("efficiency_table.txt");
    fout << "# y_low y_high pt_low pt_high efficiency error\n";
    fout << std::fixed << std::setprecision(6);

    for (int ix = 1; ix <= nYBins; ++ix) {
        double y_low  = yBins[ix-1];
        double y_high = yBins[ix];
        for (int iy = 1; iy <= nPtBins; ++iy) {
            double pt_low  = ptBins[iy-1];
            double pt_high = ptBins[iy];
            double num = hNum_ptr->GetBinContent(ix, iy);
            double den = hDen_ptr->GetBinContent(ix, iy);
            double eff = hEff->GetBinContent(ix, iy);
            double err = (den > 0) ? sqrt(eff * (1 - eff) / den) : 0.0;
            fout << y_low << " " << y_high << " "
                 << pt_low << " " << pt_high << " "
                 << eff << " " << err << "\n";
        }
    }
    fout.close();
    std::cout << "效率表格已写入 efficiency_table.txt" << std::endl;

    auto draw = [](TH2D* h, TString name) {
        TCanvas* c = new TCanvas("c", name, 900, 700);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.18);
        gPad->SetLogy();
        h->SetStats(0);
        h->SetMinimum(0);
        h->SetMaximum(0.5);
        h->GetXaxis()->SetTitle("|y|");
        h->GetYaxis()->SetTitle("p_{T} [GeV/c]");
        h->GetYaxis()->SetMoreLogLabels(false);
        h->GetYaxis()->SetNoExponent(false);
        h->Draw("COLZ TEXT");
        c->SaveAs(name + ".png");
        delete c;
    };

    draw(hEff, "H2D_efficiency");
}
