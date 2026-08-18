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

// ============================================================
// 原始筛选函数
// ============================================================
RVec<ROOT::Math::PxPyPzMVector> get_mumu_candidates(
    const RVec<float>& muPx,
    const RVec<float>& muPy,
    const RVec<float>& muPz,
    const RVec<float>& muCharge,
    const RVec<float>& mumuonlymu1Idx,
    const RVec<float>& mumuonlymu2Idx,
    unsigned int nMu,
    const RVec<float>& mumuonlyVtxCL,
    const RVec<int>& muIsPatSoftMuon,
    const RVec<std::string>& TrigNames,
    const RVec<unsigned int>& TrigRes,
    const RVec<int>& muIsJpsiTrigMatch,
    const RVec<float>& mumuonlyctau
) {
    RVec<ROOT::Math::PxPyPzMVector> candidates;

    size_t nPairs = mumuonlymu1Idx.size();
    if (nPairs != mumuonlymu2Idx.size()) {
        std::cerr << "Warning: mumuonlymu1Idx.size() != mumuonlymu2Idx.size()" << std::endl;
        nPairs = std::min(nPairs, mumuonlymu2Idx.size());
    }

    for (size_t i = 0; i < nPairs; ++i) {
        bool TrigDiMuon_4_3 = false;
        for (size_t t = 0; t < TrigNames.size(); ++t) {
            if (TrigNames[t].find("HLT_DoubleMu4_3_LowMass_v") != std::string::npos && TrigRes[t] == 1) {
                TrigDiMuon_4_3 = true;
                break;
            }
        }
        if (!TrigDiMuon_4_3) continue;

        int idx1 = static_cast<int>(mumuonlymu1Idx[i]);
        int idx2 = static_cast<int>(mumuonlymu2Idx[i]);
        if (idx1 < 0 || idx1 >= (int)muPx.size() || idx2 < 0 || idx2 >= (int)muPx.size()) continue;

        ROOT::Math::PxPyPzMVector mu1(muPx[idx1], muPy[idx1], muPz[idx1], MUON_MASS);
        ROOT::Math::PxPyPzMVector mu2(muPx[idx2], muPy[idx2], muPz[idx2], MUON_MASS);

        float pt1 = mu1.pt(), eta1 = mu1.eta();
        float pt2 = mu2.pt(), eta2 = mu2.eta();

        bool pass1 = (pt1 > 5.0 && fabs(eta1) < 2.4);
        bool pass2 = (pt2 > 5.0 && fabs(eta2) < 2.4);

        bool oppositeCharge = (muCharge[idx1] * muCharge[idx2] < 0);

        ROOT::Math::PxPyPzMVector dimu = mu1 + mu2;
        float mumu_mass = dimu.M();

        bool passMass = (mumu_mass >= 2.5 && mumu_mass <= 3.5);

        bool passSoft = (muIsPatSoftMuon[idx1] == 1 && muIsPatSoftMuon[idx2] == 1);
        bool passVtx = (mumuonlyVtxCL[i] > 0.01);

        float dphi = mu1.phi() - mu2.phi();
        if (dphi > M_PI) dphi -= 2 * M_PI;
        if (dphi < -M_PI) dphi += 2 * M_PI;
        bool passDphi = (dphi > 0);

        bool passMatch = (muIsJpsiTrigMatch[idx1] == 1 && muIsJpsiTrigMatch[idx2] == 1);

bool passCtau = (mumuonlyctau[i] >= -0.03 && mumuonlyctau[i] <= 0.12);

        bool pass_all = TrigDiMuon_4_3 && pass1 && pass2 && oppositeCharge && passMass &&
                         passSoft && passVtx && passDphi && passMatch &&passCtau;

        if (pass_all) {
            candidates.push_back(dimu);
        }
    }
    return candidates;
}

// ============================================================
// 主函数
// ============================================================
void make_mc_efficiency_2D() {
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("0.6f");
    gStyle->SetPalette(kRainBow);

    ROOT::EnableImplicitMT(8);

    TString fileName = "/eos/user/y/yingyinw/ntuple_practice/CMSSW_13_0_13/src/UserCode/MultiLepPAT/test/Run_Crab_Script/crab_Official_MC_JpsiPT5_MultiLepPAT_v4/results/mymultilep_*.root";

    ROOT::RDataFrame df("mkcands/X_data", fileName);

    // ============================================================
    // 分母：MC truth J/psi
    // ============================================================
    auto df_gen =
        df.Filter("Offical_Pythia_J_MC_xPx.size()>0")
        .Define(
            "genP4",
            [](const RVec<float>& px, const RVec<float>& py,
               const RVec<float>& pz, const RVec<float>& m) {
                RVec<ROOT::Math::PxPyPzMVector> vec;
                for (size_t i = 0; i < px.size(); ++i) {
                    vec.push_back(ROOT::Math::PxPyPzMVector(px[i], py[i], pz[i], m[i]));
                }
                return vec;
            },
            {
                "Offical_Pythia_J_MC_xPx",
                "Offical_Pythia_J_MC_xPy",
                "Offical_Pythia_J_MC_xPz",
                "Offical_Pythia_J_MC_xM"
            }
        )
        .Define(
            "genJpsiPt",
            [](const RVec<ROOT::Math::PxPyPzMVector>& v) {
                RVec<float> pts;
                for (auto& p : v) pts.push_back(static_cast<float>(p.pt()));
                return pts;
            },
            {"genP4"}
        )
        .Define(
            "genJpsiAbsY",
            [](const RVec<ROOT::Math::PxPyPzMVector>& v) {
                RVec<float> ys;
                for (auto& p : v) ys.push_back(static_cast<float>(fabs(p.Rapidity())));
                return ys;
            },
            {"genP4"}
        );

    // ============================================================
    // 分子：在分母基础上加 analysis
    auto df_ana = df_gen.Define(
        "mumuCandidates",
        get_mumu_candidates,
        {
            "muPx","muPy","muPz","muCharge",
            "mumuonlymu1Idx","mumuonlymu2Idx","nMu",
            "mumuonlyVtxCL","muIsPatSoftMuon","TrigNames","TrigRes","muIsJpsiTrigMatch","mumuonlyctau"
        }
    );

    auto df_ana_vars = df_ana
        .Define("anaPt", [](const RVec<ROOT::Math::PxPyPzMVector>& v) {
            RVec<float> pts;
            for (auto& p : v) pts.push_back(static_cast<float>(p.pt()));
            return pts;
        }, {"mumuCandidates"})
        .Define("anaAbsY", [](const RVec<ROOT::Math::PxPyPzMVector>& v) {
            RVec<float> ys;
            for (auto& p : v) ys.push_back(static_cast<float>(fabs(p.Rapidity())));
            return ys;
        }, {"mumuCandidates"});

    // ============================================================
    // 分桶
    // ============================================================
    std::vector<double> yBins = {0.0, 0.3, 0.6, 0.9, 1.2, 1.8, 2.4};
    std::vector<double> ptBins = {
        10.0,10.5,  11.0,11.5, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0,
        18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0,
        28.0, 29.0, 30.0, 32.0, 34.0, 36.0, 38.0, 42.0, 46.0, 50.0,
        60.0, 75.0, 95.0, 120.0
    };
    int nYBins = yBins.size() - 1;
    int nPtBins = ptBins.size() - 1;

    // ============================================================
    // 填充直方图
    // ============================================================
    auto hDen = df_gen.Histo2D(
        {"hDen", "MC Truth;|y|;p_{T} [GeV]",
         nYBins, yBins.data(), nPtBins, ptBins.data()},
        "genJpsiAbsY", "genJpsiPt"
    );

    auto hNum = df_ana_vars.Histo2D(
        {"hNum", "Analysis;|y|;p_{T} [GeV]",
         nYBins, yBins.data(), nPtBins, ptBins.data()},
        "anaAbsY", "anaPt"
    );

    // ============================================================
    // 效率：显式逐 bin 计算
    // ============================================================
    TH2D* hDen_ptr = (TH2D*)hDen.GetPtr();
    TH2D* hNum_ptr = (TH2D*)hNum.GetPtr();

    std::cout << "\nMC Truth entries: " << hDen_ptr->GetEntries() << std::endl;
    std::cout << "Analysis entries:  " << hNum_ptr->GetEntries() << std::endl;

//    if (hDen_ptr->GetEntries() == 0) {
//        std::cerr << "ERROR: MC Truth histogram is empty!" << std::endl;
//        return;
//    }

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

    // ============================================================
    // 输出效率表格到 txt（修正循环变量对应关系）
    // ============================================================
    std::ofstream fout("efficiency_table.txt");
    fout << "# y_low y_high pt_low pt_high efficiency error\n";
    fout << std::fixed << std::setprecision(6);

    // ix 对应 |y| 轴（X轴），iy 对应 pT 轴（Y轴）
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

    // ============================================================
    // 绘图
    // ============================================================
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

  //  draw(hDen_ptr, "H2D_mctruth");
  //  draw(hNum_ptr, "H2D_analysis");
    draw(hEff,     "H2D_efficiency");
}
