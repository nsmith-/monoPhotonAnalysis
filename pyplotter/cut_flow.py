#!/usr/bin/env python
import ROOT
import CMS_lumi, tdrstyle
import argparse
import plot_functions as plotter
import operator

def main():
    plot_info = getPlotArgs()
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
    hist_opts = "hist"
    line_color = ROOT.kBlack
    fill_color = ROOT.kWhite
    cut_flow = { "Medium Photon" : 500, "pfMET > 90 GeV" : 281,
                 "#Delta#phi > 2 " : 273, "electron veto" : 273,
                 "muon veto" : 271 }
    list_cut_flow = sorted(cut_flow.items(), key=operator.itemgetter(1), reverse=True)
    hist = getCutFlowHist(list_cut_flow)
    plotter.setHistAttributes(hist, plot_info, line_color, fill_color)
    plotter.makePlot(hist, hist_opts, plot_info)  
def getCutFlowHist(cut_flow):
    hist = ROOT.TH1F("hist", "Cut Flow", 1, 0, 1)
    for cut in cut_flow:
        hist.Fill(cut[0], cut[1])
    return hist
def getPlotArgs():
    parser = plotter.getBasicParser()
    return vars(parser.parse_args())
    
if __name__ == "__main__":
    main()
