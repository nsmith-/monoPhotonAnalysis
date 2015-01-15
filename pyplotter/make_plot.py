#!/usr/bin/env python
import ROOT
import CMS_lumi, tdrstyle
import argparse
import plot_functions as plotter

def main():
    plot_info = getPlotArgs()
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
    hist = plotter.getHistFromFile(plot_info)
    plotter.setHistAttributes(hist, plot_info, ROOT.kRed+4, ROOT.kOrange-8)
    if type(hist) == "<class '__main__.TH2F'>":
        hist_opts = "colz"
    else:
        hist_opts = "hist"
    plotter.makePlot(hist, hist_opts, plot_info)
def getPlotArgs():
    parser = plotter.getBasicParser()
    parser.add_argument("-n", "--file_name", type=str, required=True,
                        help="Name of root file where plots are stored")
    return vars(parser.parse_args())
    
if __name__ == "__main__":
    main()
