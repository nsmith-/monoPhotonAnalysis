#!/bin/env python
import ROOT
import CMS_lumi, tdrstyle
import argparse
import plot_functions as plotter
import operator

def main():
    plot_info = getPlotArgs()

    cut_flow_vars = ["No cuts", "Medium Photon ID", "pfMET > 90 GeV",  "#Delta#phi > 2",
                    "e veto", "#mu veto"]
    data = [4519553, 500, 281, 273, 273, 271]
    wjets = [5941270, 74, 18, 8, 8, 7]
    qcd = [487072, 211, 13, 7, 7, 7]
    znunu = [50000, 19269, 18821, 18643, 18643, 18635]

    samples = [data, wjets, qcd, znunu]
    normalizations = [271, 5., 13., 193.]
    
    data_hist = getCutFlowHist(cut_flow_vars, data, 1, plot_info, ROOT.kBlack)
    stacked_hist = ROOT.THStack()
    stacked_hist.Add(data_hist, "hist")
    plot_info["ylabel"] = "Events"
    i = 0
    colors = [ROOT.kBlack, ROOT.kRed-4, ROOT.kBlue-4, ROOT.kGreen-5, ROOT.kYellow+8]
    canvas = plotter.getCanvas()
    plotter.setTDRStyle(canvas, 1, 13, plot_info["printCMS"]) 
    
    names = ["Data", "W+jets", "QCD", "Z(#nu#nu)#gamma"] 
    
    legend = ROOT.TLegend(0.50, 0.6, 0.90, 0.80)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    hist = {}
    for i in range(0, len(samples)):
        wgt = normalizations[i]/samples[i][-1]
        print "weight is %s" % wgt
        print "normalized %s to %s" % (samples[i][-1], normalizations[i])
        hist.update({names[i] 
            : getCutFlowHist(cut_flow_vars, samples[i], wgt, plot_info, colors[i])})
        stacked_hist.Add(hist[names[i]])
        legend.AddEntry(hist[names[i]], names[i], "l")
        i += 1
    stacked_hist.Draw()
    for i in range(0, len(cut_flow_vars)):
        stacked_hist.GetXaxis().SetBinLabel(i+1, cut_flow_vars[i])
       #draw the lumi text on the canvas
    if plot_info["logy"]:
        canvas.SetLogy()
    if plot_info["logx"]:
        canvas.SetLogx()
    stacked_hist.Draw("nostack")
    legend.Draw()
    #draw the lumi text on the canvas
    plotter.setTDRStyle(canvas, 1, 13, plot_info["printCMS"]) 
    stacked_hist.GetXaxis().SetTitle(plot_info["xlabel"])
    if plot_info["ylabel"] == "":
        plot_info["ylabel"] = "Events / %s GeV" % int(hist.GetBinWidth(1))
    stacked_hist.GetYaxis().SetTitle(plot_info["ylabel"])
    #hist.SetTitleOffset(1.3, "y")
    #hist.SetTitleOffset(1.1, "x")
    canvas.cd()
    canvas.Update()
    canvas.RedrawAxis()
    #frame = canvas.GetFrame()
    #frame.Draw()
    #legend.SetFillColor(ROOT.kWhite)
    #legend.AddEntry(hist, legendName)

    #legend.Draw("same")
    canvas.Print(plot_info["output_file"]) 

def getCutFlowHist(names, values, wgt, plot_info, color):
    hist = ROOT.TH1F("hist", "Cut Flow", len(names), 0, len(names))
    i = 1
    for value in values:
        hist.SetBinContent(i, value*wgt)
        i += 1
    fill_color = ROOT.kWhite    
    plotter.setHistAttributes(hist, plot_info, color, fill_color)
    return hist
def getPlotArgs():
    parser = plotter.getBasicParser()
    return vars(parser.parse_args())
    
if __name__ == "__main__":
    main()
