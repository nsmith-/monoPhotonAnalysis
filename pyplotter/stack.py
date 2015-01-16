import ROOT
import numpy
import CMS_lumi, tdrstyle
import argparse

# what to draw
#label = "#slash{E}_{T} [GeV]"
#drawString = "pfMET"
label = "E_{T}^{#gamma} [GeV]"
drawString = "phoEt[selectedPhoton]"
cutString = ""
xrange = [100, 1100]
yrange = [0.1, 100]
fill = True
normalize = False
# qcd znunu wJet sig
colors = [ROOT.kGreen-5, ROOT.kBlue-6, ROOT.kGray, ROOT.kRed]

tdrstyle.setTDRStyle()

fdata = ROOT.TFile("../data_selection.root")
fqcd = ROOT.TFile("../qcd_selection.root")
fznunu = ROOT.TFile("../znunu_selection.root")
fwJets = ROOT.TFile("../wJets_selection.root")
fsig = ROOT.TFile("../signal_selection.root")

data = fdata.Get("EventTree")
qcd = fqcd.Get("EventTree")
znunu = fznunu.Get("EventTree")
wJets = fwJets.Get("EventTree")
signal = fsig.Get("EventTree")

canvas = ROOT.TCanvas("canvas", "canvas")
stack = ROOT.THStack("stack", "")

hdata = ROOT.TH1F("data", "Data", 20, xrange[0], xrange[1])
hqcd = ROOT.TH1F("qcd", "QCD", 20, xrange[0], xrange[1])
hznunu = ROOT.TH1F("znunu", "Z#nu#nu", 20, xrange[0], xrange[1])
hwJets = ROOT.TH1F("wJets", "W+Jets", 20, xrange[0], xrange[1])
hsignal = ROOT.TH1F("signal", "Signal", 20, xrange[0], xrange[1])

data.Draw(drawString+">>data", cutString, "goff")
hdata.SetMarkerStyle(ROOT.kFullCircle)

qcd.Draw(drawString+">>qcd", cutString, "goff")
if ( fill ) : hqcd.SetFillColor(colors[0])
hqcd.SetLineColor(colors[0])

znunu.Draw(drawString+">>znunu", cutString, "goff")
if ( fill ) : hznunu.SetFillColor(colors[1])
hznunu.SetLineColor(colors[1])

wJets.Draw(drawString+">>wJets", cutString, "goff")
if ( fill ) : hwJets.SetFillColor(colors[2])
hwJets.SetLineColor(colors[2])

signal.Draw(drawString+">>signal", cutString, "goff")
hsignal.SetLineColor(colors[3])
hsignal.SetLineStyle(ROOT.kDashed)

# rescale MC
hqcd.Sumw2()
hqcd.Scale(13/hqcd.Integral())
hznunu.Sumw2()
hznunu.Scale(193/hznunu.Integral())
hwJets.Sumw2()
hwJets.Scale(5/hwJets.Integral())
hsignal.Scale(70/hsignal.Integral())
bgTotal = [-1]
for i in range(1,21) :
  # first calculate uncertainty
  bgTotal.append(hqcd.GetBinError(i)+hznunu.GetBinError(i)+hwJets.GetBinError(i))

# background uncertainty without affect stack total
hbgUncertain = ROOT.TH1F("bgUncertain", "Background Uncertainty", 20, xrange[0], xrange[1])
hbgUncertain.SetFillStyle(3013)
hbgUncertain.SetFillColor(ROOT.kBlack)
hbgUncertain.SetLineColor(ROOT.kWhite)
for i in range(1,21) :
  print i
  hbgUncertain.SetBinContent(i, numpy.sqrt(bgTotal[i]))
  hsignal.SetBinContent(i, hsignal.GetBinContent(i)-hbgUncertain.GetBinContent(i))

bgDraw = "hist"
stack.Add(hqcd, bgDraw)
stack.Add(hwJets, bgDraw)
stack.Add(hznunu, bgDraw)
stack.Add(hbgUncertain, bgDraw)
stack.Add(hsignal, bgDraw)

stack.SetMinimum(yrange[0])
stack.SetMaximum(yrange[1])
#stack.Draw("nostack")
stack.Draw()
hdata.Draw("E same")

# Now, make pretty
stack.GetXaxis().SetTitle(label)
stack.GetYaxis().SetTitle("Events / 50 GeV")
canvas.SetLogy(True)
legend = ROOT.TLegend(0.65, 0.7, 0.95, 0.85)
legend.AddEntry(hdata, "Data", "ep")
legend.AddEntry(hqcd, "QCD", "f")
legend.AddEntry(hznunu, "Z#nu#nu", "f")
legend.AddEntry(hwJets, "W+Jets", "f")
legend.AddEntry(hbgUncertain, "Bkg. Stat.", "f")
legend.AddEntry(hsignal, "SM+ADD(n_{ED}=2, M_{D}=3)", "l")
legend.SetBorderSize(0)
legend.SetFillColor(0)

plot_info = vars(getBasicParser()) 
#legend.Draw()
makePlot(stack, hdata, "hist", legend, plot_info)
def makePlot (hist_stack, data_hist, hist_opts, legend, plot_info):
    #legend = ROOT.TLegend(.5 ,.65 ,.885 ,.875)
    canvas = getCanvas()
    if plot_info["logy"]:
        canvas.SetLogy()
    if plot_info["logx"]:
        canvas.SetLogx()
    #draw the lumi text on the canvas
    hist_stack.Draw(hist_opts)
    hist_stack.GetXaxis().SetTitle(plot_info["xlabel"])
    if plot_info["ylabel"] == "":
        plot_info["ylabel"] = "Events / %s GeV" % int(hist.GetBinWidth(1))
    hist_stack.GetYaxis().SetTitle(plot_info["ylabel"])
    hist_stack.SetTitleOffset(1.3, "y")
    hist_stack.SetTitleOffset(1.1, "x")
    data_hist.Draw("same e1")
    setTDRStyle(canvas, 1, 13, plot_info["printCMS"]) 
    canvas.cd()
    canvas.Update()
    canvas.RedrawAxis()
    #frame = canvas.GetFrame()
    #frame.Draw()
    #legend.SetFillColor(ROOT.kWhite)
    #legend.AddEntry(hist, legendName)

    #legend.Draw("same")
    canvas.Print(plot_info["output_file"]) 
def setTDRStyle(canvas, luminosity, energy, printCMS):
    tdrstyle.setTDRStyle() 
    if printCMS == "right" or printCMS == "left":
        if energy == 13:
            CMS_lumi.lumi_13TeV = "%s fb^{-1}" % luminosity
            if printCMS == "left":
                iPos = 11
            else:
                iPos = 13
            CMS_lumi.writeExtraText = 1
            CMS_lumi.extraText = "Preliminary"
            CMS_lumi.CMS_lumi(canvas, 4, iPos)
def getCanvas():
    H_ref = 600; 
    W_ref = 800; 
    W = W_ref
    H  = H_ref

    T = 0.08*H_ref
    B = 0.12*H_ref 
    L = 0.12*W_ref
    R = 0.04*W_ref

    canvas = ROOT.TCanvas("c2","c2",50,50,W,H)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W )
    canvas.SetRightMargin( R/W )
    canvas.SetTopMargin( T/H )
    canvas.SetBottomMargin( B/H )
    canvas.SetTickx(0)
    canvas.SetTicky(0) 
    return canvas

def getBasicParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_file', type=str, required=True,
                        help="Name produced plot file (type pdf/png/jpg etc.).")
    parser.add_argument('-t', '--tree_name', type=str, required=False, default="EventTree",
                        help="Plot group (folder in root file)")  
    parser.add_argument('-v', '--tree_var', type=str, required=False,
                        help="Variable name in root tree")  
    parser.add_argument('--xlabel', type=str, required=False, default="", 
                        help="x axis label")
    parser.add_argument('--ylabel', type=str, required=False, default="", 
                        help="y axis label")
    parser.add_argument('--xmin', type=float, required=False, default=0, 
                        help="minimum x value")
    parser.add_argument('--xmax', type=float, required=False, default=0, 
                        help="maximum x value")   
    parser.add_argument('--ymin', type=float, required=False, default=0, 
                        help="minimum y value")
    parser.add_argument('--ymax', type=float, required=False, default=0, 
                        help="maximum y value")
    parser.add_argument('--rebin', type=int, required=False, default=0, 
                        help="Number of bins to group together (1D only)")
    parser.add_argument('--logy', action='store_true',
                        help="Set y axis to logarithmic scale")
    parser.add_argument('--logx', action='store_true', 
                        help="Set x axis to logarithmic scale")
    parser.add_argument('--printCMS', type=str, default="left",required=False,
                        choices=["left","right"], help="""print 'CMS preliminary' 
                        in left (or right) upper corner""")
    return parser
