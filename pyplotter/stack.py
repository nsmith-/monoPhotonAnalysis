import ROOT
import numpy
import CMS_lumi, tdrstyle

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
legend.Draw()

CMS_lumi.lumi_13TeV = "%s fb^{-1}" % 1
iPos = 13
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.CMS_lumi(canvas, 4, iPos)
