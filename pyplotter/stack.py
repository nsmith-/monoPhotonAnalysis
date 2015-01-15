import ROOT
import CMS_lumi, tdrstyle

tdrstyle.setTDRStyle()

fdata = ROOT.TFile("../data_selection.root")
fqcd = ROOT.TFile("../qcd_selection.root")
fznunu = ROOT.TFile("../znunu_selection.root")
fwJets = ROOT.TFile("../wJets_selection.root")

data = fdata.Get("EventTree")
qcd = fqcd.Get("EventTree")
znunu = fznunu.Get("EventTree")
wJets = fwJets.Get("EventTree")

canvas = ROOT.TCanvas("canvas", "canvas")
stack = ROOT.THStack("stack", "")

hdata = ROOT.TH1F("data", "Data", 20, 100, 1200)
hqcd = ROOT.TH1F("qcd", "QCD", 20, 100, 1200)
hznunu = ROOT.TH1F("znunu", "Z#nu#nu", 20, 100, 1200)
hwJets = ROOT.TH1F("wJets", "W+Jets", 20, 100, 1200)

# what to draw
drawString = "phoEt[selectedPhoton]"
cutString = ""

data.Draw(drawString+">>data", cutString, "goff")
hdata.SetMarkerStyle(ROOT.kFullCircle)

qcd.Draw(drawString+">>qcd", cutString, "goff")
hqcd.SetFillColor(ROOT.kBlue)

znunu.Draw(drawString+">>znunu", cutString, "goff")
hznunu.SetFillColor(ROOT.kRed)

wJets.Draw(drawString+">>wJets", cutString, "goff")
hwJets.SetFillColor(ROOT.kGreen)

# rescale MC
hqcd.Scale(13/hqcd.Integral())
hznunu.Scale(193/hznunu.Integral())
hwJets.Scale(5/hwJets.Integral())

bgDraw = ""
stack.Add(hqcd, bgDraw)
stack.Add(hwJets, bgDraw)
stack.Add(hznunu, bgDraw)

stack.SetMaximum(120)
stack.Draw()
hdata.Draw("E0 same")
CMS_lumi.lumi_13TeV = "%s fb^{-1}" % 1
iPos = 13
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.CMS_lumi(canvas, 4, iPos)
