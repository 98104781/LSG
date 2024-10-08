import Lipids.GenerateLipids as GL
import copy

''' EXAMPLE EXAMPLE EXAMPLE

class LipidClass(GL.Glycerolipid):

  List of adducts to generate spectra for = {

  "Adduct":{
    Fragment:Intensity,
    Fragment:Intensity,
    Fragment:Intensity},

  def __init__(self, sn1, sn2, sn3):
    sn3 = ['Head', Headgroup mass]
    super().__init__(LipidClass.adducts, sn1, sn2, sn3)

EXAMPLE EXAMPLE EXAMPLE '''

# https://lipidmaps.org/resources/standard?lipid_category=GL Lots of good example spectra here

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class MG(GL.Glycerolipid):

  ##### Metabolites 2016, 6(3), 25; https://doi.org/10.3390/metabo6030025
  ##### LipidBlast DOI:10.1038/nmeth.2551
  ##### http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
  tooltip = 'Monoacylglycerol'
  tailOrganisation = ['A']

  adducts = {  # adduct:{spectra}
  "[M+H]+":{
    GL.MA      :100,
    GL.MA_s_H2O :50,
    GL.FAkH     :50},

  "[M+H-H2O]+":{
    GL.MA      :100,
    GL.FAkH     :50},

  "[M+Na]+":{    
    GL.MA      :100,
    GL.MA_s_H2O :50,
    GL.FAkA     :10,
    GL.FAkH     :50},

  "[M+NH4]+":{ # Matches LipidBlast    
    GL.MA        :1,
    GL.MH       :10,
    GL.MH_s_H2O:100,
    GL.FAH       :1,
    GL.FAkH     :50},
  }

  def __init__(self, sn1):
    super().__init__(MG.adducts, sn1=sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class DG(GL.Glycerolipid):

  ##### Metabolites 2016, 6(3), 25; https://doi.org/10.3390/metabo6030025
  ##### LipidBlast DOI:10.1038/nmeth.2551
  ##### https://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=586&LM_ID=LMGL02010009&TRACK_ID=76 
  ##### http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
  tooltip = 'Diacylglycerol'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M+H]+":{
    GL.MA       :25,
    GL.MA_s_H2O:100,
    GL.MH_s_FA  :25,
    GL.FAkH      :1},

  "[M+H-H2O]+":{
    GL.MA      :100,
    GL.MH_s_FA :100,
    GL.FAkH      :1},

  "[M+Na]+":{    
    GL.MA       :50,
    GL.MA_s_H2O:100,
    GL.MA_s_FA  :20,
    GL.MH_s_FA  :50,
    GL.FAkA     :10,
    GL.FAkH      :1},

  "[M+NH4]+":{ # Matches LipidBlast    
    GL.MA       :10, # Molecular ion
    GL.MH        :5, # Present in lipidBlast
    GL.MH_s_H2O:100, # Present in lipidBlast
    GL.MH_s_FA :100, # Present in lipidBlast
    GL.FAkH      :1}
  }

  def __init__(self, sn1, sn2):
    super().__init__(DG.adducts, sn1=sn1, sn2=sn2)

# ~ # ~ # ~ # ~ # ~ # ~ #

class TG(GL.Glycerolipid):

  ##### Metabolites 2016, 6(3), 25; https://doi.org/10.3390/metabo6030025
  ##### LipidBlast DOI:10.1038/nmeth.2551
  tooltip = 'Triacylglycerol'
  tailOrganisation = ['AAA']

  adducts = {  # adduct:{spectra}
  "[M+Na]+":{ # Approximate to LipidBlast   
    GL.MA       :25,# Molecular ion
    GL.MA_s_FA  :10,# Present in lipidBlast
    GL.MH_s_FA :100,
    GL.FAkA      :5,
    GL.FAkH     :10},

  "[M+NH4]+":{ # Approximate to LipidBlast
    GL.MA       :25, # Molecular ion
    GL.MH        :5, # Present in lipidBlast
    GL.MH_s_FA :100, # Present in lipidBlast
    GL.FAkH     :10}
  }

  def __init__(self, sn1, sn2, sn3):
    super().__init__(TG.adducts, sn1, sn2, sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

class oTG(GL.Glycerolipid):

  tooltip = 'Ether/Plasmanyl Triacylglycerol'
  tailOrganisation = ['OAA']
  givenName = 'Ether-TG'

  adducts = {
  #"[M+Na]+":{
  #  GL.MA       :25,
  #  GL.MA_s_FA  :10,
  #  GL.MH_s_FA :100,
  #  GL.FAkA      :5,
  #  GL.FAkH     :10},

  "[M+NH4]+":{# http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
    GL.MA       :25,
    GL.MH        :5,
    GL.MH_s_FA :100,
    GL.FAkH     :10}
  }

  def __init__(self, sn1, sn2, sn3):
    super().__init__(oTG.adducts, sn1, sn2, sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

class DGGA(GL.Glycerolipid):

  tooltip = 'Diacylglyceryl glucuronide'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    # Taken from http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html 
    # 10.1021/ac502511a
    GL.MA        :100,
    GL.MA_s_FAk    :5,
    GL.MA_s_FA     :1,
    GL.FAH        :50,
    GL.C9H13O8     :5},

  "[M+NH4]+":{ 
    # Taken from http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html  
    GL.MA          :5,
    GL.MA_s_H2O    :0,
    GL.MA_s_2H2O   :0,
    GL.MH_s_HG     :5,
    GL.MH_s_HG_H2O:60,
    GL.MA_s_FA     :0,
    GL.MH_s_HG_FA:100,
    GL.FAkH        :0}
    }

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=194.042653, type='Headgroup', chnops={'C':6, 'H':10, 'O':7},
    smiles='OC1C(O)C(O)C(C(=O)O)OC1')
    super().__init__(DGGA.adducts, sn1=sn1, sn2=sn2, sn3=headgroup)

class ADGGA(GL.Glycerolipid):

  tooltip = 'Acyl-diacylglyceryl glucuronide'
  tailOrganisation = ['AAA']
  specificTailOrganisation = ['A','AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Taken from http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
    GL.MA       :100,
    GL.FAH       :50},

  "[M+NH4]+":{ # Taken from http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html 
    GL.MA         :5,
    GL.MH_s_HG    :0,
    GL.MH_s_HG_H2O:5,
    GL.HGH_s_H2O:100,
    GL.MA_s_FA    :0,
    GL.MH_s_HG_FA:40,
    GL.FAkH      :50}
    }

  def __init__(self, sn1, sn2, sn3):
    headgroup = GL.sn(mass=194.042653, type='Headgroup', chnops={'C':6, 'H':10, 'O':7},
    smiles='O('+sn1.smiles+')C1C(O)C(O)C(C(=O)O)OC1', hgtails = [sn1])
    super().__init__(ADGGA.adducts, sn1=sn2, sn2=sn3, sn3=headgroup)
    self.specificname=f"{self.lipid_class} {sn1.name}/{sn2.name}_{sn3.name}"

# ~ # ~ # ~ # ~ # ~ # ~ #

class MGDG(GL.Glycerolipid): # Monogalactosyl diacylglycerol

  tooltip = 'Monoglycosyl diacylglycerol'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M+Na]+":{ # https://doi.org/10.1111/j.1440-1835.2010.00582.x   
    GL.MA          :5, # This spectra agrees with lipidblast and the attached
    GL.MA_s_FA   :100},# reference. Ion ratios are off, but depend on isomerism.

  "[M+NH4]+":{ # https://doi.org/10.1002/pld3.183 Figure S4     
    GL.MA          :5, # LipidBlast includes fragment 'MA_s_FA' but
    GL.MH_s_HG     :5, # the fragment is not present in the reference.
    GL.MH_s_HG_H2O:60, # Included here at 0 intensity so it can be modified
    GL.MA_s_FA     :0, # in the program if needed, otherwise it will be
    GL.MH_s_HG_FA:100} # ignored on spectra generation.
    }

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=180.063388, type='Headgroup', chnops={'C':6, 'H':12, 'O':6},
    smiles='OC1C(O)C(O)C(CO)OC1')
    super().__init__(MGDG.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class oMGDG(GL.Glycerolipid): # Monogalactosyl diacylglycerol

  tooltip = 'Ether/Plasmanyl Monoglycosyl diacylglycerol'
  tailOrganisation = ['OA']
  givenName = 'Ether-MGDG'

  adducts = {  # adduct:{spectra}

  "[M+NH4]+":{# http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
    GL.MA          :5,
    GL.MH_s_HG     :5,
    GL.MH_s_HG_H2O:60,
    GL.MA_s_FA     :0,
    GL.MH_s_HG_OA:100,
    GL.FAkH        :2}
    }

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=180.063388, type='Headgroup', chnops={'C':6, 'H':12, 'O':6},
    smiles='OC1C(O)C(O)C(CO)OC1')
    super().__init__(oMGDG.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class SQDG(GL.Glycerolipid): # Sulphoquinovosyl diacylglycerol

  tooltip = 'Sulfoquinovosyl diacylglycerol'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Closely matches LipidBlast  
    GL.MA         :40, # Molecular ion
    GL.MA_s_FA    :10, # Present in lipidBlast
    GL.FAH        :10, # Present in lipidBlast
    GL.C6H9O7S   :100, # Present in lipidBlast
    GL.O3SH       :20},

  "[M+NH4]+":{ # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
    GL.MA              :5,
    GL.MH_s_HG         :5,
    GL.MH_s_HG_H2O    :15,
    GL.MH_s_FA        :10,
    GL.MH_s_HG_FA    :100,
    GL.FAkH            :5}
    }

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=244.025287, type='Headgroup', chnops={'C':6, 'H':12, 'O':8, 'S':1},
    smiles='OC1C(O)C(O)C(CS(=O)(=O)O)OC1')
    super().__init__(SQDG.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class SGDG(GL.Glycerolipid): # Sulphoglycerol diacylglycerol

  tooltip = 'Sulfoglycerol diacylglycerol'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    GL.MA    :100,
    GL.MA_s_HG :0,
    GL.MA_s_FAk:0,
    GL.MA_s_FA:10,
    GL.FAH     :3,
    GL.C6H9O8S :3,
    GL.O4SH   :20}
    }

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=260.020201964, type='Headgroup', chnops={'C':6, 'H':12, 'O':9, 'S':1},
    smiles='OC1C(OS(=O)(=O)O)C(O)C(CO)OC1')
    super().__init__(SGDG.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class oSGDG(GL.Glycerolipid): # Sulphoglycerol diacylglycerol

  tooltip = 'Ether/Plasmanyl Sulfoglycerol diacylglycerol'
  tailOrganisation = ['OA']
  givenName = 'Ether-SGDG'

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    GL.MH       :100,
    GL.MH_s_HG    :0,
    GL.MH_s_FAk   :0,
    GL.MH_s_FAcyl:10,
    GL.FAH        :4,
    GL.C6H9O8S    :4,
    GL.O4SH      :10}
    }

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=260.020201964, type='Headgroup', chnops={'C':6, 'H':12, 'O':9, 'S':1},
    smiles='OC1C(OS(=O)(=O)O)C(O)C(CO)OC1')
    super().__init__(oSGDG.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class DGDG(GL.Glycerolipid): # Digalactosyl diacylglycerol

  tooltip = 'Diglycosyl diacylglycerol'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{  
    GL.MA              :5,
    GL.C15H27O13      :15,
    GL.C15H25O12      :40,
    GL.C15H23O11      :15,
    GL.C9H15O7        :15,
    GL.C6H5O3         :15,
    GL.FAH           :100
    },

  "[M+Na]+":{ # https://doi.org/10.1111/j.1440-1835.2010.00582.x    
    GL.MA              :5, # Molecular ion
    GL.MA_H2O_s_Gal    :5,
    GL.MA_s_FA       :100, # Present in lipidBlast
    GL.MA_s_FA_Gal_H2O:50, # Present in lipidBlast
    },

  "[M+NH4]+":{ # https://doi.org/10.1002/pld3.183 Figure S4   
    GL.MA          :5,
    GL.MH_s_Gal_H2O:5,
    GL.MH_s_HG    :70,
    GL.MH_s_HG_H2O:70,
    GL.MA_s_FA     :0,
    GL.MH_s_HG_FA:100}
    }

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=342.116212, type='Headgroup', chnops={'C':12, 'H':22, 'O':11},
    smiles='OC1C(O)C(O)C(CO(C2OC(CO)C(O)C(O)C2(O)))OC1')
    super().__init__(DGDG.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class DGCC(GL.Glycerolipid): # Diacylglyceryl-3-O-carboxyhydroxymethylcholine

  tooltip = 'Diacylglyceryl-3-O-carboxyhydroxymethylcholine' 
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M+H]+":{ #  http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
    GL.MA        :100,
    GL.MA_s_FAk   :70,
    GL.MA_s_FA    :40,
    GL.C6H14NO2   :10,
    GL.C5H14NO    :10}}

  #"[M+NH4]+":{    
  #  GL.MA          :5}
  #  }

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=177.100108, type='Headgroup', chnops={'C':7, 'H':15, 'N':1, 'O':4},
    smiles='C[N+](C)(C)CCOC(C(=O)([O-]))')
    # headgroup mass has -H to maintain neutral charge
    super().__init__(DGCC.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class DGTS(GL.Glycerolipid): # N-trimethylhomoserine diacylglycerol

  tooltip = 'N-trimethylhomoserine diacylglycerol' 
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M+H]+":{ #  https://doi.org/10.1002/rcm.7847
    GL.MA        :100,
    GL.MA_s_FAk   :70,
    GL.MA_s_FA    :40,
    GL.MA_s_allFAk:20,
    GL.C7H14N1O2   :5}}

  #"[M+NH4]+":{    
  #  GL.MA          :5}
  #  }

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=161.105193, type='Headgroup', chnops={'C':7, 'H':15, 'N':1, 'O':3},
    smiles='[O-]C(=O)C([N+](C)(C)(C))CC')
    # headgroup mass has -H to maintain neutral charge
    super().__init__(DGTS.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class AC2PIM1(GL.Glycerolipid):

  tooltip = 'Phosphatidylinositol-mannoside diacylglycerol' 
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Matches LipidBlast
    GL.MA            :50, # Molecular ion
    GL.MA_H2O_s_Gal  :20, # Present in lipidBlast
    GL.MA_s_FA      :100, # Present in lipidBlast
    GL.MH_PO4_s_HG   :60, # Present in lipidBlast
    GL.MH_PO4_s_HG_FA:40, # Present in lipidBlast
    GL.FAH           :40}}# Present in lipidBlast

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=422.082541, type='Headgroup', chnops={'C':12, 'H':23, 'O':14, 'P':1},
    smiles='OC1C(O)C(O)C(O)C(OC2OC(CO)C(O)C(O)C2O)C1OP(O)(=O)')
    super().__init__(AC2PIM1.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class AC2PIM2(GL.Glycerolipid):

  tooltip = 'Phosphatidylinositol-dimannoside diacylglycerol'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Matches LipidBlast
    GL.MA            :10, # Molecular ion
    GL.MA_H2O_s_Gal   :5, # Present in lipidBlast
    GL.MA_s_FA      :100, # Present in lipidBlast
    GL.MA_s_allFA    :20, # Present in lipidBlast
    GL.FAH           :40}}# Present in lipidBlast

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=584.135365, type='Headgroup', chnops={'C':18, 'H':33, 'O':19, 'P':1},
    smiles='O(C3OC(CO)C(O)C(O)C3O)C1C(O)C(O)C(O)C(OC2OC(CO)C(O)C(O)C2O)C1OP(O)(=O)')
    super().__init__(AC2PIM2.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class AC3PIM2(GL.Glycerolipid):

  tooltip = 'Acyl-phosphatidylinositol-dimannoside diacylglycerol'
  tailOrganisation = ['AAA']
  specificTailOrganisation = ['A','AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Matches LipidBlast
    GL.MA            :20, # Molecular ion
    GL.MA_H2O_s_Gal   :5, # Present in lipidBlast
    GL.MA_s_FAk      :10, # Present in lipidBlast
    GL.MA_s_FA      :100, # Present in lipidBlast
    GL.MA_s_allFA    :20, # Present in lipidBlast
    GL.HGA_s_H2O     :40, # Present in lipidBlast
    GL.FAH           :40}}# Present in lipidBlast

  def __init__(self, sn1, sn2, sn3):
    headgroup = GL.sn(mass=584.135365, type='Headgroup', chnops={'C':18, 'H':33, 'O':19, 'P':1},
    smiles='O(C3OC(CO)C(O)C(O)C3O)C1C(O)C(O)C(O)C(OC2OC(CO('+sn1.smiles+'))C(O)C(O)C2O)C1OP(O)(=O)', hgtails = [sn1])
    super().__init__(AC3PIM2.adducts, sn1=sn2, sn2=sn3, sn3=headgroup)
    self.specificname=f"{self.lipid_class} {sn1.name}/{sn2.name}_{sn3.name}"

# ~ # ~ # ~ # ~ # ~ # ~ #

class AC4PIM2(GL.Glycerolipid):

  tooltip = 'Diacyl-phosphatidylinositol-dimannoside diacylglycerol'
  tailOrganisation = ['AAAA']
  specificTailOrganisation = ['AA','AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Matches LipidBlast
    GL.MA            :20, # Molecular ion
    GL.MA_s_FA      :100, # Present in lipidBlast
    GL.MA_s_allFA    :10, # Present in lipidBlast
    GL.HGA_s_H2O     :30, # Present in lipidBlast
    GL.FAH           :40}}# Present in lipidBlast

  def __init__(self, sn1, sn2, sn3, sn4):
    headgroup = GL.sn(mass=584.135365, type='Headgroup', chnops={'C':18, 'H':33, 'O':19, 'P':1},
    smiles='O(C3OC(CO)C(O)C(O)C3O)C1C(O)C(O)C(O('+sn2.smiles+'))C(OC2OC(CO('+sn1.smiles+'))C(O)C(O)C2O)C1OP(O)(=O)', hgtails = [sn1, sn2])
    super().__init__(AC4PIM2.adducts, sn1=sn3, sn2=sn4, sn3=headgroup)
    self.specificname=f"{self.lipid_class} {sn1.name}_{sn2.name}/{sn3.name}_{sn4.name}"

# ~ # ~ # ~ # ~ # ~ # ~ #

class BMP(GL.Glycerolipid):

  ##### http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
  tooltip = 'Bismonoacylglycerophosphate'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M+NH4]+":{
    GL.MA              :5,
    GL.MH_s_C3H9O6P    :5,
    GL.HGH_s_PO4     :100,
    GL.MH_s_HG_H2O   :100}
    }

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=172.013674, type='Headgroup', chnops={'C':3, 'H':9, 'O':6, 'P':1},
    smiles=sn1.inverseSmiles+'OCC(O)COP(=O)(O)', hgtails = [sn1])
    super().__init__(BMP.adducts, sn1=sn2, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class HBMP(GL.Glycerolipid):

  ##### http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
  tooltip = 'Hemibismonoacylglycerophosphate'
  tailOrganisation = ['AAA']
  specificTailOrganisation = ['A','AA']

  adducts = {  # adduct:{spectra}
  "[M+NH4]+":{
    GL.MA              :5,
    GL.HGH_s_PO4      :20,
    GL.MH_s_HG_H2O   :100,
    GL.MH_s_HG_FA     :20},

  "[M-H]-":{
    GL.MA            :100,
    GL.MA_s_FAk        :1,
    GL.MA_s_FA         :1,
    GL.HGA_s_H2O       :1,
    GL.FAH            :25}
    }

  def __init__(self, sn1, sn2, sn3):
    headgroup = GL.sn(mass=172.013674, type='Headgroup', chnops={'C':3, 'H':9, 'O':6, 'P':1},
    smiles=sn1.inverseSmiles+'OCC(O)COP(=O)(O)', hgtails = [sn1])
    super().__init__(HBMP.adducts, sn1=sn2, sn2=sn3, sn3=headgroup)
    self.specificname=f"{self.lipid_class} {sn1.name}/{sn2.name}_{sn3.name}"

# ~ # ~ # ~ # ~ # ~ # ~ #

class DLCL(GL.Glycerolipid):

  tooltip = 'Dilysocardiolipin'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
  "[M-H]-":{
    GL.MA            :100,
    GL.FAH            :10,
    GL.MA_C3H8O8P2_s_HG:0,
    GL.MA_PO4_s_HG    :20,
    GL.C3H7O5P_FA     :10,
    GL.C3H5O4P_FA     :10,
    GL.C3H6O5P        :50},

  "[M-2H]2-":{
    GL.MA              :0,
    GL.FAH           :100,
    GL.MA_C3H8O8P2_s_HG:0,
    GL.MH_PO4_s_HG     :0,
    GL.C3H7O5P_FA      :0,
    GL.C3H5O4P_FA      :0,
    GL.C3H6O5P       :100,
    GL.H2O4P           :0,
    GL.O3P             :5}
    }

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=326.016783, type='Headgroup', chnops={'C':6, 'H':16, 'O':11, 'P':2},
    smiles=sn1.inverseSmiles+'OCC(O)COP(=O)(O)OCC(O)COP(O)(=O)', hgtails = [sn1])
    super().__init__(DLCL.adducts, sn1=sn2, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class MLCL(GL.Glycerolipid):

  tooltip = 'Lysocardiolipin'
  tailOrganisation = ['AAA']
  specificTailOrganisation = ['A','AA']

  adducts = {  # adduct:{spectra}
  # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
  "[M-H]-":{
    GL.MA            :100,
    GL.FAH            :10,
    GL.MA_C3H8O8P2_s_HG:0,
    GL.MA_PO4_s_HG    :20,
    GL.C3H7O5P_FA     :10,
    GL.C3H5O4P_FA     :10,
    GL.C3H6O5P        :50},
    }

  def __init__(self, sn1, sn2, sn3):
    headgroup = GL.sn(mass=326.016783, type='Headgroup', chnops={'C':6, 'H':16, 'O':11, 'P':2},
    smiles=sn1.inverseSmiles+'OCC(O)COP(=O)(O)OCC(O)COP(O)(=O)', hgtails = [sn1])
    super().__init__(MLCL.adducts, sn1=sn2, sn2=sn3, sn3=headgroup)
    self.specificname=f"{self.lipid_class} {sn1.name}/{sn2.name}_{sn3.name}"

# ~ # ~ # ~ # ~ # ~ # ~ #

class CL(GL.Glycerolipid):

  ##### http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
  tooltip = 'Cardiolipin'
  tailOrganisation = ['AAAA']
  specificTailOrganisation = ['AA','AA']

  adducts = {  # adduct:{spectra}

  "[M-H]-":{
    GL.MA             :20,
    GL.FAH           :100,
    GL.MA_C3H8O8P2_s_HG:5,
    GL.MA_PO4_s_HG    :20,
    GL.C3H7O5P_FA     :10,
    GL.C3H5O4P_FA     :10,
    GL.C3H6O5P        :50},

  "[M-2H]2-":{
    GL.MA            :20,
    GL.FAH          :100,
    GL.C3H5O4P_FA    :10,
    GL.C3H6O5P       :10}}

  def __init__(self, sn1, sn2, sn3, sn4):
    headgroup = GL.sn(mass=326.016783, type='Headgroup', chnops={'C':6, 'H':16, 'O':11, 'P':2},
    smiles=sn1.inverseSmiles+'OCC(O'+sn2.smiles+')COP(=O)(O)OCC(O)COP(O)(=O)', hgtails = [sn1, sn2])
    super().__init__(CL.adducts, sn1=sn3, sn2=sn4, sn3=headgroup)
    self.specificname=f"{self.lipid_class} {sn1.name}_{sn2.name}/{sn3.name}_{sn4.name}"

# ~ # ~ # ~ # ~ # ~ # ~ #

class PA(GL.Glycerolipid):

  tooltip = 'Phosphatidic acid'
  tailOrganisation = ['AA']

  adducts = {  
  "[M-H]-":{ # "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
    GL.MA      :15, # Molecular ion
    GL.MH_s_FAk:15, # Present in lipidBlast
    GL.MH_s_FA :30, # Present in lipidBlast
    GL.FAH    :100, # Present in lipidBlast
    GL.C3H8O6P  :5,
    GL.C3H6O5P :30,
    GL.H2O4P    :5,
    GL.O3P      :5},

  "[M+H]+":{
    GL.MA          :90,
    GL.MA_s_HG_H2O:100,
    GL.MA_s_HG_FA  :10,
    GL.FAkH        :10},

  "[M+Na]+":{ # Matches LipidBlast
    GL.MA           :5, # Molecular ion
    GL.MA_s_H2O     :5, # Present in lipidblast
    GL.MA_s_HG_H2O:100, # Present in lipidblast
    GL.MH_s_HG_H2O :30, # Present in lipidblast
    GL.MA_s_FAk     :2, # Present in lipidblast
    GL.MA_s_FA      :2}} # Present in lipidblast

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=97.976894576, type='Headgroup', chnops={'H':3, 'O':4, 'P':1},
    smiles='O=P(O)(O)')
    super().__init__(PA.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class oPA(GL.Glycerolipid):

  tooltip = 'Ether/Plasmanyl Phosphatidic acid'
  tailOrganisation = ['OA']
  givenName = 'Ether-PA'

  adducts = {  
  "[M-H]-":{ # "[M-H]-" https://doi.org/10.1194/jlr.D003715 
    GL.MA         :15,
    GL.MA_s_H2O    :2,
    GL.MH_s_FAk   :15,
    GL.MH_s_FAcyl:100,
    GL.FAH         :2,
    GL.C3H8O6P     :0,
    GL.C3H6O5P     :0,
    GL.H2O4P       :0,
    GL.O3P         :0}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=97.976894576, type='Headgroup', chnops={'H':3, 'O':4, 'P':1},
    smiles='O=P(O)(O)')
    super().__init__(oPA.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PMe(GL.Glycerolipid):

  tooltip = 'Phosphatidylmethanol'
  tailOrganisation = ['AA']

  adducts = {  # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
  "[M-H]-":{
    GL.MA      :15,
    GL.MH_s_FAk:15,
    GL.MH_s_FA :30,
    GL.FAH    :100,
    GL.C4H10O6P :5,
    GL.C4H8O5P :30,
    GL.HGA      :5,
    GL.O3P      :5}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=111.992545, type='Headgroup', chnops={'C':1, 'H':5, 'O':4, 'P':1},
    smiles='O=P(OC)(O)')
    super().__init__(PMe.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PEt(GL.Glycerolipid):

  tooltip = 'Phosphatidylethanol'
  tailOrganisation = ['AA']

  adducts = { # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
              # https://link.springer.com/article/10.1007/s00216-010-3458-5 
  "[M-H]-":{
    GL.MA      :15,
    GL.MH_s_FAk:15,
    GL.MH_s_FA :30,
    GL.FAH    :100,
    GL.C4H10O6P :5,
    GL.C4H8O5P :30,
    GL.HGA      :5,
    GL.O3P      :5}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=126.008195, type='Headgroup', chnops={'C':2, 'H':7, 'O':4, 'P':1},
    smiles='O=P(OCC)(O)')
    super().__init__(PEt.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PC(GL.Glycerolipid):

  ##### "[M+H]+", "[M+Na]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  ##### "[M+H]+", "[M+Na]+"  10.1016/S1044-0305(03)00064-3
  ##### "[M+H]+" https://doi.org/10.1039/C5AY00776C 
  tooltip = 'Phosphatidylcholine'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}

  "[M+Hac-H]-":{ # https://doi.org/10.1016/j.ijms.2006.04.001
    GL.MA           :10,
    GL.M_s_CH3     :100,
    GL.M_s_CH3_FA   :2,
    GL.M_s_CH3_FAk  :2,    
    GL.FAH         :100,
    GL.C3H6O5P      :10,
    GL.C2O2H3       :50},

  "[M+CO3H]-":{ # https://doi.org/10.1038/s41467-021-23161-5, https://doi.org/10.1039/C9SC03521D
    GL.MA            :0,
    GL.MH          :100,
    GL.MA_s_TMA_H2O :20,
    GL.M_s_TMA      :20,
    GL.MH_s_TMA    :100,
    GL.C3H5O4P_FA    :0,
    GL.C3H7O5P_FA    :0,
    GL.M_s_FAk_TMAb  :0,
    GL.M_s_FA_TMAb :100,
    GL.MH_PO4_s_HG  :50,
    GL.FAH          :70,
    GL.C3H6O5P       :1},

  "[M+H]+":{ # Matches literature but not lipidblast
    GL.MA          :10, # Lipidblast seems to provide a completely
    GL.MA_s_H2O     :0, # different spectra which is not consistant
    GL.MA_s_TMA     :0, # with literature. Perhaps a result of
    GL.MA_s_HG_H2O  :0, # instrumentation.
    GL.MA_s_FA_TMA  :0,
    GL.MH_s_FAk     :1,
    GL.MH_s_FA      :1,
    GL.C5H15NO4P  :100,
    GL.C5H12N       :2},

  "[M+Na]+":{ # Looks Good, Matches lipidblast
    GL.MA         :100, # Molecular ion
    GL.MA_s_TMA    :30, # Present in lipidblast
    GL.MA_s_HG_H2O :50, # Present in lipidblast
    GL.MH_s_HG_H2O :10,
    GL.MA_s_FA_TMA  :2, # Present in lipidblast
    GL.MA_s_FAk     :1, # Present in lipidblast
    GL.MA_s_FA      :1, # Present in lipidblast
    GL.MH_s_FA      :1,
    GL.FAkH         :1,
    GL.C5H15NO4P    :5,
    GL.C2H5NaO4P   :50}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=183.066044, type='Headgroup', chnops={'C':5, 'H':14, 'N':1, 'O':4, 'P':1},
    smiles='C[N+](C)(C)CCOP([O-])(=O)')
    # headgroup mass has -H to maintain neutral charge
    super().__init__(PC.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class oPC(GL.Glycerolipid):

  ##### http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
  tooltip = 'Ether/Plasmanyl Phosphatidylcholine'
  tailOrganisation = ['OA']
  givenName = 'Ether-PC'

  adducts = {  # adduct:{spectra}

  "[M+Hac-H]-":{ # https://doi.org/10.1016/j.ijms.2006.04.001
    GL.MA           :10,
    GL.M_s_CH3     :100,
    GL.M_s_CH3_FA   :2,
    GL.M_s_CH3_FAk  :2,    
    GL.FAH         :100,
    GL.C3H6O5P      :10,
    GL.C2O2H3       :50},

  "[M+H]+":{# http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html 
    GL.MA          :10,
    GL.MA_s_HG_H2O  :1,
    GL.C5H15NO4P  :100,
    GL.C5H12N       :2},

  "[M+Na]+":{# https://doi.org/10.1002/jms.491
    GL.MA           :100,
    GL.MA_s_TMA      :50,
    GL.MA_s_HG_H2O   :20,
    GL.MH_s_HG_H2O    :5,
    GL.MH_s_HG_FA_H2O :1,
    GL.MH_s_FAk       :1,
    GL.MH_s_FA        :1,
    GL.FAkH           :1,
    GL.C5H15NO4P      :5,
    GL.C2H5NaO4P     :25}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=183.066044, type='Headgroup', chnops={'C':5, 'H':14, 'N':1, 'O':4, 'P':1},
    smiles='C[N+](C)(C)CCOP([O-])(=O)')
    # headgroup mass has -H to maintain neutral charge
    super().__init__(oPC.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class pPC(GL.Glycerolipid):

  ##### https://doi.org/10.1002/jms.491
  tooltip = 'Plasmenyl Phosphatidylcholine'
  tailOrganisation = ['PA']
  givenName = 'Plasmenyl-PC'

  adducts = {  # adduct:{spectra}
  "[M+H]+":{# https://doi.org/10.1002/jms.491
    GL.MA          :10,
    GL.MA_s_FAk     :2,
    GL.MA_s_FA      :2,
    GL.C5H15NO4P  :100,
    GL.C5H12N       :2},

  "[M+Na]+":{# https://doi.org/10.1002/jms.491
    GL.MA           :100,
    GL.MA_s_TMA      :50,
    GL.MA_s_HG_H2O    :5,
    GL.MH_s_HG_H2O    :5,
    GL.MH_s_HG_FA_H2O :5,
    GL.MH_s_FAk       :0,
    GL.MH_s_FA        :0,
    GL.FAkH           :0,
    GL.C5H15NO4P      :5,
    GL.C2H5NaO4P     :50}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=183.066044, type='Headgroup', chnops={'C':5, 'H':14, 'N':1, 'O':4, 'P':1},
    smiles='C[N+](C)(C)CCOP([O-])(=O)')
    # headgroup mass has -H to maintain neutral charge
    super().__init__(pPC.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PE(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  ##### "[M+H]+", "[M-H]-" https://doi.org/10.1039/C5AY00776C
  tooltip = 'Phosphatidylethanolamine'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Looks Good, Matches lipidblast
    GL.MA       :15, # Molecular ion
    GL.MH_s_FAk  :5, # Present in lipidblast
    GL.MH_s_FA   :2, # Present in lipidblast
    GL.FAH     :100, # Present in lipidblast
    GL.FAH_PUFA :20, # Present in lipidex
    GL.C5H11NO5P :5,
    GL.C3H6O5P  :20,
    GL.C2H7NO4P  :3,
    GL.H2O4P     :5,
    GL.O3P       :5},
  
  "[M+H]+":{ # Looks Good, Matches lipidblast
    GL.MA          :10, # Molecular ion
    GL.MA_s_H2O     :1, # Present in lipidblast
    GL.MA_s_HG_H2O:100, # Present in lipidblast
    GL.MH_s_FAk     :1, # Present in lipidblast
    GL.MH_s_FA      :1, # Present in lipidblast
    GL.MA_s_HG_FA   :1,
    GL.FAkH         :1},# Present in lipidblast
    
  "[M+Na]+":{ # Looks Good, Matches lipidblast
    GL.MA          :10, # Molecular ion
    GL.MA_s_AZD   :100, # Present in lipidblast
    GL.MH_s_HG      :0,
    GL.MA_s_HG_H2O :40, # Present in lipidblast
    GL.MH_s_HG_H2O  :0,
    GL.MA_s_FA_AZD  :1, # Present in lipidblast
    GL.MH_s_FA      :0,
    GL.MA_s_FA      :0,
    GL.FAkH         :0}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1},
    smiles='NCCOP(O)(=O)')
    super().__init__(PE.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class oPE(GL.Glycerolipid):

  tooltip = 'Ether/Plasmanyl Phosphatidylethanolamine'
  tailOrganisation = ['OA']
  givenName = 'Ether-PE'

  adducts = { 
  "[M-H]-":{ # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
    GL.MA       :25,
    GL.MH_s_FAk :10,
    GL.MH_s_FA   :2,
    GL.FAH     :100,
    GL.FAH_PUFA :20, # Present in lipidex
    GL.C5H11NO5P:20,
    GL.C3H6O5P   :0,
    GL.C2H7NO4P  :1,
    GL.H2O4P     :0,
    GL.O3P       :0},
  
  "[M+H]+":{ # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
    GL.MA              :2,
    GL.MA_s_HG_H2O   :100,
    GL.FAkH            :2}}#,
    
  #"[M+Na]+":{}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1},
    smiles='NCCOP(O)(=O)')
    super().__init__(oPE.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class pPE(GL.Glycerolipid):

  tooltip = 'Plasmenyl Phosphatidylethanolamine'
  tailOrganisation = ['PA']
  givenName = 'Plasmenyl-PE'

  adducts = { 
  # https://doi.org/10.1016/j.aca.2012.05.035
  # https://doi.org/10.1016/j.plipres.2021.101111
  "[M-H]-":{
    GL.MA       :25,
    GL.MH_s_FAk :10,
    GL.MH_s_FA  :20,
    GL.FAH     :100,
    GL.FAH_PUFA :20, # Present in lipidex
    GL.C5H11NO5P:20,
    GL.C3H6O5P   :5,
    GL.C2H7NO4P  :5,
    GL.H2O4P     :5,
    GL.O3P       :5},
  
  # https://doi.org/10.1016/j.jasms.2004.07.009
  # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
  "[M+H]+":{
    GL.MA              :15,
    GL.MA_s_HG_H2O     :10,
    GL.MA_s_FA          :1,
    GL.HGA_FA_s_H2O    :50,
    GL.C3H7O2_FA_s_H2O:100,
    GL.HGA_FA_s_H2O_PO4 :5,
    GL.FAkA             :1
    }}#,
    
  #"[M+Na]+":{}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1},
    smiles='NCCOP(O)(=O)')
    super().__init__(pPE.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class MMPE(GL.Glycerolipid):

  tooltip = 'N-methyl-phosphatidylethanolamine'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    GL.MA          :15,
    GL.FAH        :100,
    GL.MH_s_FAk    :10,
    GL.MH_s_FA     :10,
    GL.C3H6O5P     :10,
    GL.H2O4P        :5,
    GL.O3P          :5},
  
  "[M+H]+":{
    GL.MA           :1,
    GL.MA_s_HG_H2O:100,
    GL.MA_s_HG_FA  :20,
    GL.C3H11NO4P    :5}}
    
  #"[M+Na]+":{
  #  GL.MA          :10}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=155.034744, type='Headgroup', chnops={'C':3, 'H':10, 'N':1, 'O':4, 'P':1},
    smiles='N(C)CCOP(O)(=O)')
    super().__init__(MMPE.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class DMPE(GL.Glycerolipid):

  tooltip = 'N,N-dimethyl-phosphatidylethanolamine'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    GL.MA          :15,
    GL.FAH        :100,
    GL.MH_s_FAk    :10,
    GL.MH_s_FA     :10,
    GL.C4H11NO4P    :1,
    GL.C3H6O5P     :10,
    GL.H2O4P        :5,
    GL.O3P          :5},
  
  "[M+H]+":{
    GL.MA           :1,
    GL.MA_s_HG_H2O:100,
    GL.MA_s_HG_FA  :20,
    GL.FAkH         :5,
    GL.C4H13NO4P   :50}}
    
  #"[M+Na]+":{
  #  GL.MA          :10}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=169.050394, type='Headgroup', chnops={'C':4, 'H':12, 'N':1, 'O':4, 'P':1},
    smiles='N(C)(C)CCOP(O)(=O)')
    super().__init__(DMPE.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class NAPE(GL.Glycerolipid):

  tooltip = 'N-Acyl-phosphatidylethanolamine'
  tailOrganisation = ['A','AA']
  givenName = 'PE-N(FA)'

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    GL.MA            :15,
    GL.MH_s_FAk       :5,
    GL.MH_s_FA       :15,
    GL.MH_PO4_s_HG_FA :1,
    GL.MA_s_FA_FAk    :1,
    GL.FAH_b        :100,
    GL.HGA            :2,
    GL.HGA_s_H2O      :1,
    GL.C5H11NO5P      :0,
    GL.C3H6O5P       :10,
    GL.C2H7NO4P       :0,
    GL.H2O4P          :5,
    GL.O3P            :5}
    }

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1, sn2, sn3):
    headgroup = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1},
    smiles=sn1.inverseSmiles+'NCCOP(O)(=O)', hgtails=[sn1])
    super().__init__(NAPE.adducts, sn1=sn2, sn2=sn3, sn3=headgroup)
    self.name=f"PE-N({sn1.name}) {sn2.name}_{sn3.name}"

# ~ # ~ # ~ # ~ # ~ # ~ #

class PG(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  tooltip = 'Phosphatidylglycerol'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Looks Good, Matches lipidblast
    GL.MA            :15, # Molecular ion
    GL.MH_s_FAk      :20, # Present in lipidblast
    GL.MH_s_FA       :20, # Present in lipidblast
    GL.FAH          :100, # Present in lipidblast
    GL.MH_PO4_s_HG_FA:20, # Present in lipidblast
    GL.C6H12O7P       :2,
    GL.C3H8O6P        :2,
    GL.C3H6O5P       :10,
    GL.H2O4P          :5,
    GL.O3P            :5},

  "[M+H]+":{ # Matches literature
    GL.MA          :10,
    GL.MA_s_HG_H2O:100,
    GL.MA_s_HG_FA  :10,
    GL.FAkH        :10}}#,
    
  #"[M+Na]+":{
  #  GL.MA          :10,
  #  GL.HG_NL_H2O_A:100,
  #  GL.HG_FA_NL_A  :10,
  #  GL.FAkH        :10,
  #  GL.HGA         :80}
  # }

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=172.013674, type='Headgroup', chnops={'C':3, 'H':9, 'O':6, 'P':1},
    smiles='OCC(O)COP(O)(=O)')
    super().__init__(PG.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PI(GL.Glycerolipid):

  ##### "[M-H]-" https://doi.org/10.1039/C5AY00776C
  tooltip = 'Phosphatidylinositol'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Looks Good, Matches lipidblast
    GL.MA            :2, # Molecular ion
    GL.MH_s_FAk      :1, # Present in lipidblast
    GL.MH_s_FA       :2, # Present in lipidblast
    GL.MH_PO4_s_HG_FA:1, # Present in lipidblast
    GL.FAH         :100, # Present in lipidblast 
    GL.C9H16O10P     :5, 
    GL.C9H14O9P      :5,
    GL.C6H12O9P      :5, 
    GL.C6H10O8P     :40,
    GL.C6H8O7P      :15, 
    GL.C3H6O5P      :20,
    GL.H2O4P         :5,
    GL.O3P           :5},

  "[M+Na-2H]-":{
    GL.MA         :0,
    GL.MH_s_FA   :35,
    GL.MA_s_2FAk :35,
    GL.FAH      :100,
    GL.C3H5O4P_FA :1,
    GL.C9H16O10P:100,
    GL.HGA_s_H2O  :5,
    GL.C9H14O9P  :35, 
    GL.C6H10O8P   :5,
    GL.C6H8O7P   :80,
    GL.C3H7NaO6P :35,
    GL.C3H6O5P  :100,
    GL.O3P       :30}, 

  "[M+H]+":{
    GL.MA           :5,
    GL.MH_s_HG_H2O:100,
    GL.MH_s_HG_FA  :25},

  "[M+Na]+":{ # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
    GL.MA           :5,
    GL.MH           :0,
    GL.MH_s_HG_H2O:100,
    GL.HGA        :100},

  "[M+NH4]+":{ # https://doi.org/10.1002/pld3.183 Figure S4   
    GL.MA           :5,
    GL.MH           :5,
    GL.MH_s_HG_H2O:100,
    GL.MH_s_HG_FA  : 5}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=260.029718, type='Headgroup', chnops={'C':6, 'H':13, 'O':9, 'P':1},
    smiles='OC1C(O)C(O)C(O)C(O)C1OP(O)(=O)')
    super().__init__(PI.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class oPI(GL.Glycerolipid):

  ##### "[M-H]-" https://doi.org/10.1194/jlr.D003715
  tooltip = 'Ether/Plasmanyl Phosphatidylinositol'
  tailOrganisation = ['OA']
  givenName = 'Ether-PI'

  adducts = {  # adduct:{spectra}

  "[M-H]-":{
    GL.MA           :100,
    GL.MA_PO4_s_HG    :1,
    GL.MH_s_FAk      :10,
    GL.MH_s_FAcyl    :25,
    GL.MH_PO4_s_HG_FAk:2,
    GL.MH_PO4_s_HG_FA:50,
    GL.FAH           :10,
    GL.C9H16O10P      :0, 
    GL.C9H14O9P       :0,
    GL.C6H12O9P       :4, 
    GL.C6H10O8P       :5,
    GL.C6H8O7P        :1, 
    GL.C3H6O5P        :0,
    GL.H2O4P          :0,
    GL.O3P            :0}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=260.029718, type='Headgroup', chnops={'C':6, 'H':13, 'O':9, 'P':1},
    smiles='OC1C(O)C(O)C(O)C(O)C1OP(O)(=O)')
    super().__init__(oPI.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PIP(GL.Glycerolipid):

  tooltip = 'Phosphatidylinositol phosphate'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  #"[M+Na-2H]-":{
  #  GL.MA            :2,
  #  GL.MH_s_PO3      :2,
  #  GL.MH_s_FA       :1,
  #  GL.MH_s_FA_H2O   :1,
  #  GL.MH_s_FAk_PO3  :1,
  #  GL.MH_s_FA_PO3   :1,
  #  GL.HG_FA_NL_B    :5,
  #  GL.FAH          :50,
  #  GL.C6H10NaO11P2:100, 
  #  GL.C6H8NaO10P2  :10,
  #  GL.C6H12O9P      :5, 
  #  GL.C6H10O8P     :50,
  #  GL.C6H8O7P      :20, 
  #  GL.C3H6O5P      :10},

  "[M-H]-":{ # Looks Good 10.1016/S1044-0305(00)00172-0 
    GL.MH             :20,
    GL.MH_s_H2O       :10,
    GL.MH_s_PO3        :2,
    GL.MH_s_PO4       :10,
    GL.MH_P2O6_s_HG    :2,
    GL.MH_s_FAk        :1,
    GL.MH_s_FA         :4,
    GL.MH_s_FA_H2O     :1,
    GL.MH_s_FAk_PO3    :4,
    GL.MH_s_FA_PO3     :4,
    GL.MH_PO4_s_HG_FA :10,
    GL.MH_PO4_s_HG_FAk:10,
    GL.FAH           :100,
    GL.C6H11O11P2     :75, 
    GL.C6H9O10P2      :15,
    GL.C6H12O9P        :5, 
    GL.C6H10O8P       :75,
    GL.C6H8O7P        :75, 
    GL.C3H6O5P        :40,
    GL.H2O4P           :5,
    GL.O3P             :5},

  "[M-2H]2-":{ # Looks Good 10.1016/S1044-0305(00)00172-0 
    GL.MA             :50,
    GL.MH_s_H2O        :0,
    GL.MH_s_PO3        :8,
    GL.MH_P2O6_s_HG    :0,
    GL.MA_s_FAk       :20,
    GL.MH_s_FA        :20,
    GL.MH_s_FA_H2O     :0,
    GL.MH_s_FAk_PO3    :2,
    GL.MH_s_FA_PO3     :1,
    GL.MH_PO4_s_HG_FA  :0,
    GL.MH_PO4_s_HG_FAk:10,
    GL.FAH           :100,
    GL.C6H11O11P2      :5, 
    GL.C6H9O10P2       :5,
    GL.C6H12O9P        :5, 
    GL.C6H10O8P       :75,
    GL.C6H8O7P        :20, 
    GL.C3H6O5P        :10,
    GL.H2O4P           :5,
    GL.O3P             :5}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=339.996048, type='Headgroup', chnops={'C':6, 'H':14, 'O':12, 'P':2},
    smiles='OC1C(O)C(O)C(OP(=O)(O)O)C(O)C1OP(O)(=O)')
    super().__init__(PIP.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PIP2(GL.Glycerolipid):

  tooltip = 'Phosphatidylinositol diphosphate'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra} 
  #"[M+Na-2H]-":{
  #},

  "[M-H]-":{ # Looks Good 10.1016/S1044-0305(00)00172-0 
    GL.MH             :10,
    GL.MH_s_H2O       :10,
    GL.MH_s_PO3       :10,
    GL.MH_s_PO4       :75,
    GL.MH_s_PO4_H2O   :10,
    GL.MH_P2O6_s_HG   :10,
    GL.MH_s_FAk        :0,
    GL.MH_s_FA         :0,
    GL.MH_s_FA_H2O     :1,
    GL.MH_s_FAk_PO3    :0,
    GL.MH_s_FA_PO3     :0,
    GL.MH_PO4_s_HG_FA :10,
    GL.MH_PO4_s_HG_FAk:10,
    GL.FAH           :100,
    GL.C6H11O11P2     :10, 
    GL.C6H9O10P2      :15,
    GL.C6H12O9P        :5, 
    GL.C6H10O8P       :75,
    GL.C6H8O7P        :75, 
    GL.C3H6O5P        :40,
    GL.H2O4P           :5,
    GL.O3P             :5},

  "[M-2H]2-":{ # Looks Good 10.1016/S1044-0305(00)00172-0 
    GL.MA             :50,
    GL.MH_s_H2O        :0,
    GL.MH_s_PO3       :50,
    GL.MH_P2O6_s_HG    :0,
    GL.MA_s_FAk       :20,
    GL.MH_s_FA        :20,
    GL.MH_s_FA_H2O     :0,
    GL.MH_s_FAk_PO3    :2,
    GL.MH_s_FA_PO3     :1,
    GL.MH_PO4_s_HG_FAk:10,
    GL.FAH           :100,
    GL.C6H11O11P2      :5, 
    GL.C6H9O10P2       :5,
    GL.C6H12O9P        :5, 
    GL.C6H10O8P       :75,
    GL.C6H8O7P        :20, 
    GL.C3H6O5P        :10,
    GL.H2O4P           :5,
    GL.O3P             :5}}

  #"[M-3H]3-":{
  #}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=419.962378, type='Headgroup', chnops={'C':6, 'H':15, 'O':15, 'P':3},
    smiles='OC1C(OP(=O)(O)O)C(O)C(OP(=O)(O)O)C(O)C1OP(O)(=O)')
    super().__init__(PIP2.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PIPMe3(GL.Glycerolipid):

  tooltip = 'Methylated Phosphatidylinositol phosphate'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra} 
  # https://doi.org/10.1194/jlr.D069989
  "[M+H]+":{
    GL.MA             :10,
    GL.MA_s_HG_H2O   :100,
    GL.HGH             :0,
    GL.MA_s_HG_FA     :10},

  "[M+NH4]+":{
    GL.MA             :10,
    GL.MH              :5,
    GL.MH_s_HG_H2O   :100,
    GL.HGA             :0,
    GL.HGH             :0,
    GL.MA_s_HG_FA      :0,
    GL.MH_s_HG_FA     :10}
    }

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=382.042998, type='Headgroup', chnops={'C':9, 'H':20, 'O':12, 'P':2},
    smiles='OC1C(O)C(O)C(OP(=O)(OC)OC)C(O)C1OP(OC)(=O)')
    super().__init__(PIPMe3.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PI2PMe5(GL.Glycerolipid):

  tooltip = 'Methylated Phosphatidylinositol diphosphate'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra} 
  # https://doi.org/10.1194/jlr.D069989
  "[M+H]+":{
    GL.MA             :10,
    GL.MA_s_HG_H2O   :100,
    GL.HGH             :0,
    GL.MA_s_HG_FA     :10},

  "[M+NH4]+":{
    GL.MA             :10,
    GL.MH              :5,
    GL.MH_s_HG_H2O   :100,
    GL.HGA             :0,
    GL.HGH             :0,
    GL.MA_s_HG_FA      :0,
    GL.MH_s_HG_FA     :10}
    }

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=490.040628, type='Headgroup', chnops={'C':11, 'H':25, 'O':15, 'P':3},
    smiles='OC1C(OP(=O)(OC)OC)C(O)C(OP(=O)(OC)OC)C(O)C1OP(OC)(=O)')
    super().__init__(PI2PMe5.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PS(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  ##### "[M+H]+" https://doi.org/10.1039/C5AY00776C
  tooltip = 'Phosphatidylserine'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Looks Good, Matches lipidblast
    GL.MA             :15, # Molecular ion
    GL.MH_PO4_s_HG   :100, # Present in lipidblast
    GL.MH_PO4_s_HG_FAk:25, # Present in lipidblast
    GL.MH_PO4_s_HG_FA :50, # Present in lipidblast
    GL.FAH            :50, # Present in lipidblast
    GL.C3H6O5P        :30,
    GL.H2O4P           :5,
    GL.O3P             :5},

  "[M+H]+":{ # Looks Good, Matches lipidblast
    GL.MA          :10, # Molecular ion
    GL.MH_PO4_s_HG  :1, # Present in lipidblast
    GL.MA_s_HG_H2O:100, # Present in lipidblast
    GL.MA_s_FAk     :1, # Present in lipidblast
    GL.MA_s_FA      :1, # Present in lipidblast
    GL.MA_s_HG_FA   :1,
    GL.FAkH         :1},
    
  "[M+Na]+":{ # Looks Good, Matches lipidblast
    GL.MA          :10, # Molecular ion
    GL.MA_PO4_s_HG:100, # Present in lipidblast
    GL.MA_s_PO4    :30, # Present in lipidblast
    GL.MA_s_HG_H2O:100, # Present in lipidblast
    GL.MH_s_HG_H2O :50, # Present in lipidblast
    GL.MA_s_HG_FA   :0,
    GL.FAkH         :0,
    GL.HGA          :0}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=185.008923, type='Headgroup', chnops={'C':3, 'H':8, 'N':1, 'O':6, 'P':1},
    smiles='[O-]C(=O)C([NH3+])COP(O)(=O)')
    super().__init__(PS.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class oPS(GL.Glycerolipid):

  tooltip = 'Ether/Plasmanyl Phosphatidylserine'
  tailOrganisation = ['OA']
  givenName = 'Ether-PS'

  adducts = {  # adduct:{spectra} 
  "[M-H]-":{# https://doi.org/10.1016/j.aca.2012.05.035
    GL.MA             :10,
    GL.MH_PO4_s_HG    :10,
    GL.MH_PO4_s_HG_FA:100,
    GL.MH_PO4_s_HG_FAk:25,
    GL.FAH             :1,
    GL.C3H6O5P        :10,
    GL.H2O4P           :5,
    GL.O3P             :5
  }}#,

  #"[M+H]+":{},
    
  #"[M+Na]+":{}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=185.008923, type='Headgroup', chnops={'C':3, 'H':8, 'N':1, 'O':6, 'P':1},
    smiles='[O-]C(=O)C([NH3+])COP(O)(=O)')
    super().__init__(oPS.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class pPS(GL.Glycerolipid):

  tooltip = 'Plasmenyl Phosphatidylserine'
  tailOrganisation = ['PA']
  givenName = 'Plasmenyl-PS'

  adducts = {  # adduct:{spectra} 
  "[M-H]-":{# https://doi.org/10.1016/j.aca.2012.05.035
    GL.MA             :10,
    GL.MH_PO4_s_HG    :10,
    GL.MH_PO4_s_HG_FA:100,
    GL.MH_PO4_s_HG_FAk:25,
    GL.FAH             :1,
    GL.C3H6O5P        :10,
    GL.H2O4P           :5,
    GL.O3P             :5
  }}#,

  #"[M+H]+":{},
    
  #"[M+Na]+":{}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=185.008923, type='Headgroup', chnops={'C':3, 'H':8, 'N':1, 'O':6, 'P':1},
    smiles='[O-]C(=O)C([NH3+])COP(O)(=O)')
    super().__init__(pPS.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PT(GL.Glycerolipid):

  ##### "[M-H]-" https://doi.org/10.1194/jlr.D003715
  tooltip = 'Phosphatidylthreonine'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    GL.MA             :15,
    GL.MH_PO4_s_HG   :100,
    GL.MH_PO4_s_HG_FAk:15,
    GL.MH_PO4_s_HG_FA :15,
    GL.FAH            :50,
    GL.C3H6O5P        :30,
    GL.H2O4P           :5,
    GL.O3P             :5}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=199.024573, type='Headgroup', chnops={'C':4, 'H':10, 'N':1, 'O':6, 'P':1},
    smiles='[O-]C(=O)C([NH3+])C(C)OP(O)(=O)')
    super().__init__(PT.adducts, sn1, sn2, headgroup)

# # ~ # ~ # ~ # ~ # ~ # ~ #

# class PAsp(GL.Glycerolipid):

#   ##### "[M-H]-" https://doi.org/10.1016/j.jchromb.2011.04.033
#   tooltip = 'Phosphatidyl-L-Aspartate? (Not confirmed)'
#   tailOrganisation = ['AA']

#   adducts = {  # adduct:{spectra}
#   "[M-H]-":{
#     GL.MA             :1,
#     GL.MH_PO4_s_HG    :1,
#     GL.MH_PO4_s_HG_FAk:1,
#     GL.MH_PO4_s_HG_FA :5,
#     GL.FAH           :100,
#     GL.C3H6O5P        :50,
#     GL.H2O4P           :1,
#     GL.O3P            :30}}

#   # sn3 = headgroup
#   def __init__(self, sn1, sn2):
#     headgroup = GL.sn(mass=213.003838, type='Headgroup', chnops={'C':4, 'H':8, 'N':1, 'O':7, 'P':1},
#     smiles='[O-]C(=O)CC([NH3+])C(=O)OP(O)(=O)')
#     super().__init__(PAsp.adducts, sn1, sn2, headgroup)

# # ~ # ~ # ~ # ~ # ~ # ~ #

# class PHpa(GL.Glycerolipid):

#   ##### "[M-H]-" 
#   tooltip = 'Phosphatidyl Hydroxypyruvic acid? (Not confirmed)'
#   tailOrganisation = ['AA']

#   adducts = {  # adduct:{spectra}
#   "[M-H]-":{
#     GL.MA             :1,
#     GL.MH_PO4_s_HG    :1,
#     GL.MH_PO4_s_HG_FAk:1,
#     GL.MH_PO4_s_HG_FA :5,
#     GL.FAH           :100,
#     GL.C3H6O5P        :50,
#     GL.H2O4P           :1,
#     GL.O3P            :30}}

#   # sn3 = headgroup
#   def __init__(self, sn1, sn2):
#     headgroup = GL.sn(mass=183.977289, type='Headgroup', chnops={'C':3, 'H':5, 'O':7, 'P':1},
#     smiles='OC(=O)CC(=O)OP(O)(=O)')
#     super().__init__(PHpa.adducts, sn1, sn2, headgroup)

# # ~ # ~ # ~ # ~ # ~ # ~ #

class PPA(GL.Glycerolipid):

  tooltip = 'Pyrophosphatidic acid'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # https://www.lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=779&LM_ID=LMGP11010002&TRACK_ID=235
    GL.MA     :100,
    GL.MH_s_H2O:40,
    GL.MH_s_FA  :2,
    GL.MH_s_FAk :2,
    GL.FAH      :5,
    GL.HO6P2   :50}}#,

  #"[M-2H]2-":{
  #  GL.MA     :10
  # },

  #"[M-3H]3-":{
  #  GL.MA     :10}
  # }#,
 
  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=177.943224, type='Headgroup', chnops={'H':4, 'O':7, 'P':2},
    smiles='O=P(O)(O)OP(=O)(O)')
    super().__init__(PPA.adducts, sn1, sn2, headgroup)

class CDPDG(GL.Glycerolipid):

  tooltip = "Cytidine-5'-diphosphate diacylglycerol"
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    # https://www.lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=952&LM_ID=LMGP13010003&TRACK_ID=265
    # https://www.lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=952&LM_ID=LMGP13010003&TRACK_ID=267
    # https://www.lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=952&LM_ID=LMGP13010003&TRACK_ID=266
    GL.MA            :100,
    GL.MH_s_H2O        :0,
    GL.MH_s_HCNO       :2,
    GL.MH_s_C4H7N3O    :2,
    GL.MA_P2O6_s_HG   :50,
    GL.MA_s_FA         :5,
    GL.MA_PO4_s_HG     :5,
    GL.MA_P2O6_s_HG_FAk:5,
    GL.MA_P2O6_s_HG_FA :5,
    GL.HGA_s_H2O      :25,
    GL.HGH_s_PO3       :5,
    GL.HGH_s_PO4      :10,
    GL.C5H7O9P2       :10,
    GL.FAH             :5,
    GL.C3H5O7P2        :1,
    GL.HO6P2           :1}
    }#,

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=403.01818, type='Headgroup', chnops={'C':9, 'H':15, 'N':3, 'O':11, 'P':2},
    smiles='OC1C(O)C(N2C=CC(N)=NC2=O)OC1COP(=O)(O)OP(=O)(O)')
    super().__init__(CDPDG.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class LDGCC(GL.Glycerolipid): # Diacylglyceryl-3-O-carboxyhydroxymethylcholine

  tooltip = 'Lyso-3-O-carboxyhydroxymethylcholine' 
  tailOrganisation = ['A']

  adducts = {  # adduct:{spectra}
  "[M+H]+":{ #  http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
    GL.MA        :100,
    GL.C6H14NO2   :10,
    GL.C5H14NO    :10}}

  #"[M+NH4]+":{    
  #  GL.MA          :5}
  #  }

  def __init__(self, sn1):
    headgroup = GL.sn(mass=177.100108, type='Headgroup', chnops={'C':7, 'H':15, 'N':1, 'O':4},
    smiles='C[N+](C)(C)CCOC(C(=O)([O-]))')
    # headgroup mass has -H to maintain neutral charge
    super().__init__(LDGCC.adducts, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class LDGTS(GL.Glycerolipid): # LysoN-trimethylhomoserine

  tooltip = 'Lyso-N-trimethylhomoserine' 
  tailOrganisation = ['A']

  adducts = {  # adduct:{spectra}
  "[M+H]+":{ #  https://doi.org/10.1002/rcm.7847
    GL.MA          :50,
    GL.MA_s_H2O    :20,
    GL.MA_s_allFAk:100,
    GL.C7H14N1O2   :5}}

  #"[M+NH4]+":{    
  #  GL.MA          :5}
  #  }

  def __init__(self, sn1):
    headgroup = GL.sn(mass=161.105193, type='Headgroup', chnops={'C':7, 'H':15, 'N':1, 'O':3},
    smiles='[O-]C(=O)C([N+](C)(C)(C))CC')
    # headgroup mass has -H to maintain neutral charge
    super().__init__(LDGTS.adducts, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class LPA(GL.Glycerolipid):

  ##### [M-H]- = https://doi.org/10.1002/lipd.12172, Not a fragmentation study!
  ##### [M-H]- = https://doi.org/10.1016/j.jchromb.2010.03.030, Neither...
  tooltip = 'Lyso-phosphatidic acid'
  tailOrganisation = ['A']
  
  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Matches lipidblast
    GL.MA      :10, # Present in lipidblast
    GL.FAH     :10, 
    GL.MH_s_FAk :5,
    GL.C3H6O5P:100, # Present in lipidblast
    GL.H2O4P    :2,
    GL.O3P     :10} # Present in lipidblast
    }

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    headgroup = GL.sn(mass=97.976895, type='Headgroup', chnops={'H':3, 'O':4, 'P':1},
    smiles='O=P(O)(O)')
    super().__init__(LPA.adducts, sn1=sn1, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class LPC(GL.Glycerolipid):

  ##### https://doi.org/10.1016/S1044-0305(99)00158-0
  ##### LipidBlast DOI:10.1038/nmeth.2551
  tooltip = 'Lyso-phosphatidylcholine'
  tailOrganisation = ['A']

  adducts = {  # adduct:{spectra}

  "[M+Hac-H]-":{ # https://doi.org/10.1016/j.ijms.2006.04.001
    GL.MA          :10,
    GL.M_s_CH3     :50,
    GL.FAH        :100,
    GL.C3H6O5P      :5,
    GL.C2O2H3      :10},

  "[M+H]+":{ # Matches lipidblast
    GL.MA        :80, # Molecular ion
    GL.MA_s_H2O  :20, # Present in lipidblast
    GL.MH_s_FAk   :1,
    GL.MH_s_FA    :5, # Present in lipidblast
    GL.FAkH       :1, # Present in lipidblast
    GL.C5H15NO4P:100, # Present in lipidblast
    GL.C5H14NO   :80,
    GL.C5H12N     :2}}#,

  #"[M+Na]+":{ # Assuming [M+Na] fragments like [M+Li]
  #  GL.MA          :10,
  #  GL.MA_s_TMA     :5,
  #  GL.HG_NL_H2O_A :50,
  #  GL.HG_NL_H2O_C:100,
  #  GL.MA_s_FA_TMA :10,
  #  GL.MH_s_FA     :10,
  #  GL.MA_s_FA     :10,
  #  GL.FAkH        :5}
  # }

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    headgroup = GL.sn(mass=183.066044, type='Headgroup', chnops={'C':5, 'H':14, 'N':1, 'O':4, 'P':1},
    smiles='C[N+](C)(C)CCOP([O-])(=O)')
    # headgroup mass has -H to maintain neutral charge
    super().__init__(LPC.adducts, sn1=sn1, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class oLPC(GL.Glycerolipid):

  tooltip = 'Ether/Plasmanyl Lyso Phosphatidylcholine'
  tailOrganisation = ['O']
  givenName = 'Ether-LPC'

  adducts = {  # adduct:{spectra}

  "[M+H]+":{
    GL.MA        :80, # Molecular ion
    GL.MA_s_H2O   :0,
    GL.C5H15NO4P :10,
    GL.C5H14NO  :100,
    GL.C5H12N     :2}}#,

  # sn3 = headgroup
  def __init__(self, sn1):
    headgroup = GL.sn(mass=183.066044, type='Headgroup', chnops={'C':5, 'H':14, 'N':1, 'O':4, 'P':1},
    smiles='C[N+](C)(C)CCOP([O-])(=O)')
    # headgroup mass has -H to maintain neutral charge
    super().__init__(oLPC.adducts, sn1=sn1, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class pLPC(GL.Glycerolipid):

  tooltip = 'Plasmenyl Lyso Phosphatidylcholine'
  tailOrganisation = ['P']
  givenName = 'Plasmenyl-LPC'

  adducts = {  # adduct:{spectra}

  "[M+H]+":{
    GL.MA        :80, # Molecular ion
    GL.MA_s_H2O   :0,
    GL.C5H15NO4P :10,
    GL.C5H14NO  :100,
    GL.C5H12N     :2}}#,

  # sn3 = headgroup
  def __init__(self, sn1):
    headgroup = GL.sn(mass=183.066044, type='Headgroup', chnops={'C':5, 'H':14, 'N':1, 'O':4, 'P':1},
    smiles='C[N+](C)(C)CCOP([O-])(=O)')
    # headgroup mass has -H to maintain neutral charge
    super().__init__(pLPC.adducts, sn1=sn1, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class LPE(GL.Glycerolipid):

  tooltip = 'Lyso-phosphatidylethanolamine'
  tailOrganisation = ['A']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Matches lipidblast
    GL.MA            :15, # Molecular ion
    GL.MH_PO4_s_HG_H2O:1, # Present in lipidblast
    GL.MH_s_FAk       :5, # Present in lipidblast
    GL.MH_s_FA       :15, # Present in lipidblast
    GL.FAH          :100, # Present in lipidblast
    GL.FAH_PUFA      :20, # Present in lipidex
    GL.C5H11NO5P      :5, # Present in lipidblast
    GL.C3H6O5P       :10,
    GL.C2H7NO4P       :3,
    GL.H2O4P          :5,
    GL.O3P            :5},
  
  "[M+H]+":{ # Matches lipidblast
    GL.MA              :50,  # Molecular ion
    GL.MA_s_H2O        :25,  # Present in lipidblast
    GL.MA_s_HG_H2O     :100, # Present in lipidblast
    GL.MH_PO4_s_HG_H2O :10,  # Present in lipidblast
    GL.MH_s_FAk        :1,   # Present in lipidblast
    GL.MH_s_FA         :1,   # Present in lipidblast
    GL.FAkH            :0,
    GL.C2H8NO         :10}
    }

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    headgroup = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1},
    smiles='NCCOP(O)(=O)')
    super().__init__(LPE.adducts, sn1=sn1, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class oLPE(GL.Glycerolipid):

  tooltip = 'Ether/Plasmanyl Lyso-phosphatidylethanolamine'
  tailOrganisation = ['O']
  givenName = 'Ether-LPE'

  adducts = {  # adduct:{spectra}

  "[M+H]+":{ # Matches lipidblast
    GL.MA          :100, # Molecular ion
    GL.MA_s_H2O     :20,
    GL.MA_PO4_s_HG  :10,
    GL.C2H8NO       :10}
    }

  # sn3 = headgroup
  def __init__(self, sn1):
    headgroup = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1},
    smiles='NCCOP(O)(=O)')
    super().__init__(oLPE.adducts, sn1=sn1, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class pLPE(GL.Glycerolipid):

  tooltip = 'Plasmenyl Lyso-phosphatidylethanolamine'
  tailOrganisation = ['P']
  givenName = 'Plasmenyl-LPE'

  adducts = {  # adduct:{spectra}

  "[M+H]+":{ # Matches lipidblast
    GL.MA               :100, # Molecular ion
    GL.MA_s_H2O          :20,
    GL.MA_PO4_s_HG        :5,
    GL.HGA_FA_s_H2O       :1,
    GL.HGA_FA_s_PO4      :20,
    GL.HGA_FA_s_H2O_PO4  :40}}

  # sn3 = headgroup
  def __init__(self, sn1):
    headgroup = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1},
    smiles='NCCOP(O)(=O)')
    super().__init__(pLPE.adducts, sn1=sn1, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class LNAPE(GL.Glycerolipid):

  tooltip = 'Lyso-N-Acyl-phosphatidylethanolamine'
  tailOrganisation = ['A','A']
  givenName = 'LPE-N(FA)'

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    GL.MA            :15,
    GL.MH_PO4_s_HG_H2O:1,
    GL.MH_s_FAk       :5,
    GL.MH_s_FA       :15,
    GL.FAH_b        :100,
    GL.HGA            :0,
    GL.C5H11NO5P      :0,
    GL.C3H6O5P       :10,
    GL.C2H7NO4P       :0,
    GL.H2O4P          :5,
    GL.O3P            :5}
    }

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1},
    smiles=sn1.inverseSmiles+'NCCOP(O)(=O)', hgtails=[sn1])
    super().__init__(LNAPE.adducts, sn1=sn2, sn3=headgroup)
    self.name=f"LPE-N({sn1.name}) {sn2.name}"

# ~ # ~ # ~ # ~ # ~ # ~ #

class GNAPE(GL.Glycerolipid):

  tooltip = 'Glycero-N-Acyl-phosphatidylethanolamine'
  tailOrganisation = ['A']
  givenName = 'DLPE-N(FA)'

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    GL.MA             :1,
    GL.MH_PO4_s_HG_H2O:1,
    GL.MH_s_FAk       :5,
    GL.MH_s_FA       :15,
    GL.HGA            :0,
    GL.C5H11NO5P      :1,
    GL.C3H8O6P       :75,
    GL.C3H6O5P      :100,
    GL.C2H7NO4P       :0,
    GL.H2O4P          :5,
    GL.O3P           :50}
    }

  # Should have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    headgroup = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1},
    smiles=sn1.inverseSmiles+'NCCOP(O)(=O)', hgtails=[sn1])
    super().__init__(GNAPE.adducts, sn3=headgroup)
    self.name=f"LPE-N({sn1.name})"

# ~ # ~ # ~ # ~ # ~ # ~ #

class LPG(GL.Glycerolipid):

  tooltip = 'Lyso-phosphatidylglycerol'
  tailOrganisation = ['A']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    GL.MH           :100,
    GL.MH_PO4_s_HG    :1,
    GL.MH_PO4_s_HG_H2O:1,
    GL.FAH          :100,
    GL.MH_s_FAk       :0,
    GL.C6H12O7P       :0,
    GL.C3H6O5P       :50,
    GL.H2O4P          :1,
    GL.O3P            :1},

  "[M+H]+":{
    GL.MA          :10,
    GL.MA_s_HG_H2O:100,
    GL.MA_s_HG_FA  :10,
    GL.FAkH        :10}
    }

  # sn3 = headgroup
  def __init__(self, sn1):
    headgroup = GL.sn(mass=172.013674, type='Headgroup', chnops={'C':3, 'H':9, 'O':6, 'P':1},
    smiles='OCC(O)COP(O)(=O)')
    super().__init__(LPG.adducts, sn1=sn1, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class LPI(GL.Glycerolipid):

  tooltip = 'Lyso-phosphatidylinositol'
  tailOrganisation = ['A']

  adducts = {  # adduct:{spectra}
  #"[M+Na-2H]-":{
  #  GL.C9H16O10P:50, # Not confident what the mz 192.988  is,
  #  GL.C6H10O8P :10, # perhaps as a [Phosphoglycerol+Na-2H]-.
  #  GL.C6H8O7P  :10, # Appears in "[M+Na-2H]-" spectra.
  #  Gl.C3H7NaO6P:80, # Additional mz 355.04  not included.
  #  GL.C3H6O5P  :80},

  "[M-H]-":{ # Needs Validation
    GL.MA             :2,
    GL.MH_PO4_s_HG    :0,
    GL.MH_PO4_s_HG_H2O:1,
    GL.FAH          :100, 
    GL.C9H16O10P     :10, 
    GL.C6H10O8P      :45,  # mz 233.001 can be a major -
    GL.C6H8O7P       :10,  # or minor fragment depending -
    GL.C3H6O5P       :45,  # on if its a sn1 or sn2 lyso lipid.
    GL.H2O4P          :1,
    GL.O3P            :5}
    } 

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    headgroup = GL.sn(mass=260.029718, type='Headgroup', chnops={'C':6, 'H':13,'O':9, 'P':1},
    smiles='OC1C(O)C(O)C(O)C(O)C1OP(O)(=O)')
    super().__init__(LPI.adducts, sn1=sn1, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class LPS(GL.Glycerolipid):

  tooltip = 'Lyso-phosphatidylserine'
  tailOrganisation = ['A']

  adducts = {  # adduct:{spectra}
  "[M+H]+":{
    GL.MA          :10,
    GL.MA_s_HG_H2O:100,
    GL.MA_s_HG_FA  :10,
    GL.FAkH        :10},

    "[M-H]-":{
    GL.MA           :0,
    GL.MA_PO4_s_HG  :1,
    GL.MA_s_HG_H2O  :0,
    GL.MA_s_FAk     :0,
    GL.MA_s_FA      :0,
    GL.FAH          :1,
    GL.C3H8O6P      :5,
    GL.C3H6O5P    :100,
    GL.H2O4P        :1,
    GL.O3P          :5},

    "[M+Na-2H]-":{
    GL.MA           :1,
    GL.MA_PO4_s_HG  :1,
    GL.MA_s_HG_H2O  :0,
    GL.MA_s_FAk    :10,
    GL.MA_s_FA     :10,
    GL.FAH          :5,
    GL.C3H7NaO6P    :1,
    GL.C3H8O6P      :0,
    GL.C3H6O5P    :100,
    GL.HGH_s_PO3    :1,
    GL.H2O4P        :0,
    GL.O3P          :5}

    }#,




  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    headgroup = GL.sn(mass=185.008923, type='Headgroup', chnops={'C':3, 'H':8, 'N':1, 'O':6, 'P':1},
    smiles='[O-]C(=O)C([NH3+])COP(O)(=O)')
    super().__init__(LPS.adducts, sn1=sn1, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class LNAPS(GL.Glycerolipid):

  tooltip = 'Lyso-N-Acyl-phosphatidylserine'
  tailOrganisation = ['A','A']
  givenName = 'LPS-N(FA)'

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    GL.MA          :10,
    GL.MA_PO4_s_HG:100,
    GL.MA_s_HG_H2O  :0,
    GL.MA_s_FAk     :0,
    GL.MA_s_FA      :0,
    GL.C3H6O5P     :50,}
    }#,

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=185.008923, type='Headgroup', chnops={'C':3, 'H':8, 'N':1, 'O':6, 'P':1},
    smiles='OC(=O)C(N'+sn1.smiles+')COP(O)(=O)', hgtails=[sn1])
    super().__init__(LNAPS.adducts, sn1=sn2, sn3=headgroup)
    self.name=f"LPS-N({sn1.name}) {sn2.name}"

# ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class Acyldeoxysphinganine(GL.Sphingolipid):

  tooltip = 'N-Acyl-deoxysphinganine'
  base_types = ['Dihydrodeoxysphinganine', 'Deoxysphinganine'] # 18:0;O, 18:1;O
  tailOrganisation = ['B','A']
  givenName = 'Cer[NMS]' # Cannot find proper name. Perhaps as NDS is Nonhydo FA, Dihydroceramide,  NMS could mean nonhydro FA, Monohydroceramide?

  adducts = {
  "[M+H]+":{
    GL.MA         :100,
    GL.MA_s_H2O    :40,
    GL.Cer_B       :15,
    GL.MA_s_FA    :100,
    GL.MA_s_FA_H2O  :0,
    GL.Cer_U       :30}
    }

  #"[M+H-H2O]+":{ #
  #  GL.MA         :100
  #  }}#,

  def __init__(self, base, sn1):
    super().__init__(Acyldeoxysphinganine.adducts, base, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class Acylsphinganine(GL.Sphingolipid):

  tooltip = 'N-Acyl-sphinganine'
  base_types = ['Sphinganine'] # 18:0;O2
  tailOrganisation = ['B','A']
  givenName = 'Cer[NDS]'

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # https://doi.org/10.1016/j.biochi.2016.07.012
    GL.MA           :30,
    GL.MA_s_H2O      :2,
    GL.MA_s_MeOH    :10,
    GL.MA_s_CH2O_H2O:10,
    GL.Cer_B        :10,
    GL.Cer_C         :2,
    GL.Cer_P        :10,
    GL.Cer_R        :25,
    GL.Cer_S        :20,
    GL.Cer_T       :100,
    GL.Cer_U         :5,
    GL.FAkH         :10},
    
  "[M+H]+":{
    GL.MA          :75,
    GL.MA_s_H2O   :100,
    GL.MA_s_2H2O    :5,
    GL.Cer_B       :10,
    GL.MA_s_FA     :30,
    GL.MA_s_FA_H2O :20,
    GL.Cer_U       :20,
    GL.Cer_D       :1}
    }

  #"[M+H-H2O]+":{ #
  #  GL.MA         :100
  #  }}#,

  def __init__(self, base, sn1):
    super().__init__(Acylsphinganine.adducts, base, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class Acylsphingosine(GL.Sphingolipid):

  tooltip = 'N-Acyl-sphingosine'
  base_types = ['Sphingosine', 'Sphingadiene'] # 18:1;O2, 18:2;O2
  tailOrganisation = ['B','A']
  givenName = 'Cer[NS]'

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ 
    # https://doi.org/10.1006/abio.2001.5536, 
    # https://doi.org/10.1016/S1044-0305(02)00358-6, 
    # https://doi.org/10.1016/j.biochi.2016.07.012, 
    # https://doi.org/10.1002/rcm.878 
    # https://www.lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=562&LM_ID=LMSP02010003&TRACK_ID=94
    GL.MA           :30,
    GL.MA_s_H2O      :2,
    GL.MA_s_CH2O    :20, # This fragment not present with sphinganine
    GL.MA_s_MeOH    :10,
    GL.MA_s_CH2O_H2O:10,
    GL.Cer_P        :10,
    GL.Cer_R        :25,
    GL.Cer_S        :20,
    GL.Cer_T       :100,
    GL.Cer_U         :5,
    GL.FAH          :25,
    GL.FAkH         :10},

  "[M+Hac-H]-":{ 
    GL.MA           :10,
    GL.MH           :30,
    GL.MH_s_H2O      :2,
    GL.MH_s_CH2O    :20,
    GL.MH_s_MeOH    :10,
    GL.MH_s_CH2O_H2O:10,
    GL.Cer_P        :10,
    GL.Cer_R        :25,
    GL.Cer_S        :20,
    GL.Cer_T       :100,
    GL.Cer_U         :5,
    GL.FAH          :25,
    GL.FAkH         :10,
    GL.C2O2H3       :50},

  "[M+H]+":{ 
    # https://doi.org/10.1002/bmc.4790, 
    # https://www.lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=482&LM_ID=LMSP02010002&TRACK_ID=326
    GL.MA            :5,
    GL.MA_s_H2O     :60,
    GL.MA_s_2H2O     :5,
    GL.MA_s_CH2O     :0,
    GL.MA_s_CH2O_H2O :2,
    GL.MA_s_FA      :15,
    GL.MA_s_FA_H2O: 100,
    GL.Cer_U        :10,
    GL.Cer_D        :5},

  "[M+H-H2O]+":{ #
    GL.MA           :60,
    GL.MA_s_H2O     :10,
    GL.MA_s_CH2O     :1,
    GL.Cer_R       :100,
    GL.Cer_Bb       :15,
    GL.Cer_U         :2,
    GL.Cer_D         :5}}

  def __init__(self, base, sn1):
    super().__init__(Acylsphingosine.adducts, base, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class Acylphytosphingosine(GL.Sphingolipid):

  tooltip = 'N-Acyl-phytosphingosine'
  base_types = ['Phytosphingosine'] # 18:0;O3
  tailOrganisation = ['B','A']
  givenName = 'Cer[NP]'

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # https://doi.org/10.1002/rcm.878
    GL.MA       :10,
    GL.MA_s_H2O :40,
    GL.MA_s_MeOH:50,
    GL.MA_s_2H2O:50,    
    GL.Cer_P     :5,
    GL.Cer_Q     :5,    
    GL.Cer_Rb   :50,
    GL.Cer_S     :5,
    GL.Cer_T     :1,
    GL.Cer_U     :5,
    GL.Cer_W    :75,
    GL.Cer_X   :100,
    GL.FAH      :20,
    GL.FAkH      :1},
    
  "[M+H]+":{
    GL.MA        :100,
    GL.MA_s_H2O   :50,
    GL.MA_s_2H2O   :5,
    GL.Cer_B       :5,
    GL.Cer_Bb     :15,
    GL.Cer_R      :20,
    GL.Cer_R_s_H2O:15
    }}#,

  #"[M+H-H2O]+":{ #
  #  GL.MA         :100
  #  }}#,

  def __init__(self, base, sn1):
    super().__init__(Acylphytosphingosine.adducts, base, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class CerP(GL.Sphingolipid):

  tooltip = 'N-Acyl-ceramide-1-phosphate'
  base_types = ['Sphinganine', 'Sphingosine', 'Phytosphingosine']  # 18:0;O2, 18:1;O2
  tailOrganisation = ['B','A']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Matches lipidblast, 
    # https://doi.org/10.1016/j.biochi.2016.07.012, 
    # https://www.lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=560&LM_ID=LMSP02050001&TRACK_ID=335
    GL.MA         :10, # Molecular ion
    GL.MA_s_H2O   :10, # Present in lipidblast
    GL.MA_s_FAk   :20, # Present in lipidblast
    GL.MA_s_FA    :20, # Present in lipidblast
    GL.H2O4P     :100, # Present in lipidblast
    GL.H2O2P       :0,
    GL.O3P       :100},# Present in lipidblast
    
  "[M+H]+":{
    # https://www.lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=562&LM_ID=LMSP02050001&TRACK_ID=331
    GL.MA             :10,
    GL.MA_s_H2O       :10,
    GL.MA_s_HG        :10,
    GL.MA_s_HG_H2O    :10,
    GL.MA_s_HG_2H2O   :10,
    GL.MA_s_HG_FA_H2O:100}
    }

  #"[M+H-H2O]+":{ #
  #  GL.MA         :100
  #  }}#,

  def __init__(self, base, sn1):
    headgroup = GL.sn(mass=97.976895, type='Headgroup', chnops={'H':3, 'O':4, 'P':1},
    smiles='O=P(O)(O)')
    super().__init__(CerP.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #    
'''
class SpPA(GL.Sphingolipid): # While this class works, disabled as ceramide varability isn't yet built in.

  tooltip = 'Sphingosine-1-phosphate'
  base_types = ['Sphingosine']  # 18:0;O2, 18:1;O2
  No_Tails = 0
  tailOrganisation = ['B']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    GL.MA         :10,
    GL.H2O4P     :100,
    GL.O3P       :100},
    
  "[M+H]+":{
    GL.MA           :10,
    GL.MA_s_HG_2H2O:100}
    }

  #"[M+H-H2O]+":{ #
  #  GL.MA         :100
  #  }}#,

  def __init__(self, base):
    headgroup = GL.sn(mass=97.976895, type='Headgroup', chnops={'H':3, 'O':4, 'P':1},
    smiles='O=P(O)(O)')
    super().__init__(SpPA.adducts, base, headgroup)
'''
# ~ # ~ # ~ # ~ # ~ # ~ #

class CerPC(GL.Sphingolipid):

  tooltip = 'N-Acyl-ceramide-1-phosphocholine'
  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  tailOrganisation = ['B','A']
  givenName = 'PC-Cer'

  adducts = {  # adduct:{spectra}
  # https://www.lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=647&LM_ID=LMSP03010002&TRACK_ID=320
  # https://www.lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=731&LM_ID=LMSP03010001&TRACK_ID=82
  "[M+H]+":{ # Looks Good, Matches lipidblast, 
    GL.MA            :10, # Molecular ion
    GL.MA_s_H2O       :1, # Present in lipidblast
    GL.MA_s_TMA       :1, # Present in lipidblast
    GL.MA_s_TMA_H2O   :1, # Present in lipidblast
    GL.MA_s_HG_H2O    :1, # Present in lipidblast
    GL.MA_s_HG_FA_H2O :5, 
    GL.C5H15NO4P    :100,
    GL.C5H12N       :5},

  "[M+Li]+":{
    GL.MA           :100,
    GL.MA_s_H2O       :0,
    GL.MA_s_TMA      :60,
    GL.MA_s_HG_H2O  :100,
    GL.MH_s_HG_H2O    :3,
    GL.MH_s_HG_2H2O   :5,
    GL.FA_C2H3N       :3,
    GL.Cer_R          :1},

    "[M+Hac-H]-":{ # http://prime.psc.riken.jp/compms/static/images/figure/lipid/SM_[M+CH3COO]-.png
    GL.MA           :10,
    GL.M_s_CH3     :100,
    GL.M_s_CH3_FAk   :2,    
    GL.C4H11NO4P    :25,
    GL.C2H4O4P       :2,
    GL.H2O4P         :5,
    GL.O3P           :5,
    GL.C2O2H3        :1},

    "[M-CH3]-":{ # In source fragment
    GL.MA          :100,
    GL.M_s_CH3_FAk   :0,    
    GL.C4H11NO4P    :40,
    GL.C2H4O4P       :0,
    GL.H2O4P         :0,
    GL.O3P          :80}}

  def __init__(self, base, sn1):
    headgroup = GL.sn(mass=183.066044, type='Headgroup', chnops={'C':5, 'H':14, 'N':1, 'O':4, 'P':1},
    smiles='C[N+](C)(C)CCOP([O-])(=O)')
    # headgroup mass has -H to maintain neutral charge
    super().__init__(CerPC.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class CerPE(GL.Sphingolipid):

  tooltip = 'N-Acyl-ceramide-1-phosphoethanolamine'
  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  tailOrganisation = ['B','A']
  givenName = 'PE-Cer'

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ 
    # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
    GL.MA        :100,
    GL.MA_s_H2O    :1,
    GL.MA_s_FAk    :5,
    GL.MA_s_FA     :1,
    GL.C2H7NO4P   :25,
    GL.H2O4P      :25,
    GL.H2O2P       :0,
    GL.O3P        :25},

  "[M+H]+":{
    GL.MA            :50, # Molecular ion
    GL.MA_s_H2O       :0,
    GL.MA_s_HG_H2O  :100,
    GL.MA_s_HG_2H2O  :50,
    GL.MA_s_HG_FA_H2O:40,
    GL.FA_C2H3N      :10}
    }

  def __init__(self, base, sn1):
    headgroup = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1},
    smiles='NCCOP(O)(=O)')
    super().__init__(CerPE.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class CerPI(GL.Sphingolipid):

  tooltip = 'N-Acyl-ceramide-1-phosphoinositol'
  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  tailOrganisation = ['B','A']
  givenName = 'PI-Cer'

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ 
    # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
    GL.MA        :100,
    GL.MA_s_H2O    :1,
    GL.MA_PO4_s_HG :2,
    GL.MA_s_Gal    :0,
    GL.MH_PO4_s_HG_FA:3,
    GL.MA_s_FAk    :5,
    GL.MA_s_FA     :1,
    GL.C6H12O9P    :0,
    GL.C6H10O8P   :40,
    GL.C6H8O7P    :15, 
    GL.H2O4P      :25,
    GL.O3P        :25}
    }

  def __init__(self, base, sn1):
    headgroup = GL.sn(mass=260.029718, type='Headgroup', chnops={'C':6, 'H':13, 'O':9, 'P':1},
    smiles='OC1C(O)C(O)C(O)C(O)C1OP(O)(=O)')
    super().__init__(CerPI.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class ACer(GL.Sphingolipid):

  tooltip = 'N-Acyl-ceramide'
  base_types = ['Sphinganine']  # 18:0;O2
  tailOrganisation = ['B','AA']
  specificTailOrganisation = ['B','A','A']
  givenName = 'A(FA)-Cer'

  adducts = {  # adduct:{spectra}
  "[M+H]+":{
    GL.MA             :25,
    GL.MA_s_H2O        :5,
    GL.MA_s_HG        :25,
    GL.MA_s_HG_H2O   :100,
    GL.MA_s_HG_2H2O   :10,
    GL.MA_s_2FA       :10,
    GL.MA_s_2FAk       :5,
    GL.MA_s_HG_FA     :15,
    GL.MA_s_HG_FAk    :10,
    GL.MA_s_HG_FA_H2O :10,
    GL.Cer_U           :0,
    GL.Cer_D           :0,
    GL.Cer_P           :0}
    }

  def __init__(self, base, sn1, sn2):
    headgroup = GL.sn(mass=GL.masses['H2O'], type='Headgroup', chnops={'H':2, 'O':1},
    smiles=sn1.inverseSmiles, hgtails=[sn1])
    super().__init__(ACer.adducts, base, sn1=sn2, headgroup=headgroup)
    self.name=f"A(FA)-Cer{base.lipidSuffix if base.lipidSuffix not in self.lipid_class else ''} {base.name}/{sn1.name}_{sn2.name}"
    self.specificname=f"A({sn1.name})-Cer{base.lipidSuffix if base.lipidSuffix not in self.lipid_class else ''} {base.name}/{sn2.name}"

# ~ # ~ # ~ # ~ # ~ # ~ #

class HexCer(GL.Sphingolipid):

  tooltip = 'N-Acyl-ceramide-1-hexose'
  base_types = ['Sphinganine', 'Sphingosine', 'Sphingadiene']  # 18:0;O2, 18:1;O2, 18:2;O2
  tailOrganisation = ['B','A']
  givenName = 'Hex-Cer'

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ 
    # https://www.lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=642&LM_ID=LMSP0501AA01&TRACK_ID=330
    GL.MA           :30,
    GL.MA_s_H2O      :0,
    GL.MA_s_2H2O     :0,
    GL.MA_s_HG     :100,
    GL.MA_s_HG_H2O   :5,
    GL.Cer_P         :1,
    GL.Cer_R         :1,
    GL.Cer_S         :1,
    GL.Cer_T         :1,
    GL.Cer_U         :1,
    GL.FAH           :1,
    GL.FAkH          :1},

  "[M+H]+":{ 
    # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html 
    # https://www.lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=644&LM_ID=LMSP0501AA01&TRACK_ID=323
    GL.MA              :5,
    GL.MA_s_H2O       :25,
    GL.MA_s_HG         :5,
    GL.MA_s_HG_H2O    :25,
    GL.MA_s_HG_2H2O   :10,
    GL.MA_s_HG_FA     :25,
    GL.MA_s_HG_FAk     :0,
    GL.MA_s_HG_FA_H2O:100,
    GL.Cer_U           :5,
    GL.Cer_D          :10}
    }

  def __init__(self, base, sn1):
    headgroup = GL.sn(mass=180.063388, type='Headgroup', chnops={'C':6, 'H':12, 'O':6},
    smiles='OC1C(O)C(O)C(CO)OC1')
    super().__init__(HexCer.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class AHexCer(GL.Sphingolipid):

  tooltip = 'N-Acyl-ceramide-1-Acylhexose'
  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  tailOrganisation = ['B','AA']
  specificTailOrganisation = ['B','A','A']
  givenName = 'Hex(FA)-Cer'

  adducts = {  # adduct:{spectra}
  "[M+H]+":{ # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
    GL.MA              :5,
    GL.MA_s_H2O       :25,
    GL.MA_s_HG         :5,
    GL.MA_s_HG_H2O    :25,
    GL.MA_s_HG_2H2O   :10,
    GL.MA_s_HG_FA     :25,
    GL.MA_s_HG_FAk     :0,
    GL.MA_s_HG_FA_H2O:100,
    GL.Cer_U           :1,
    GL.Cer_D          :10,
    GL.Cer_P          :10}
    }

  def __init__(self, base, sn1, sn2):
    headgroup = GL.sn(mass=180.063388, type='Headgroup', chnops={'C':6, 'H':12, 'O':6},
    smiles='OC1C(O)C(O)C(CO'+sn1.smiles+')OC1', hgtails=[sn1])
    super().__init__(AHexCer.adducts, base, sn1=sn2, headgroup=headgroup)
    self.name=f"Hex(FA)-Cer{base.lipidSuffix if base.lipidSuffix not in self.lipid_class else ''} {base.name}/{sn1.name}_{sn2.name}"
    self.specificname=f"Hex({sn1.name})-Cer{base.lipidSuffix if base.lipidSuffix not in self.lipid_class else ''} {base.name}/{sn2.name}"

# ~ # ~ # ~ # ~ # ~ # ~ #

class Hex2Cer(GL.Sphingolipid):

  tooltip = 'N-Acyl-ceramide-1-dihexose'
  base_types = ['Sphinganine', 'Sphingosine', 'Sphingadiene']  # 18:0;O2, 18:1;O2, 18:2;O2
  tailOrganisation = ['B','A']
  givenName = 'Hex2-Cer'

  adducts = {  # adduct:{spectra}
  "[M+H]+":{ 
    # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html 
    # https://www.lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=806&LM_ID=LMSP0501AB02&TRACK_ID=322
    GL.MA              :5,
    GL.MA_s_H2O       :25,
    GL.MA_s_Gal        :5,  
    GL.MA_H2O_s_Gal   :25,
    GL.MA_s_HG         :5,
    GL.MA_s_HG_H2O    :25,
    GL.MA_s_HG_2H2O   :10,
    GL.MA_s_HG_FA     :25,
    GL.MA_s_HG_FAk     :0,
    GL.MA_s_HG_FA_H2O:100,
    GL.Cer_U           :5,
    GL.Cer_D          :10}
    }

  def __init__(self, base, sn1):
    headgroup = GL.sn(mass=342.116212, type='Headgroup', chnops={'C':12, 'H':22, 'O':11},
    smiles='OC1C(O)C(O(C2OC(CO)C(O)C(O)C2(O)))C(CO)OC1')
    super().__init__(Hex2Cer.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class Hex3Cer(GL.Sphingolipid):

  tooltip = 'N-Acyl-ceramide-1-trihexose'
  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  tailOrganisation = ['B','A']
  givenName = 'Hex3-Cer'

  adducts = {  # adduct:{spectra}
  "[M+H]+":{ # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
    GL.MA              :5,
    GL.MA_s_H2O       :25,
    GL.MA_H2O_s_Gal    :5,
    GL.MA_s_Gal       :25,
    GL.MA_s_Gal_H2O    :5,
    GL.MA_2H2O_s_2Gal  :5,  
    GL.MA_H2O_s_2Gal  :25,
    GL.MA_s_2Gal       :5,
    GL.MA_s_HG         :5,
    GL.MA_s_HG_H2O   :100,
    GL.MA_s_HG_2H2O   :10,
    GL.MA_s_HG_FA      :5,
    GL.MA_s_HG_FAk     :0,
    GL.MA_s_HG_FA_H2O :25,
    GL.Cer_U           :5,
    GL.Cer_D          :10}
    }

  def __init__(self, base, sn1):
    headgroup = GL.sn(mass=504.169035, type='Headgroup', chnops={'C':18, 'H':32, 'O':16},
    smiles='OC1C(O)C(O(C2OC(CO)C(O(C3OC(CO)C(O)C(O)C3(O)))C(O)C2(O)))C(CO)OC1')
    super().__init__(Hex3Cer.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class SHexCer(GL.Sphingolipid):

  tooltip = 'Sulfatide'
  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  tailOrganisation = ['B','A']
  givenName = 'SHex-Cer'

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ #
    GL.MA      :10, # Molecular ion
    GL.MA_s_H2O:10,
    GL.Cer_T    :5,
    GL.MA_s_HG  :5,
    GL.MA_s_FAk:20,
    GL.MA_s_FA :20,
    GL.C6H9O8S:100,
    GL.O4SH   :100}
    }

  def __init__(self, base, sn1):
    headgroup = GL.sn(mass=260.020201964, type='Headgroup', chnops={'C':6, 'H':12, 'O':9, 'S':1},
    smiles='OC1C(OS(=O)(=O)O)C(O)C(CO)OC1')
    super().__init__(SHexCer.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class GM1a(GL.Sphingolipid):

  tooltip = 'GM1a'
  base_types = ['Sphinganine']
  tailOrganisation = ['B','A']
  givenName = 'GM1a'
  gm = GL.gm1a

  adducts = {  # adduct:{spectra}
  "[M+H]+":{
    GL.MA       :20,
    GL.MA_s_H2O :5,

    GL.GMAY0     :3,
    GL.GMAY2     :50,
    GL.GMAY4     :30,

    GL.GMAZ4     :5,
    GL.GMAY5     :30,
    GL.GMAZ5     :5,

    GL.GMAB2     :100,
    GL.GMAB3     :12,
    GL.GMAB4     :1,
    GL.GMAB5     :2,

    GL.MA_s_HG     :5,
    GL.MA_s_HG_H2O :20,
    GL.HGA_s_H2O   :25,
  }}

  def __init__(self, base, sn1):
    headgroup = GL.gm1a[0]
    super().__init__(GM1a.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ # 

class GM2(GL.Sphingolipid):

  tooltip = 'GM2'
  base_types = ['Sphinganine']
  tailOrganisation = ['B','A']
  givenName = 'GM2'
  gm = GL.gm2

  adducts = {  # adduct:{spectra}
  "[M+H]+":{
    GL.MA       :50,
    GL.MA_s_H2O :50,

    GL.GMAY0     :30,
    GL.GMAZ0     :5,
    GL.GMAY1     :5,
    GL.GMAZ1     :1,
    GL.GMAY2     :1,
    GL.GMAZ2     :1,
    GL.GMAY3     :20,
    GL.GMAZ3     :5,

    GL.GMAB1     :2,
    GL.GMAB3     :15,
    GL.GMAC3     :1,   

    GL.MA_s_HG     :5,
    GL.MA_s_HG_H2O :10,
    GL.HGA_s_H2O   :100,
  }}

  def __init__(self, base, sn1):
    headgroup = GL.gm2[0]
    super().__init__(GM2.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ # 

class GM3(GL.Sphingolipid):

  tooltip = 'GM3'
  base_types = ['Sphinganine']
  tailOrganisation = ['B','A']
  givenName = 'GM3'
  gm = GL.gm3

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    GL.MA:100,
    GL.GMAY0:1,
    GL.GMAB0:25},

  "[M+H]+":{
    GL.MA       :5,
    GL.MA_s_H2O :75,
    GL.GMAY0     :30,
    GL.GMAZ0     :25,
    GL.GMAY1     :85,
    GL.GMAZ1     :45,
    GL.GMAB0     :10,
    GL.GMAB1     :10,  
    GL.MA_s_HG     :70,
    GL.MA_s_HG_H2O :100,
    GL.MA_s_HG_2H2O:10,
    GL.HGA_s_H2O   :3}}#,

  #"[M+NH4]+":{
  #  GL.MA:50,
  #  GL.MA_s_H2O:75,
  #  GL.GMY0:1,
  #  GL.MH_s_HG_H2O:100}}

  def __init__(self, base, sn1):
    headgroup = GL.gm3[0]
    super().__init__(GM3.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ # 

class GD1a(GL.Sphingolipid):

  tooltip = 'GD1a'
  base_types = ['Sphinganine']
  tailOrganisation = ['B','A']
  givenName = 'GD1a'
  gm = GL.gd1a

  adducts = {  # adduct:{spectra}
  "[M+H]+":{
    GL.MA          :40,
    GL.MA_s_H2O    :10,
    GL.MA_s_2H2O   :1,

    GL.GMAX2        :100,
    GL.GMAX4        :1,
    GL.GMAX5        :10,

    GL.GMAY1        :1,
    GL.GMAY3        :2,
    GL.GMAY6        :2,

    GL.GMAB0        :2,  
    GL.GMAB1        :2,  
    GL.GMAB3        :10,

    GL.MA_s_HG     :5,

    GL.HGA_s_H2O   :10,
  }}

  def __init__(self, base, sn1):
    headgroup = GL.gd1a[0]
    super().__init__(GD1a.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ # 

class GD1b(GL.Sphingolipid):

  tooltip = 'GD1b'
  base_types = ['Sphinganine']
  tailOrganisation = ['B','A']
  givenName = 'GD1b'
  gm = GL.gd1b

  adducts = {  # adduct:{spectra}
  "[M+H]+":{
    GL.MA          :100,
    GL.MA_s_H2O    :50,
    GL.MA_s_2H2O   :1,

    GL.GMAX0        :2,
    GL.GMAX1        :2,
    GL.GMAX3        :5,
    GL.GMAX5        :5,

    GL.GMAY2        :5,
    GL.GMAY5        :1,
    GL.GMAY7        :10,
  
    GL.GMAZ2        :2,
    GL.GMAZ7        :2,

    GL.GMAB1        :2,
    GL.GMAB2        :40,
    GL.GMAB3        :2,
    GL.GMAB4        :2,
    GL.GMAB5        :20,  

    GL.MA_s_HG     :5,

    GL.HGA_s_H2O   :10,
    GL.HGA_s_2H2O  :5,
  }}

  def __init__(self, base, sn1):
    headgroup = GL.gd1b[0]
    super().__init__(GD1b.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ # 

class GD2(GL.Sphingolipid):

  tooltip = 'GD2'
  base_types = ['Sphinganine']
  tailOrganisation = ['B','A']
  givenName = 'GD2'
  gm = GL.gd2

  adducts = {  # adduct:{spectra}
  "[M+H]+":{
    GL.MA          :100,
    GL.MA_s_H2O    :50,
    GL.MA_s_2H2O   :5,

    GL.GMAX0        :10,
    GL.GMAX1        :10,
    GL.GMAX2        :1,
    GL.GMAX4        :2,

    GL.GMAY0        :10,
    GL.GMAY2        :1,
    GL.GMAY3        :10,
    GL.GMAY5        :5,

    GL.GMAZ0        :5,
    GL.GMAZ3        :2,
    GL.GMAZ5        :1,

    GL.GMAA1        :1,   

    GL.GMAB1        :2,  
    GL.GMAB3        :10,

    GL.MA_s_HG     :2,
    GL.MA_s_HG_H2O :5,

    GL.HGA_s_H2O   :30,
  }}

  def __init__(self, base, sn1):
    headgroup = GL.gd2[0]
    super().__init__(GD2.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ # 

class GD3(GL.Sphingolipid):

  tooltip = 'GD3'
  base_types = ['Sphinganine']
  tailOrganisation = ['B','A']
  givenName = 'GD3'
  gm = GL.gd3

  adducts = {  # adduct:{spectra}
  "[M+H]+":{
    GL.MA       :40,
    GL.MA_s_H2O :40,

    GL.GMAX1     :50,

    GL.GMAY0     :20,
    GL.GMAZ0     :10,
    GL.GMAY1     :75,
    GL.GMAZ1     :20,
    GL.GMAY2     :50,
    GL.GMAZ2     :20,

    GL.GMAA0     :10,

    GL.GMAB0     :30,
    GL.GMAB1     :60,  
    GL.GMAB2     :100,

    GL.MA_s_HG     :30,
    GL.MA_s_HG_H2O :20,

    GL.HGA_s_H2O   :30,
  }}

  def __init__(self, base, sn1):
    headgroup = GL.gd3[0]
    super().__init__(GD3.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ # 

class GT1a(GL.Sphingolipid):

  tooltip = 'GT1a'
  base_types = ['Sphinganine']
  tailOrganisation = ['B','A']
  givenName = 'GT1a'
  gm = GL.gt1a

  adducts = {  # adduct:{spectra}
  "[M-2H]2-":{ #  https://doi.org/10.1002/jssc.202001248
    GL.GMHY0 :50,
    GL.GMHY1 :75,
    GL.GMHY3  :1,
    GL.GMAY0 :25,
    GL.GMHB1:100
  }}

  def __init__(self, base, sn1):
    headgroup = GL.gt1a[0]
    super().__init__(GT1a.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ # 

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class FA(GL.OtherLipid):

  tooltip = 'Free Fatty acid'
  tailOrganisation = ['A']

  adducts = {
    "[M-H]-":{
    GL.MH            :10,
    GL.MH_s_H2O     :100
    }}

  def __init__(self, sn1): # Mass and formula of body is water because body == tail
    body = GL.Other(name='FA ', mass=18.010564684, chnops={'H':2, 'O':1}, smiles='O'+sn1.smiles)
    super().__init__(FA.adducts, body, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #    
'''
class FAHFA(GL.OtherLipid):

  tooltip = 'Fatty Acid Hydroxyl Fatty Acid\nEnsure Hydroxy Fatty Acids are available'
  tailOrganisation = ['A','A']

  adducts = {
    # 10.1016/j.cell.2014.09.035
    "[M-H]-":{
    GL.MH             :10,
    GL.MH_s_H2O        :0,
    GL.MH_s_FAk       :10,
    GL.MH_s_FA        :10,
    GL.FAH           :100
    }}

  def __init__(self, sn1, sn2):
    smiles = sn1.smiles
    smiles = smiles.replace('(O)', '(O'+sn2.smiles+')', 1)
    body = GL.Other(name='FAHFA '+sn1.name, mass=sn1.mass, chnops=sn1.formula, smiles='O'+smiles)
    if sn1.oh > 0: super().__init__(FAHFA.adducts, body, sn1=sn2)
    else: pass # sn1 = GL.sn(c=sn1.c, d=sn1.d, type=sn1.type, me=sn1.me, oh=sn1.oh+1, dt=sn1.dt)
'''
# ~ # ~ # ~ # ~ # ~ # ~ #

class ZE(GL.OtherLipid):

  tooltip = 'Zymosteryl Ester'
  tailOrganisation = ['A']

  adducts = {  # adduct:{spectra}
  "[M+NH4]+":{ #
    GL.MA      :10, # Molecular ion
    GL.MA_s_FAk :0,
    GL.MA_s_FA  :0,
    GL.MH_s_FAk :0,
    GL.MH_s_FA:100,
    GL.FAkA     :0,
    GL.FAkH     :0,
    GL.C13H19  :20, # Present in lipidblast
    GL.C12H17  :20, # Present in lipidblast
    GL.C11H15  :20, # Present in lipidblast
    GL.C10H15  :20} # Present in lipidblast
    }

  def __init__(self, sn1):
    body = GL.Other(name='Zymosteryl', mass=384.339216, chnops={'C':27, 'H':44, 'O':1},
    smiles='CC(C)=CCCC(C)C1CCC2C3CC=C4CC(CCC4(C)C3CCC12C)O'+sn1.smiles)
    super().__init__(ZE.adducts, body, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class CE(GL.OtherLipid):

  tooltip = 'Cholesteryl Ester'
  tailOrganisation = ['A']

  adducts = {  # adduct:{spectra}
  "[M+NH4]+":{ #
    GL.MA      :10, # Molecular ion
    GL.MA_s_FAk :0,
    GL.MA_s_FA  :0,
    GL.MH_s_FAk :0,
    GL.MH_s_FA:100,
    GL.FAkA     :0,
    GL.FAkH     :0,
    GL.C13H19  :20, # Present in lipidblast
    GL.C12H17  :20, # Present in lipidblast
    GL.C11H15  :20, # Present in lipidblast
    GL.C10H15  :20} # Present in lipidblast
    }

  def __init__(self, sn1):
    body = GL.Other(name='Cholesteryl', mass=386.354866092, chnops={'C':27, 'H':46, 'O':1},
    smiles='CC(C)CCCC(C)C1CCC2C3CC=C4CC(CCC4(C)C3CCC12C)O'+sn1.smiles)
    super().__init__(CE.adducts, body, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class BRSE(GL.OtherLipid):

  tooltip = 'Brassicasterol Ester'
  tailOrganisation = ['A']

  adducts = {  # adduct:{spectra}
  "[M+NH4]+":{ #
    GL.MA      :10, # Molecular ion
    GL.MA_s_FAk :0,
    GL.MA_s_FA  :0,
    GL.MH_s_FAk :0,
    GL.MH_s_FA:100,
    GL.FAkA     :0,
    GL.FAkH     :0,
    GL.C13H19  :20,
    GL.C12H17  :20,
    GL.C11H15  :20,
    GL.C10H15  :20}
    }

  def __init__(self, sn1):
    body = GL.Other(name='Brassicasterol', mass=398.354866, chnops={'C':28, 'H':46, 'O':1},
    smiles='CC(C)C(C)C=CC(C)C1CCC2C3CC=C4CC(CCC4(C)C3CCC12C)O'+sn1.smiles)
    super().__init__(BRSE.adducts, body, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class CASE(GL.OtherLipid):

  tooltip = 'Camposterol Ester'
  tailOrganisation = ['A']

  adducts = {  # adduct:{spectra}
  "[M+NH4]+":{ #
    GL.MA      :10, # Molecular ion
    GL.MA_s_FAk :0,
    GL.MA_s_FA  :0,
    GL.MH_s_FAk :0,
    GL.MH_s_FA:100,
    GL.FAkA     :0,
    GL.FAkH     :0,
    GL.C13H19  :20,
    GL.C12H17  :20,
    GL.C11H15  :20,
    GL.C10H15  :20}
    }

  def __init__(self, sn1):
    body = GL.Other(name='Camposterol', mass=400.370516, chnops={'C':28, 'H':48, 'O':1},
    smiles='CC(C)C(C)CCC(C)C1CCC2C3CC=C4CC(CCC4(C)C3CCC12C)O'+sn1.smiles)
    super().__init__(CASE.adducts, body, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class SISE(GL.OtherLipid):

  tooltip = 'Sitosterol Ester'
  tailOrganisation = ['A']

  adducts = {  # adduct:{spectra}
  "[M+NH4]+":{ #
    GL.MA      :10, # Molecular ion
    GL.MA_s_FAk :0,
    GL.MA_s_FA  :0,
    GL.MH_s_FAk :0,
    GL.MH_s_FA:100,
    GL.FAkA     :0,
    GL.FAkH     :0,
    GL.C13H19  :20,
    GL.C12H17  :20,
    GL.C11H15  :20,
    GL.C10H15  :20}
    }

  def __init__(self, sn1):
    body = GL.Other(name='Sitosterol', mass=414.386166, chnops={'C':29, 'H':50, 'O':1},
    smiles='CC(C)C(CC)CCC(C)C1CCC2C3CC=C4CC(CCC4(C)C3CCC12C)O'+sn1.smiles)
    super().__init__(SISE.adducts, body, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class STSE(GL.OtherLipid):

  tooltip = 'Sigmasterol Ester'
  tailOrganisation = ['A']

  adducts = {  # adduct:{spectra}
  "[M+NH4]+":{ #
    GL.MA      :10, # Molecular ion
    GL.MA_s_FAk :0,
    GL.MA_s_FA  :0,
    GL.MH_s_FAk :0,
    GL.MH_s_FA:100,
    GL.FAkA     :0,
    GL.FAkH     :0,
    GL.C13H19  :20,
    GL.C12H17  :20,
    GL.C11H15  :20,
    GL.C10H15  :20}
    }

  def __init__(self, sn1):
    body = GL.Other(name='Sigmasterol', mass=412.370516, chnops={'C':29, 'H':48, 'O':1},
    smiles='CC(C)C(CC)C=CC(C)C1CCC2C3CC=C4CC(CCC4(C)C3CCC12C)O'+sn1.smiles)
    super().__init__(STSE.adducts, body, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class LE(GL.OtherLipid):

  tooltip = 'Lanosteryl Ester'
  tailOrganisation = ['A']

  adducts = {  # adduct:{spectra}
  "[M+NH4]+":{ #
    GL.MA      :10, # Molecular ion
    GL.MA_s_FAk :0,
    GL.MA_s_FA  :0,
    GL.MH_s_FAk :0,
    GL.MH_s_FA:100,
    GL.FAkA     :0,
    GL.FAkH     :0,
    GL.C13H19  :20,
    GL.C12H17  :20,
    GL.C11H15  :20,
    GL.C10H15  :20}
    }

  def __init__(self, sn1):
    body = GL.Other(name='Lanosteryl', mass=426.386166, chnops={'C':30, 'H':50, 'O':1},
    smiles='C/C(C)=C\CCC(C)C1CCC2(C)C=3CCC4C(C)(C)C(CCC4(C)C3CCC12C)O'+sn1.smiles)
    super().__init__(LE.adducts, body, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class AC(GL.OtherLipid):

  tooltip = 'Acyl-carnitine'
  tailOrganisation = ['A']

  adducts = {
    "[M+H]+":{ # Taken from LipidMatch
    GL.MA          :10, # Molecular ion
    GL.MA_s_TMA    :10,
    GL.MA_s_FA     :10,
    GL.MA_s_FA_TMA:100,
    GL.FAkH        :10,
    GL.C3H10N      :10
    }}

  def __init__(self, sn1):
    body = GL.Other(name='AC', mass=161.105193, chnops={'C':7, 'H':15, 'N':1, 'O':3},
    smiles='C[N+](C)(C)CC(O'+sn1.smiles+')CC(=O)[O-]')
    super().__init__(AC.adducts, body, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class AcylCoA(GL.OtherLipid):

  tooltip = 'Acyl-CoA'
  tailOrganisation = ['A']
  
  adducts = {
    "[M+H]+":{ # https://pubs.acs.org/doi/pdf/10.1021/acs.analchem.6b03623
    GL.MA             :10,
    GL.C11H19N2O2S_FA:100,
    GL.C10H16N5O10P2  :40
    }}

  def __init__(self, sn1):
    body = GL.Other(name='AcylCoA', mass=767.115206, chnops={'C':21, 'H':36, 'N':7, 'O':16, 'P':3, 'S':1},
    smiles='OP(=O)(O)OC1C(O)C(n2cnc3c(N)ncnc32)OC1COP(=O)(O)OP(=O)(O)OC(C)(C)CC(O)C(=O)NCCC(=O)NCCS'+sn1.smiles)
    super().__init__(AcylCoA.adducts, body, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class AsFA(GL.OtherLipid):

  tooltip = 'Arsenolipid Fatty acid'
  tailOrganisation = ['A']

  adducts = {
    "[M+H]+":{ # https://doi.org/10.1016/j.jpba.2023.115628, https://doi.org/10.1021/om4011092 
    GL.MH            :80,
    GL.MH_s_H2O      :50,
    GL.AsC2H8O      :100,
    GL.AsC2H6        :40,
    GL.AsC2H4         :1
    }}

  def __init__(self, sn1):
  # Mass and formula +H2O due to -H2O in Otherlipid class
    body = GL.Other(name='AsFA ', mass=137.966201, chnops={'C':2, 'H':7, 'O':2, 'As':1}, smiles='O'+sn1.smiles+'[As](=O)(C)C')
    super().__init__(AsFA.adducts, body, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class AsPL(GL.Glycerolipid):

  tooltip = 'Arsenosugar phospholipid'
  tailOrganisation = ['AA']

  adducts = {  
  "[M-H]-":{ # https://doi.org/10.1007/s11745-017-4266-x, https://doi.org/10.1016/j.jpba.2023.115628
    GL.MA            :50,
    GL.MH_s_FAk      :20,
    GL.MH_s_FA       :20,
    GL.MH_PO4_s_HG   :10,
    GL.MH_PG_s_HG    :20,
    GL.MH_PO4_s_HG_FA:40, 
    GL.HGA_s_H2O    :100,
    GL.FAH           :40,
    GL.C3H6O5P       :20},

 "[M+H]+":{ # http://dx.doi.org/10.1071/EN14069
    GL.MA         :1,
    GL.MA_s_HG_H2O:5,
    GL.MA_s_HG_FA :0,
    GL.HGA      :100,
    GL.HGA_s_PG  :50,
    GL.HGA_s_H2O :20,
    GL.HGA_s_PO3 :15,
    GL.FAkH       :0,
    GL.C5H5O2    :50
    }}#
  
  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=408.016655, type='Headgroup', chnops={'C':10, 'H':22, 'As':1, 'O':10, 'P':1},
    smiles='OC(COC1OC(C[As](C)(C)=O)C(O)C1O)COP(=O)(O)')
    super().__init__(AsPL.adducts, sn1, sn2, headgroup)

class AsPC(GL.Glycerolipid):

  tooltip = 'Arseno Phosphatidylcholine'
  tailOrganisation = ['A','A']
  arsenotails = {}

  adducts = {  # adduct:{spectra}

  "[M+H]+":{ # https://doi.org/10.1002/ange.201512031
    GL.MA           :1,
    GL.MA_s_H2O     :0,
    GL.MA_s_HG_H2O :30,
    GL.AsFAkH      :10,
    GL.AsFAH       :50,
    GL.C5H15NO4P  :100,
    GL.C2H6O4P     :20,
    GL.C5H14NO     :60,
    GL.C5H12N      :80}}
    
    # "[M+2H]2+":{ #
    # GL.MA          :10}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):

    try: 
      arsenotail = AsPC.arsenotails[sn1.name]
    except:
      arsenotail = copy.deepcopy(sn1)
      arsenotail.formula.update({'C':2, 'H':5, 'O':1, 'As':1})
      arsenotail.mass += 119.955637
      arsenotail.smiles = sn1.smiles+'[As](=O)(C)C'
      arsenotail.name += ':As'
      AsPC.arsenotails[sn1.name] = arsenotail
      AsPC.arsenotails[arsenotail.name] = arsenotail

    headgroup = GL.sn(mass=183.066044, type='Headgroup', chnops={'C':5, 'H':14, 'N':1, 'O':4, 'P':1},
    smiles='C[N+](C)(C)CCOP([O-])(=O)')
    # headgroup mass has -H to maintain neutral charge
    super().__init__(AsPC.adducts, arsenotail, sn2, headgroup)

# class AsPE(GL.Glycerolipid):

#   tooltip = 'Arseno Phosphatidylethanolamine'
#   tailOrganisation = ['A','A']
#   arsenotails = {}

#   adducts = {  # adduct:{spectra}

#   "[M+H]+":{ # https://doi.org/10.1002/ange.201512031
#     GL.MA           :0,
#     GL.MA_s_H2O     :0,
#     GL.MA_s_HG_H2O :50,
#     GL.AsFAkH     :100,
#     GL.AsFAH       :40},
    
#     "[M+2H]2+":{ #
#     GL.MA          :10}}

#   # sn3 = headgroup
#   def __init__(self, sn1, sn2):

#     try: 
#       arsenotail = AsPC.arsenotails[sn1.name]
#     except:
#       arsenotail = copy.deepcopy(sn1)
#       arsenotail.formula.update({'C':2, 'H':5, 'O':1, 'As':1})
#       arsenotail.mass += 119.955637
#       arsenotail.smiles = sn1.smiles+'[As](=O)(C)C'
#       arsenotail.name += ':As'
#       AsPC.arsenotails[sn1.name] = arsenotail
#       AsPC.arsenotails[arsenotail.name] = arsenotail

#     headgroup = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1},
#     smiles='NCCOP(O)(=O)')
#     super().__init__(AsPE.adducts, arsenotail, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PHEG(GL.Glycerolipid):

  tooltip = 'Phosphatidyl-O-[N-(2-hydroxyethyl) glycine]'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # https://doi.org/10.1007/s11745-017-4266-x
    GL.MA            :100,
    GL.MA_PO4_s_HG    :25,
    GL.MH_s_FAk        :2,
    GL.MH_s_FA         :2,
    GL.C3H7O5P_FA      :1,
    GL.MH_PO4_s_HG_FA:100,
    GL.FAH            :75,
    GL.FAH_PUFA       :20,
    GL.C3H6O5P        :25,
    GL.H2O4P           :0,
    GL.O3P             :0}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=199.024573, type='Headgroup', chnops={'C':4, 'H':9, 'N':1, 'O':6, 'P':1},
    smiles='OC(=O)CNCCOP(O)(=O)')
    super().__init__(PHEG.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PECM(GL.Glycerolipid):

  tooltip = 'Phosphatidylethanolamine-carboxymethyl'
  tailOrganisation = ['AA']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # https://doi.org/10.1007/s11745-017-4266-x
    GL.MA          :100,
    GL.MA_s_MeOH    :15,
    GL.MH_s_FAk      :2,
    GL.MH_s_FA       :2,
    GL.C3H7O5P_FA    :1,
    GL.MH_PO4_s_HG_FA:2,
    GL.FAH          :75,
    GL.FAH_PUFA     :20,
    GL.C3H6O5P       :2,
    GL.H2O4P         :0,
    GL.O3P           :0}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=199.024573, type='Headgroup', chnops={'C':4, 'H':9, 'N':1, 'O':6, 'P':1},
    smiles='COC(=O)NCCOP(O)(=O)')
    super().__init__(PECM.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #