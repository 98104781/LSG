import GenerateLipids as GL

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

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class MAG(GL.Glycerolipid):

  ##### Metabolites 2016, 6(3), 25; https://doi.org/10.3390/metabo6030025 # Needs work ?
  tooltip = 'Monoacylglycerol'
  tailOrganisation = ['T']
  No_Tails = 1

  adducts = {  # adduct:{spectra}
  "[M+H]+":{ # Needs Validation
    GL.MA      :100,
    GL.MA_s_H2O :50,
    GL.FAkH     :50},

  "[M+H-H2O]+":{ # Needs Validation
    GL.MA      :100,
    GL.FAkH     :50},

  "[M+Na]+":{ # Needs Validation    
    GL.MA      :100,
    GL.MA_s_H2O :50,
    GL.FAkA     :10,
    GL.FAkH     :50},

  "[M+NH4]+":{ # Matches LipidBlast    
    GL.MA        :1, # Fragmentation pattern matches
    GL.MH       :10, # LipidBlast, but the masses in LB are
    GL.MH_s_H2O:100, # consistantly off those predicted here.
    GL.FAH       :1, # Double checked the formula mass, appears
    GL.FAkH     :50},# to be correctly predicted here.
  }

  def __init__(self, sn1):
    super().__init__(MAG.adducts, sn1=sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class DAG(GL.Glycerolipid):

  tooltip = 'Diacylglycerol'
  No_Tails = 2
  tailOrganisation = ['TT']

  adducts = {  # adduct:{spectra}
  "[M+H]+":{ # Needs Validation
    GL.MA       :25,
    GL.MA_s_H2O:100,
    GL.MH_s_FA  :25,
    GL.FAkH      :1},

  "[M+H-H2O]+":{ # Needs Validation
    GL.MA      :100,
    GL.MH_s_FA :100,
    GL.FAkH      :1},

  "[M+Na]+":{ # Needs Validation    
    GL.MA       :50,
    GL.MA_s_H2O:100,
    GL.MA_s_FA  :20,
    GL.MH_s_FA  :50,
    GL.FAkA     :10,
    GL.FAkH      :1},

  "[M+NH4]+":{ # Matches LipidBlast    
    GL.MA       :25, # Molecular ion
    GL.MH       :10, # Present in lipidBlast
    GL.MH_s_H2O:100, # Present in lipidBlast
    GL.MH_s_FA :100, # Present in lipidBlast
    GL.FAkH      :1}
  }

  def __init__(self, sn1, sn2):
    super().__init__(DAG.adducts, sn1=sn1, sn2=sn2)

# ~ # ~ # ~ # ~ # ~ # ~ #

class TAG(GL.Glycerolipid):

  tooltip = 'Triacylglycerol'
  No_Tails = 3
  tailOrganisation = ['TTT']

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
    super().__init__(TAG.adducts, sn1, sn2, sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

class DGGA(GL.Glycerolipid):

  tooltip = 'Diacylglyceryl glucuronide'
  No_Tails = 2
  tailOrganisation = ['TT']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Taken from http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
    GL.MA        :100,
    GL.MA_s_FAk    :5,
    GL.MA_s_FA     :1,
    GL.FAH        :50},

  "[M+NH4]+":{ # Taken from http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html  
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
  No_Tails = 3
  tailOrganisation = ['T','TT']

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

# ~ # ~ # ~ # ~ # ~ # ~ #

class MGDG(GL.Glycerolipid): # Monogalactosyl diacylglycerol

  tooltip = 'Monoglycosyl diacylglycerol'
  No_Tails = 2
  tailOrganisation = ['TT']

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

class SQDG(GL.Glycerolipid): # Sulphoquinovosyl diacylglycerol

  tooltip = 'Sulfoquinovosyl diacylglycerol'
  No_Tails = 2
  tailOrganisation = ['TT']

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

class DGDG(GL.Glycerolipid): # Digalactosyl diacylglycerol

  tooltip = 'Diglycosyl diacylglycerol'
  No_Tails = 2
  tailOrganisation = ['TT']

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
  No_Tails = 2
  tailOrganisation = ['TT']

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
  No_Tails = 2
  tailOrganisation = ['TT']

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
  No_Tails = 2
  tailOrganisation = ['TT']

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
  No_Tails = 2
  tailOrganisation = ['TT']

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
  No_Tails = 3
  tailOrganisation = ['T','TT']

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

# ~ # ~ # ~ # ~ # ~ # ~ #

class AC4PIM2(GL.Glycerolipid):

  tooltip = 'Diacyl-phosphatidylinositol-dimannoside diacylglycerol'
  No_Tails = 4
  tailOrganisation = ['TT','TT']

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

# ~ # ~ # ~ # ~ # ~ # ~ #

class BMP(GL.Glycerolipid):

  tooltip = 'Bismonoacylglycerophosphate'
  No_Tails = 2
  tailOrganisation = ['TT']

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

  tooltip = 'Hemibismonoacylglycerophosphate'
  No_Tails = 3
  tailOrganisation = ['T','TT']

  adducts = {  # adduct:{spectra}
  "[M+NH4]+":{
    GL.MA              :5,
    GL.HGH_s_PO4      :20,
    GL.MH_s_HG_H2O   :100,
    GL.MH_s_HG_FA     :20},

  "[M-H]-":{
    GL.MA              :5}
    }

  def __init__(self, sn1, sn2, sn3):
    headgroup = GL.sn(mass=172.013674, type='Headgroup', chnops={'C':3, 'H':9, 'O':6, 'P':1},
    smiles=sn1.inverseSmiles+'OCC(O)COP(=O)(O)', hgtails = [sn1])
    super().__init__(HBMP.adducts, sn1=sn2, sn2=sn3, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class CL(GL.Glycerolipid):

  tooltip = 'Cardiolipin'
  No_Tails = 4
  tailOrganisation = ['TT','TT']

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

# ~ # ~ # ~ # ~ # ~ # ~ #

class PA(GL.Glycerolipid):

  tooltip = 'Phosphatidic acid'
  No_Tails = 2
  tailOrganisation = ['TT']

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

  "[M+H]+":{ # Needs Validation
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
    headgroup = GL.sn(mass=97.976895, type='Headgroup', chnops={'H':3, 'O':4, 'P':1},
    smiles='O=P(O)(O)')
    super().__init__(PA.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PMe(GL.Glycerolipid):

  tooltip = 'Phosphatidylmethanol'
  No_Tails = 2
  tailOrganisation = ['TT']

  adducts = {  
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
  No_Tails = 2
  tailOrganisation = ['TT']

  adducts = {  
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
  No_Tails = 2
  tailOrganisation = ['TT']

  adducts = {  # adduct:{spectra}
  "[M+H]+":{ # Matches literature but not lipidblast
    GL.MA          :10, # Lipidblast seems to provide a completely
    GL.MA_s_H2O     :0, # different spectra which is not consistant
    GL.MA_s_TMA     :0, # with literature. Perhaps a result of
    GL.MA_s_HG_H2O  :0, # instrumentation.
    GL.MA_s_FA_TMA  :0,
    GL.MH_s_FAk     :1,
    GL.MH_s_FA      :1,
    GL.C5H15NO4P  :100},

  "[M+Na]+":{ # Looks Good, Matches lipidblast
    GL.MA          :10, # Molecular ion
    GL.MA_s_TMA   :100, # Present in lipidblast
    GL.MA_s_HG_H2O :60, # Present in lipidblast
    GL.MH_s_HG_H2O  :0,
    GL.MA_s_FA_TMA  :4, # Present in lipidblast
    GL.MA_s_FAk     :2, # Present in lipidblast
    GL.MA_s_FA      :2, # Present in lipidblast
    GL.MH_s_FA      :0,
    GL.FAkH         :0,
    GL.C2H5NaO4P    :0}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=183.066044, type='Headgroup', chnops={'C':5, 'H':14, 'N':1, 'O':4, 'P':1},
    smiles='C[N+](C)(C)CCOP([O-])(=O)')
    # headgroup mass has -H to maintain neutral charge
    super().__init__(PC.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PE(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  ##### "[M+H]+", "[M-H]-" https://doi.org/10.1039/C5AY00776C
  tooltip = 'Phosphatidylethanolamine'
  No_Tails = 2
  tailOrganisation = ['TT']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Looks Good, Matches lipidblast
    GL.MA       :15, # Molecular ion
    GL.MH_s_FAk  :5, # Present in lipidblast
    GL.MH_s_FA   :2, # Present in lipidblast
    GL.FAH     :100, # Present in lipidblast
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
    GL.MA_s_HG_FA   :0,
    GL.FAkH         :1},# Present in lipidblast
    
  "[M+Na]+":{ # Looks Good, Matches lipidblast
    GL.MA          :10, # Molecular ion
    GL.MA_s_AZD   :100, # Present in lipidblast
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

class MMPE(GL.Glycerolipid):

  tooltip = 'N-methyl-phosphatidylethanolamine'
  No_Tails = 2
  tailOrganisation = ['TT']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    GL.MA          :15,
    GL.FAH        :100,
    GL.MH_s_FAk    :10,
    GL.MH_s_FA     :10,
    GL.C3H6O5P     :10,
    GL.H2O4P        :5,
    GL.O3P          :5}}
  
  #"[M+H]+":{
  #  GL.MA          :10,
  #  GL.MA_s_HG_H2O :10}}
    
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
  No_Tails = 2
  tailOrganisation = ['TT']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    GL.MA          :15,
    GL.FAH        :100,
    GL.MH_s_FAk    :10,
    GL.MH_s_FA     :10,
    GL.C3H6O5P     :10,
    GL.H2O4P        :5,
    GL.O3P          :5}}#,
  
  #"[M+H]+":{
  #  GL.MA          :10},
    
  #"[M+Na]+":{
  #  GL.MA          :10}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=169.050394, type='Headgroup', chnops={'C':4, 'H':12, 'N':1, 'O':4, 'P':1},
    smiles='N(C)(C)CCOP(O)(=O)')
    super().__init__(DMPE.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PG(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  tooltip = 'Phosphatidylglycerol'
  No_Tails = 2
  tailOrganisation = ['TT']

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
  No_Tails = 2
  tailOrganisation = ['TT']

  adducts = {  # adduct:{spectra}
  #"[M+Na-2H]-":{
  #  GL.MH        :2,
  #  GL.MH_s_FA   :4,
  #  GL.HG_FA_NL_B:1,
  #  GL.FAH     :100,
  #  GL.C9H16O10P:40,
  #  GL.C9H14O9P :10, 
  #  GL.C6H10O8P :10,
  #  GL.C6H8O7P  :40,
  #  Gl.C3H7NaO6P:50,
  #  GL.C3H6O5P  :60}, 

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

  "[M+NH4]+":{ # https://doi.org/10.1002/pld3.183 Figure S4   
    GL.MA           :5,
    GL.MH           :5,
    GL.MH_s_HG_H2O:100}}

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=260.029718, type='Headgroup', chnops={'C':6, 'H':13, 'O':9, 'P':1},
    smiles='OC1C(O)C(O)C(O)C(O)C1OP(O)(=O)')
    super().__init__(PI.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PIP(GL.Glycerolipid):

  tooltip = 'Phosphatidylinositol phosphate'
  No_Tails = 2
  tailOrganisation = ['TT']

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

class PI2P(GL.Glycerolipid):

  tooltip = 'Phosphatidylinositol diphosphate'
  No_Tails = 2
  tailOrganisation = ['TT']

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
    super().__init__(PI2P.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PS(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  ##### "[M+H]+" https://doi.org/10.1039/C5AY00776C
  tooltip = 'Phosphatidylserine'
  No_Tails = 2
  tailOrganisation = ['TT']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Looks Good, Matches lipidblast
    GL.MA            :15, # Molecular ion
    GL.MH_PO4_s_HG  :100, # Present in lipidblast
    GL.HG_FAk_NL_B   :25, # Present in lipidblast
    GL.MH_PO4_s_HG_FA:50, # Present in lipidblast
    GL.FAH           :50, # Present in lipidblast
    GL.C3H6O5P       :30,
    GL.H2O4P          :5,
    GL.O3P            :5},

  "[M+H]+":{ # Looks Good, Matches lipidblast
    GL.MA          :10, # Molecular ion
    GL.MH_PO4_s_HG :40, # Present in lipidblast
    GL.MA_s_HG_H2O:100, # Present in lipidblast
    GL.MA_s_FAk    :20, # Present in lipidblast
    GL.MA_s_FA     :20, # Present in lipidblast
    GL.MA_s_HG_FA   :0,
    GL.FAkH         :0},
    
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

class PPA(GL.Glycerolipid):

  tooltip = 'Pyrophosphatidic acid'
  No_Tails = 2
  tailOrganisation = ['TT']

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

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyDGCC(GL.Glycerolipid): # Diacylglyceryl-3-O-carboxyhydroxymethylcholine

  tooltip = 'Lyso-3-O-carboxyhydroxymethylcholine' 
  No_Tails = 1
  tailOrganisation = ['T']

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
    super().__init__(lyDGCC.adducts, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyDGTS(GL.Glycerolipid): # LysoN-trimethylhomoserine

  tooltip = 'Lyso-N-trimethylhomoserine' 
  No_Tails = 1
  tailOrganisation = ['T']

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
    super().__init__(lyDGTS.adducts, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPA(GL.Glycerolipid):

  ##### [M-H]- = https://doi.org/10.1002/lipd.12172, Not a fragmentation study!
  ##### [M-H]- = https://doi.org/10.1016/j.jchromb.2010.03.030, Neither...
  tooltip = 'Lyso-phosphatidic acid'
  No_Tails = 1
  tailOrganisation = ['T']
  
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
    super().__init__(lyPA.adducts, sn1=sn1, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPC(GL.Glycerolipid):

  tooltip = 'Lyso-phosphatidylcholine'
  No_Tails = 1
  tailOrganisation = ['T']

  adducts = {  # adduct:{spectra}
  "[M+H]+":{ # Matches lipidblast
    GL.MA         :5, # Molecular ion
    GL.MA_s_H2O  :20, # Present in lipidblast
    GL.MH_s_FAk  :10,
    GL.MH_s_FA    :5, # Present in lipidblast
    GL.FAkH       :1, # Present in lipidblast
    GL.C5H15NO4P:100} # Present in lipidblast
    }#,

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
    super().__init__(lyPC.adducts, sn1=sn1, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPE(GL.Glycerolipid):

  tooltip = 'Lyso-phosphatidylethanolamine'
  No_Tails = 1
  tailOrganisation = ['T']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Matches lipidblast
    GL.MA            :15, # Molecular ion
    GL.MH_PO4_s_HG_H2O:1, # Present in lipidblast
    GL.MH_s_FAk       :5, # Present in lipidblast
    GL.MH_s_FA       :15, # Present in lipidblast
    GL.FAH          :100, # Present in lipidblast
    GL.C5H11NO5P      :5,# Present in lipidblast
    GL.C3H6O5P       :10,
    GL.C2H7NO4P       :3,
    GL.H2O4P          :5,
    GL.O3P            :5},
  
  "[M+H]+":{ # Matches lipidblast
    GL.MA          :10, # Molecular ion
    GL.MA_s_H2O    :20, # Present in lipidblast
    GL.MA_s_HG_H2O:100, # Present in lipidblast
    GL.MH_PO4_s_HG_H2O :20, # Present in lipidblast
    GL.MH_s_FAk     :1, # Present in lipidblast
    GL.MH_s_FA      :1, # Present in lipidblast
    GL.FAkH         :0}
    }

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    headgroup = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1},
    smiles='NCCOP(O)(=O)')
    super().__init__(lyPE.adducts, sn1=sn1, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyNAPE(GL.Glycerolipid):

  tooltip = 'Lyso-N-Acyl-phosphatidylethanolamine'
  No_Tails = 2
  tailOrganisation = ['T','T']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{
    GL.MA            :15,
    GL.MH_PO4_s_HG_H2O:1,
    GL.MH_s_FAk       :5,
    GL.MH_s_FA       :15,
    GL.FAH          :100,
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
    super().__init__(lyNAPE.adducts, sn1=sn2, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPG(GL.Glycerolipid):

  #####  Requires reference! Needs work!
  tooltip = 'Lyso-phosphatidylglycerol'
  No_Tails = 1
  tailOrganisation = ['T']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Needs Validation
    GL.MH           :100,
    GL.MH_PO4_s_HG    :1,
    GL.MH_PO4_s_HG_H2O:1,
    GL.FAH          :100,
    GL.C3H6O5P       :50,
    GL.H2O4P          :1,
    GL.O3P            :1},

  "[M+H]+":{ # Needs Validation
    GL.MA          :10,
    GL.MA_s_HG_H2O:100,
    GL.MA_s_HG_FA  :10,
    GL.FAkH        :10}
    }

  # sn3 = headgroup
  def __init__(self, sn1):
    headgroup = GL.sn(mass=172.013674, type='Headgroup', chnops={'C':3, 'H':9, 'O':6, 'P':1},
    smiles='OCC(O)COP(O)(=O)')
    super().__init__(lyPG.adducts, sn1=sn1, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPI(GL.Glycerolipid):

  #####  Requires reference!
  tooltip = 'Lyso-phosphatidylinositol'
  No_Tails = 1
  tailOrganisation = ['T']

  adducts = {  # adduct:{spectra}
  #"[M+Na-2H]-":{
  #  GL.C9H16O10P:50, # Not confident what the mz 192.988  is,
  #  GL.C6H10O8P :10, # perhaps as a [Phosphoglycerol+Na-2H]-.
  #  GL.C6H8O7P  :10, # Appears in "[M+Na-2H]-" spectra.
  #  Gl.C3H7NaO6P:80, # Additional mz 355.04  not included.
  #  GL.C3H6O5P  :80},

  "[M-H]-":{ # Needs Validation
    GL.MA             :2,
    GL.MH_PO4_s_HG    :1,
    GL.MH_PO4_s_HG_H2O:1,
    GL.FAH          :100, 
    GL.C9H16O10P     :10, 
    GL.C6H10O8P      :45,  # mz 233.001 can be a major -
    GL.C6H8O7P       :10,  # or minor fragment depending -
    GL.C3H6O5P       :45}  # on if its a sn1 or sn2 lyso lipid.
    } 

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    headgroup = GL.sn(mass=260.029718, type='Headgroup', chnops={'C':6, 'H':13,'O':9, 'P':1},
    smiles='OC1C(O)C(O)C(O)C(O)C1OP(O)(=O)')
    super().__init__(lyPI.adducts, sn1=sn1, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPS(GL.Glycerolipid):

  #####  Requires reference! Needs work!
  tooltip = 'Lyso-phosphatidylserine'
  No_Tails = 1
  tailOrganisation = ['T']

  adducts = {  # adduct:{spectra}
  "[M+H]+":{ # Needs Validation
    GL.MA          :10,
    GL.MA_s_HG_H2O:100,
    GL.MA_s_HG_FA  :10,
    GL.FAkH        :10}
    }#,

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    headgroup = GL.sn(mass=185.008923, type='Headgroup', chnops={'C':3, 'H':8, 'N':1, 'O':6, 'P':1},
    smiles='[O-]C(=O)C([NH3+])COP(O)(=O)')
    super().__init__(lyPS.adducts, sn1=sn1, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyNAPS(GL.Glycerolipid):

  #####  Requires reference! Needs work!
  tooltip = 'Lyso-N-Acyl-phosphatidylserine'
  No_Tails = 2
  tailOrganisation = ['T','T']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Needs Validation
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
    super().__init__(lyNAPS.adducts, sn1=sn2, sn3=headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class Acylsphinganine(GL.Sphingolipid):

  tooltip = 'N-Acyl-sphinganine'
  base_types = ['Sphinganine'] # 18:0;O2
  No_Tails = 1
  tailOrganisation = ['B','T']

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
    GL.MA           :5,
    GL.MA_s_H2O     :5,
    GL.MA_s_2H2O    :5,
    GL.MA_s_FA     :10,
    GL.MA_s_FA_H2O:100,
    GL.Cer_U       :10,
    GL.Cer_D       :10}
    }

  #"[M+H-H2O]+":{ #
  #  GL.MA         :100
  #  }}#,

  def __init__(self, base, sn1):
    super().__init__(Acylsphinganine.adducts, base, sn1)

class Acylsphingosine(GL.Sphingolipid):

  tooltip = 'N-Acyl-sphingosine'
  base_types = ['Sphingosine'] # 18:1;O2
  No_Tails = 1
  tailOrganisation = ['B','T']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # https://doi.org/10.1006/abio.2001.5536, https://doi.org/10.1016/S1044-0305(02)00358-6, https://doi.org/10.1016/j.biochi.2016.07.012, https://doi.org/10.1002/rcm.878
    GL.MA           :30,
    GL.MA_s_H2O      :2,
    GL.MA_s_CH2O    :20, # This fragment not present with sphinganine
    GL.MA_s_MeOH    :10,
    GL.MA_s_CH2O_H2O:10,
    GL.Cer_P        :10,
    GL.Cer_R        :26,
    GL.Cer_S        :20,
    GL.Cer_T       :100,
    GL.Cer_U         :5,
    GL.FAH          :25,
    GL.FAkH         :10},
    
  "[M+H]+":{ # https://doi.org/10.1002/bmc.4790
    GL.MA          :20,
    GL.MA_s_H2O    :15,
    GL.MA_s_FA     :10,
    GL.MA_s_FA_H2O:100,
    GL.Cer_U       :10}
    }#,

  #"[M+H-H2O]+":{ #
  #  GL.MA         :100
  #  }}#,

  def __init__(self, base, sn1):
    super().__init__(Acylsphingosine.adducts, base, sn1)

class Acylphytosphingosine(GL.Sphingolipid):

  tooltip = 'N-Acyl-phytosphingosine'
  base_types = ['Phytosphingosine'] # 18:0;O3
  No_Tails = 1
  tailOrganisation = ['B','T']

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
    GL.FAkH      :1}
    }#,
    
  #"[M+H]+":{
  #  GL.MA          :20
  #  }}#,

  #"[M+H-H2O]+":{ #
  #  GL.MA         :100
  #  }}#,

  def __init__(self, base, sn1):
    super().__init__(Acylphytosphingosine.adducts, base, sn1)

class CerP(GL.Sphingolipid):

  tooltip = 'N-Acyl-ceramide-1-phosphate'
  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  No_Tails = 1
  tailOrganisation = ['B','T']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # Matches lipidblast, https://doi.org/10.1016/j.biochi.2016.07.012
    GL.MA         :10, # Molecular ion
    GL.MA_s_H2O   :10, # Present in lipidblast
    GL.MA_s_FAk   :20, # Present in lipidblast
    GL.MA_s_FA    :20, # Present in lipidblast
    GL.H2O4P     :100, # Present in lipidblast
    GL.H2O2P       :0,
    GL.O3P       :100},# Present in lipidblast
    
  "[M+H]+":{
    GL.MA             :10,
    GL.MA_s_H2O       :10,
    GL.MA_s_HG_H2O    :10,
    GL.MA_s_HG_2H2O   :10,
    GL.MA_s_HG_FA_H2O:100
    }}

  #"[M+H-H2O]+":{ #
  #  GL.MA         :100
  #  }}#,

  def __init__(self, base, sn1):
    headgroup = GL.sn(mass=97.976895, type='Headgroup', chnops={'H':3, 'O':4, 'P':1},
    smiles='O=P(O)(O)')
    super().__init__(CerP.adducts, base, sn1, headgroup)

class CerPC(GL.Sphingolipid):

  tooltip = 'N-Acyl-ceramide-1-phosphocholine'
  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  No_Tails = 1
  tailOrganisation = ['B','T']

  adducts = {  # adduct:{spectra}
  "[M+H]+":{ # Looks Good, Matches lipidblast
    GL.MA          :10, # Molecular ion
    GL.MA_s_H2O   :100, # Present in lipidblast
    GL.MA_s_TMA     :1, # Present in lipidblast
    GL.MA_s_TMA_H2O :5, # Present in lipidblast
    GL.MA_s_HG_H2O  :2, # Present in lipidblast
    GL.C5H15NO4P  :100}
    }

  def __init__(self, base, sn1):
    headgroup = GL.sn(mass=183.066044, type='Headgroup', chnops={'C':5, 'H':14, 'N':1, 'O':4, 'P':1},
    smiles='C[N+](C)(C)CCOP([O-])(=O)')
    # headgroup mass has -H to maintain neutral charge
    super().__init__(CerPC.adducts, base, sn1, headgroup)

class CerPE(GL.Sphingolipid):

  tooltip = 'N-Acyl-ceramide-1-phosphoethanolamine'
  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  No_Tails = 1
  tailOrganisation = ['B','T']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
    GL.MA        :100,
    GL.MA_s_H2O    :1,
    GL.MA_s_FAk    :5,
    GL.MA_s_FA     :1,
    GL.C2H7NO4P   :25,
    GL.H2O4P      :25,
    GL.H2O2P       :0,
    GL.O3P        :25}
    }

  def __init__(self, base, sn1):
    headgroup = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1},
    smiles='NCCOP(O)(=O)')
    super().__init__(CerPE.adducts, base, sn1, headgroup)

class CerPI(GL.Sphingolipid):

  tooltip = 'N-Acyl-ceramide-1-phosphoinositol'
  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  No_Tails = 1
  tailOrganisation = ['B','T']

  adducts = {  # adduct:{spectra}
  "[M-H]-":{ # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
    GL.MA        :100,
    GL.MA_s_H2O    :1,
    GL.MA_PO4_s_HG :2,
    GL.MH_PO4_s_HG_FA:3,
    GL.MA_s_FAk    :5,
    GL.MA_s_FA     :1,
    GL.C6H10O8P   :40,
    GL.C6H8O7P    :15, 
    GL.H2O4P      :25,
    GL.O3P        :25}
    }

  def __init__(self, base, sn1):
    headgroup = GL.sn(mass=260.029718, type='Headgroup', chnops={'C':6, 'H':13, 'O':9, 'P':1},
    smiles='OC1C(O)C(O)C(O)C(O)C1OP(O)(=O)')
    super().__init__(CerPI.adducts, base, sn1, headgroup)

class HexCer(GL.Sphingolipid):

  tooltip = 'N-Acyl-ceramide-1-hexose'
  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  No_Tails = 1
  tailOrganisation = ['B','T']

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
    GL.Cer_U           :5,
    GL.Cer_D          :10}
    }

  def __init__(self, base, sn1):
    headgroup = GL.sn(mass=180.063388, type='Headgroup', chnops={'C':6, 'H':12, 'O':6},
    smiles='OC1C(O)C(O)C(CO)OC1')
    super().__init__(HexCer.adducts, base, sn1, headgroup)

class AHexCer(GL.Sphingolipid): # doesnt work with how generation is set up...

  tooltip = 'N-Acyl-ceramide-1-Acylhexose'
  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  No_Tails = 2
  tailOrganisation = ['B','T','T']

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

class Hex2Cer(GL.Sphingolipid):

  tooltip = 'N-Acyl-ceramide-1-dihexose'
  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  No_Tails = 1
  tailOrganisation = ['B','T']

  adducts = {  # adduct:{spectra}
  "[M+H]+":{ # http://prime.psc.riken.jp/compms/msdial/lipidnomenclature.html
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

class Hex3Cer(GL.Sphingolipid):

  tooltip = 'N-Acyl-ceramide-1-trihexose'
  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  No_Tails = 1
  tailOrganisation = ['B','T']

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

class Sulfatide(GL.Sphingolipid):

  tooltip = 'Sulfatide'
  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  No_Tails = 1
  tailOrganisation = ['B','T']

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
    super().__init__(Sulfatide.adducts, base, sn1, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class FA(GL.OtherLipid):

  tooltip = 'Free Fatty acid'
  No_Tails = 1
  tailOrganisation = ['T']

  adducts = {
    "[M-H]-":{
    GL.MH            :100,
    GL.MH_s_H2O        :0
    }}

  def __init__(self, sn1):
    body = GL.Other(name='FA', smiles='O'+sn1.smiles)
    super().__init__(FA.adducts, body, sn1)

class CE(GL.OtherLipid):

  tooltip = 'Cholesteryl Ester'
  No_Tails = 1
  tailOrganisation = ['T']

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

class BRSE(GL.OtherLipid):

  tooltip = 'Brassicasterol Ester'
  No_Tails = 1
  tailOrganisation = ['T']

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

class CASE(GL.OtherLipid):

  tooltip = 'Camposterol Ester'
  No_Tails = 1
  tailOrganisation = ['T']

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

class SISE(GL.OtherLipid):

  tooltip = 'Sitosterol Ester'
  No_Tails = 1
  tailOrganisation = ['T']

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

class STSE(GL.OtherLipid):

  tooltip = 'Sigmasterol Ester'
  No_Tails = 1
  tailOrganisation = ['T']

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

class LE(GL.OtherLipid):

  tooltip = 'Lanosteryl Ester'
  No_Tails = 1
  tailOrganisation = ['T']

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

class AC(GL.OtherLipid):

  tooltip = 'Acyl-carnitine'
  No_Tails = 1
  tailOrganisation = ['T']

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

class AcylCoA(GL.OtherLipid):

  tooltip = 'Acyl-CoA'
  No_Tails = 1
  tailOrganisation = ['T']
  
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