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

  No_Tails = 2
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

  No_Tails = 3
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

class AcylGlcADG(GL.Glycerolipid):

  No_Tails = 3
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Taken from LipidMatch 
    GL.MA       :10,
    GL.FAH      :10}
  }

  def __init__(self, sn1, sn2, sn3):
    headgroup = GL.sn(mass=194.042653, type='Headgroup', chnops={'C':6, 'H':10, 'O':7},
    smiles='O('+sn1.smiles+')C1C(O)C(O)C(C(=O)O)OC1', hgtails = [sn1])
    super().__init__(AcylGlcADG.adducts, sn1=sn2, sn2=sn3, sn3=headgroup)

class MGDG(GL.Glycerolipid): # Monogalactosyl diacylglycerol

  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M+Na]+":{ # https://doi.org/10.1111/j.1440-1835.2010.00582.x   
    GL.MA          :5, # This spectra agrees with lipidblast and the attached
    GL.MA_s_FA   :100},# reference. Ion ratios are off, but depend on isomerism.

  "[M+NH4]+":{ # https://doi.org/10.1002/pld3.183 Figure S4     
    GL.MA          :5, # LipidBlast includes fragment 'MA_s_FA' but
    GL.MH_s_HG     :5, # the fragment is not present in the reference.
    GL.MH_s_HG_H2O:60, # Included here at 0 intensity so it can be modified
    GL.MA_s_FA     :0, # in the program if needed, otherwise it will be
    GL.HG_FA_NL_C:100} # ignored on spectra generation.
    }

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=180.063388, type='Headgroup', chnops={'C':6, 'H':12, 'O':6},
    smiles='OC1C(O)C(O)C(CO)OC1')
    super().__init__(MGDG.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class SQDG(GL.Glycerolipid): # Sulphoquinovosyl diacylglycerol

  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Closely matches LipidBlast  
    GL.MA         :40, # Molecular ion
    GL.MA_s_FA    :10, # Present in lipidBlast
    GL.FAH        :10, # Present in lipidBlast
    GL.C6H9O7S   :100, # Present in lipidBlast
    GL.O3SH       :20}}

  #"[M+NH4]+":{    
  #  GL.MA          :5,
  #  GL.
  #  }

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=244.025287, type='Headgroup', chnops={'C':6, 'H':12, 'O':8, 'S':1},
    smiles='OC1C(O)C(O)C(CS(=O)(=O)O)OC1')
    super().__init__(SQDG.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class DGDG(GL.Glycerolipid): # Digalactosyl diacylglycerol

  No_Tails = 2
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
    GL.MA_s_Gal_H2O    :5,
    GL.MA_s_FA       :100, # Present in lipidBlast
    GL.MA_s_FA_Gal_H2O:50, # Present in lipidBlast
    },

  "[M+NH4]+":{ # https://doi.org/10.1002/pld3.183 Figure S4   
    GL.MA          :5,
    GL.MH_s_Gal_H2O:5,
    GL.MH_s_HG    :70,
    GL.MH_s_HG_H2O:70,
    GL.MA_s_FA     :0,
    GL.HG_FA_NL_C:100}
    }

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=342.116212, type='Headgroup', chnops={'C':12, 'H':22, 'O':11},
    smiles='OC1C(O)C(O)C(CO(C2OC(CO)C(O)C(O)C2(O)))OC1')
    super().__init__(DGDG.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class DGTS(GL.Glycerolipid): # N-trimethylhomoserine diacylglycerol

  No_Tails = 2
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

  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Matches LipidBlast
    GL.MA            :50, # Molecular ion
    GL.MA_s_Gal_H2O  :20, # Present in lipidBlast
    GL.MA_s_FA      :100, # Present in lipidBlast
    GL.MH_PO4_s_HG   :60, # Present in lipidBlast
    GL.MH_PO4_s_HG_FA:40, # Present in lipidBlast
    GL.FAH           :40}}# Present in lipidBlast

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=422.082541, type='Headgroup', chnops={'C':12, 'H':23, 'O':14, 'P':1},
    smiles='OC1C(O)C(O)C(O)C(OC2OC(CO)C(O)C(O)C2O)C1OP(O)(=O)')
    super().__init__(AC2PIM1.adducts, sn1, sn2, headgroup)

class AC2PIM2(GL.Glycerolipid):

  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Matches LipidBlast
    GL.MA            :10, # Molecular ion
    GL.MA_s_Gal_H2O   :5, # Present in lipidBlast
    GL.MA_s_FA      :100, # Present in lipidBlast
    GL.MA_s_allFA    :20, # Present in lipidBlast
    GL.FAH           :40}}# Present in lipidBlast

  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=584.135365, type='Headgroup', chnops={'C':18, 'H':33, 'O':19, 'P':1},
    smiles='O(C3OC(CO)C(O)C(O)C3O)C1C(O)C(O)C(O)C(OC2OC(CO)C(O)C(O)C2O)C1OP(O)(=O)')
    super().__init__(AC2PIM2.adducts, sn1, sn2, headgroup)

class AC3PIM2(GL.Glycerolipid):

  No_Tails = 3
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Matches LipidBlast
    GL.MA            :20, # Molecular ion
    GL.MA_s_Gal_H2O   :5, # Present in lipidBlast
    GL.MA_s_FAk      :10, # Present in lipidBlast
    GL.MA_s_FA      :100, # Present in lipidBlast
    GL.MA_s_allFA    :20, # Present in lipidBlast
    GL.HGA_s_H2O     :40, # Present in lipidBlast
    GL.FAH           :40}}# Present in lipidBlast

  def __init__(self, sn1, sn2, sn3):
    headgroup = GL.sn(mass=584.135365, type='Headgroup', chnops={'C':18, 'H':33, 'O':19, 'P':1},
    smiles='O(C3OC(CO)C(O)C(O)C3O)C1C(O)C(O)C(O)C(OC2OC(CO('+sn1.smiles+'))C(O)C(O)C2O)C1OP(O)(=O)', hgtails = [sn1])
    super().__init__(AC3PIM2.adducts, sn1=sn2, sn2=sn3, sn3=headgroup)

class AC4PIM2(GL.Glycerolipid):

  No_Tails = 4
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

class Cardiolipin(GL.Glycerolipid):

  No_Tails = 4
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
    super().__init__(Cardiolipin.adducts, sn1=sn3, sn2=sn4, sn3=headgroup)

class PA(GL.Glycerolipid):

  No_Tails = 2
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

class PC(GL.Glycerolipid):

  ##### "[M+H]+", "[M+Na]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  ##### "[M+H]+", "[M+Na]+"  10.1016/S1044-0305(03)00064-3
  ##### "[M+H]+" https://doi.org/10.1039/C5AY00776C 
  No_Tails = 2
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
  No_Tails = 2
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

  No_Tails = 2
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
    headgroup = GL.sn(mass=155.034744, type='Headgroup', chnops={'C':3, 'H':8, 'N':1, 'O':4, 'P':1},
    smiles='N(C)CCOP(O)(=O)')
    super().__init__(MMPE.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class DMPE(GL.Glycerolipid):

  No_Tails = 2
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
    headgroup = GL.sn(mass=169.050394, type='Headgroup', chnops={'C':4, 'H':8, 'N':1, 'O':4, 'P':1},
    smiles='N(C)(C)CCOP(O)(=O)')
    super().__init__(DMPE.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PG(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  No_Tails = 2
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
  No_Tails = 2
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

  #####  Requires reference!
  No_Tails = 2
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

  "[M-H]-":{ # Looks Good
    GL.MH            :1,
    GL.MH_s_H2O      :2,
    GL.MH_s_PO3      :2,
    GL.MH_P2O6_s_HG  :2,
    GL.MH_s_FA       :4,
    GL.MH_s_FA_H2O   :1,
    GL.MH_s_FAk_PO3  :4,
    GL.MH_s_FA_PO3   :4,
    GL.MH_PO4_s_HG_FA:5,
    GL.FAH         :100,
    GL.C6H11O11P2   :10, 
    GL.C6H9O10P2    :15,
    GL.C6H12O9P      :5, 
    GL.C6H10O8P     :90,
    GL.C6H8O7P      :90, 
    GL.C3H6O5P      :40}}#,

  #"[M-2H]2-":{ # Looks Good
  #  GL.MA          :5,
  #  GL.FAH       :100,
  #  GL.C6H11O11P2  :2, 
  #  GL.C6H9O10P2   :1,
  #  GL.C6H10O8P   :40,
  #  GL.C6H8O7P     :1, 
  #  GL.C3H6O5P    :10}
  # }

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=339.996048, type='Headgroup', chnops={'C':6, 'H':14, 'O':12, 'P':2},
    smiles='OC1C(O)C(O)C(OP(=O)(O)O)C(O)C1OP(O)(=O)')
    super().__init__(PIP.adducts, sn1, sn2, headgroup)

# ~ # ~ # ~ # ~ # ~ # ~ #
'''
class PI2P(GL.Glycerolipid):

  #####  Requires reference! Needs work!
  No_Tails = 2
  adducts = {  # adduct:{spectra} 

  #"[M+Na-2H]-":{
  #},

  #"[M-H]-":{
  #},

  #"[M-2H]2-":{
  #},

  #"[M-3H]3-":{
  #}
  } 

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    headgroup = GL.sn(mass=419.962378, type='Headgroup', chnops={'C':6, 'H':15, 'O':15, 'P':3},
    smiles='OC1C(OP(=O)(O)O)C(O)C(OP(=O)(O)O)C(O)C1OP(O)(=O)')
    super().__init__(PI2P.adducts, sn1, sn2, headgroup)
'''
# ~ # ~ # ~ # ~ # ~ # ~ #

class PS(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  ##### "[M+H]+" https://doi.org/10.1039/C5AY00776C 
  No_Tails = 2
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

  No_Tails = 2
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

class lyPA(GL.Glycerolipid):

  ##### [M-H]- = https://doi.org/10.1002/lipd.12172, Not a fragmentation study!
  ##### [M-H]- = https://doi.org/10.1016/j.jchromb.2010.03.030, Neither...
  No_Tails = 1
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

  No_Tails = 1
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

  No_Tails = 1
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Matches lipidblast
    GL.MA            :15, # Molecular ion
    GL.MH_PO4_s_HG_H2O:1, # Present in lipidblast
    GL.MH_s_FAk       :5, # Present in lipidblast
    GL.MH_s_FA       :15, # Present in lipidblast
    GL.FAH          :100, # Present in lipidblast
    #GL.C5H11NO5P     :5,# Present in lipidblast
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

class lyPG(GL.Glycerolipid):

  #####  Requires reference! Needs work!
  No_Tails = 1
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
  No_Tails = 1
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
  No_Tails = 1
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

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class Acylsphinganine(GL.Sphingolipid):

  base_types = ['Sphinganine'] # 18:0;O2
  No_Tails = 1
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
    GL.FAkH         :10}
    }
    
  #"[M+H]+":{
  #  GL.MA         :100
  #  }}

  #"[M+H-H2O]+":{ #
  #  GL.MA         :100
  #  }}#,

  def __init__(self, base, sn1):
    super().__init__(Acylsphinganine.adducts, base, sn1)

class Acylsphingosine(GL.Sphingolipid):

  base_types = ['Sphingosine'] # 18:1;O2
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # https://doi.org/10.1006/abio.2001.5536, https://doi.org/10.1016/S1044-0305(02)00358-6, https://doi.org/10.1016/j.biochi.2016.07.012, https://doi.org/10.1002/rcm.878
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
    
  "[M+H]+":{ # https://doi.org/10.1002/bmc.4790
    GL.MA       :20,
    GL.MA_s_H2O :15,
    GL.MA_s_FA  :10,
    GL.Cer_R   :100,
    GL.Cer_U    :10}
    }#,

  #"[M+H-H2O]+":{ #
  #  GL.MA         :100
  #  }}#,

  def __init__(self, base, sn1):
    super().__init__(Acylsphingosine.adducts, base, sn1)

class Acylphytosphingosine(GL.Sphingolipid):

  base_types = ['Phytosphingosine'] # 18:0;O3
  No_Tails = 1
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

  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  No_Tails = 1
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
    GL.MA          :10,
    GL.MA_s_H2O    :10,
    GL.MA_s_HG_H2O :10,
    GL.MA_s_HG_2H2O:10,
    GL.Cer_R      :100,
    }}

  #"[M+H-H2O]+":{ #
  #  GL.MA         :100
  #  }}#,

  def __init__(self, base, sn1):
    headgroup = GL.sn(mass=97.976895, type='Headgroup', chnops={'H':3, 'O':4, 'P':1},
    smiles='O=P(O)(O)')
    super().__init__(CerP.adducts, base, sn1, headgroup)

class CerPC(GL.Sphingolipid):

  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  No_Tails = 1
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

class Sulfatide(GL.Sphingolipid):

  base_types = ['Sphinganine', 'Sphingosine']  # 18:0;O2, 18:1;O2
  No_Tails = 1
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

class Cholesteryl_Ester(GL.OtherLipid):

  No_Tails = 1
  adducts = {  # adduct:{spectra}

  "[M+NH4]+":{ #
    GL.MA      :10, # Molecular ion
    GL.C27H45 :100, # Present in lipidblast
    GL.C13H19  :20, # Present in lipidblast
    GL.C12H17  :20, # Present in lipidblast
    GL.C11H15  :20, # Present in lipidblast
    GL.C10H15  :20}  # Present in lipidblast
    }
    
  def __init__(self, sn1):
    body = GL.Other(name='Cholesteryl', mass=386.354866092, chnops={'C':27, 'H':46, 'O':1},
    smiles='CC(C)CCCC(C)C1CCC2C3CC=C4CC(CCC4(C)C3CCC12C)OX') # smiles has an 'X' in it to
    super().__init__(Cholesteryl_Ester.adducts, body, sn1) # indicate where tail is added

class AcylCarnitine(GL.OtherLipid):

  No_Tails = 1
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
    body = GL.Other(name='AcylCarnitine', mass=161.105193, chnops={'C':7, 'H':15, 'N':1, 'O':3},
    smiles='C[N+](C)(C)CC(OX)CC(=O)[O-]')
    super().__init__(AcylCarnitine.adducts, body, sn1)
