from re import A
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

  #####  Requires reference! # Needs work ?
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

  "[M+NH4]+":{ # Needs Validation    
    GL.MA       :25,
    GL.MH       :10,
    GL.MH_s_H2O:100,
    GL.MH_s_FA :100,
    GL.FAkH      :1}
  }

  def __init__(self, sn1, sn2):
    super().__init__(DAG.adducts, sn1=sn1, sn2=sn2)

# ~ # ~ # ~ # ~ # ~ # ~ #

class TAG(GL.Glycerolipid):

  #####  Requires reference! # Needs work ?
  No_Tails = 3
  adducts = {  # adduct:{spectra}

  "[M+Na]+":{ # Needs Validation   
    GL.MA       :25,
    GL.MA_s_FA  :10,
    GL.MH_s_FA :100,
    GL.FAkA      :5,
    GL.FAkH     :10},

  "[M+NH4]+":{ # Needs Validation    
    GL.MA       :25,
    GL.MH_s_FA :100,
    GL.FAkH     :10}
  }

  def __init__(self, sn1, sn2, sn3):
    super().__init__(TAG.adducts, sn1, sn2, sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

class MGDG(GL.Glycerolipid): # Monogalactosyl diacylglycerol

  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M+Na]+":{ # https://doi.org/10.1111/j.1440-1835.2010.00582.x   
    GL.MA          :5, # This spectra agrees with lipidblast and the attached
    GL.MA_s_FA   :100},# reference. Ion ratios are off, but depend on isomerism.

  "[M+NH4]+":{ # https://doi.org/10.1002/pld3.183 Figure S4     
    GL.MA          :5, # LipidBlast includes fragment 'MA_s_FA' but
    GL.HG_NL_C     :5, # the fragment is not present in the reference.
    GL.HG_NL_H2O_C:60, # Included here at 0 intensity so it can be modified
    GL.MA_s_FA     :0, # in the program if needed, otherwise it will be
    GL.HG_FA_NL_C:100} # ignored on spectra generation.
    }

  def __init__(self, sn1, sn2):
    sn3 = GL.sn(mass=180.063388, type='Headgroup', chnops={'C':6, 'H':12, 'O':6})
    super().__init__(MGDG.adducts, sn1, sn2, sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

class SQDG(GL.Glycerolipid): # Sulphoquinovosyl diacylglycerol

  #####  Requires reference! # Needs work ?
  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M-H]-":{  
    GL.MA         :40,
    GL.MA_s_FA    :10,
    GL.FAH        :10,
    GL.C6H9O7S   :100,
    GL.O3S        :20}}

  #"[M+NH4]+":{    
  #  GL.MA          :5
  #  }

  def __init__(self, sn1, sn2):
    sn3 = GL.sn(mass=244.025287, type='Headgroup', chnops={'C':6, 'H':12, 'O':8, 'S':1})
    super().__init__(SQDG.adducts, sn1, sn2, sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

class DGDG(GL.Glycerolipid): # Digalactosyl diacylglycerol

  #####  Requires reference! # Needs work ?
  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M+Na]+":{ # https://doi.org/10.1111/j.1440-1835.2010.00582.x    
    GL.MA          :5,
    GL.MA_s_Gal_H2O:5,
    GL.MA_s_FA   :100,
    GL.MA_s_FA_Gal:50,
    },

  "[M+NH4]+":{ # https://doi.org/10.1002/pld3.183 Figure S4   
    GL.MA          :5,
    GL.MH_s_Gal_H2O:5,
    GL.HG_NL_C    :70,
    GL.HG_NL_H2O_C:70,
    GL.MA_s_FA     :0,
    GL.HG_FA_NL_C:100}
    }

  def __init__(self, sn1, sn2):
    sn3 = GL.sn(mass=342.116212, type='Headgroup', chnops={'C':12, 'H':22, 'O':11})
    super().__init__(DGDG.adducts, sn1, sn2, sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

class DGTS(GL.Glycerolipid): # N-trimethylhomoserine diacylglycerol

  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M+H]+":{ #  https://doi.org/10.1002/rcm.7847
    GL.MA        :100, # Based off a single spectra...
    GL.MA_s_FAk   :70,
    GL.MA_s_FA    :40,
    GL.MA_s_allFAk:20,
    GL.C7H14N1O2   :5}
  }

  #"[M+NH4]+":{    
  #  GL.MA          :5}
  #  }

  def __init__(self, sn1, sn2):
    sn3 = GL.sn(mass=161.105193, type='Headgroup', chnops={'C':7, 'H':15, 'N':1, 'O':3})
    # headgroup mass has -H to maintain neutral charge
    super().__init__(DGTS.adducts, sn1, sn2, sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PA(GL.Glycerolipid):

  #####  Is it just me or is there a lack of fragmentation studies for PA / LyPA...
  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033 Thank you Mr. Hsu and Mr. Turk !!
  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Looks Good
    GL.MA      :15,
    GL.MH_s_FAk:15,
    GL.MH_s_FA :30,
    GL.FAH    :100,
    GL.C3H8O6P  :5, 
    GL.C3H6O5P :30,
    GL.H2O4P    :5,
    GL.O3P      :5},

  "[M+H]+":{ # Needs Validation
    GL.MA          :90,
    GL.HG_NL_H2O_A:100,
    GL.HG_FA_NL_A  :10,
    GL.FAkH        :10}
    }

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    sn3 = GL.sn(mass=97.976895, type='Headgroup', chnops={'H':3, 'O':4, 'P':1})
    super().__init__(PA.adducts, sn1, sn2, sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PC(GL.Glycerolipid):

  ##### "[M+H]+", "[M+Na]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  ##### "[M+H]+", "[M+Na]+"  10.1016/S1044-0305(03)00064-3
  ##### "[M+H]+" https://doi.org/10.1039/C5AY00776C 
  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M+H]+":{ # Looks Good
    GL.MA          :10,
    GL.MH_s_FAk     :1,
    GL.MH_s_FA      :1,
    GL.C5H15NO4P  :100}
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

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    sn3 = GL.sn(mass=183.066044, type='Headgroup', chnops={'C':5, 'H':14, 'N':1, 'O':4, 'P':1})
    # headgroup mass has -H to maintain neutral charge
    super().__init__(PC.adducts, sn1, sn2, sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PE(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  ##### "[M+H]+", "[M-H]-" https://doi.org/10.1039/C5AY00776C
  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Looks Good
    GL.MA       :15,
    GL.MH_s_FAk  :5,
    GL.MH_s_FA   :2,
    GL.FAH     :100, 
    GL.C5H11NO5P :5,
    GL.C3H6O5P  :20,
    GL.C2H7NO4P  :3,
    GL.H2O4P     :5,
    GL.O3P       :5},
  
  "[M+H]+":{ # Looks Good
    GL.MA          :10,
    GL.HG_NL_H2O_A:100,
    GL.HG_FA_NL_A  :10,
    GL.FAkH        :10}
    }#,
    
  #"[M+Na]+":{ # Assuming [M+Na] fragments like [M+Li]
  #  GL.MA          :10,
  #  GL.MA_s_AZD     :5,
  #  GL.HG_NL_H2O_A :50,
  #  GL.HG_NL_H2O_C:100,
  #  GL.MA_s_FA_AZD :10,
  #  GL.MH_s_FA     :10,
  #  GL.MA_s_FA     :10,
  #  GL.FAkH        :5}
  # }

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    sn3 = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1})
    super().__init__(PE.adducts, sn1, sn2, sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PG(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Looks Good
    GL.MA       :15,
    GL.MH_s_FAk  :5,
    GL.MH_s_FA   :2,
    GL.FAH     :100,
    GL.C6H12O7P  :2,
    GL.C3H8O6P   :2,
    GL.C3H6O5P  :10,
    GL.H2O4P     :5,
    GL.O3P       :5},

  "[M+H]+":{ # Looks Good
    GL.MA          :10,
    GL.HG_NL_H2O_A:100,
    GL.HG_FA_NL_A  :10,
    GL.FAkH        :10}
    }#,
    
  #"[M+Na]+":{
  #  GL.MA          :10,
  #  GL.HG_NL_H2O_A:100,
  #  GL.HG_FA_NL_A  :10,
  #  GL.FAkH        :10,
  #  GL.HGA         :80}
  # }

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    sn3 = GL.sn(mass=172.013674, type='Headgroup', chnops={'C':3, 'H':9, 'O':6, 'P':1})
    super().__init__(PG.adducts, sn1, sn2, sn3)

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

  "[M-H]-":{ # Looks Good
    GL.MA        :2,
    GL.MH_s_FAk  :1,
    GL.MH_s_FA   :2,
    GL.HG_FA_NL_B:1,
    GL.FAH     :100, 
    GL.C9H16O10P :5, 
    GL.C9H14O9P  :5,
    GL.C6H12O9P  :5, 
    GL.C6H10O8P :40,
    GL.C6H8O7P  :15, 
    GL.C3H6O5P  :20},

  "[M+NH4]+":{ # https://doi.org/10.1002/pld3.183 Figure S4   
    GL.MA           :5,
    GL.MH           :5,
    GL.HG_NL_H2O_C:100}
    }

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    sn3 = GL.sn(mass=260.029718, type='Headgroup', chnops={'C':6, 'H':13, 'O':9, 'P':1})
    super().__init__(PI.adducts, sn1, sn2, sn3)

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
    GL.MH          :1,
    GL.MH_s_H2O    :2,
    GL.MH_s_PO3    :2,
    GL.HG_NL_2B    :2,
    GL.MH_s_FA     :4,
    GL.MH_s_FA_H2O :1,
    GL.MH_s_FAk_PO3:4,
    GL.MH_s_FA_PO3 :4,
    GL.HG_FA_NL_B  :5,
    GL.FAH       :100,
    GL.C6H11O11P2 :10, 
    GL.C6H9O10P2  :15,
    GL.C6H12O9P    :5, 
    GL.C6H10O8P   :90,
    GL.C6H8O7P    :90, 
    GL.C3H6O5P    :40}
    }#,

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
    sn3 = GL.sn(mass=339.996048, type='Headgroup', chnops={'C':6, 'H':14, 'O':12, 'P':2})
    super().__init__(PIP.adducts, sn1, sn2, sn3)

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
    sn3 = ['Head', 419.962378]
    sn3 = GL.sn(mass=419.962378, type='Headgroup', chnops={'C':6, 'H':15, 'O':15, 'P':3})
    super().__init__(PI2P.adducts, sn1, sn2, sn3)
'''
# ~ # ~ # ~ # ~ # ~ # ~ #

class PS(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  ##### "[M+H]+" https://doi.org/10.1039/C5AY00776C 
  No_Tails = 2
  adducts = {  # adduct:{spectra}
  
  "[M-H]-":{ # Looks Good
    GL.MA         :15,
    GL.HG_NL_B   :100,
    GL.HG_FAk_NL_B:25,
    GL.HG_FA_NL_B :50,
    GL.FAH        :50,
    GL.C3H6O5P    :30,
    GL.H2O4P       :5,
    GL.O3P         :5},

  "[M+H]+":{ # Looks Good
    GL.MA          :10,
    GL.HG_NL_H2O_A:100,
    GL.HG_FA_NL_A  :10,
    GL.FAkH        :10}
    }#,
    
  #"[M+Na]+":{ # Assuming [M+Na] fragments like [M+Li]
  #  GL.MA          :10,
  #  GL.HG_NL_H2O_A:100,
  #  GL.HG_FA_NL_A  :10,
  #  GL.FAkH        :10,
  #  GL.HGA         :80}
  # }

  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    sn3 = GL.sn(mass=185.008923, type='Headgroup', chnops={'C':3, 'H':8, 'N':1, 'O':6, 'P':1})
    super().__init__(PS.adducts, sn1, sn2, sn3)

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
    GL.HO6P2   :50}
    }#,

  #"[M-2H]2-":{
  #  GL.MA     :10
  # },

  #"[M-3H]3-":{
  #  GL.MA     :10}
  # }#,
    
  # sn3 = headgroup
  def __init__(self, sn1, sn2):
    sn3 = GL.sn(mass=177.943224, type='Headgroup', chnops={'H':4, 'O':7, 'P':2})
    super().__init__(PPA.adducts, sn1, sn2, sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPA(GL.Glycerolipid):

  ##### [M-H]- = https://doi.org/10.1002/lipd.12172, Not a fragmentation study!
  ##### [M-H]- = https://doi.org/10.1016/j.jchromb.2010.03.030, Neither...
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Needs Validation
    GL.MA      :10,
    GL.FAH     :10, 
    GL.MH_s_FAk :5,
    GL.C3H6O5P:100,
    GL.H2O4P    :2,
    GL.O3P     :10}
    }

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    sn3 = GL.sn(mass=97.976895, type='Headgroup', chnops={'H':3, 'O':4, 'P':1})
    super().__init__(lyPA.adducts, sn1=sn1, sn3=sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPC(GL.Glycerolipid):

  #####  Requires reference!
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  "[M+H]+":{ # Needs Validation
    GL.MA         :5,
    GL.MH_s_FAk  :10,
    GL.MH_s_FA    :5,
    GL.FAkH       :1,
    GL.C5H15NO4P:100}
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
    sn3 = GL.sn(mass=183.066044, type='Headgroup', chnops={'C':5, 'H':14, 'N':1, 'O':4, 'P':1})
    # headgroup mass has -H to maintain neutral charge
    super().__init__(lyPC.adducts, sn1=sn1, sn3=sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPE(GL.Glycerolipid):

  #####  Requires reference! Needs work!
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Needs Validation
    GL.MA      :15,
    GL.MH_s_FAk :5,
    GL.MH_s_FA :15,
    GL.FAH    :100, 
    GL.C5H11NO5P:5,
    GL.C3H6O5P :10,
    GL.C2H7NO4P :3,
    GL.H2O4P    :5,
    GL.O3P      :5},
  
  "[M+H]+":{ # Needs Validation
    GL.MA          :10,
    GL.MA_s_H2O    :20,
    GL.HG_NL_H2O_A:100,
    GL.FAkH        :10}
    }

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    sn3 = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1})
    super().__init__(lyPE.adducts, sn1=sn1, sn3=sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPG(GL.Glycerolipid):

  #####  Requires reference! Needs work!
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Needs Validation
    GL.MH         :100,
    GL.HG_NL_B      :1,
    GL.HG_NL_H2O_B  :1,
    GL.FAH        :100,
    GL.C3H6O5P     :50,
    GL.H2O4P        :1,
    GL.O3P          :1},

  "[M+H]+":{ # Needs Validation
    GL.MA          :10,
    GL.HG_NL_H2O_A:100,
    GL.HG_FA_NL_A  :10,
    GL.FAkH        :10}
    }

  # sn3 = headgroup
  def __init__(self, sn1):
    sn3 = GL.sn(mass=172.013674, type='Headgroup', chnops={'C':3, 'H':9, 'O':6, 'P':1})
    super().__init__(lyPG.adducts, sn1=sn1, sn3=sn3)

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
    GL.MA         :2,
    GL.HG_NL_B    :1,
    GL.HG_NL_H2O_B:1,
    GL.FAH      :100, 
    GL.C9H16O10P :10, 
    GL.C6H10O8P  :45,  # mz 233.001 can be a major -
    GL.C6H8O7P   :10,  # or minor fragment depending -
    GL.C3H6O5P   :45}  # on if its a sn1 or sn2 lyso lipid.
    } 

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    sn3 = GL.sn(mass=260.029718, type='Headgroup', chnops={'C':6, 'H':13,'O':9, 'P':1})
    super().__init__(lyPI.adducts, sn1=sn1, sn3=sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPS(GL.Glycerolipid):

  #####  Requires reference! Needs work!
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  "[M+H]+":{ # Needs Validation
    GL.MA          :10,
    GL.HG_NL_H2O_A:100,
    GL.HG_FA_NL_A  :10,
    GL.FAkH        :10}
    }#,

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    sn3 = GL.sn(mass=185.008923, type='Headgroup', chnops={'C':3, 'H':8, 'N':1, 'O':6, 'P':1})
    super().__init__(lyPS.adducts, sn1=sn1, sn3=sn3)

# ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class Cer_d18(GL.Sphingolipid):

  #####  Requires reference! Needs work!
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # https://doi.org/10.1006/abio.2001.5536, https://doi.org/10.1016/S1044-0305(02)00358-6, https://doi.org/10.1016/j.biochi.2016.07.012
    GL.MA        :30,
    GL.MA_s_H2O   :2,
    GL.MA_s_Al   :20,
    GL.MA_s_MeOH :10,
    GL.MA_s_CH2O :10,
    GL.FA_C2H5NO :20,
    GL.FA_C2H5N :100,
    GL.FAH       :20,
    GL.FA_N      :20,
    GL.B_LC      :30,
    GL.FAkH      :20,
    GL.B_s_C2H6O2:10
    },
    
  "[M+H]+":{ # https://doi.org/10.1002/bmc.4790
    GL.MA          :20,
    GL.MA_s_H2O    :15,
    GL.MA_s_FA     :10,
    GL.FA_C2H3N    :10,
    GL.FA_N        :10,
    GL.B_LC       :100
    },

  "[M+H-H2O]+":{ #
    GL.MA         :100
    }}#,

  # sn3 = headgroup
  def __init__(self, sn1):
    base = GL.base(18, 1, type='Base', oh=2)
    super().__init__(Cer_d18.adducts, base=base, sn1=sn1)