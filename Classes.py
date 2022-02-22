import GenerateLipids as GL

''' EXAMPLE EXAMPLE EXAMPLE

class LipidClass(GL.Glycerolipid):

  List of adducts to generate spectra for = {

  "Adduct":{
    Fragment:Intensity,
    Fragment:Intensity,
    Fragment:Intensity},

  def __init__(self, sn3, sn2, sn1):
    sn3 = ['Head', Headgroup mass]
    super().__init__(LipidClass.adducts, sn3, sn2, sn1)

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
    GL.FAkH     :50,
    GL.FAkA     :10},

  "[M+NH4]+":{ # Needs Validation    
    GL.MA      :100,
    GL.MH      :100,
    GL.MA_s_H2O :50,
    GL.FAkH     :50,
    GL.FAkA     :10},
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
    GL.FAkH      :1,
    GL.MH_s_FA  :25},

  "[M+H-H2O]+":{ # Needs Validation
    GL.MA      :100,
    GL.FAkH      :1,
    GL.MH_s_FA :100},

  "[M+Na]+":{ # Needs Validation    
    GL.MA       :50,
    GL.MA_s_H2O:100,
    GL.FAkH      :1,
    GL.MH_s_FA  :50,
    GL.FAkA     :10,
    GL.MA_s_FA  :20},

  "[M+NH4]+":{ # Needs Validation    
    GL.MA       :25,
    GL.MH       :10,
    GL.MA_s_H2O:100,
    GL.FAkH      :1,
    GL.MH_s_FA  :25}
  }

  def __init__(self, sn2, sn1):
    super().__init__(DAG.adducts, sn2=sn2, sn1=sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class TAG(GL.Glycerolipid):

  #####  Requires reference! # Needs work ?
  No_Tails = 3
  adducts = {  # adduct:{spectra}

  "[M+Na]+":{ # Needs Validation   
    GL.MA       :25,
    GL.MH_s_FA :100,
    GL.MA_s_FA  :10,
    GL.FAkH     :10,
    GL.FAkA      :5},

  "[M+NH4]+":{ # Needs Validation    
    GL.MA       :25,
    GL.MH_s_FA :100,
    GL.FAkH     :10}
  }

  def __init__(self, sn3, sn2, sn1):
    super().__init__(TAG.adducts, sn3, sn2, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPA(GL.Glycerolipid):

  ##### [M-H]- = https://doi.org/10.1002/lipd.12172, Not a fragmentation study!
  ##### [M-H]- = https://doi.org/10.1016/j.jchromb.2010.03.030, Neither...
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Needs Validation
    GL.MA      :10,
    GL.MH_s_FAk :5,
    GL.FAH     :10, 
    GL.C3H6O5P:100,
    GL.H2O4P    :5,
    GL.O3P     :20}}

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    sn3 = GL.sn(mass=97.976895, type='Headgroup', chnops={'H':3, 'O':4, 'P':1})
    super().__init__(lyPA.adducts, sn3, sn1=sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PA(GL.Glycerolipid):

  #####  Is it just me or is there a lack of fragmentation studies for PA / LyPA...
  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033 Thank you Mr. Hsu and Mr. Turk !!
  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Looks Good
    GL.MA      :15,
    GL.MH_s_FA :30,
    GL.MH_s_FAk:15,
    GL.FAH    :100,
    GL.C3H8O6P  :5, 
    GL.C3H6O5P :30,
    GL.H2O4P    :5,
    GL.O3P      :5},

  "[M+H]+":{ # Needs Validation
    GL.MA          :90,
    GL.HG_NL_H2O_A:100,
    GL.HG_FA_NL_A  :10,
    GL.FAkH        :10}}

  # sn3 = headgroup
  def __init__(self, sn2, sn1):
    sn3 = GL.sn(mass=97.976895, type='Headgroup', chnops={'H':3, 'O':4, 'P':1})
    super().__init__(PA.adducts, sn3, sn2, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPC(GL.Glycerolipid):

  #####  Requires reference!
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  "[M+H]+":{ # Needs Validation
    GL.MA         :5,
    GL.MH_s_FA    :5,
    GL.MH_s_FAk  :10,
    GL.FAkH       :1,
    GL.C5H15NO4P:100}}#,

  #"[M+Na]+":{ # Assuming [M+Na] fragments like [M+Li]
  #  GL.MA          :10,
  #  GL.MA_s_TMA     :5,
  #  GL.HG_NL_H2O_A :50,
  #  GL.HG_NL_H2O_C:100,
  #  GL.MA_s_FA_TMA :10,
  #  GL.MH_s_FA     :10,
  #  GL.MA_s_FA     :10,
  #  GL.FAkH        :5}}

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    sn3 = GL.sn(mass=183.066044, type='Headgroup', chnops={'C':5, 'H':14, 'N':1, 'O':4, 'P':1})
    # headgroup mass has -H to maintain neutral charge
    super().__init__(lyPC.adducts, sn3, sn1=sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PC(GL.Glycerolipid):

  ##### "[M+H]+", "[M+Na]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  ##### "[M+H]+", "[M+Na]+"  10.1016/S1044-0305(03)00064-3
  ##### "[M+H]+" https://doi.org/10.1039/C5AY00776C 
  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M+H]+":{ # Looks Good
    GL.MA           :5,
    GL.MH_s_FA      :1,
    GL.MH_s_FAk     :5,
    GL.C5H15NO4P  :100}}#,

  #"[M+Na]+":{ # Assuming [M+Na] fragments like [M+Li]
  #  GL.MA          :10,
  #  GL.MA_s_TMA     :5,
  #  GL.HG_NL_H2O_A :50,
  #  GL.HG_NL_H2O_C:100,
  #  GL.MA_s_FA_TMA :10,
  #  GL.MH_s_FA     :10,
  #  GL.MA_s_FA     :10,
  #  GL.FAkH        :5}}

  # sn3 = headgroup
  def __init__(self, sn2, sn1):
    sn3 = GL.sn(mass=183.066044, type='Headgroup', chnops={'C':5, 'H':14, 'N':1, 'O':4, 'P':1})
    # headgroup mass has -H to maintain neutral charge
    super().__init__(PC.adducts, sn3, sn2, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPE(GL.Glycerolipid):

  #####  Requires reference! Needs work!
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Needs Validation
    GL.MA      :15,
    GL.MH_s_FA :15,
    GL.MH_s_FAk :5,
    GL.FAH    :100, 
    GL.C5H11NO5P:5,
    GL.C3H6O5P :10,
    GL.C2H7NO4P :3,
    GL.H2O4P    :5,
    GL.O3P      :5},
  
  "[M+H]+":{ # Needs Validation
    GL.MA          :10,
    GL.HG_NL_H2O_A:100,
    GL.HG_FA_NL_A  :10,
    GL.FAkH        :10}}

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    sn3 = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1})
    super().__init__(lyPE.adducts, sn3, sn1=sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PE(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  ##### "[M+H]+", "[M-H]-" https://doi.org/10.1039/C5AY00776C
  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Looks Good
    GL.MA       :15,
    GL.MH_s_FA   :2,
    GL.MH_s_FAk  :5,
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
    GL.FAkH        :10}}#,
    
  #"[M+Na]+":{ # Assuming [M+Na] fragments like [M+Li]
  #  GL.MA          :10,
  #  GL.MA_s_AZD     :5,
  #  GL.HG_NL_H2O_A :50,
  #  GL.HG_NL_H2O_C:100,
  #  GL.MA_s_FA_AZD :10,
  #  GL.MH_s_FA     :10,
  #  GL.MA_s_FA     :10,
  #  GL.FAkH        :5}}

  # sn3 = headgroup
  def __init__(self, sn2, sn1):
    sn3 = GL.sn(mass=141.019094, type='Headgroup', chnops={'C':2, 'H':8, 'N':1, 'O':4, 'P':1})
    super().__init__(PE.adducts, sn3, sn2, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPG(GL.Glycerolipid):

  #####  Requires reference! Needs work!
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  "[M+H]+":{ # Needs Validation
    GL.MA          :10,
    GL.HG_NL_H2O_A:100,
    GL.HG_FA_NL_A  :10,
    GL.FAkH        :10}}#,
  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    sn3 = GL.sn(mass=172.013674, type='Headgroup', chnops={'C':3, 'H':9, 'O':6, 'P':1})
    super().__init__(lyPG.adducts, sn3, sn1=sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PG(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Looks Good
    GL.MA       :15,
    GL.MH_s_FA   :2,
    GL.MH_s_FAk  :5,
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
    GL.FAkH        :10}}#,
    
  #"[M+Na]+":{
  #  GL.MA          :10,
  #  GL.HG_NL_H2O_A:100,
  #  GL.HG_FA_NL_A  :10,
  #  GL.FAkH        :10,
  #  GL.HGA         :80}}

  # sn3 = headgroup
  def __init__(self, sn2, sn1):
    sn3 = GL.sn(mass=172.013674, type='Headgroup', chnops={'C':3, 'H':9, 'O':6, 'P':1})
    super().__init__(PG.adducts, sn3, sn2, sn1)

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
    GL.HG_NL_H2O_B:1,
    GL.FAH      :100, 
    GL.C9H16O10P :10, 
    GL.C6H10O8P  :45,  # mz 233.001 can be a major -
    GL.C6H8O7P   :10,  # or minor fragment depending -
    GL.C3H6O5P   :45}} # on if its a sn1 or sn2 lyso lipid.

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    sn3 = GL.sn(mass=260.029718, type='Headgroup', chnops={'C':6, 'H':13,'O':9, 'P':1})
    super().__init__(lyPI.adducts, sn3, sn1=sn1)

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
    GL.MH_s_FA   :2,
    GL.MH_s_FAk  :1,
    GL.HG_FA_NL_B:1,
    GL.FAH     :100, 
    GL.C9H16O10P :5, 
    GL.C9H14O9P  :5,
    GL.C6H12O9P  :5, 
    GL.C6H10O8P :40,
    GL.C6H8O7P  :15, 
    GL.C3H6O5P  :20}}

  # sn3 = headgroup
  def __init__(self, sn2, sn1):
    sn3 = GL.sn(mass=260.029718, type='Headgroup', chnops={'C':6, 'H':13, 'O':9, 'P':1})
    super().__init__(PI.adducts, sn3, sn2, sn1)

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
    GL.C3H6O5P    :40}}#,

  #"[M-2H]2-":{ # Looks Good
  #  GL.MA          :5,
  #  GL.FAH       :100,
  #  GL.C6H11O11P2  :2, 
  #  GL.C6H9O10P2   :1,
  #  GL.C6H10O8P   :40,
  #  GL.C6H8O7P     :1, 
  #  GL.C3H6O5P    :10}}

  # sn3 = headgroup
  def __init__(self, sn2, sn1):
    sn3 = GL.sn(mass=339.996048, type='Headgroup', chnops={'C':6, 'H':14, 'O':12, 'P':2})
    super().__init__(PIP.adducts, sn3, sn2, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #
'''
class PI2P(GL.Glycerolipid):

  #####  Requires reference! Needs work!
  No_Tails = 2
  adducts = {  # adduct:{spectra} 

  "[M+Na-2H]-":{
  },

  "[M-H]-":{
  },

  "[M-2H]2-":{
  },

  "[M-3H]3-":{
  }} 

  # sn3 = headgroup
  def __init__(self, sn2, sn1):
    sn3 = ['Head', 419.962378]
    sn3 = GL.sn(mass=419.962378, type='Headgroup', chnops={'C':6, 'H':15, 'O':15, 'P':3})
    super().__init__(PI2P.adducts, sn3, sn2, sn1)
'''
# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPS(GL.Glycerolipid):

  #####  Requires reference! Needs work!
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  "[M+H]+":{ # Needs Validation
    GL.MA          :10,
    GL.HG_NL_H2O_A:100,
    GL.HG_FA_NL_A  :10,
    GL.FAkH        :10}}#,

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    sn3 = GL.sn(mass=185.008923, type='Headgroup', chnops={'C':3, 'H':8, 'N':1, 'O':6, 'P':1})
    super().__init__(lyPS.adducts, sn3, sn1=sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PS(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  ##### "[M+H]+" https://doi.org/10.1039/C5AY00776C 
  No_Tails = 2
  adducts = {  # adduct:{spectra}
  
  "[M-H]-":{ # Looks Good
    GL.MA         :15,
    GL.HG_NL_B   :100,
    GL.HG_FA_NL_B :50,
    GL.HG_FAk_NL_B:25,
    GL.FAH        :50,
    GL.C3H6O5P    :30,
    GL.H2O4P       :5,
    GL.O3P         :5},

  "[M+H]+":{ # Looks Good
    GL.MA          :10,
    GL.HG_NL_H2O_A:100,
    GL.HG_FA_NL_A  :10,
    GL.FAkH        :10}}#,
    
  #"[M+Na]+":{ # Assuming [M+Na] fragments like [M+Li]
  #  GL.MA          :10,
  #  GL.HG_NL_H2O_A:100,
  #  GL.HG_FA_NL_A  :10,
  #  GL.FAkH        :10,
  #  GL.HGA         :80}}

  # sn3 = headgroup
  def __init__(self, sn2, sn1):
    sn3 = GL.sn(mass=185.008923, type='Headgroup', chnops={'C':3, 'H':8, 'N':1, 'O':6, 'P':1})
    super().__init__(PS.adducts, sn3, sn2, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #