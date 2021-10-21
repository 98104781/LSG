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
    super().__init__(MAG.adducts, sn1)

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
    super().__init__(DAG.adducts, sn2, sn1)

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
    GL.MA       :10,
    GL.MH_s_FAk  :5,
    GL.FAH      :10, 
    152.995833 :100,
    96.969618    :5,
    78.959053   :20}}

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    sn3 = ['Head', 97.976895]
    super().__init__(lyPA.adducts, sn3, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PA(GL.Glycerolipid):

  #####  Is it just me or is there a lack of fragmentation studies for PA / LyPA...
  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033 Thank you Mr. Hsu and Mr. Turk !!
  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Looks Good
    GL.MA          :15,
    GL.MH_s_FA     :30,
    GL.MH_s_FAk    :15,
    GL.FAH        :100,
    171.006398      :5, 
    152.995833     :30,
    96.969618       :5,
    78.959053       :5},

  "[M+H]+":{ # Needs Validation
    GL.MA          :90,
    GL.HG_NL_H2O_A:100,
    GL.HG_FA_NL_A  :10,
    GL.FAkH        :10}}

  # sn3 = headgroup
  def __init__(self, sn2, sn1):
    sn3 = ['Head', 97.976895]
    super().__init__(PA.adducts, sn3, sn2, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPC(GL.Glycerolipid):

  #####  Requires reference!
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  "[M+H]+":{ # Needs Validation
    GL.MA       :5,
    GL.MH_s_FA  :5,
    GL.MH_s_FAk:10,
    GL.FAkH     :1,
    184.07332 :100}}#,

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
    sn3 = ['Head', 183.066044]
    # headgroup mass has -H to maintain neutral charge
    super().__init__(lyPC.adducts, sn3, sn1)

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
    184.07332     :100}}#,

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
    sn3 = ['Head', 183.066044]
    # headgroup mass has -H to maintain neutral charge
    super().__init__(PC.adducts, sn3, sn2, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPE(GL.Glycerolipid):

  #####  Requires reference! Needs work!
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Needs Validation
    GL.MA       :15,
    GL.MH_s_FA  :15,
    GL.MH_s_FAk  :5,
    GL.FAH     :100, 
    196.038032   :5,
    152.995833  :10,
    140.011817   :3,
    96.969618    :5,
    78.959053    :5},
  
  "[M+H]+":{ # Needs Validation
    GL.MA          :10,
    GL.HG_NL_H2O_A:100,
    GL.HG_FA_NL_A  :10,
    GL.FAkH        :10}}

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    sn3 = ['Head', 141.019094]
    super().__init__(lyPE.adducts, sn3, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PE(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  ##### "[M+H]+", "[M-H]-" https://doi.org/10.1039/C5AY00776C
  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Looks Good
    GL.MA          :15,
    GL.MH_s_FA      :2,
    GL.MH_s_FAk     :5,
    GL.FAH        :100, 
    196.038032      :5,
    152.995833     :20,
    140.011817      :3,
    96.969618       :5,
    78.959053       :5},
  
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
    sn3 = ['Head', 141.019094]
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
    sn3 = ['Head', 172.013674]
    super().__init__(lyPG.adducts, sn3, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PG(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  No_Tails = 2
  adducts = {  # adduct:{spectra}

  "[M-H]-":{ # Looks Good
    GL.MA          :15,
    GL.MH_s_FA      :2,
    GL.MH_s_FAk     :5,
    GL.FAH        :100,
    227.032612      :2,
    171.006398      :2,
    152.995833     :10,
    96.969618       :5,
    78.959053       :5},

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
    sn3 = ['Head', 172.013674]
    super().__init__(PG.adducts, sn3, sn2, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class lyPI(GL.Glycerolipid):

  #####  Requires reference!
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  #"[M+Na-2H]-":{
  #  315.048656   :50, # Not confident what the mz 192.988  is,
  #  241.011877   :10, # perhaps as a [Phosphoglycerol+Na-2H]-.
  #  223.001312   :10, # Appears in "[M+Na-2H]-" spectra.
  #  192.988      :80, # Additional mz 355.04  not included.
  #  152.995833   :80},

  "[M-H]-":{ # Needs Validation
    GL.MA         :2,
    GL.HG_NL_H2O_B:1,
    GL.FAH      :100, 
    315.048656   :10, 
    241.011877   :45,  # mz 233.001 can be a major -
    223.001312   :10,  # or minor fragment depending -
    152.995833   :45}} # on if its a sn1 or sn2 lyso lipid.

  # Should Lyso GPLs have [M+H-H2O]+ ?

  # sn3 = headgroup
  def __init__(self, sn1):
    sn3 = ['Head', 260.029718]
    super().__init__(lyPI.adducts, sn3, sn1)

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
  #  315.048656  :40,
  #  297.038092  :10, # Not confident what the mz 192.988  is, 
  #  241.011877  :10, # perhaps as a [Phosphoglycerol+Na-2H]-.
  #  223.001312  :40, # Appears in "[M+Na-2H]-" spectra.
  #  192.988     :50, # Additional mz 355.04  not included.
  #  152.995833  :60}, 

  "[M-H]-":{ # Looks Good
    GL.MA        :2,
    GL.MH_s_FA   :2,
    GL.MH_s_FAk  :1,
    GL.HG_FA_NL_B:1,
    GL.FAH     :100, 
    315.048656   :5, 
    297.038092   :5,
    259.022442   :5, 
    241.011877  :40,
    223.001312  :15, 
    152.995833  :20}}

  # sn3 = headgroup
  def __init__(self, sn2, sn1):
    sn3 = ['Head', 260.029718]
    super().__init__(PI.adducts, sn3, sn2, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PIP(GL.Glycerolipid):

  #####  Requires reference!
  No_Tails = 2
  adducts = {  # adduct:{spectra}

  #"[M+Na-2H]-":{
  #  GL.MA          :2,
  #  GL.MH_s_PO3    :2,
  #  GL.MH_s_FA     :1,
  #  GL.MH_s_FA_H2O :1,
  #  GL.MH_s_FAk_PO3:1,
  #  GL.MH_s_FA_PO3 :1,
  #  GL.HG_FA_NL_B  :5,
  #  GL.FAH        :50,
  #  342.960152   :100, 
  #  324.949587    :10,
  #  259.022442     :5, 
  #  241.011877    :50,
  #  223.001312    :20, 
  #  152.995833    :10},

  "[M-H]-":{ # Looks Good
    GL.MH          :2,
    GL.MH_s_PO3    :2,
    GL.MH_s_FA     :4,
    GL.MH_s_FA_H2O :1,
    GL.MH_s_FAk_PO3:4,
    GL.MH_s_FA_PO3 :4,
    GL.HG_FA_NL_B  :5,
    GL.FAH       :100,
    320.978207    :10, 
    302.967642    :10,
    259.022442     :5, 
    241.011877    :60,
    223.001312    :30, 
    152.995833    :30},

  "[M-2H]2-":{ # Looks Good
    GL.MA          :5,
    GL.FAH       :100,
    320.978207     :2, 
    302.967642     :1,
    241.011877    :40,
    223.001312     :1, 
    152.995833    :10}}

  # sn3 = headgroup
  def __init__(self, sn2, sn1):
    sn3 = ['Head', 339.996048]
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
    sn3 = ['Head', 185.008923]
    super().__init__(lyPS.adducts, sn3, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PS(GL.Glycerolipid):

  ##### "[M+H]+" https://doi.org/10.1016/j.jchromb.2009.02.033
  ##### "[M+H]+" https://doi.org/10.1039/C5AY00776C 
  No_Tails = 2
  adducts = {  # adduct:{spectra}
  
  "[M-H]-":{ # Looks Good
    GL.MA          :15,
    GL.HG_NL_B    :100,
    GL.HG_FA_NL_B  :50,
    GL.HG_FAk_NL_B :25,
    GL.FAH         :50,
    152.995833     :30,
    96.969618       :5,
    78.959053       :5},

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
    sn3 = ['Head', 185.008923]
    super().__init__(PS.adducts, sn3, sn2, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #