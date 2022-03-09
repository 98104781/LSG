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

class sn1_lysoPI(GL.Glycerolipid):

  #####  Requires reference!
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  #"[M+Na-2H]-":{
  #  GL.C9H16O10P :50,
  #  GL.C6H10O8P  :10,
  #  GL.C6H8O7P   :10,
  #  GL.C3H7NaO6P :80,
  #  GL.C3H6O5P   :80},

  "[M-H]-":{
    GL.MA         :2,
    GL.MH_PO4_s_HG_H2O:1,
    GL.sn1      :100, 
    GL.C9H16O10P :15, 
    GL.C6H10O8P  :40,
    GL.C6H8O7P    :4,
    GL.C3H6O5P   :50}}

  # sn3 = headgroup
  def __init__(self, sn1):
    sn3 = GL.sn(mass=260.029718, type='Headgroup', chnops={'C':6, 'H':13,'O':9, 'P':1})
    super().__init__(sn1_lysoPI.adducts, sn3=sn3, sn1=sn1)
    self.name = f"{self.lipid_class} {'/'.join(snx.name for snx in self.tails if snx.name != 'Headgroup')}"

class sn2_lysoPI(GL.Glycerolipid):

  #####  Requires reference!
  No_Tails = 1
  adducts = {  # adduct:{spectra}

  #"[M+Na-2H]-":{
  #  GL.C9H16O10P :50,
  #  GL.C6H10O8P  :10,
  #  GL.C6H8O7P   :10,
  #  GL.C3H7NaO6P :80,
  #  GL.C3H6O5P   :80},

  "[M-H]-":{
    GL.MA         :2,
    GL.sn2      :100, 
    GL.C9H16O10P :10, 
    GL.C6H10O8P  :30,
    GL.C6H8O7P   :15,
    GL.C3H6O5P   :80}}

  # sn3 = headgroup
  def __init__(self, sn2):
    sn3 = GL.sn(mass=260.029718, type='Headgroup', chnops={'C':6, 'H':13,'O':9, 'P':1})
    super().__init__(sn2_lysoPI.adducts, sn3=sn3, sn2=sn2)
    self.name = f"{self.lipid_class} {'/'.join(snx.name for snx in self.tails if snx.name != 'Headgroup')}"

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
  #  GL.C3H7NaO6P:50,
  #  GL.C3H6O5P  :60}, 

  "[M-H]-":{
    GL.MA          :2,
    GL.MH_s_sn1    :1,
    GL.MH_s_sn2k   :2,
    GL.MH_s_sn2    :5,
    GL.HG_sn1_NL_B :1,
    GL.HG_sn2k_NL_B:1,
    GL.HG_sn2_NL_B :5,
    GL.sn1        :95,
    GL.sn2       :100, 
    GL.C9H16O10P   :5, 
    GL.C9H14O9P    :5,
    GL.C6H12O9P    :5, 
    GL.C6H10O8P   :40,
    GL.C6H8O7P    :15, 
    GL.C3H6O5P    :20}}

  # sn3 = headgroup
  def __init__(self, sn2, sn1):
    sn3 = GL.sn(mass=260.029718, type='Headgroup', chnops={'C':6, 'H':13,'O':9, 'P':1})
    super().__init__(PI.adducts, sn3=sn3, sn2=sn2, sn1=sn1)
    self.name = f"{self.lipid_class} {'/'.join(snx.name for snx in self.tails if snx.name != 'Headgroup')}"
    

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

  "[M-H]-":{
    GL.MH          :1,
    GL.MH_s_H2O    :2,
    GL.MH_s_PO3    :5,
    GL.MH_P2O6_s_HG    :5,
    GL.MH_s_sn2_H2O:1,
    GL.MH_s_sn2_PO3:4,
    GL.HG_sn2_NL_B :5,
    GL.sn1        :95,
    GL.sn2       :100,
    GL.C6H11O11P2 :10, 
    GL.C6H9O10P2  :15,
    GL.C6H12O9P    :5, 
    GL.C6H10O8P   :90,
    GL.C6H8O7P    :90, 
    GL.C3H6O5P    :40}}#,

  #"[M-2H]2-":{
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
    super().__init__(PIP.adducts, sn3=sn3, sn2=sn2, sn1=sn1)
    self.name = f"{self.lipid_class} {'/'.join(snx.name for snx in self.tails if snx.name != 'Headgroup')}"

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
    super().__init__(PI2P.adducts, sn3=sn3, sn2=sn2, sn1=sn1)
'''

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #