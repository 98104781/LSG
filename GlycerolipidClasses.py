from itertools import combinations_with_replacement as cwr
import GenerateLipids as GL

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

'''
class LipidClass(GL.Glycerolipid):

  List of adducts to generate spectra for = {

  "Adduct":{
    Fragment:Intensity,
    Fragment:Intensity,
    Fragment:Intensity},

    sn3 = ['Head', Headgroup mass]
'''


class lyPI(GL.Glycerolipid):

  created = []
  adducts = {  # adduct:{spectra}

  "[M+Na-2H]-":{
    315.048656:50, # Not confident what the mz 192.988  is,
    241.011877:10, # perhaps as a [Phosphoglycerol+Na-2H]-.
    223.001312:10, # Appears in "[M+Na-2H]-" spectra.
    192.988   :80, # Additional mz 355.04  not included.
    152.995833:80},

  "[M-H]-":{
    GL.MA       :2,
    GL.HG_NL_H2O:1,
    GL.FAH    :100, 
    315.048656: 10, 
    241.011877 :45,  # mz 233.001 can be a major -
    223.001312: 10,  # or minor fragment depending -
    152.995833 :45}} # on if its a sn1 or sn2 lyso lipid.

  # sn3 = headgroup
  def __init__(self, sn2, sn1):
    sn3 = ['Head', 260.029718]
    if sn1 in lyPI.created: pass
    else:
      super().__init__(lyPI.adducts, sn3, sn1)
      lyPI.created.append(sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PI(GL.Glycerolipid):

  adducts = {  # adduct:{spectra}

  "[M+Na-2H]-":{
    GL.MH      :2,
    GL.MH_s_FA :4,
    GL.HG_FA_NL:1,
    GL.FAH   :100,
    315.048656:40,
    297.038092:10, # Not confident what the mz 192.988  is, 
    241.011877:10, # perhaps as a [Phosphoglycerol+Na-2H]-.
    223.001312:40, # Appears in "[M+Na-2H]-" spectra.
    192.988   :50, # Additional mz 355.04  not included.
    152.995833:60}, 

  "[M-H]-":{
    GL.MA      :2,
    GL.MH_s_FA :2,
    GL.MH_s_FAk:1,
    GL.HG_FA_NL:1,
    GL.FAH   :100, 
    315.048656 :5, 
    297.038092 :5,
    259.022442 :5, 
    241.011877:40,
    223.001312:15, 
    152.995833:20}}

  # sn3 = headgroup
  def __init__(self, sn2, sn1):
    sn3 = ['Head', 260.029718]
    super().__init__(PI.adducts, sn3, sn2, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PIP(GL.Glycerolipid):

  adducts = {  # adduct:{spectra}

  "[M+Na-2H]-":{
    GL.MA          :2,
    GL.MH_s_PO3    :2,
    GL.MH_s_FA     :1,
    GL.MH_s_FA_H2O :1,
    GL.MH_s_FAk_PO3:1,
    GL.MH_s_FA_PO3 :1,
    GL.HG_FA_NL    :5,
    GL.FAH        :50,
    342.960152   :100, 
    324.949587    :10,
    259.022442     :5, 
    241.011877    :50,
    223.001312    :20, 
    152.995833    :10},

  "[M-H]-":{
    GL.MH          :2,
    GL.MH_s_PO3    :2,
    GL.MH_s_FA     :4,
    GL.MH_s_FA_H2O :1,
    GL.MH_s_FAk_PO3:4,
    GL.MH_s_FA_PO3 :4,
    GL.HG_FA_NL    :5,
    GL.FAH       :100,
    320.978207    :10, 
    302.967642    :10,
    259.022442     :5, 
    241.011877    :60,
    223.001312    :30, 
    152.995833    :30},

  "[M-2H]2-":{
    GL.MA         :5,
    GL.FAH      :100,
    320.978207    :2, 
    302.967642    :1,
    241.011877    :40,
    223.001312    :1, 
    152.995833    :10}}

  # sn3 = headgroup
  def __init__(self, sn2, sn1):
    sn3 = ['Head', 339.996048]
    super().__init__(PIP.adducts, sn3, sn2, sn1)

# ~ # ~ # ~ # ~ # ~ # ~ #

class PI2P(GL.Glycerolipid):

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

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

def Generate_Lipids():
  classes = GL.Glycerolipid.__subclasses__()
  for cls in classes:
    for sn1, sn2 in cwr(GL.tails, 2):
      cls(sn2, sn1)
  return GL.Glycerolipid.instances

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #