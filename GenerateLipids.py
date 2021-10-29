from collections import Counter

Masses = {

  #[Delta_mz,   Polarity,  z]
  "[M-H]-":     
   [-1.007276, 'Negative', 1],
  "[M-2H]2-":   
   [-2.014552, 'Negative', 2],
  "[M-3H]3-":   
   [-3.021828, 'Negative', 3],
  "[M+Cl]-":    
   [34.969401, 'Negative', 1],
  "[M+Na-2H]-":    
   [20.974668, 'Negative', 1],
  "[M+2Na-3H]-":    
   [42.956613, 'Negative', 1],
  "[M+H]+":     
   [1.007276,  'Positive', 1],
  "[M+H-H2O]+": 
   [-17.003289,'Positive', 1],
  "[M+Na]+":    
   [22.989221, 'Positive', 1],
  "[M+Li]+":    
   [7.015455, 'Positive', 1],
  "[M+NH4]+":   
   [18.033826, 'Positive', 1],
   
  #  Handy masses 
  "H":           
   1.007276,
  "H2O":         
   18.010565,
  "Glycerol":
   92.047344,
   "PO3":
   79.96633,
   "PO4":
   97.976895,
   "TMA": # Trimethylamine, fragment for PC+Na/Li
   59.073499,
   "AZD": # Aziridine, fragment for PE+Na/Li
   43.042199}


# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class sn:
  def __init__(self, c=0, d=0, m=0, mass=None, chnops={}, type=None):

    self.type = type
    
    if self.type == 'Acyl':
      self.name = f"{c}:{d}"
      #           O2 mass     + c*CH2 mass    - d*H2 mass
      self.mass = 31.98982926 + c*14.01565007 - d*2.01565007
      self.formula = {'C':c, 'H':(2*c-2*d),'O':2}

    elif self.type == 'Ether':
      self.name = f"{c}:{d};O"
      #           O mass    + c*CH2 mass    - d*H2 mass
      self.mass = 18.010565 + c*14.01565007 - d*2.01565007
      self.formula = {'C':c, 'H':(2*c-2*d),'O':1}

    elif self.type == 'Methyl': # Methyl branched Acyl kind
      self.name = f"{c}:{d};{m}M" # Perhaps exclude? Mass identical to acyl c+m
      #           O2 mass     + c*CH2 mass    - d*H2 mass
      self.mass = 31.98982926 + (c+m)*14.01565007 - d*2.01565007
      self.formula = {'C':c+m, 'H':(2*(c+m)-2*d),'O':2}

    elif self.type == 'Headgroup':
      self.name = 'Headgroup'
      self.mass = mass
      self.formula = chnops

    else:
      self.name = '0:0'
      self.mass = Masses["H2O"]
      self.formula = {'H':2,'O':1}

def generate_tails(n):

    tail_list = []
    
    [tail_list.append(sn(c, d, type='Acyl'))
     for c in range(n[0], n[1] + 1)
      for d in range(n[2], n[3] + 1)
       if d <= (c-1)/2]
    
    return tail_list

class Glycerolipid:

  instances = []  # All created GPLs stored here
  default_tail = sn() # Keeps tail as OH

  def __init__(self, adduct_spectra_dict, sn3=default_tail, sn2=default_tail, sn1=default_tail):

    self.adduct_spectra_dict = adduct_spectra_dict  # contains adduct and info to generate spectra
    self.tails = [sn1, sn2, sn3]
    self.lipid_class = type(self).__name__  # Takes name from class which generated it
    self.mass = round(Masses["Glycerol"] + sum([snx.mass - Masses["H2O"] for snx in self.tails]), 6) 
    self.name = f"{self.lipid_class} {'_'.join(snx.name for snx in self.tails if snx.name != 'Headgroup')}"
    self.spectra = {}  # spectra for all adducts will be stored here

    # Works out CHNOPS for Glycerolipid
    formula = Counter({'C':3, 'H':8, 'O':3}) # Glycerol
    for snx in self.tails:
      formula.subtract({'H':2,'O':1}) # -H2O for bonds
      formula += snx.formula
    self.formula = dict(formula)

    Glycerolipid.instances.append(self) # Finish by appending to list of all GPLs!

  def resolve_frag(self, adduct, input={}):

    try: self.spectra[adduct]
    except KeyError: self.spectra[adduct] = []

    for mz in input:

      if input[mz] == 0:
        continue
      else: pass

      try: # Test if float, int pair
        frag = [float(mz), input[mz]]
        return [frag] # [] for below

      except: # If not float:
        try:  # Test if function, int pair
          x = mz(self, Masses[adduct])
          frag = [round(x,6), input[mz]]
          if frag not in self.spectra[adduct]:
            return [frag] # [] for below
          else: continue

        except: # If not function:
          try:  # May be list of functions
            frag = [] # [] for this
            for fragment in x:
                y = [round(fragment,6), input[mz]]
                if y not in frag:
                  frag.append(y)
                else: continue  
            return frag # [[],[]]

          except: 
            print(f"Error, fragment: {mz}")
            continue

  def resolve_spectra(self, adduct_spectra):

    for adduct in adduct_spectra:
      try: self.spectra[adduct]
      except KeyError: self.spectra[adduct] = []

      for frag in adduct_spectra[adduct]:
        self.spectra[adduct].extend(self.resolve_frag(adduct,{frag:adduct_spectra[adduct][frag]}))


# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

#  Fragments which may be generated for the spectra using the following functions
#  Isomer-ambiguous as well as isomer-specific fragments are provided.
#  Fragment list is nowhere near exhaustive, though hopefully sufficient!

# ~ # ~ # ~ # [M +/- adduct]

def MA(lipid, adduct):  # [M+-A]+-
  return (lipid.mass + adduct[0])/adduct[2]

def MA_s_H2O(lipid, adduct):  # [M+-A-H2O]+-
  return MA(lipid, adduct) - (Masses['H2O']/adduct[2])

def MA_s_2H2O(lipid, adduct):  # [M+-A-2H2O]+-
  return MA(lipid, adduct) - (2*Masses['H2O']/adduct[2])

def MA_s_PO3(lipid, adduct):  # [M+-A-PO3]+-
  return MA(lipid, adduct) - (Masses['PO3']/adduct[2])

def MA_s_PO4(lipid, adduct):  # [M+-A-PO4]+-
  return MA(lipid, adduct) - (Masses['PO4']/adduct[2])

# ~ # Fragments for PC+Na/Li

def MA_s_TMA(lipid, adduct):  # [M+-A-TMA]+-
  return MA(lipid, adduct) - (Masses['TMA']/adduct[2])

def MA_s_AZD(lipid, adduct):  # [M+-H-AZD]+-
  return MA(lipid, adduct) - (Masses['AZD']/adduct[2])

# ~ # ~ # ~ # [M +/- adduct] - fatty acid

def MA_s_FA(lipid, adduct):  # [M+-H-FA]+-
  for tail in lipid.tails:   # FA = fatty acid
    if tail.type == 'Acyl':
      yield MA(lipid, adduct) - tail.mass

def MA_s_sn1(lipid, adduct):  # [M+-H-FA]+-
  tail = lipid.tails[0]   # FA = fatty acid
  return MA(lipid, adduct) - tail.mass

def MA_s_sn2(lipid, adduct):  # [M+-H-FA]+-
  tail = lipid.tails[1]   # FA = fatty acid
  return MA(lipid, adduct) - tail.mass

# ~ #

def MA_s_FA_H2O(lipid, adduct):  # [M+-H-FA-H2O]+-
  for tail in lipid.tails:       # FA = fatty acid
    if tail.type == 'Acyl':
      yield MA_s_H2O(lipid, adduct) - tail.mass

def MA_s_sn1_H2O(lipid, adduct):  # [M+-H-FA-H2O]+-
  tail = lipid.tails[0]       # FA = fatty acid
  return MA_s_H2O(lipid, adduct) - tail.mass

def MA_s_sn2_H2O(lipid, adduct):  # [M+-H-FA-H2O]+-
  tail = lipid.tails[1]       # FA = fatty acid
  return MA_s_H2O(lipid, adduct) - tail.mass

# ~ #

def MA_s_FAk(lipid, adduct):  # [M+-H-FAk]+-
  for tail in lipid.tails:    # FAk = fatty acid (ketone)
    if tail.type == 'Acyl':
      yield MA(lipid, adduct) + Masses['H2O'] - tail.mass

def MA_s_sn1k(lipid, adduct):  # [M+-H-FAk]+-
  tail = lipid.tails[0]     # FAk = fatty acid (ketone)
  return MA(lipid, adduct) + Masses['H2O'] - tail.mass

def MA_s_sn2k(lipid, adduct):  # [M+-H-FAk]+-
  tail = lipid.tails[1]     # FAk = fatty acid (ketone)
  return MA(lipid, adduct) + Masses['H2O'] - tail.mass

# ~ #

def MA_s_FA_PO3(lipid, adduct):  # [M+-H-FA-PO3]+-
  for tail in lipid.tails:       # FA = fatty acid
    if tail.type == 'Acyl':
      yield MA_s_PO3(lipid, adduct) - tail.mass

def MA_s_sn1_PO3(lipid, adduct):  # [M+-H-FA-PO3]+-
  tail = lipid.tails[0]       # FA = fatty acid
  return MA_s_PO3(lipid, adduct) - tail.mass

def MA_s_sn2_PO3(lipid, adduct):  # [M+-H-FA-PO3]+-
  tail = lipid.tails[1]       # FA = fatty acid
  return MA_s_PO3(lipid, adduct) - tail.mass

# ~ # 

def MA_s_FAk_PO3(lipid, adduct):  # [M+-H-FAk-PO3]+-
  for tail in lipid.tails:    # FAk = fatty acid (ketone)
    if tail.type == 'Acyl':
      yield MA_s_PO3(lipid, adduct) + Masses['H2O'] - tail.mass

def MA_s_sn1k_PO3(lipid, adduct):  # [M+-H-FAk-PO3]+-
  tail = lipid.tails[0]    # FAk = fatty acid (ketone)
  return MA_s_PO3(lipid, adduct) + Masses['H2O'] - tail.mass

def MA_s_sn2k_PO3(lipid, adduct):  # [M+-H-FAk-PO3]+-
  tail = lipid.tails[1]    # FAk = fatty acid (ketone)
  return MA_s_PO3(lipid, adduct) + Masses['H2O'] - tail.mass

# ~ # Fragments for PE / PC+Na/Li

def MA_s_FA_TMA(lipid, adduct):  # [M+-H-FA-TMA]+-
  for tail in lipid.tails:       # FA = fatty acid
    if tail.type == 'Acyl':
      yield MA_s_TMA(lipid, adduct) - tail.mass

def MA_s_sn1_TMA(lipid, adduct):  # [M+-H-FA-TMA]+-
  tail = lipid.tails[0]       # FA = fatty acid
  return MA_s_TMA(lipid, adduct) - tail.mass

def MA_s_sn2_TMA(lipid, adduct):  # [M+-H-FA-TMA]+-
  tail = lipid.tails[1]       # FA = fatty acid
  return MA_s_TMA(lipid, adduct) - tail.mass

# ~ # 

def MA_s_FAk_TMA(lipid, adduct):  # [M+-H-FAk-TMA]+-
  for tail in lipid.tails:    # FAk = fatty acid (ketone)
    if tail.type == 'Acyl':
      yield MA_s_TMA(lipid, adduct) + Masses['H2O'] - tail.mass

def MA_s_sn1k_TMA(lipid, adduct):  # [M+-H-FAk-TMA]+-
  tail = lipid.tails[0]    # FAk = fatty acid (ketone)
  return MA_s_TMA(lipid, adduct) + Masses['H2O'] - tail.mass

def MA_s_sn2k_TMA(lipid, adduct):  # [M+-H-FAk-TMA]+-
  tail = lipid.tails[1]    # FAk = fatty acid (ketone)
  return MA_s_TMA(lipid, adduct) + Masses['H2O'] - tail.mass

# ~ # 

def MA_s_FA_AZD(lipid, adduct):  # [M+-H-FA-AZD]+-
  for tail in lipid.tails:       # FA = fatty acid
    if tail.type == 'Acyl':
      yield MA_s_AZD(lipid, adduct) - tail.mass

def MA_s_sn1_AZD(lipid, adduct):  # [M+-H-FA-AZD]+-
  tail = lipid.tails[0]       # FA = fatty acid
  return MA_s_AZD(lipid, adduct) - tail.mass

def MA_s_sn2_AZD(lipid, adduct):  # [M+-H-FA-AZD]+-
  tail = lipid.tails[1]       # FA = fatty acid
  return MA_s_AZD(lipid, adduct) - tail.mass

# ~ # 

def MA_s_FAk_AZD(lipid, adduct):  # [M+-H-FAk-AZD]+-
  for tail in lipid.tails:    # FAk = fatty acid (ketone)
    if tail.type == 'Acyl':
      yield MA_s_AZD(lipid, adduct) + Masses['H2O'] - tail.mass

def MA_s_sn1k_AZD(lipid, adduct):  # [M+-H-FAk-AZD]+-
  tail = lipid.tails[0]    # FAk = fatty acid (ketone)
  return MA_s_AZD(lipid, adduct) + Masses['H2O'] - tail.mass

def MA_s_sn2k_AZD(lipid, adduct):  # [M+-H-FAk-AZD]+-
  tail = lipid.tails[1]    # FAk = fatty acid (ketone)
  return MA_s_AZD(lipid, adduct) + Masses['H2O'] - tail.mass

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # [M +/- H]

def MH(lipid, adduct):  # [M+-H]+-
  if adduct[1] == 'Positive':
    return (lipid.mass + Masses['H'])
  else: return (lipid.mass - Masses['H'])

def MH_s_H2O(lipid, adduct):  # [M+-H-H2O]+-
  return MH(lipid, adduct) - Masses['H2O']
    
def MH_s_2H2O(lipid, adduct):  # [M+-H-2H2O]+-
  return MH(lipid, adduct) - 2*Masses['H2O']

def MH_s_PO3(lipid, adduct):  # [M+-H-PO3]+-
  return MH(lipid, adduct) - Masses['PO3']

def MH_s_PO4(lipid, adduct):  # [M+-H-PO4]+-
  return MH(lipid, adduct) - Masses['PO4']

# ~ # Fragments for PE / PC+Na/Li

def MH_s_TMA(lipid, adduct):  # [M+-H-TMA]+-
  return MH(lipid, adduct) - Masses['TMA']

def MH_s_AZD(lipid, adduct):  # [M+-H-AZD]+-
  return MH(lipid, adduct) - Masses['AZD']

# ~ # ~ # ~ # [M +/- H] - fatty acid

def MH_s_FA(lipid, adduct):  # [M+-H-FA]+-
  for tail in lipid.tails:   # FA = fatty acid
    if tail.type == 'Acyl':
      yield MH(lipid, adduct) - tail.mass

def MH_s_sn1(lipid, adduct):  # [M+-H-FA]+-
  tail = lipid.tails[0]   # FA = fatty acid
  return MH(lipid, adduct) - tail.mass

def MH_s_sn2(lipid, adduct):  # [M+-H-FA]+-
  tail = lipid.tails[1]   # FA = fatty acid
  return MH(lipid, adduct) - tail.mass

# ~ # 

def MH_s_FA_H2O(lipid, adduct):  # [M+-H-FA-H2O]+-
  for tail in lipid.tails:       # FA = fatty acid
    if tail.type == 'Acyl':
      yield MH_s_H2O(lipid, adduct) - tail.mass

def MH_s_sn1_H2O(lipid, adduct):  # [M+-H-FA-H2O]+-
  tail = lipid.tails[0]       # FA = fatty acid
  return MH_s_H2O(lipid, adduct) - tail.mass

def MH_s_sn2_H2O(lipid, adduct):  # [M+-H-FA-H2O]+-
  tail = lipid.tails[1]       # FA = fatty acid
  return MH_s_H2O(lipid, adduct) - tail.mass

# ~ # 

def MH_s_FAk(lipid, adduct):  # [M+-H-FAk]+-
  for tail in lipid.tails:    # FAk = fatty acid (ketone)
    if tail.type == 'Acyl':
      yield MH(lipid, adduct) + Masses['H2O'] - tail.mass

def MH_s_sn1k(lipid, adduct):  # [M+-H-FAk]+-
  tail = lipid.tails[0]    # FAk = fatty acid (ketone)
  return MH(lipid, adduct) + Masses['H2O'] - tail.mass

def MH_s_sn2k(lipid, adduct):  # [M+-H-FAk]+-
  tail = lipid.tails[1]    # FAk = fatty acid (ketone)
  return MH(lipid, adduct) + Masses['H2O'] - tail.mass

# ~ # 

def MH_s_FA_PO3(lipid, adduct):  # [M+-H-FA-PO3]+-
  for tail in lipid.tails:       # FA = fatty acid
    if tail.type == 'Acyl':
      yield MH_s_PO3(lipid, adduct) - tail.mass

def MH_s_sn1_PO3(lipid, adduct):  # [M+-H-FA-PO3]+-
  tail = lipid.tails[0]       # FA = fatty acid
  return MH_s_PO3(lipid, adduct) - tail.mass

def MH_s_sn2_PO3(lipid, adduct):  # [M+-H-FA-PO3]+-
  tail = lipid.tails[1]       # FA = fatty acid
  return MH_s_PO3(lipid, adduct) - tail.mass

# ~ # 

def MH_s_FAk_PO3(lipid, adduct):  # [M+-H-FAk]+-
  for tail in lipid.tails:    # FAk = fatty acid (ketone)
    if tail.type == 'Acyl':
      yield MH_s_PO3(lipid, adduct) + Masses['H2O'] - tail.mass

def MH_s_sn1k_PO3(lipid, adduct):  # [M+-H-FAk]+-
  tail = lipid.tails[0]    # FAk = fatty acid (ketone)
  return MH_s_PO3(lipid, adduct) + Masses['H2O'] - tail.mass

def MH_s_sn2k_PO3(lipid, adduct):  # [M+-H-FAk]+-
  tail = lipid.tails[1]    # FAk = fatty acid (ketone)
  return MH_s_PO3(lipid, adduct) + Masses['H2O'] - tail.mass

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # [M +/- 2H]

def M2H(lipid, adduct):  # [M+-2H]2+-
  if adduct[1] == 'Positive':
    return (lipid.mass + 2*Masses['H'])/2
  else: return (lipid.mass - 2*Masses['H'])/2

def M2H_s_H2O(lipid, adduct):  # [M+-2H-H2O]2+-
  return M2H(lipid, adduct) - (Masses['H2O']/2)

def M2H_s_2H2O(lipid, adduct):  # [M+-2H-2H2O]2+-
  return M2H(lipid, adduct) - Masses['H2O']

# ~ # ~ # ~ # [M +/- 2H] - fatty acid

def M2H_s_FA(lipid, adduct):   # [M+-2H-FA]+-
  for tail in lipid.tails:     # FA = fatty acid
    if tail.type == 'Acyl':
      yield M2H(lipid, adduct) - tail.mass/2

def M2H_s_sn1(lipid, adduct):   # [M+-2H-FA]+-
  tail = lipid.tails[0]   # FA = fatty acid
  return M2H(lipid, adduct) - tail.mass/2

def M2H_s_sn2(lipid, adduct):   # [M+-2H-FA]+-
  tail = lipid.tails[1]   # FA = fatty acid
  return M2H(lipid, adduct) - tail.mass/2

# ~ # 

def M2H_s_FAk(lipid, adduct):  # [M+-2H-FAk]+-
  for tail in lipid.tails:       # FAk = fatty acid (ketone)
    if tail.type == 'Acyl':
      yield M2H(lipid, adduct) + Masses['H2O'] - tail.mass/2

def M2H_s_sn1k(lipid, adduct):  # [M+-2H-FAk]+-
  tail = lipid.tails[0]       # FAk = fatty acid (ketone)
  return M2H(lipid, adduct) + Masses['H2O'] - tail.mass/2

def M2H_s_sn2k(lipid, adduct):  # [M+-2H-FAk]+-
  tail = lipid.tails[1]       # FAk = fatty acid (ketone)
  return M2H(lipid, adduct) + Masses['H2O'] - tail.mass/2

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # Free fatty acids and fatty acid ketones

def FAH(lipid, adduct):  # Free fatty acid
  for tail in lipid.tails:
    if tail.type == 'Acyl':
    #if tail[0] not in ['Head', '0:0', ';O']:
      yield tail.mass - Masses['H']

def sn1(lipid, adduct):  # Free fatty acid
  tail = lipid.tails[0]
  return tail.mass - Masses['H']

def sn2(lipid, adduct):  # Free fatty acid
  tail = lipid.tails[1]
  return tail.mass - Masses['H']

# ~ # 

def FAkH(lipid, adduct):  # Free fatty acid (ketone)
  for tail in lipid.tails:
    if tail.type == 'Acyl':
      if adduct[1] == 'Positive':
        yield tail.mass - Masses['H2O'] + Masses['H']
      else: yield tail.mass - Masses['H2O'] - Masses['H']

def sn1k(lipid, adduct):  # Free fatty acid (ketone)
  tail = lipid.tails[0]
  if adduct[1] == 'Positive':
    return tail.mass - Masses['H2O'] + Masses['H']
  else: return tail.mass - Masses['H2O'] - Masses['H']

def sn2k(lipid, adduct):  # Free fatty acid (ketone)
  tail = lipid.tails[1]
  if adduct[1] == 'Positive':
    return tail.mass - Masses['H2O'] + Masses['H']
  else: return tail.mass - Masses['H2O'] - Masses['H']

# ~ # 

def FAkA(lipid, adduct):
  for tail in lipid.tails:
    if tail.type == 'Acyl':
      yield tail.mass - Masses['H2O'] + adduct[0]

def sn1kA(lipid, adduct):
  tail = lipid.tails[0]
  return tail.mass - Masses['H2O'] + adduct[0]

def sn2kA(lipid, adduct):
  tail = lipid.tails[1]
  return tail.mass - Masses['H2O'] + adduct[0]

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # Headgroup neutral losses, variants
# A = Headgroup NL, with phosphate
# B = Headgroup NL, without phosphate
# C = Headgroup NL, MA -> MH

# ~ # Sometimes the headgroup takes the phosphate with it

def HG_NL_A(lipid, adduct):  # Headgroup neutral loss
  for tail in lipid.tails:
    if tail.type == 'Headgroup':
      hg = tail  # Fragments common to some GPLs.
      return MA(lipid, adduct) - hg.mass + Masses['H2O']

def HG_NL_H2O_A(lipid, adduct):
  return HG_NL_A(lipid, adduct) - Masses['H2O']

# ~ # 

def HG_FA_NL_A(lipid, adduct):
  for tail in lipid.tails:
    if tail.type == 'Acyl':
      yield HG_NL_A(lipid, adduct) - tail.mass

def HG_sn1_NL_A(lipid, adduct):
  tail = lipid.tails[0]
  return HG_NL_A(lipid, adduct) - tail.mass

def HG_sn2_NL_A(lipid, adduct):
  tail = lipid.tails[1]
  return HG_NL_A(lipid, adduct) - tail.mass

# ~ # 

def HG_FAk_NL_A(lipid, adduct):
  for tail in lipid.tails:
    if tail.type == 'Acyl':
      yield HG_NL_A(lipid, adduct) + Masses['H2O'] - tail.mass

def HG_sn1k_NL_A(lipid, adduct):
  tail = lipid.tails[0]
  return HG_NL_A(lipid, adduct) + Masses['H2O'] - tail.mass

def HG_sn2k_NL_A(lipid, adduct):
  tail = lipid.tails[1]
  return HG_NL_A(lipid, adduct) + Masses['H2O'] - tail.mass

# ~ # 

def HG_FA_NL_H2O_A(lipid, adduct):
  for tail in lipid.tails:
    if tail.type == 'Acyl':
      yield HG_NL_H2O_A(lipid, adduct) - tail.mass

def HG_sn1_NL_H2O_A(lipid, adduct):
  tail = lipid.tails[0]
  return HG_NL_H2O_A(lipid, adduct) - tail.mass

def HG_sn2_NL_H2O_A(lipid, adduct):
  tail = lipid.tails[1]
  return HG_NL_H2O_A(lipid, adduct) - tail.mass

# ~ #  Sometimes the headgroup leaves the phosphate behind

def HG_NL_B(lipid, adduct):  # Headgroup neutral loss
  for tail in lipid.tails:
    if tail.type == 'Headgroup':
      hg = tail  # Fragments common to some GPLs.
      return MH(lipid, adduct) - (hg.mass - Masses['PO4'])

def HG_NL_2B(lipid, adduct):  # Headgroup neutral loss
  for tail in lipid.tails:
    if tail.type == 'Headgroup':
      hg = tail  # Special case for PIPs
      return MH(lipid, adduct) - (hg.mass - 2*Masses['PO3'])

def HG_NL_H2O_B(lipid, adduct):
  return HG_NL_B(lipid, adduct) - Masses['H2O']

# ~ # 

def HG_FA_NL_B(lipid, adduct):
  for tail in lipid.tails:
    if tail.type == 'Acyl':
      yield HG_NL_B(lipid, adduct) - tail.mass

def HG_sn1_NL_B(lipid, adduct):
  tail = lipid.tails[0]
  return HG_NL_B(lipid, adduct) - tail.mass

def HG_sn2_NL_B(lipid, adduct):
  tail = lipid.tails[1]
  return HG_NL_B(lipid, adduct) - tail.mass

# ~ # 

def HG_FAk_NL_B(lipid, adduct):
  for tail in lipid.tails:
    if tail.type == 'Acyl':
      yield HG_NL_B(lipid, adduct) + Masses['H2O'] - tail.mass

def HG_sn1k_NL_B(lipid, adduct):
  tail = lipid.tails[0]
  return HG_NL_B(lipid, adduct) + Masses['H2O'] - tail.mass

def HG_sn2k_NL_B(lipid, adduct):
  tail = lipid.tails[1]
  return HG_NL_B(lipid, adduct) + Masses['H2O'] - tail.mass

# ~ # 

def HG_FA_NL_H2O_B(lipid, adduct):
  for tail in lipid.tails:
    if tail.type == 'Acyl':
      yield HG_NL_H2O_B(lipid, adduct) - tail.mass

def HG_sn1_NL_H2O_B(lipid, adduct):
  tail = lipid.tails[0]
  return HG_NL_H2O_B(lipid, adduct) - tail.mass

def HG_sn2_NL_H2O_B(lipid, adduct):
  tail = lipid.tails[1]
  return HG_NL_H2O_B(lipid, adduct) - tail.mass    

# ~ #  Sometimes the headgroup leaves with the adduct, leaving M+/-H

def HG_NL_C(lipid, adduct):  # Headgroup neutral loss
  for tail in lipid.tails:
    if tail.type == 'Headgroup':
      hg = tail  # Fragments common to some GPLs.
      return MH(lipid, adduct) - hg.mass + Masses['H2O']

def HG_NL_H2O_C(lipid, adduct):
  return HG_NL_C(lipid, adduct) - Masses['H2O']

# ~ # 

def HG_FA_NL_C(lipid, adduct):
  for tail in lipid.tails:
    if tail.type == 'Acyl':
      yield HG_NL_C(lipid, adduct) - tail.mass

def HG_sn1_NL_C(lipid, adduct):
  tail = lipid.tails[0]
  return HG_NL_C(lipid, adduct) - tail.mass

def HG_sn2_NL_C(lipid, adduct):
  tail = lipid.tails[1]
  return HG_NL_C(lipid, adduct) - tail.mass

# ~ # 

def HG_FAk_NL_C(lipid, adduct):
  for tail in lipid.tails:
    if tail.type == 'Acyl':
      yield HG_NL_C(lipid, adduct) + Masses['H2O'] - tail.mass

def HG_sn1k_NL_C(lipid, adduct):
  tail = lipid.tails[0]
  return HG_NL_C(lipid, adduct) + Masses['H2O'] - tail.mass

def HG_sn2k_NL_C(lipid, adduct):
  tail = lipid.tails[1]
  return HG_NL_C(lipid, adduct) + Masses['H2O'] - tail.mass

# ~ # 

def HG_FA_NL_H2O_C(lipid, adduct):
  for tail in lipid.tails:
    if tail.type == 'Acyl':
      yield HG_NL_H2O_A(lipid, adduct) - tail.mass

def HG_sn1_NL_H2O_C(lipid, adduct):
  tail = lipid.tails[0]
  return HG_NL_H2O_A(lipid, adduct) - tail.mass

def HG_sn2_NL_H2O_C(lipid, adduct):
  tail = lipid.tails[1]
  return HG_NL_H2O_A(lipid, adduct) - tail.mass

# ~ #  Sometimes the headgroup flies off the the adduct!

def HGA(lipid, adduct):  # Headgroup + Adduct
  for tail in lipid.tails:
    if tail.type == 'Headgroup':
      hg = tail  # Fragments common to some GPLs.
      return (hg.mass + adduct[0])/adduct[2]

# ~ # Goodness gracious that's a lot of fragments...