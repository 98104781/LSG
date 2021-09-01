
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
   97.976895}

#  [Cmin, Cmax, Dmin, Dmax]
Tails_to_Generate = [8, 26, 0, 6]

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class Adduct_Spectra:

  def __init__(self, lipid, adduct, mz_list):

    self.spectrum = []

    for mz in mz_list:
      fragment = mz(lipid, Masses[adduct])
      intensity = mz_list[mz]  

      # Try to call function 'mz' to generate fragment:
      try: self.spectrum.append([round(fragment, 6), intensity])

      except: # If it fails (probably not a callable function)
        try:  # It could be a list of fragments and iterate through:
          for a in fragment:
            ion = [round(a, 6), intensity]
            if ion not in self.spectrum:
              self.spectrum.append([round(a, 6), intensity])

        except: # If that fails, it could just be a float!
          try: self.spectrum.append([float(mz), intensity])

          # If that fails, I don't know anymore...
          except: print(f"Error assigning fragment: {mz}")


class Glycerolipid:

  instances = []  # All created GPLs stored here
  default_tail = ['0:0', Masses["H2O"]] # Keeps tail as OH


  def __init__(self, adducts, sn3=default_tail, sn2=default_tail, sn1=default_tail):

    self.adducts = adducts  # contains adduct and info to generate spectra
    self.tails = [sn3, sn2, sn1]
    self.lipid_class = type(self).__name__  # Takes name from class which generated it
    self.mass = round(Masses["Glycerol"] + sum([sn[1] - Masses["H2O"] for sn in self.tails]), 6) 
    self.name = f"{self.lipid_class} {'_'.join(sn[0] for sn in self.tails if sn[0] != 'Head')}"
    self.spectra = {}  # spectra for all adducts will be stored here

    Glycerolipid.instances.append(self) # Finish by appending to list of all GPLs!


  def generate_spectra(self):

    for adduct in self.adducts:
      self.spectra.update(
        {adduct: Adduct_Spectra(self, adduct, self.adducts[adduct]).spectrum})


def generate_tails(n):

    tail_list = []
    #                 [Name      , O2 mass     + c*CH2 mass    - d*H2 mass   ]
    [tail_list.append([f"{c}:{d}", 31.98982926 + c*14.01565007 - d*2.01565007])
     for c in range(n[0], n[1] + 1)
      for d in range(n[2], n[3] + 1)
       if d <= (c-1)/2]

    return tail_list

tails = generate_tails(Tails_to_Generate)

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

#  Fragments which may be generated for the spectra:

def MA(lipid, adduct):  # [M+-A]+-
  return (lipid.mass + adduct[0])/adduct[2]

def MA_s_H2O(lipid, adduct):  # [M+-A-H2O]+-
    return MA(lipid, adduct) - (Masses['H2O']/adduct[2])

def MA_s_2H2O(lipid, adduct):  # [M+-A-2H2O]+-
    return MA(lipid, adduct) - (2*Masses['H2O']/adduct[2])

def MA_s_PO3(lipid, adduct):  # [M+-A-PO3]+-
    return MA(lipid, adduct) - (Masses['PO3']/adduct[2])

# ~ # ~ # ~ #

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

# ~ # ~ # ~ #

def M2H(lipid, adduct):  # [M+-2H]2+-
  if adduct[1] == 'Positive':
    return (lipid.mass + 2*Masses['H'])/2
  else: return (lipid.mass - 2*Masses['H'])/2

def M2H_s_H2O(lipid, adduct):  # [M+-2H-H2O]2+-
  return M2H(lipid, adduct) - (Masses['H2O']/2)

def M2H_s_2H2O(lipid, adduct):  # [M+-2H-2H2O]2+-
  return M2H(lipid, adduct) - Masses['H2O']

def M2H_s_FA(lipid, adduct):   # [M+-2H-FA]+-
  for tail in lipid.tails:     # FA = fatty acid
    yield M2H(lipid, adduct) - tail[1]/2

def M2H_sub_FAk(lipid, adduct):  # [M+-2H-FAk]+-
  for tail in lipid.tails:       # FAk = fatty acid (ketone)
    yield M2H(lipid, adduct) + Masses['H2O'] - tail[1]/2

# ~ # ~ # ~ #

def MH_s_FA(lipid, adduct):  # [M+-H-FA]+-
  for tail in lipid.tails:   # FA = fatty acid
    if tail[0] not in ['Head', '0:0']:
      yield MH(lipid, adduct) - tail[1]

def MH_s_FA_H2O(lipid, adduct):  # [M+-H-FA-H2O]+-
  for tail in lipid.tails:       # FA = fatty acid
    if tail[0] not in ['Head', '0:0']:
      yield MH_s_H2O(lipid, adduct) - tail[1]

def MH_s_FAk(lipid, adduct):  # [M+-H-FAk]+-
  for tail in lipid.tails:    # FAk = fatty acid (ketone)
    if tail[0] not in ['Head', '0:0']:
      yield MH(lipid, adduct) + Masses['H2O'] - tail[1]

def MH_s_FA_PO3(lipid, adduct):  # [M+-H-FA-PO3]+-
  for tail in lipid.tails:       # FA = fatty acid
    if tail[0] not in ['Head', '0:0']:
      yield MH_s_PO3(lipid, adduct) - tail[1]

def MH_s_FAk_PO3(lipid, adduct):  # [M+-H-FAk]+-
  for tail in lipid.tails:    # FAk = fatty acid (ketone)
    if tail[0] not in ['Head', '0:0']:
      yield MH_s_PO3(lipid, adduct) + Masses['H2O'] - tail[1]

# ~ # ~ # ~ #

def FAH(lipid, adduct):  # Free fatty acid
  for tail in lipid.tails:
    if tail[0] not in ['Head', '0:0']:
      yield tail[1] - Masses['H']

def FAkH(lipid, adduct):  # Free farty acid (ketone)
  for tail in lipid.tails:
    if tail[0] not in ['Head', '0:0']:
      if adduct[1] == 'Positive':
        yield tail[1] - Masses['H2O'] + Masses['H']
      else: yield tail[1] - Masses['H2O'] - Masses['H']

def FAkA(lipid, adduct):
  for tail in lipid.tails:
    if tail[0] not in ['Head', '0:0']:
      yield tail[1] - Masses['H2O'] + adduct[0]

# ~ # ~ # ~ #

def HG_NL(lipid, adduct):  # Headgroup neutral loss
  for tail in lipid.tails:
    if tail[0] == 'Head':
      HG = tail  # Fragments common to some GPLs.
      if adduct[1] == 'Positive':
        return MA(lipid, adduct) - HG[1]
      if adduct[1] == 'Negative':
        return MH(lipid, adduct) - (HG[1] - Masses['PO4'])

def HG_NL_H2O(lipid, adduct):
  return HG_NL(lipid, adduct) - Masses['H2O']

def HG_FA_NL(lipid, adduct):
  for tail in lipid.tails:
    if tail[0] not in ['Head', '0:0']:
      yield HG_NL(lipid, adduct) - tail[1]

def HG_FA_NL_H2O(lipid, adduct):
  for tail in lipid.tails:
    if tail[0] not in ['Head', '0:0']:
      yield HG_NL_H2O(lipid, adduct) - tail[1]