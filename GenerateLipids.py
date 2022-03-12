from collections import Counter
import copy

Masses = {

  #[Delta_mz,       Polarity,  z,   chnops, smiles]
  "[M-H]-":     
   [ -1.007276467, 'Negative', -1, {'H':-1},         ''],

  "[M-2H]2-":   
   [ -2.014552934, 'Negative', -2, {'H':-2},         ''],

  "[M-3H]3-":   
   [ -3.021829401, 'Negative', -3, {'H':-3},         ''],

  "[M+Cl]-": # Cl 35, as is major isotope    
   [ 34.96830412,  'Negative', -1, {'Cl':1},         '[Cl-].'],
   # Calculated by Mass of CL[35] - Mass of e-

  "[M+Na-2H]-":    
   [ 20.974668486, 'Negative', -1, {'Na':1, 'H':-2}, '[Na+].'],

  "[M+2Na-3H]-":    
   [ 42.956613439, 'Negative', -1, {'Na':2, 'H':-3}, '[Na+].[Na+].'],

  "[M+H]+":     
   [  1.007276467, 'Positive',  1, {'H':1},          '[H+].'],

  "[M+H-H2O]+": 
   [-17.003288217, 'Positive',  1, {'H':-1, 'O':-1}, '[H+].'],

  "[M+Na]+":    
   [ 22.98922142,  'Positive',  1, {'Na':1},         '[Na+].'],
   # Calculated by Mass of Na[23] - Mass of e-

  "[M+Li]+": # Li 7 as is major isotope    
   [  7.01545542,  'Positive' , 1, {'Li':1},         '[Li+].'],
   # Calculated by Mass of Li[7] - Mass of e-

  "[M+NH4]+":   
   [ 18.033825568, 'Positive',  1, {'N':1, 'H':4},   '[NH4+].'],
   
  #  Handy masses 
  "e-":           
    0.000548580,

  "H":           
    1.007825032,

  "H+":           
    1.007276467,

  "OH-": # Calculated by H2O - H+           
   17.003288217, 

  "H2O":         
   18.010564684,

  "Glycerol":
   92.047344116,

  "NH3":
   17.026549101,

  "PO3H":
   79.966329892,

  "PO4H3":
   97.976894576,

  "TMA": # Trimethylamine, fragment for PC+Na/Li
   59.073499293,

  "AZD": # Aziridine, fragment for PE+Na/Li
   43.042199165}

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class sn:
  '''Moiety to substitute onto glycerol backbone\n
  e.g. sn1, sn2, sn3 fatty acids on a triacylglycerol\n
  c = carbon number of fatty acid\n
  d = desaturation number of fatty acid\n
  m = number of branching methyl-groups (currently unused)\n
  mass = mass reset group is provided (e.g. a headgroup)\n
  chnops = formula if preset group is provided (e.g. a headgroup)\n
  type = 'Acyl', 'Ether', 'Methyl', 'Headgroup', anything else will return water\n
  providing no parameters also returns a water (-OH), which does not modify backbone.
  '''
  def __init__(self, c=0, d=0, mass=None, chnops={}, smiles='', type=None, hgtails=[], me=0, oh=0, dt=0):

    self.type = type # Tail list is shared. Deepcopy made as to not overwrite
    self.hgtails = copy.deepcopy(hgtails) # the tails in the tail list.
    string = [] # string is created in list so that constituent order can be reversed

    self.c  = c
    self.d  = d
    self.me = me
    self.oh = oh
    self.dt = dt
    
    if self.type == 'Acyl':
      self.name = f"{self.c}:{self.d}"
      #           O2 mass     + c*CH2 mass    - d*H2 mass
      self.mass = 31.98982924 + self.c*14.01565007 - self.d*2.01565007
      self.formula = Counter({'C':self.c, 'H':(2*self.c-2*self.d),'O':2})
      string = ['C(=O)',self.d*'C=C',self.oh*'C(O)',(self.c-1-2*self.d-self.oh)*'C']
      # string is used to generate smiles for the tails. Made as a list first
      # in case the order needs to be reversed later on.
      self.smiles = ''.join(string) # smiles order can be reversed.
      self.inverseSmiles = ''.join(string[::-1])

    elif self.type == 'Ether':
      self.name = f"O-{self.c}:{self.d}"
      #           H2O mass      + c*CH2 mass    - d*H2 mass
      self.mass = Masses['H2O'] + self.c*14.01565007 - self.d*2.01565007
      self.formula = Counter({'C':self.c, 'H':(2*self.c-2*(self.d-1)),'O':1})
      string = ['C',self.d*'C=C',self.oh*'C(O)',(self.c-1-2*self.d-self.oh)*'C']
      self.smiles = ''.join(string) # Repeated in each as to not overwrite
      self.inverseSmiles = ''.join(string[::-1]) # self.smiles set for headgroup

    elif self.type == 'Vinyl':
      self.name = f"P-{self.c}:{self.d}"
      #           H2O mass      + c*CH2 mass    - d*H2 mass
      self.mass = Masses['H2O'] + self.c*14.01565007 - (self.d+1)*2.01565007
      self.formula = Counter({'C':self.c, 'H':(2*self.c-2*self.d),'O':1})
      string = ['\C=C/',self.d*'C=C',self.oh*'C(O)',(self.c-2-2*self.d-self.oh)*'C']
      self.smiles = ''.join(string) # Repeated in each as to not overwrite
      self.inverseSmiles = ''.join(string[::-1]) # self.smiles set for headgroup

    elif self.type == 'Headgroup':
      self.name = 'Headgroup'
      self.mass = mass
      self.formula = Counter(chnops)
      self.smiles = smiles
      for tail in self.hgtails: # Headgroup can be acyl-functionalised
        tail.type = 'HeadTail'
        self.mass += (tail.mass-Masses['H2O'])
        self.formula.update(tail.formula)
        self.formula.subtract({'H':2,'O':1})

    else: # If nothing, just give it values for water
      self.name = '0:0'
      self.mass = Masses['H2O']
      self.formula = Counter({'H':2,'O':1})
      self.smiles = ''

    # Perhaps exclude?  Identical to fatty acid of c = c+m
    if self.me > 0: # Methyl branching of fatty acid
      self.name += f";{self.me}-M" 
      self.mass += self.me*14.01565007
      self.formula += {'C':self.me, 'H':2*self.me}
    if self.oh > 0: # Hydroxy functionalisation of fatty acid
      self.name += f";O{self.oh}"
      self.mass += self.oh*15.99491462
      self.formula += {'O':self.oh}
    if self.dt > 0: # deuterium labelled fatty acids
      self.name += f"(D{self.dt})" # Deuterium doesn't update smiles currently.
      self.mass += self.dt*1.006276746 # Calculated by D - H
      self.formula += {'H':-self.dt, 'D':self.dt}

  def __hash__(self):
    return hash(('name', self.name))
  def __eq__(self, other):
    return self.name == other.name
  def __lt__(self, other):
    if self.c  < other.c: return True
    elif self.c == other.c and self.d < other.d: return True
    else: return False

def generate_tails(n, type):
  tail_list = []
  [tail_list.append(sn(c, d, type=type, oh=oh))
    for c in range(n[0], n[1] + 1)
    for d in range(n[2], n[3] + 1)
      if d <= (c-1)/2
      for oh in range(0, n[4]+1)
      if d+oh <= (c-1)/2]
  return tail_list

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class base:
  '''Moiety to use as backbone for sphingoids\n
  c = carbon number of base\n
  type = determines oh groups and desaturation [Sphinganine, Sphingosine, Phytosphingosine]\n
  providing no parameters also returns a water (-OH), which does not modify backbone.
  '''
  def __init__(self, c=3, type=None, dt=0):

    self.type = type
    # The bases, (ceramide bodies) can fall into one of several categories.
    # The categories are listed here (instead of just allowing the loop to
    # generate a range) because the fragmentation patterns change depending on
    # the positioning of functional groups.
    if self.type == 'Sphinganine':
      self.name = f"{c}:0;O2"
      #           H2O mass      + c*CH2 mass    + 2* O mass
      self.mass = Masses['H2O'] + c*14.01565007 + 2*15.99491462
      self.formula = Counter({'C':c, 'H':(2+2*c),'O':3})
      self.smiles = 'OCC(N)C(O)'
      self.smiles += (c-3)*'C'

    elif self.type == 'Sphingosine':
      self.name = f"{c}:1;O2"
      #           H2O mass      + c*CH2 mass    - H2 mass    + 2* O mass
      self.mass = Masses['H2O'] + c*14.01565007 - 2.01565007 + 2*15.99491462
      self.formula = Counter({'C':c, 'H':(2*c),'O':3})
      self.smiles = 'OCC(N)C(O)/C=C/'
      self.smiles += (c-5)*'C'

    elif self.type == 'Phytosphingosine':
      self.name = f"{c}:0;O3"
      #           H2O mass      + c*CH2 mass    + 3* O mass
      self.mass = Masses['H2O'] + c*14.01565007 + 3*15.99491462
      self.formula = Counter({'C':c, 'H':(2+2*c),'O':4})
      self.smiles = 'OCC(N)C(O)C(O)'
      self.smiles += (c-4)*'C'
      # Despite being based around an ammonia group, H2O is still
      # used as the default mass. These groups are attached to an
      # ammonia group in the 'Sphingolipid' class like how fatty 
      # acids are attached to a glycerol in the 'Glycerolipid' class.
      # The water is removed when the bonds are formed. This is to be
      # consistant with the headgroups which are attached via an ether /
      # ester on the base backbone and not the ammonia group.
    else: # If nothing, just give it values for water
      self.name = '0:0'
      self.mass = Masses['H2O']
      self.formula = Counter({'H':2,'O':1})
      self.smiles = ''

    if dt > 0: # deuterium labels
      self.name += f"(D{dt})"
      self.mass += dt*1.006276746 # Calculated by D - H
      self.formula += {'H':-dt, 'D':dt}

def generate_base_tails(n):
  base_dict = {}
  for type in n[2]:
    base_list = []
    [base_list.append(base(c, type))
      for c in range(n[0], n[1] + 1)
      if c > 6]
    base_dict[type] = base_list
  return base_dict

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class Other:
  def __init__(self, name='H2O', mass=Masses['H2O'], chnops={'H':2,'O':1}, smiles='', dt=0):
    
    self.name = name
    self.mass = mass
    self.formula = Counter(chnops)
    self.smiles=smiles

    if dt > 0: # deuterium labels
      self.name += f"(D{dt})"
      self.mass += dt*1.006276746 # Calculated by D - H
      self.formula += {'H':-dt, 'D':dt}

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #
class Lipid:
  def __init__(self, adducts):

    self.adducts = adducts
    self.spectra = {}
    self.smiles = ''

  def resolve_spectra(self, adduct, spectra={}):
    x = []
    for fragment, intensity in spectra.items():
      fgmt = fragment(self, adduct, intensity)
      try: x.extend(fgmt)
      except:
        try: x.append(fgmt)
        except:
          print('Error assigning', fragment)
    self.spectra[adduct] = sorted(set(x), reverse=True)

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class Glycerolipid(Lipid):
  def __init__(self, adducts, sn1=sn(), sn2=sn(), sn3=sn()):
    super().__init__(adducts)

    self.tails = []
    if sn3.type in ['Headgroup']: 
      self.tails.extend([tail for tail in sn3.hgtails])
    self.tails.extend([sn1, sn2, sn3])

    self.lipid_class = type(self).__name__  # Takes name from class which generated it
    self.name = f"{self.lipid_class} "
    self.name += f"{'_'.join(snx.name for snx in self.tails if snx.type != 'Headgroup')}" # Headtail mass factored into headgroup
    self.mass = round(Masses['Glycerol'] + sum([snx.mass-Masses['H2O'] for snx in self.tails if snx.type != 'HeadTail']), 6)

    if sn3.type in ['Acyl', 'Ether', 'Vinyl']: string = sn3.inverseSmiles
    else: string = sn3.smiles # TAGs need the first tail reversed.
    self.smiles = f"{string}OCC(O{sn2.smiles})CO{sn1.smiles}"

    self.formula = Counter({'C':3, 'H':8, 'O':3})  # Glycerol
    for snx in self.tails:  # Works out CHNOPS for lipid
      if snx.type != 'HeadTail': # Headtails factored into headgroup
        self.formula.update(snx.formula)
        self.formula.subtract({'H':2,'O':1}) # -H2O for bonding

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class Sphingolipid(Lipid):
  def __init__(self, adducts, base=base(), sn1=sn(), headgroup=sn()):
    super().__init__(adducts)

    self.tails = [base, sn1, headgroup] # Base and headgroup included to be consistant with fragment generation 
    self.lipid_class = type(self).__name__ # Takes name from class which generated it
    self.name = f"{self.lipid_class} {'_'.join(snx.name for snx in self.tails if snx.name not in ['Headgroup', '0:0'])}"
    self.mass = round(Masses['NH3'] + sum([snx.mass-Masses['H2O'] for snx in self.tails]), 6)
    self.smiles= f'{headgroup.smiles}{base.smiles[:5]}{sn1.smiles}{base.smiles[5:]}'

    self.formula = Counter({'H':3,'N':1}) # Unlike GPLs, sphingoids built around the base
    for snx in self.tails: # Works out CHNOPS for lipid
      self.formula.update(snx.formula)
      self.formula.subtract({'H':2,'O':1}) # -H2O for bonding)

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class OtherLipid(Lipid): # For cholsterols?
  def __init__(self, adducts, body=Other(), sn1=sn()):
    super().__init__(adducts)

    self.tails = [sn1]
    self.lipid_class = type(self).__name__
    self.name = f"{body.name} {(sn1.name if sn1.name not in ['Headgroup', '0:0'] else '')}"
    self.mass = round(body.mass + sn1.mass - Masses['H2O'], 6)
    self.smiles = body.smiles

    self.formula = body.formula
    self.formula.update(sn1.formula)
    self.formula.subtract({'H':2,'O':1})

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #
#  Fragments which may be generated for the spectra using the following functions
#  Isomer-ambiguous as well as isomer-specific fragments are provided.
#  Fragment list is nowhere near exhaustive, though hopefully sufficient!

class Fragment:
  ''' Fragments are used to generate mass spectra for lipids'''
  def __init__(self, lipid, adduct, intensity, fragmentType = None):

    self.lipid = lipid
    self.adduct = adduct
    self.intensity = intensity
    self.mass = round(self.MZ(), 6)

    if fragmentType is not None:
      self.fragmentType = fragmentType
    else: self.fragmentType = type(self)

  def MZ(self):
    return None
  def Formula(self):
    return None
  def Charge(self):
    return (Masses[self.adduct][2]/
    abs(Masses[self.adduct][2]))
  def Comment(self):
    return ""
  def Smiles(self):
    return ""

  def __hash__(self):
    return hash(('mass', self.mass))
  def __eq__(self, other):
    return self.mass == other.mass
  def __lt__(self, other):
    return self.mass < other.mass

# ~ # ~ # ~ # [M +/- adduct]

class MA(Fragment):
  '''[ MA ]\n
  Fragment for molecular ion with adduct'''
  def MZ(self):
    return (self.lipid.mass + Masses[self.adduct][0])/abs(Masses[self.adduct][2])
  def Formula(self):
    formula = Counter(self.lipid.formula)
    formula.update(Masses[self.adduct][3])
    return formula
  def Charge(self):
    return Masses[self.adduct][2]
  def Comment(self):
    return self.adduct

class MA_s_H2O(MA):
  '''[ MA - H2O ]\n
  Fragment for adducted molecular ion, with loss of water'''
  def MZ(self):
    return super().MZ() - (Masses['H2O']/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'H':2,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-H2O]')
    return comment

class MA_s_2H2O(MA):
  '''[ MA - H4O2 ]\n
  Fragment for adducted molecular ion, with loss of two waters'''
  def MZ(self):
    return super().MZ() - (2*Masses['H2O']/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'H':4,'O':2})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-2H2O]')
    return comment

class MA_s_PO3(MA):
  '''[ MA - PO3 ]\n
  Fragment for adducted molecular ion, with loss of phosphite'''
  def MZ(self):
    return super().MZ() - (Masses['PO3H']/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'P':1,'O':3, 'H':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-HPO3]')
    return comment


class MA_s_PO4(MA):
  '''[ MA - PO4 ]\n
  Fragment for adducted molecular ion, with loss of phosphate'''
  def MZ(self):
    return super().MZ() - (Masses['PO4H3']/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'P':1, 'O':4, 'H':3})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-H3PO4]')
    return comment

# ~ # Fragments for DGDG, AcPIMs

class MA_s_Gal(MA):
  '''[ MA - Galactose - H2O ]\n
  Fragment for adducted molecular ion, with loss of galactose, leaving -OH\n
  galactose NL common to DGDG lipids'''
  def MZ(self):
    return super().MZ() - (180.063388/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'C':6, 'H':12,'O':6})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-C6H12O6]')
    return comment

class MA_s_Gal_H2O(MA):
  '''[ MA - Galactose - H2O ]\n
  Fragment for adducted molecular ion, with loss of galactose, leaving -OH\n
  galactose NL common to DGDG lipids'''
  def MZ(self):
    return super().MZ() - ((180.063388+Masses['H2O'])/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'C':6, 'H':14,'O':7})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-C6H14O7]')
    return comment

class MA_H2O_s_Gal(MA):
  '''[ MA - Galactose - H2O ]\n
  Fragment for adducted molecular ion, with loss of galactose, leaving -OH\n
  galactose NL common to DGDG lipids'''
  def MZ(self):
    return super().MZ() - (162.052823/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'C':6, 'H':10,'O':5})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-C6H10O5]')
    return comment

class MA_s_2Gal(MA_s_Gal):
  '''[ MA - Galactose - H2O ]\n
  Fragment for adducted molecular ion, with loss of galactose, leaving -OH\n
  galactose NL common to DGDG lipids'''
  def MZ(self):
    return super().MZ() - (180.063388/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'C':6, 'H':12,'O':6})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-C6H12O6]')
    return comment

class MA_H2O_s_2Gal(MA_s_Gal):
  '''[ MA - Galactose - H2O ]\n
  Fragment for adducted molecular ion, with loss of galactose, leaving -OH\n
  galactose NL common to DGDG lipids'''
  def MZ(self):
    return super().MZ() - (162.052823/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'C':6, 'H':10,'O':5})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-C6H10O5]')
    return comment

class MA_2H2O_s_2Gal(MA_H2O_s_Gal):
  '''[ MA - Galactose - H2O ]\n
  Fragment for adducted molecular ion, with loss of galactose, leaving -OH\n
  galactose NL common to DGDG lipids'''
  def MZ(self):
    return super().MZ() - (162.052823/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'C':6, 'H':10,'O':5})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-C6H10O5]')
    return comment

# ~ # Fragments for PC+Na/Li

class MA_s_TMA(MA):
  '''[ MA - C3H9N ]\n
  Fragment for adducted molecular ion, with loss of trimethylamine\n
  Common for Phosphatidylcholines'''
  def MZ(self):
    return super().MZ() - (Masses['TMA']/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'C':3, 'H':9 ,'N':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-C3H9N]')
    return comment

class MA_s_AZD(MA):
  '''[ MA - C2H5N ]\n
  Fragment for adducted molecular ion, with loss of aziridine\n
  Common for Phosphatidylcholines'''
  def MZ(self):
    return super().MZ() - (Masses['AZD']/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'C':2, 'H':5 ,'N':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-C2H5N]')
    return comment

class MA_s_TMA_H2O(MA_s_TMA):
  '''[ MA - C3H9N ]\n
  Fragment for adducted molecular ion, with loss of trimethylamine\n
  Common for Phosphatidylcholines'''
  def MZ(self):
    return super().MZ() - (Masses['H2O']/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-H2O]')
    return comment

class MA_s_AZD_H2O(MA_s_AZD):
  '''[ MA - C2H5N ]\n
  Fragment for adducted molecular ion, with loss of aziridine\n
  Common for Phosphatidylcholines'''
  def MZ(self):
    return super().MZ() - (Masses['H2O']/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-H2O]')
    return comment

# ~ # Fragments for Ceramides

class MA_s_MeOH(MA):
  '''[ MA - CH4O ]\n
  Fragment for adducted molecular ion, with loss of methanol\n
  Common for ceramides in -ve ESI'''
  def MZ(self):         #    MeOH mass  
    return super().MZ() - (32.026214784/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'C':1, 'H':4 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-CH4O]')
    return comment

class MA_s_CH2O(MA):
  '''[ MA - CH2O ]\n
  Fragment for adducted molecular ion, with loss of formaldehyde\n
  Common for ceramides in -ve ESI'''
  def MZ(self):         #  Formaldehyde mass  
    return super().MZ() - (30.010564684/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'C':1, 'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-CH2O]')
    return comment

class MA_s_CH2O_H2O(MA_s_CH2O):
  '''[ MA - CH2O ]\n
  Fragment for adducted molecular ion, with loss of formaldehyde\n
  Common for ceramides in -ve ESI'''
  def MZ(self):          
    return super().MZ() - (Masses['H2O']/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'H':2 ,'O':1})
    return formula  
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-H2O]')
    return comment  

# ~ # ~ # ~ # [M +/- adduct] - fatty acid

def MA_s_FA(lipid, adduct, intensity):
  '''[ MA - (ROOH) ]\n
  Fragment for adducted molecular ion, with loss of a free-fatty acid\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl', 'Ether', 'Vinyl']:
      yield MA_s_FAx(lipid, adduct, intensity, MA_s_FA, tail)

class MA_s_FAx(MA):
  '''[ MA - (ROOH) ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)
  def MZ(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = 1.006277
    else: x = 0 # Rearrangement on tail loss exchanges a hydrogen
    # https://doi.org/10.1016/j.algal.2016.05.016
    return super().MZ() - ((self.tail.mass+x)/abs(Masses[self.adduct][2]))
  def Formula(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = {'H':1, 'D':-1}
    else: x = {} # Rearrangement on tail loss exchanges a hydrogen
    formula = super().Formula()
    formula.update(x)
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn1(MA):
  '''[ MA - (ROOH) ] (sn1)\n
  Fragment for adducted molecular ion, with loss of sn1 as free-fatty acid'''
  def MZ(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = 1.006277
    else: x = 0 # Rearrangement on tail loss exchanges a hydrogen
    return super().MZ() - ((self.lipid.tails[0].mass+x)/abs(Masses[self.adduct][2]))
  def Formula(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = {'H':1, 'D':-1}
    else: x = {} # Rearrangement on tail loss exchanges a hydrogen
    formula = super().Formula()
    formula.update(x)
    formula.subtract(self.lipid.tails[0].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn2(MA):
  '''[ MA - (ROOH) ] (sn2)\n
  Fragment for adducted molecular ion, with loss of sn2 as free-fatty acid'''
  def MZ(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = 1.006277
    else: x = 0 # Rearrangement on tail loss exchanges a hydrogen
    return super().MZ() - ((self.lipid.tails[1].mass+x)/abs(Masses[self.adduct][2]))
  def Formula(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = {'H':1, 'D':-1}
    else: x = {} # Rearrangement on tail loss exchanges a hydrogen
    formula = super().Formula()
    formula.update(x)
    formula.subtract(self.lipid.tails[1].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

# ~ #

class MA_s_allFA(MA):
  '''[ MA - (ROOH) ] (ALL)\n
  Fragment for adducted molecular ion, with loss of ALL fatty acid tails'''
  def MZ(self):
    mass = super().MZ()
    for tail in self.lipid.tails:
      if tail.type not in ['Headgroup', 'HeadTail']:
        mass -= (tail.mass/abs(Masses[self.adduct][2]))
    return mass
  def Formula(self):
    formula = super().Formula()
    for tail in self.lipid.tails:
      if tail.type not in ['Headgroup', 'HeadTail']:
        formula.subtract(tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)-(ROOH)')
    comment += ' ('
    comment += ', '.join(tail.name for tail in self.lipid.tails if tail.type not in ['Headgroup', 'HeadTail'])
    comment += ')'
    return comment  

# ~ #

def MA_s_FA_H2O(lipid, adduct, intensity):
  '''[ MA - (ROOH) - H2O ]\n
  Fragment for adducted molecular ion with loss of a free fatty acid AND water\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MA_s_FAx_H2O(lipid, adduct, intensity, MA_s_FA_H2O, tail)

class MA_s_FAx_H2O(MA_s_H2O):
  '''[ MA - (ROOH) - H2O ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - (self.tail.mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn1_H2O(MA_s_H2O):
  '''[ MA - (ROOH) - H2O ] (sn1)\n
  Fragment for adducted molecular ion with loss of sn1 as a free fatty acid AND water'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[0].mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn2_H2O(MA_s_H2O):
  '''[ MA - (ROOH) - H2O ] (sn2)\n
  Fragment for adducted molecular ion with loss of sn2 as a free fatty acid AND water'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[1].mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

# ~ #

def MA_s_FAk(lipid, adduct, intensity):
  '''[ MA - (R=O) ]\n
  Fragment for adducted molecular ion with loss of a fatty acid ketene\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MA_s_FAkx(lipid, adduct, intensity, MA_s_FAk, tail)

class MA_s_FAkx(MA):
  '''[ MA - (R=O) ]\n
  Fatty acid Ketene\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    mz = super().MZ()
    mz += (Masses['H2O']/abs(Masses[self.adduct][2])) 
    mz -= (self.tail.mass/abs(Masses[self.adduct][2]))
    return mz
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn1k(MA):
  '''[ MA - (R=O) ] (sn1)\n
  Fragment for adducted molecular ion with loss of a fatty acid ketene'''
  def MZ(self):
    return super().MZ() - ((self.lipid.tails[0].mass-Masses['H2O'])/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn2k(MA):
  '''[ MA - (R=O) ] (sn2)\n
  Fragment for adducted molecular ion with loss of a fatty acid ketene'''
  def MZ(self):
    return super().MZ() - ((self.lipid.tails[1].mass-Masses['H2O'])/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment  

# ~ #

class MA_s_allFAk(MA):
  '''[ MA - (R=O) ] (ALL)\n
  Fragment for adducted molecular ion, with loss of ALL fatty acid tails'''
  def MZ(self):
    mass = super().MZ()
    for tail in self.lipid.tails:
      if tail.type not in ['Headgroup', 'HeadTail']:
        mass -= (tail.mass-Masses['H2O']/abs(Masses[self.adduct][2]))
    return mass
  def Formula(self):
    formula = super().Formula()
    for tail in self.lipid.tails:
      if tail.type not in ['Headgroup', 'HeadTail']:
        formula.subtract(tail.formula)
        formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)-(R=O)')
    comment += ' ('+', '.join(tail.name for tail in self.lipid.tails if tail.type not in ['Headgroup', 'HeadTail'])+')'
    return comment

# ~ #

def MA_s_FA_PO3(lipid, adduct, intensity):
  '''[ MA - (ROOH) - PO3 ]\n
  Fragment for adducted molecular ion with loss of a free fatty acid AND phosphite\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MA_s_FA_PO3x(lipid, adduct, intensity, MA_s_FA_PO3, tail)

class MA_s_FA_PO3x(MA_s_PO3):
  '''[ MA - (ROOH) - PO3 ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - (self.tail.mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn1_PO3(MA_s_PO3):
  '''[ MA - (ROOH) - PO3 ] (sn1)\n
  Fragment for adducted molecular ion with loss of sn1 as a free fatty acid AND phosphite'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[0].mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn2_PO3(MA_s_PO3):
  '''[ MA - (ROOH) - PO3 ] (sn2)\n
  Fragment for adducted molecular ion with loss of sn1 as a free fatty acid AND phosphite'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[1].mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

# ~ #

def MA_s_FAk_PO3(lipid, adduct, intensity):
  '''[ MA - (R=O) - PO3 ]\n
  Fragment for adducted molecular ion with loss of a fatty acid ketene AND phosphite\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MA_s_FAk_PO3x(lipid, adduct, intensity, MA_s_FAk_PO3, tail)

class MA_s_FAk_PO3x(MA_s_PO3):
  '''[ MA - (R=O) - PO3 ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - ((self.tail.mass-Masses['H2O'])/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn1k_PO3(MA_s_PO3):
  '''[ MA - (R=O) - PO3 ] (sn1)\n
  Fragment for adducted molecular ion with loss of sn1 as a fatty acid ketene AND phosphite'''
  def MZ(self):
    return super().MZ() - ((self.lipid.tails[0].mass-Masses['H2O'])/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn2k_PO3(MA_s_PO3):
  '''[ MA - (R=O) - PO3 ] (sn2)\n
  Fragment for adducted molecular ion with loss of sn1 as a fatty acid ketene AND phosphite'''
  def MZ(self):
    return super().MZ() - ((self.lipid.tails[1].mass-Masses['H2O'])/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment  

# ~ # Fragments for PE / PC+Na/Li

def MA_s_FA_TMA(lipid, adduct, intensity):
  '''[ MA - (ROOH) - C3H9N ]\n
  Fragment for adducted molecular ion with loss of a free fatty acid AND trimethylamine\n
  Common for Phosphatidylcholines\n
  Method used to generate multiple objects'''  
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MA_s_FA_TMAx(lipid, adduct, intensity, MA_s_FA_TMA, tail)

class MA_s_FA_TMAx(MA_s_TMA):
  '''[ MA - (ROOH) - C3H9N ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - (self.tail.mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn1_TMA(MA_s_TMA):
  '''[ MA - (ROOH) - C3H9N ]\n
  Fragment for adducted molecular ion with loss of a free fatty acid AND trimethylamine\n
  Common for Phosphatidylcholines'''  
  def MZ(self):
    return super().MZ() - (self.lipid.tails[0].mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn2_TMA(MA_s_TMA):
  '''[ MA - (ROOH) - C3H9N ]\n
  Fragment for adducted molecular ion with loss of a free fatty acid AND trimethylamine\n
  Common for Phosphatidylcholines'''  
  def MZ(self):
    return super().MZ() - (self.lipid.tails[1].mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

# ~ #

def MA_s_FAk_TMA(lipid, adduct, intensity):
  '''[ MA - (R=O) - C3H9N ]\n
  Fragment for adducted molecular ion with loss of a fatty acid ketene AND trimethylamine\n
  Common for Phosphatidylcholines\n
  Method used to generate multiple objects'''  
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MA_s_FAk_TMAx(lipid, adduct, intensity, MA_s_FAk_TMA, tail)

class MA_s_FAk_TMAx(MA_s_TMA):
  '''[ MA - (R=O) - C3H9N ]\n
  Do not use this class, intended for use in loop''' 
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - ((self.tail.mass-Masses['H2O'])/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn1k_TMA(MA_s_TMA):  # [M+-H-FAk-TMA]+-
  '''[ MA - (R=O) - C3H9N ] (sn1)\n
  Fragment for adducted molecular ion with loss of a fatty acid ketene AND trimethylamine\n
  Common for Phosphatidylcholines'''  
  def MZ(self):
    return super().MZ() - ((self.lipid.tails[0].mass-Masses['H2O'])/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn2k_TMA(MA_s_TMA):  # [M+-H-FAk-TMA]+-
  '''[ MA - (R=O) - C3H9N ] (sn2)\n
  Fragment for adducted molecular ion with loss of a fatty acid ketene AND trimethylamine\n
  Common for Phosphatidylcholines'''  
  def MZ(self):
    return super().MZ() - ((self.lipid.tails[1].mass-Masses['H2O'])/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment  
    
# ~ #

def MA_s_FA_AZD(lipid, adduct, intensity):  # [M+-H-FA-AZD]+-
  '''[ MA - (ROOH) - C2H5N ]\n
  Fragment for adducted molecular ion with loss of a free fatty acid AND aziridine\n
  Common for Phosphatidylcholines\n
  Method used to generate multiple objects'''  
  for tail in lipid.tails:       # FA = fatty acid
    if tail.type in ['Acyl']:
      yield MA_s_FA_AZDx(lipid, adduct, intensity, MA_s_FA_AZD, tail)

class MA_s_FA_AZDx(MA_s_AZD):  # [M+-H-FA-AZD]+-
  '''[ MA - (ROOH) - C2H5N ]\n
  Do not use this class, intended for use in loop''' 
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - (self.tail.mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn1_AZD(MA_s_AZD):  # [M+-H-FA-AZD]+-
  '''[ MA - (ROOH) - C2H5N ]\n
  Fragment for adducted molecular ion with loss of a free fatty acid AND aziridine\n
  Common for Phosphatidylcholines'''  
  def MZ(self):
    return super().MZ() - (self.lipid.tails[0].mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn2_AZD(MA_s_AZD):  # [M+-H-FA-AZD]+-
  '''[ MA - (ROOH) - C2H5N ]\n
  Fragment for adducted molecular ion with loss of a free fatty acid AND aziridine\n
  Common for Phosphatidylcholines'''  
  def MZ(self):
    return super().MZ() - (self.lipid.tails[1].mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

# ~ #

def MA_s_FAk_AZD(lipid, adduct, intensity):  # [M+-H-FAk-AZD]+-
  '''[ MA - (R=O) - C2H5N ]\n
  Fragment for adducted molecular ion with loss of a fatty acid ketene AND aziridine\n
  Common for Phosphatidylcholines\n
  Method used to generate multiple objects'''  
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MA_s_FAk_AZDx(lipid, adduct, intensity, MA_s_FAk_AZD, tail)

class MA_s_FAk_AZDx(MA_s_AZD):  # [M+-H-FAk-AZD]+-
  '''[ MA - (R=O) - C2H5N ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - ((self.tail.mass-Masses['H2O'])/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn1k_AZD(MA_s_AZD):  # [M+-H-FAk-AZD]+-
  '''[ MA - (R=O) - C2H5N ] (sn1)\n
  Fragment for adducted molecular ion with loss of a fatty acid ketene AND aziridine\n
  Common for Phosphatidylcholines'''  
  def MZ(self):
    return super().MZ() - ((self.lipid.tails[0].mass-Masses['H2O'])/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn2k_AZD(MA_s_AZD):  # [M+-H-FAk-AZD]+-
  '''[ MA - (R=O) - C2H5N ] (sn2)\n
  Fragment for adducted molecular ion with loss of a fatty acid ketene AND aziridine\n
  Common for Phosphatidylcholines'''  
  def MZ(self):
    return super().MZ() - ((self.lipid.tails[1].mass-Masses['H2O'])/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment  

# ~ #

def MA_s_FA_Gal(lipid, adduct, intensity):
  '''[ MA - (ROOH) ]\n
  Fragment for adducted molecular ion, with loss of a free-fatty acid\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MA_s_FA_Galx(lipid, adduct, intensity, MA_s_FA_Gal, tail)

class MA_s_FA_Galx(MA_s_Gal):
  '''[ MA - (ROOH) ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)
  def MZ(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = 1.006277
    else: x = 0 # Rearrangement on tail loss exchanges a hydrogen
    # https://doi.org/10.1016/j.algal.2016.05.016
    return super().MZ() - ((self.tail.mass+x)/abs(Masses[self.adduct][2]))
  def Formula(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = {'H':1, 'D':-1}
    else: x = {} # Rearrangement on tail loss exchanges a hydrogen
    formula = super().Formula()
    formula.update(x)
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn1_Gal(MA_s_Gal):
  '''[ MA - (ROOH) ] (sn1)\n
  Fragment for adducted molecular ion, with loss of sn1 as free-fatty acid'''
  def MZ(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = 1.006277
    else: x = 0 # Rearrangement on tail loss exchanges a hydrogen
    return super().MZ() - ((self.lipid.tails[0].mass+x)/abs(Masses[self.adduct][2]))
  def Formula(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = {'H':1, 'D':-1}
    else: x = {} # Rearrangement on tail loss exchanges a hydrogen
    formula = super().Formula()
    formula.update(x)
    formula.subtract(self.lipid.tails[0].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn2_Gal(MA_s_Gal):
  '''[ MA - (ROOH) ] (sn2)\n
  Fragment for adducted molecular ion, with loss of sn2 as free-fatty acid'''
  def MZ(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = 1.006277
    else: x = 0 # Rearrangement on tail loss exchanges a hydrogen
    return super().MZ() - ((self.lipid.tails[1].mass+x)/abs(Masses[self.adduct][2]))
  def Formula(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = {'H':1, 'D':-1}
    else: x = {} # Rearrangement on tail loss exchanges a hydrogen
    formula = super().Formula()
    formula.update(x)
    formula.subtract(self.lipid.tails[1].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

def MA_s_FA_Gal_H2O(lipid, adduct, intensity):
  '''[ MA - (ROOH) ]\n
  Fragment for adducted molecular ion, with loss of a free-fatty acid\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MA_s_FA_Gal_H2Ox(lipid, adduct, intensity, MA_s_FA_Gal_H2O, tail)

class MA_s_FA_Gal_H2Ox(MA_H2O_s_Gal):
  '''[ MA - (ROOH) ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)
  def MZ(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = 1.006277
    else: x = 0 # Rearrangement on tail loss exchanges a hydrogen
    # https://doi.org/10.1016/j.algal.2016.05.016
    return super().MZ() - ((self.tail.mass+x)/abs(Masses[self.adduct][2]))
  def Formula(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = {'H':1, 'D':-1}
    else: x = {} # Rearrangement on tail loss exchanges a hydrogen
    formula = super().Formula()
    formula.update(x)
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn1_Gal_H2O(MA_H2O_s_Gal):
  '''[ MA - (ROOH) ] (sn1)\n
  Fragment for adducted molecular ion, with loss of sn1 as free-fatty acid'''
  def MZ(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = 1.006277
    else: x = 0 # Rearrangement on tail loss exchanges a hydrogen
    return super().MZ() - ((self.lipid.tails[0].mass+x)/abs(Masses[self.adduct][2]))
  def Formula(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = {'H':1, 'D':-1}
    else: x = {} # Rearrangement on tail loss exchanges a hydrogen
    formula = super().Formula()
    formula.update(x)
    formula.subtract(self.lipid.tails[0].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_sn2_Gal_H2O(MA_H2O_s_Gal):
  '''[ MA - (ROOH) ] (sn2)\n
  Fragment for adducted molecular ion, with loss of sn2 as free-fatty acid'''
  def MZ(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = 1.006277
    else: x = 0 # Rearrangement on tail loss exchanges a hydrogen
    return super().MZ() - ((self.lipid.tails[1].mass+x)/abs(Masses[self.adduct][2]))
  def Formula(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = {'H':1, 'D':-1}
    else: x = {} # Rearrangement on tail loss exchanges a hydrogen
    formula = super().Formula()
    formula.update(x)
    formula.subtract(self.lipid.tails[1].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # [M +/- H]

class MH(Fragment):
  '''[ M-H ] or [ M+H ]\n
  Fragment for molecular ion (de)protonation'''
  def MZ(self):
    if Masses[self.adduct][1] == 'Positive':
      return self.lipid.mass + Masses['H+']
    else:
      return self.lipid.mass - Masses['H+']
  def Formula(self):
    formula = Counter(self.lipid.formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.update({'H':1})
    else:
      formula.subtract({'H':1})
    return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1
  def Comment(self):
    if Masses[self.adduct][1] == 'Positive':
      return '[M+H]+'
    else:
      return '[M-H]-'

class MH_s_H2O(MH):
  '''[ M(+/-)H - H2O ]\n
  Fragment for (de)protonated molecular ion with loss of water'''
  def MZ(self):
    return super().MZ() - Masses['H2O']
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'H':2,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-H2O]')
    return comment  

class MH_s_2H2O(MH):
  '''[ M(+/-)H - 2H2O ]\n
  Fragment for (de)protonated molecular ion with loss of two waters'''
  def MZ(self):
    return super().MZ() - 2*Masses['H2O']
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'H':4,'O':2})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-2H2O]')
    return comment  

class MH_s_PO3(MH):
  '''[ M(+/-)H - PO3 ]\n
  Fragment for (de)protonated molecular ion with loss of phosphite'''
  def MZ(self):
    return super().MZ() - Masses['PO3H']
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'H':1, 'P':1,'O':3})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-HPO3]')
    return comment  

class MH_s_PO4(MH):
  '''[ M(+/-)H - PO4 ]\n
  Fragment for (de)protonated molecular ion with loss of phosphate'''
  def MZ(self):
    return super().MZ() - Masses['PO4H3']
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'H':3, 'P':1,'O':4})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-H3PO4]')
    return comment  

class MH_s_PO4_H2O(MH_s_PO4):
  '''[ M(+/-)H - PO4 -H2O]\n
  Fragment for (de)protonated molecular ion with loss of phosphate'''
  def MZ(self):
    return super().MZ() - Masses['H2O']
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'O':1, 'H':2})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-H2O]')
    return comment  

# ~ # Fragments for DGDG

class MH_s_Gal_H2O(MH):
  '''[ MA - Galactose - H2O ]\n
  Fragment for (de)protonated molecular ion, with loss of galactose, leaving -OH\n
  galactose NL common to DGDG lipids'''
  def MZ(self):
    return super().MZ() - (162.052823/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'C':6, 'H':10,'O':5})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-C6H10O5]')
    return comment

# ~ # Fragments for PE / PC+Na/Li

class MH_s_TMA(MH):
  '''[ M(+/-)H - C3H9N ]\n
  Fragment for (de)protonated molecular ion with loss of trimethylamine'''
  def MZ(self):
      return super().MZ() - Masses['TMA']
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'C':3, 'H':9 ,'N':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-C3H9N]')
    return comment

class MH_s_AZD(MH):
  '''[ M(+/-)H - C2H5N ]\n
  Fragment for (de)protonated molecular ion with loss of aziridine'''
  def MZ(self):
      return super().MZ() - Masses['AZD']
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'C':2, 'H':5 ,'N':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-C2H5N]')
    return comment

class MH_s_C3H9O6P(MH):
  '''[ M(+/-)H - C3H9O6P ]\n
  Fragment for (de)protonated molecular ion with loss of phosphoglycerol'''
  def MZ(self):
      return super().MZ() - 172.013674
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'C':3, 'H':9 ,'O':6, 'P':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-C3H9O6P]')
    return comment

# ~ # ~ # ~ # [M +/- H] - fatty acid

def MH_s_FA(lipid, adduct, intensity):
  '''[ M(+/-)H - (ROOH) ]\n
  Fragment for (de)protonated molecular ion, with loss of a free-fatty acid\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MH_s_FAx(lipid, adduct, intensity, MH_s_FA, tail)

class MH_s_FAx(MH):
  '''[ M(+/-)H - (ROOH) ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = 1.006277
    else: x = 0 # Rearrangement on tail loss exchanges a hydrogen
    return super().MZ() - (self.tail.mass+x)
  def Formula(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = {'H':1, 'D':-1}
    else: x = {} # Rearrangement on tail loss exchanges a hydrogen
    formula = super().Formula()
    formula.update(x)
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

class MH_s_sn1(MH):
  '''[ M(+/-)H - (ROOH) ] (sn1)\n
  Fragment for (de)protonated molecular ion, with loss of sn1 as free-fatty acid'''
  def MZ(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = 1.006277
    else: x = 0 # Rearrangement on tail loss exchanges a hydrogen
    return super().MZ() - (self.lipid.tails[0].mass+x)
  def Formula(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = {'H':1, 'D':-1}
    else: x = {} # Rearrangement on tail loss exchanges a hydrogen
    formula = super().Formula()
    formula.update(x)
    formula.subtract(self.lipid.tails[0].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

class MH_s_sn2(MH):
  '''[ M(+/-)H - (ROOH) ] (sn2)\n
  Fragment for (de)protonated molecular ion, with loss of sn2 as free-fatty acid'''
  def MZ(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = 1.006277
    else: x = 0 # Rearrangement on tail loss exchanges a hydrogen
    return super().MZ() - (self.lipid.tails[1].mass+x)
  def Formula(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = {'H':1, 'D':-1}
    else: x = {} # Rearrangement on tail loss exchanges a hydrogen
    formula = super().Formula()
    formula.update(x)
    formula.subtract(self.lipid.tails[1].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

class MH_s_sn3(MH):
  '''[ M(+/-)H - (ROOH) ] (sn3)\n
  Fragment for (de)protonated molecular ion, with loss of sn2 as free-fatty acid'''
  def MZ(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = 1.006277
    else: x = 0 # Rearrangement on tail loss exchanges a hydrogen
    return super().MZ() - (self.lipid.tails[2].mass+x)
  def Formula(self):
    # If the tail is fully deuterated: special case
    if self.tail.formula['D'] and self.tail.formula['H'] == 1: x = {'H':1, 'D':-1}
    else: x = {} # Rearrangement on tail loss exchanges a hydrogen
    formula = super().Formula()
    formula.update(x)
    formula.subtract(self.lipid.tails[2].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

# ~ # 

def MH_s_FA_H2O(lipid, adduct, intensity):
  '''[ M(+/-)H - (ROOH) - H2O ]\n
  Fragment for (de)protonated molecular ion with loss of a free fatty acid AND water\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MH_s_FA_H2Ox(lipid, adduct, intensity, MH_s_FA_H2O, tail)

class MH_s_FA_H2Ox(MH_s_H2O):
  '''[ M(+/-)H - (ROOH) - H2O ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - self.tail.mass
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

class MH_s_sn1_H2O(MH_s_H2O):
  '''[ M(+/-)H - (ROOH) - H2O ] (sn1)\n
  Fragment for (de)protonated molecular ion with loss of sn1 as a free fatty acid AND water'''
  def MZ(self):
    return super().MZ() - self.lipid.tails[0].mass
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

class MH_s_sn2_H2O(MH_s_H2O):
  '''[ M(+/-)H - (ROOH) - H2O ] (sn2)\n
  Fragment for (de)protonated molecular ion with loss of sn2 as a free fatty acid AND water'''
  def MZ(self):
    return super().MZ() - self.lipid.tails[1].mass
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

# ~ # 

def MH_s_FAk(lipid, adduct, intensity):
  '''[ M(+/-)H - (R=O) ]\n
  Fragment for (de)protonated molecular ion with loss of a fatty acid ketene\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MH_s_FAkx(lipid, adduct, intensity, MH_s_FAk, tail)

class MH_s_FAkx(MH):
  '''[ M(+/-)H - (ROOH) ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - (self.tail.mass-Masses['H2O'])
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

class MH_s_sn1k(MH):
  '''[ M(+/-)H - (R=O) ] (sn1)\n
  Fragment for (de)protonated molecular ion with loss of a fatty acid ketene'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[0].mass-Masses['H2O'])
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

class MH_s_sn2k(MH):
  '''[ M(+/-)H - (R=O) ] (sn2)\n
  Fragment for (de)protonated molecular ion with loss of a fatty acid ketene'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[1].mass-Masses['H2O'])
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

class MH_s_allFAk(MH):
  '''[ M(+/-)H - (R=O) ] (ALL)\n
  Fragment for (de)protonated molecular ion, with loss of ALL fatty acids'''
  def MZ(self):
    mass = super().MZ()
    for tail in self.lipid.tails:
      if tail.type != 'Headgroup':
        mass -= (tail.mass-Masses['H2O'])
    return mass
  def Formula(self):
    formula = super().Formula()
    for tail in self.lipid.tails:
      if tail.type != 'Headgroup':
        formula.subtract(tail.formula)
        formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-(R=O)-(R=O)]')
    return comment

# ~ # 

def MH_s_FA_PO3(lipid, adduct, intensity):
  '''[ M(+/-)H - (ROOH) - PO3 ]\n
  Fragment for (de)protonated molecular ion with loss of a free fatty acid AND phosphite\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MH_s_FA_PO3x(lipid, adduct, intensity, MH_s_FA_PO3, tail)

class MH_s_FA_PO3x(MH_s_PO3):  # [M+-H-FA-PO3]+-
  '''[ M(+/-)H - (ROOH) - PO3 ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - self.tail.mass
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

class MH_s_sn1_PO3(MH_s_PO3):
  '''[ M(+/-)H - (ROOH) - PO3 ] (sn1)\n
  Fragment for (de)protonated molecular ion with loss of sn1 as a free fatty acid AND phosphite'''
  def MZ(self):
    return super().MZ() - self.lipid.tails[0].mass
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

class MH_s_sn2_PO3(MH_s_PO3):
  '''[ M(+/-)H - (ROOH) - PO3 ] (sn2)\n
  Fragment for (de)protonated molecular ion with loss of sn2 as a free fatty acid AND phosphite'''
  def MZ(self):
    return super().MZ() - self.lipid.tails[1].mass
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

# ~ # 

def MH_s_FAk_PO3(lipid, adduct, intensity):
  '''[ M(+/-)H - (R=O) - PO3 ]\n
  Fragment for (de)protonated molecular ion with loss of a fatty acid ketene AND phosphite\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MH_s_FAk_PO3x(lipid, adduct, intensity, MH_s_FAk_PO3, tail)

class MH_s_FAk_PO3x(MH_s_PO3):
  '''[ M(+/-)H - (R=O) - PO3 ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - (self.tail.mass-Masses['H2O'])
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

class MH_s_sn1k_PO3(MH_s_PO3):
  '''[ M(+/-)H - (R=O) - PO3 ] (sn1)\n
  Fragment for (de)protonated molecular ion with loss of sn1 as a fatty acid ketene AND phosphite'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[0].mass-Masses['H2O'])
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

class MH_s_sn2k_PO3(MH_s_PO3):
  '''[ M(+/-)H - (R=O) - PO3 ] (sn2)\n
  Fragment for (de)protonated molecular ion with loss of sn2 as a fatty acid ketene AND phosphite'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[1].mass-Masses['H2O'])
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # [M +/- 2H]

class M2H(Fragment):
  '''[ M-2H ] or [ M+2H ]\n
  Fragment for molecular ion double (de)protonation'''
  def MZ(self):
    if Masses[self.adduct][1] == 'Positive':
      return (self.lipid.mass + 2*Masses['H+'])/2
    else:
      return (self.lipid.mass - 2*Masses['H+'])/2
  def Formula(self):
    formula = Counter(self.lipid.formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.update({'H':2})
    else:
      formula.subtract({'H':2})
    return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 2
    else:
      return -2
  def Comment(self):
    if Masses[self.adduct][1] == 'Positive':
      return '[M+2H]2+'
    else:
      return '[M-2H]2-'

class M2H_s_H2O(M2H):
  '''[ M(+/-)2H - H2O ]\n
  Fragment for double (de)protonated molecular ion with loss of water'''
  def MZ(self):
    return super().MZ() - (Masses['H2O'])/2
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'H':2,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-H2O]')
    return comment  

class M2H_s_2H2O(M2H):
  '''[ M(+/-)2H - 2H2O ]\n
  Fragment for double (de)protonated molecular ion with loss of two waters'''
  def MZ(self):
    return super().MZ() - Masses['H2O'] # i.e. (2*Masses['H2O'])/2
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'H':4,'O':2})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-2H2O]')
    return comment  

# ~ # ~ # ~ # [M +/- 2H] - fatty acid

def M2H_s_FA(lipid, adduct, intensity):
  '''[ M(+/-)2H - (ROOH) ]\n
  Fragment for double (de)protonated molecular ion, with loss of a free-fatty acid\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield M2H_s_FAx(lipid, adduct, intensity, M2H_s_FA, tail)

class M2H_s_FAx(M2H):
  '''[ M(+/-)2H - (ROOH) ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - (self.tail.mass/2)
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

class M2H_s_sn1(M2H):
  '''[ M(+/-)2H - (ROOH) ] (sn1)\n
  Fragment for double (de)protonated molecular ion, with loss of sn1 as free-fatty acid'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[0].mass/2)
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

class M2H_s_sn2(M2H):
  '''[ M(+/-)2H - (ROOH) ] (sn2)\n
  Fragment for double (de)protonated molecular ion, with loss of sn2 as free-fatty acid'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[1].mass/2)
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

# ~ # 

def M2H_s_FAk(lipid, adduct, intensity):
  '''[ M(+/-)2H - (R=O) ]\n
  Fragment for double (de)protonated molecular ion with loss of a fatty acid ketene\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield M2H_s_FAkx(lipid, adduct, intensity, M2H_s_FAk, tail)

class M2H_s_FAkx(M2H):
  '''[ M(+/-)2H - (ROOH) ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - (self.tail.mass-Masses['H2O'])/2
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

class M2H_s_sn1k(M2H):
  '''[ M(+/-)2H - (R=O) ] (sn1)\n
  Fragment for double (de)protonated molecular ion with loss of a fatty acid ketene'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[0].mass-Masses['H2O'])/2
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

class M2H_s_sn2k(M2H):
  '''[ M(+/-)2H - (R=O) ] (sn2)\n
  Fragment for double (de)protonated molecular ion with loss of a fatty acid ketene'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[1].mass-Masses['H2O'])/2
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # Free fatty acids and fatty acid ketenes

def FAH(lipid, adduct, intensity):
  '''[ FA - H ]\n
  Fragment for a deprotonated free fatty acid\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl', 'HeadTail']:
      yield FAHx(lipid, adduct, intensity, FAH, tail)

class FAHx(Fragment):
  '''[ FA (+ / -) H ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    if Masses[self.adduct][1] == 'Positive':
      return self.tail.mass + Masses['H+']
    else:
      return self.tail.mass - Masses['H+']
  def Formula(self):
    formula = Counter(self.tail.formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.update({'H':1})
    else:
      formula.subtract({'H':1})
    return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1
  def Comment(self):
    if Masses[self.adduct][1] == 'Positive':
      return '[(ROOH)+H]+ ('+self.tail.name+')'
    else:
      return '[RCOO]- ('+self.tail.name+')'

class sn1(Fragment):
  '''[ FA - H ] (sn1)\n
  Fragment for a deprotonated free fatty acid'''
  def MZ(self):
    if Masses[self.adduct][1] == 'Positive':
      return self.lipid.tails[0].mass + Masses['H+']
    else:
      return self.lipid.tails[0].mass - Masses['H+']
  def Formula(self):
    formula = Counter(self.lipid.tails[0].formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.update({'H':1})
    else:
      formula.subtract({'H':1})
    return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1
  def Comment(self):
    if Masses[self.adduct][1] == 'Positive':
      return '[(ROOH)+H]+ ('+self.tail.name+')'
    else:
      return '[RCOO]- ('+self.tail.name+')'

class sn2(Fragment):
  '''[ FA - H ] (sn2)\n
  Fragment for a deprotonated free fatty acid'''
  def MZ(self):
    if Masses[self.adduct][1] == 'Positive':
      return self.lipid.tails[1].mass + Masses['H+']
    else:
      return self.lipid.tails[1].mass - Masses['H+']
  def Formula(self):
    formula = Counter(self.lipid.tails[1].formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.update({'H':1})
    else:
      formula.subtract({'H':1})
    return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1
  def Comment(self):
    if Masses[self.adduct][1] == 'Positive':
      return '[(ROOH)+H]+ ('+self.tail.name+')'
    else:
      return '[RCOO]- ('+self.tail.name+')'

class sn3(Fragment):
  '''[ FA - H ] (sn3)\n
  Fragment for a deprotonated free fatty acid'''
  def MZ(self):
    if Masses[self.adduct][1] == 'Positive':
      return self.lipid.tails[2].mass + Masses['H+']
    else:
      return self.lipid.tails[2].mass - Masses['H+']
  def Formula(self):
    formula = Counter(self.lipid.tails[2].formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.update({'H':1})
    else:
      formula.subtract({'H':1})
    return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1
  def Comment(self):
    if Masses[self.adduct][1] == 'Positive':
      return '[(ROOH)+H]+ ('+self.tail.name+')'
    else:
      return '[RCOO]- ('+self.tail.name+')'

# ~ # 

def C3H5O4P_FA(lipid, adduct, intensity):
  '''[ C3H5O4P + FA - H ]\n
  Fragment for a deprotonated fatty acid\n
  attached to a dehydrated phosphoglycerol\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl', 'HeadTail']:
      yield C3H5O4P_FAx(lipid, adduct, intensity, C3H5O4P_FA, tail)

class C3H5O4P_FAx(Fragment):
  '''[ C3H6O4P + FA - H ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return self.tail.mass + 135.992544 - Masses['H+']
  def Formula(self):
    formula = Counter(self.tail.formula)
    formula.update({'C':3, 'H':4, 'O':4, 'P':1})
    return formula
  def Charge(self):
    return -1
  def Comment(self):
    return '[C3H5O4P+(RCOO)]- ('+self.tail.name+')'

def C3H7O5P_FA(lipid, adduct, intensity):
  '''[ C3H7O5P + FA - H ]\n
  Fragment for a deprotonated fatty acid\n
  attached to a dehydrated phosphoglycerol\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl', 'HeadTail']:
      yield C3H7O5P_FAx(lipid, adduct, intensity, C3H7O5P_FA, tail)

class C3H7O5P_FAx(Fragment):
  '''[ C3H7O5P + FA - H ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return self.tail.mass + 154.003109 - Masses['H+']
  def Formula(self):
    formula = Counter(self.tail.formula)
    formula.update({'C':3, 'H':6, 'O':5, 'P':1})
    return formula
  def Charge(self):
    return -1
  def Comment(self):
    return '[C3H7O5P+(RCOO)]- ('+self.tail.name+')'

class C11H19N2O2S_FA(Fragment):
  '''[ C11H18N2O2S + FA + H ]\n
  AcylCoA FattyAcid Fragment'''
  def MZ(self):
    return self.lipid.tails[0].mass + 243.116175
  def Formula(self):
    formula = Counter(self.lipid.tails[0].formula)
    formula.update({'C':11, 'H':19, 'N':2, 'O':2, 'S':1})
    return formula
  def Charge(self):
    return 1
  def Comment(self):
    return '[C11H18N2O2S+(RCOO)+H]+ ('+self.lipid.tails[0].name+')'

# ~ # 

def FAkH(lipid, adduct, intensity):
  '''[ FA - H2O +/- H ]\n
  Fragment for a (de)protonated  fatty acid ketene\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield FAkHx(lipid, adduct, intensity, FAkH, tail)

class FAkHx(Fragment):
  '''[ FA - H2O +/- H ]\n
  Do not use this class, intended for use in loop''' 
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
    self.tail = tail
    super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    if Masses[self.adduct][1] == 'Positive':
      return self.tail.mass - Masses['H2O'] + Masses['H+']
    else:
      return self.tail.mass - Masses['H2O'] - Masses['H+']
  def Formula(self):
    formula = Counter(self.tail.formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.subtract({'H':1, 'O':1})
      return formula
    else:
      formula.subtract({'H':3, 'O':1})
      return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1
  def Comment(self):
    if Masses[self.adduct][1] == 'Positive':
      return '[(R=O)+H]+ ('+self.tail.name+')'
    else:
      return '[(R=O)-H]- ('+self.tail.name+')'

class sn1k(Fragment):
  '''[ FA - H2O +/- H ] (sn1)\n
  Fragment for a deprotonated free fatty acid'''
  def MZ(self):
    if Masses[self.adduct][1] == 'Positive':
      return self.lipid.tails[0].mass - Masses['H2O'] + Masses['H+']
    else:
      return self.lipid.tails[0].mass - Masses['H2O'] - Masses['H+']
  def Formula(self):
    formula = Counter(self.lipid.tails[0].formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.subtract({'H':1, 'O':1})
      return formula
    else:
      formula.subtract({'H':3, 'O':1})
      return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1
  def Comment(self):
    if Masses[self.adduct][1] == 'Positive':
      return '[(R=O)+H]+ ('+self.tail.name+')'
    else:
      return '[(R=O)-H]- ('+self.tail.name+')'

class sn2k(Fragment):
  '''[ FA - H2O +/- H ] (sn2)\n
  Fragment for a deprotonated free fatty acid'''
  def MZ(self):
    if Masses[self.adduct][1] == 'Positive':
      return self.lipid.tails[1].mass - Masses['H2O'] + Masses['H+']
    else:
      return self.lipid.tails[1].mass - Masses['H2O'] - Masses['H+']
  def Formula(self):
    formula = Counter(self.lipid.tails[1].formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.subtract({'H':1, 'O':1})
      return formula
    else:
      formula.subtract({'H':3, 'O':1})
      return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1
  def Comment(self):
    if Masses[self.adduct][1] == 'Positive':
      return '[(R=O)+H]+ ('+self.tail.name+')'
    else:
      return '[(R=O)-H]- ('+self.tail.name+')'

class sn3k(Fragment):
  '''[ FA - H2O +/- H ] (sn3)\n
  Fragment for a deprotonated free fatty acid'''
  def MZ(self):
    if Masses[self.adduct][1] == 'Positive':
      return self.lipid.tails[2].mass - Masses['H2O'] + Masses['H+']
    else:
      return self.lipid.tails[2].mass - Masses['H2O'] - Masses['H+']
  def Formula(self):
    formula = Counter(self.lipid.tails[2].formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.subtract({'H':1, 'O':1})
      return formula
    else:
      formula.subtract({'H':3, 'O':1})
      return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1
  def Comment(self):
    if Masses[self.adduct][1] == 'Positive':
      return '[(R=O)+H]+ ('+self.tail.name+')'
    else:
      return '[(R=O)-H]- ('+self.tail.name+')'

def FAkA(lipid, adduct, intensity):
  '''[ FA - H2O + Adduct ]\n
  Fragment for adducted fatty acid ketene\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield FAkAx(lipid, adduct, intensity, FAkA, tail)

class FAkAx(Fragment):
  '''[ FA - H2O + Adduct ]\n
  Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return (self.tail.mass - Masses['H2O'] + Masses[self.adduct][0])/abs(Masses[self.adduct][2])
  def Formula(self):
    formula = Counter(self.tail.formula)
    formula.subtract({'H':2, 'O':1})
    formula.update(Masses[self.adduct][3])
    return formula
  def Charge(self):
    return Masses[self.adduct][2]
  def Comment(self):
    comment = self.adduct
    comment = comment.replace('M', '(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

class sn1kA(Fragment):
  '''[ FA - H2O + Adduct ] (sn1)\n
  Fragment for adducted fatty acid ketene'''
  def MZ(self):
    return (self.lipid.tails[0].mass - Masses['H2O'] + Masses[self.adduct][0])/abs(Masses[self.adduct][2])
  def Formula(self):
    formula = Counter(self.lipid.tails[0].formula)
    formula.subtract({'H':2, 'O':1})
    formula.update(Masses[self.adduct][3])
    return formula
  def Charge(self):
    return Masses[self.adduct][2]
  def Comment(self):
    comment = self.adduct
    comment = comment.replace('M', '(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

class sn2kA(Fragment):
  '''[ FA - H2O + Adduct ] (sn2)\n
  Fragment for adducted fatty acid ketene'''
  def MZ(self):
    return (self.lipid.tails[1].mass - Masses['H2O'] + Masses[self.adduct][0])/abs(Masses[self.adduct][2])
  def Formula(self):
    formula = Counter(self.lipid.tails[1].formula)
    formula.subtract({'H':2, 'O':1})
    formula.update(Masses[self.adduct][3])
    return formula
  def Charge(self):
    return Masses[self.adduct][2]
  def Comment(self):
    comment = self.adduct
    comment = comment.replace('M', '(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

class sn3kA(Fragment):
  '''[ FA - H2O + Adduct ] (sn3)\n
  Fragment for adducted fatty acid ketene'''
  def MZ(self):
    return (self.lipid.tails[2].mass - Masses['H2O'] + Masses[self.adduct][0])/abs(Masses[self.adduct][2])
  def Formula(self):
    formula = Counter(self.lipid.tails[2].formula)
    formula.subtract({'H':2, 'O':1})
    formula.update(Masses[self.adduct][3])
    return formula
  def Charge(self):
    return Masses[self.adduct][2]
  def Comment(self):
    comment = self.adduct
    comment = comment.replace('M', '(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

# ~ # ~ # Ceramide dC fragments

class Cer_B(Fragment):
  '''Base fragment\n
  [ Base (+/-) H+ ](+/-)'''
  def MZ(self):
    if Masses[self.adduct][1] == 'Positive': # Base is created without the ammonia, so needs to be added here
      return (self.lipid.tails[0].mass - Masses['H2O'] + Masses['NH3'] + Masses['H+'])
    else:
      return (self.lipid.tails[0].mass - Masses['H2O'] + Masses['NH3'] - Masses['H+'])
  def Formula(self):
    formula = Counter(self.lipid.tails[0].formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.update({'N':1, 'H':2, 'O':-1})
    else:
      formula.update({'N':1, 'O':-1})
    return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1  
  def Comment(self):
    return 'Ceramide fragment B'

class Cer_Bb(Cer_B):
  '''Base fragment\n
  [ Base -H2O (+/-) H+ ](+/-)'''
  def MZ(self):
    return super().MZ() - Masses['H2O']
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'H':2, 'O':1})
    return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1  
  def Comment(self):
    return 'Ceramide fragment Bb'

class Cer_C(Fragment):
  '''Base fragment\n
  [ Base - MeOH (+/-) H+ ](+/-)'''
  def MZ(self):
    if Masses[self.adduct][1] == 'Positive': # Base is created without the ammonia, so needs to be added here
      return (self.lipid.tails[0].mass - 32.026214784 - Masses['H2O'] + Masses['NH3'] + Masses['H+'])
    else:
      return (self.lipid.tails[0].mass - 32.026214784 - Masses['H2O'] + Masses['NH3'] - Masses['H+'])
  def Formula(self):
    formula = Counter(self.lipid.tails[0].formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.update({'N':1, 'C':-1, 'H':-2, 'O':-2})
    else:
      formula.update({'N':1, 'C':-1, 'H':-4, 'O':-2})
    return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1  
  def Comment(self):
    return 'Ceramide fragment C'

class Cer_D(Fragment):
  '''Base fragment\n
  [ Base - Me(OH)2 (+/-) H+ ](+/-)'''
  def MZ(self):
    if Masses[self.adduct][1] == 'Positive': # Base is created without the ammonia, so needs to be added here
      return (self.lipid.tails[0].mass - 48.021129 - Masses['H2O'] + Masses['NH3'] + Masses['H+'])
    else:
      return (self.lipid.tails[0].mass - 48.021129 - Masses['H2O'] + Masses['NH3'] - Masses['H+'])
  def Formula(self):
    formula = Counter(self.lipid.tails[0].formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.update({'N':1, 'C':-1, 'H':-2, 'O':-3})
    else:
      formula.update({'N':1, 'C':-1, 'H':-4, 'O':-3})
    return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1  
  def Comment(self):
    return 'Ceramide fragment D'

class Cer_P(Fragment):
  '''[ Base - C2H6O2 - H+ ]-'''
  # Fragment 'P' https://pubs.acs.org/doi/abs/10.1021/ac00049a004
  def MZ(self):
    if Masses[self.adduct][1] == 'Positive':
      return (self.lipid.tails[0].mass - 62.036779432 + Masses['H+'])
    else:
      return (self.lipid.tails[0].mass - 62.036779432 - Masses['H+'])
  def Formula(self):
    formula = Counter(self.lipid.tails[0].formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.subtract({'C':2, 'H':5, 'O':2})
    else:
      formula.subtract({'C':2, 'H':7, 'O':2})
    return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1
  def Comment(self):
    return 'Ceramide fragment P'

class Cer_Q(Fragment):
  '''[ Base - C3H8O3 - H+ ]-'''
  # Fragment 'Q' https://doi.org/10.1002/rcm.878
  def MZ(self):
    return (self.lipid.tails[0].mass - 92.047344116 - Masses['H+'])
  def Formula(self):
    formula = Counter(self.lipid.tails[0].formula)
    formula.subtract({'C':3, 'H':8, 'O':3})
    return formula
  def Charge(self):
    return -1
  def Comment(self):
    return 'Ceramide fragment Q'

class Cer_R(Fragment):
  '''Long-Chain-Base fragment\n
  [ Base - 2H2O + H+ ]+\n or
  [ Base - NH3 - H2O - H+ ]-'''
  # Fragment 'R' https://pubs.acs.org/doi/abs/10.1021/ac00049a004
  def MZ(self):
    if Masses[self.adduct][1] == 'Positive': # Base is created without the ammonia, so needs to be added here
      return (self.lipid.tails[0].mass - 3*Masses['H2O'] + Masses['NH3'] + Masses['H+'])
    else:
      return (self.lipid.tails[0].mass - 2*Masses['H2O'] - Masses['H+'])
  def Formula(self):
    formula = Counter(self.lipid.tails[0].formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.update({'N':1, 'H':-2, 'O':-3})
    else:
      formula.subtract({'H':5, 'O':2})
    return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1  
  def Comment(self):
      return 'Ceramide fragment R'

class Cer_Rb(Fragment):
  '''[ Base - CH6O2 - H+ ]-\n
  Fragment 'R' variant for Phytosphingosine'''
  # Fragment 'R' but for Phytosphingosine
  def MZ(self):
    return (self.lipid.tails[0].mass - 50.036779432 - Masses['H+'])
  def Formula(self):
    formula = Counter(self.lipid.tails[0].formula)
    formula.subtract({'C':1, 'H':7, 'O':2})
    return formula
  def Charge(self):
    return -1
  def Comment(self):
    return 'Ceramide fragment Rb'

class Cer_S(Fragment):
  '''[ FA + C2H5NO - H2O - H+ ]-'''
  # Fragment 'S' https://pubs.acs.org/doi/abs/10.1021/ac00049a004
  #              https://doi.org/10.1002/rcm.878
  def MZ(self):
    return (self.lipid.tails[1].mass + (59.037113785-Masses['H2O']-Masses['H+']))
  def Formula(self):
    formula = Counter(self.lipid.tails[1].formula)
    formula.update({'C':2, 'H':2, 'N':1})
    return formula
  def Charge(self):
    return -1
  def Comment(self):
    return 'Ceramide fragment S'

class Cer_T(Fragment):
  '''[ FA + C2H5N - H2O - H+ ]-'''
  # Fragment 'T' https://pubs.acs.org/doi/abs/10.1021/ac00049a004
  #              https://doi.org/10.1002/rcm.878
  def MZ(self):
    return (self.lipid.tails[1].mass + (43.042199165-Masses['H2O']-Masses['H+']))
  def Formula(self):
    formula = Counter(self.lipid.tails[1].formula)
    formula.update({'C':2, 'H':2, 'N':1, 'O':-1})
    return formula
  def Charge(self):
    return -1
  def Smiles(self):
    return 'C=C[N-]'+self.lipid.tails[1].smiles
  def Comment(self):
    return 'Ceramide fragment T'

class FA_C2H3N(Fragment):
  '''[ FA + C2H3N - HO- ]+'''
  def MZ(self):
    return (self.lipid.tails[1].mass + 41.026549101 - Masses['OH-'])
  def Formula(self):
    formula = Counter(self.lipid.tails[1].formula)
    formula.update({'C':2, 'H':2, 'N':1, 'O':-1})
    return formula
  def Charge(self):
    return 1
  def Comment(self):
    chnops = self.lipid.tails[1].formula.copy()
    chnops.update({'C':2, 'H':2, 'N':1, 'O':-1})
    new_chnops = {k: v for k, v in chnops.items() if v != 0}
    tail = ''.join(''.join((key, str(val))) for (key, val) in new_chnops.items())
    comment = '['+tail+']+'
    return comment

class Cer_U(Fragment):
  '''[ FA + NH3 - H2O (+/-) H+ ]-\n
  Free fatty acid RCONH2, i.e. R-C(-OH)=NH'''
  # Fragment 'U' https://pubs.acs.org/doi/abs/10.1021/ac00049a004
  def MZ(self):
    if Masses[self.adduct][1] == 'Positive':
      return (self.lipid.tails[1].mass + (Masses['NH3']-Masses['H2O']+Masses['H+']))
    else:
      return (self.lipid.tails[1].mass + (Masses['NH3']-Masses['H2O']-Masses['H+']))
  def Formula(self):
    formula = Counter(self.lipid.tails[1].formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.update({'N':1, 'H':2, 'O':-1})
    else:  
      formula.update({'N':1, 'O':-1})
    return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1  
  def Comment(self):
    if Masses[self.adduct][1] == 'Positive':
      return '[RONH3]+'
    else:
      return 'Ceramide fragment U'

class Cer_W(Fragment):
  '''[ FA + C2H5N - H+ ]-\n
  Appears in tCers'''
  # Fragment 'W' https://pubs.acs.org/doi/abs/10.1021/ac00049a004
  def MZ(self):
    return (self.lipid.tails[1].mass + 43.042199165 - Masses['H+'])
  def Formula(self):
    formula = Counter(self.lipid.tails[1].formula)
    formula.update({'C':2, 'H':4, 'N':1})
    return formula
  def Charge(self):
    return -1
  def Comment(self):
    return 'Ceramide fragment W'

class Cer_X(Fragment):
  '''[ FA + C3H5N - H+ ]-\n
  Appears in tCers'''
  # Fragment 'X' https://pubs.acs.org/doi/abs/10.1021/ac00049a004
  def MZ(self):
    return (self.lipid.tails[1].mass + 55.042199165 - Masses['H+'])
  def Formula(self):
    formula = Counter(self.lipid.tails[1].formula)
    formula.update({'C':3, 'H':4, 'N':1})
    return formula
  def Charge(self):
    return -1
  def Comment(self):
    return 'Ceramide fragment X'

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # Headgroup neutral losses, variants

# ~ # Sometimes the headgroup takes the phosphate with it

class MA_s_HG(MA):
  def __init__(self, lipid, adduct, intensity, fragmentType=None):
      for sn in lipid.tails: # ie, sn1, sn2, or sn3
        if sn.type == 'Headgroup':
          self.headgroup = sn
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
      return super().MZ() - (self.headgroup.mass - Masses['H2O'])/abs(Masses[self.adduct][2])
  def Formula(self):
      formula = super().Formula()
      formula.subtract(self.headgroup.formula)
      formula.update({'H':2, 'O':1})
      return formula
  def Comment(self):
    comment = super().Comment()
    headgroup = ''.join(''.join((key, str(val))) for (key, val) in self.headgroup.formula.items())
    comment = comment.replace('M', 'M-'+headgroup)
    return comment

class MA_s_HG_H2O(MA_s_HG):
  def MZ(self):
      return super().MZ() - Masses['H2O']/abs(Masses[self.adduct][2])
  def Formula(self):
      formula = super().Formula()
      formula.subtract({'H':2, 'O':1})
      return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-H2O]')
    return comment  

class MA_s_HG_2H2O(MA_s_HG):
  def MZ(self):
      return super().MZ() - 2*Masses['H2O']/abs(Masses[self.adduct][2])
  def Formula(self):
      formula = super().Formula()
      formula.subtract({'H':4, 'O':2})
      return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-2H2O]')
    return comment  

# ~ # 

def MA_s_HG_FA(lipid, adduct, intensity):
  '''[ MA - Headgroup + H2O - (ROOH) ]\n
  Fragment for a 'clean' headgroup neutral loss and loss of free-fatty acid\n
  i.e. loss of headgroup, including phospate from adducted molecular ion\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MA_s_HG_FAx(lipid, adduct, intensity, MA_s_HG_FA, tail)

class MA_s_HG_FAx(MA_s_HG):
  '''Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - (self.tail.mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_HG_sn1(MA_s_HG):
  '''[ MA - Headgroup + H2O - (ROOH) ] (sn1)\n
  Fragment for a 'clean' headgroup neutral loss and loss of free-fatty acid\n
  i.e. loss of headgroup, including phospate from adducted molecular ion'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[0].mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_HG_sn2(MA_s_HG):
  '''[ MA - Headgroup + H2O - (ROOH) ] (sn2)\n
  Fragment for a 'clean' headgroup neutral loss and loss of free-fatty acid\n
  i.e. loss of headgroup, including phospate from adducted molecular ion'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[1].mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

# ~ #

def MA_s_HG_FA_H2O(lipid, adduct, intensity):
  '''[ MA - Headgroup - H2O - (ROOH) ]\n
  Fragment for a 'clean' headgroup neutral loss and loss of free-fatty acid and water\n
  i.e. loss of headgroup, including phospate and bridging -OH from adducted molecular ion\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MA_s_HG_FA_H2Ox(lipid, adduct, intensity, MA_s_HG_FA_H2O, tail)

class MA_s_HG_FA_H2Ox(MA_s_HG_H2O):
  '''Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - (self.tail.mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_HG_sn1_H2O(MA_s_HG_H2O):
  '''[ MA - Headgroup - (ROOH) ] (sn1)\n
  Fragment for a 'clean' headgroup neutral loss and loss of free-fatty acid and water\n
  i.e. loss of headgroup, including phospate and bridging -OH from adducted molecular ion'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[0].mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_HG_sn2_H2O(MA_s_HG_H2O):
  '''[ MA - Headgroup - (ROOH) ] (sn2)\n
  Fragment for a 'clean' headgroup neutral loss and loss of free-fatty acid and water\n
  i.e. loss of headgroup, including phospate and bridging -OH from adducted molecular ion'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[1].mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

# ~ # 

def MA_s_HG_FAk(lipid, adduct, intensity):
  '''[ MA - Headgroup - (R=O) ]\n
  Fragment for a 'clean' headgroup neutral loss and loss of a fatty acid ketene\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MA_s_HG_FAkx(lipid, adduct, intensity, MA_s_HG_FAk, tail)

class MA_s_HG_FAkx(MA_s_HG):
  '''Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - (self.tail.mass-Masses['H2O'])/abs(Masses[self.adduct][2])
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_HG_sn1k(MA_s_HG):
  '''[ MA - Headgroup - (R=O) ] (sn1)\n
  Fragment for a 'clean' headgroup neutral loss and loss of a fatty acid ketene'''
  def MZ(self):
    return super().MZ() - ((self.lipid.tails[0].mass-Masses['H2O'])/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment  

class MA_s_HG_sn2k(MA_s_HG):
  '''[ MA - Headgroup - (R=O) ] (sn2)\n
  Fragment for a 'clean' headgroup neutral loss and loss of a fatty acid ketene'''
  def MZ(self):
    return super().MZ() - ((self.lipid.tails[1].mass-Masses['H2O'])/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment  

# ~ #  Sometimes the headgroup leaves the phosphate behind

class MH_PO4_s_HG(MH):
  '''[ M(+/-)H - Headgroup + PO4 ]\n
  Fragment for headgroup neutral loss, excluding phosphate\n
  A = Headgroup NL, with phosphate\n
  B = Headgroup NL, without phosphate\n
  C = Headgroup NL, with phosphate MA -> MH'''
  def __init__(self, lipid, adduct, intensity, fragmentType=None):
      for sn in lipid.tails: # ie, sn1, sn2, or sn3
        if sn.type == 'Headgroup':
          self.headgroup = sn
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
      return super().MZ() - (self.headgroup.mass - Masses['PO4H3'])
  def Formula(self):
      formula = super().Formula()
      formula.subtract(self.headgroup.formula)
      formula.update({'P':1, 'O':4, 'H':3})
      return formula
  def Comment(self):
    comment = self.adduct
    chnops = self.headgroup.formula.copy()
    chnops.subtract({'P':1, 'O':4, 'H':3})
    new_chnops = {k: v for k, v in chnops.items() if v != 0}
    headgroup = ''.join(''.join((key, str(val))) for (key, val) in new_chnops.items())
    comment = comment.replace('M', 'M-'+headgroup)
    return comment

class MH_P2O6_s_HG(MH):  # Headgroup neutral loss
  '''[ M(+/-)H - Headgroup + P2O6 ]\n
  Fragment for headgroup neutral loss, excluding two phosphites\n
  This is a special case for some PIPs\n
  A = Headgroup NL, with phosphate\n
  B = Headgroup NL, without phosphate\n
  C = Headgroup NL, with phosphate MA -> MH'''
  def __init__(self, lipid, adduct, intensity, fragmentType=None):
      for sn in lipid.tails: # ie, sn1, sn2, or sn3
        if sn.type == 'Headgroup':
          self.headgroup = sn
      super().__init__(lipid, adduct, intensity, fragmentType)
  
  def MZ(self):
      return super().MZ() - (self.headgroup.mass - 2*Masses['PO3H'])
  def Formula(self):
      formula = super().Formula()
      formula.subtract(self.headgroup.formula)
      formula.update({'P':2, 'O':6, 'H':2})
      return formula
  def Comment(self):
    comment = self.adduct
    chnops = self.headgroup.formula.copy()
    chnops.subtract({'P':2, 'O':6, 'H':2})
    new_chnops = {k: v for k, v in chnops.items() if v != 0}
    headgroup = ''.join(''.join((key, str(val))) for (key, val) in new_chnops.items())
    comment = comment.replace('M', 'M-'+headgroup)
    return comment

class MA_PO4_s_HG(MA):
  '''[ MA - Headgroup + PO4 ]\n
  Fragment for headgroup neutral loss, excluding phosphate\n
  A = Headgroup NL, with phosphate\n
  B = Headgroup NL, without phosphate\n
  C = Headgroup NL, with phosphate MA -> MH'''
  def __init__(self, lipid, adduct, intensity, fragmentType=None):
      for sn in lipid.tails: # ie, sn1, sn2, or sn3
        if sn.type == 'Headgroup':
          self.headgroup = sn
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
      return super().MZ() - (self.headgroup.mass - Masses['PO4H3'])
  def Formula(self):
      formula = super().Formula()
      formula.subtract(self.headgroup.formula)
      formula.update({'P':1, 'O':4, 'H':3})
      return formula
  def Comment(self):
    comment = self.adduct
    chnops = self.headgroup.formula.copy()
    chnops.subtract({'P':1, 'O':4, 'H':3})
    new_chnops = {k: v for k, v in chnops.items() if v != 0}
    headgroup = ''.join(''.join((key, str(val))) for (key, val) in new_chnops.items())
    comment = comment.replace('M', 'M-'+headgroup)
    return comment

class MA_C3H8O8P2_s_HG(MA):
  '''[ MA - Headgroup + PO4 ]\n
  Fragment for cardiolipid'''
  def __init__(self, lipid, adduct, intensity, fragmentType=None):
      for sn in lipid.tails: # ie, sn1, sn2, or sn3
        if sn.type == 'Headgroup':
          self.headgroup = sn
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
      return super().MZ() - (self.headgroup.mass - Masses['PO4H3'] - 135.992544)
  def Formula(self):
      formula = super().Formula()
      formula.subtract(self.headgroup.formula)
      formula.update({'C':3, 'H':8, 'O':8, 'P':2})
      return formula
  def Comment(self):
    comment = self.adduct
    chnops = self.headgroup.formula.copy()
    chnops.subtract({'C':3, 'H':8, 'O':8, 'P':2})
    new_chnops = {k: v for k, v in chnops.items() if v != 0}
    headgroup = ''.join(''.join((key, str(val))) for (key, val) in new_chnops.items())
    comment = comment.replace('M', 'M-'+headgroup)
    return comment

class MH_PO4_s_HG_H2O(MH_PO4_s_HG):
  '''[ M(+/-)H - Headgroup - H2O + PO4 ]\n
  Fragment for headgroup neutral loss, excluding phosphate, but loss of water\n
  A = Headgroup NL, with phosphate\n
  B = Headgroup NL, without phosphate\n
  C = Headgroup NL, with phosphate MA -> MH'''
  def MZ(self):
      return super().MZ() - Masses['H2O']
  def Formula(self):
      formula = super().Formula()
      formula.subtract({'H':2, 'O':1})
      return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-H2O]')
    return comment  

# ~ # 

def MH_PO4_s_HG_FA(lipid, adduct, intensity):
  '''[ M(+/-)H - Headgroup + PO4 - (ROOH) ]\n
  Fragment for headgroup neutral loss, excluding phosphate and loss of free-fatty acid\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield HG_FA_NL_Bx(lipid, adduct, intensity, MH_PO4_s_HG_FA, tail)

class HG_FA_NL_Bx(MH_PO4_s_HG):
  '''Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - self.tail.mass
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

class HG_sn1_NL_B(MH_PO4_s_HG):
  '''[ M(+/-)H - Headgroup + PO4 - (ROOH) ] (sn1)\n
  Fragment for headgroup neutral loss, excluding phosphate and loss of free-fatty acid\n
  i.e. loss of headgroup, including phospate and bridging -OH from adducted molecular ion'''
  def MZ(self):
    return super().MZ() - self.lipid.tails[0].mass
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

class HG_sn2_NL_B(MH_PO4_s_HG):
  '''[ M(+/-)H - Headgroup + PO4 - (ROOH) ] (sn2)\n
  Fragment for headgroup neutral loss, excluding phosphate and loss of free-fatty acid\n
  i.e. loss of headgroup, including phospate and bridging -OH from adducted molecular ion'''
  def MZ(self):
    return super().MZ() - self.lipid.tails[1].mass
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

# ~ # 

def MH_PO4_s_HG_FAk(lipid, adduct, intensity):
  '''[ M(+/-)H - Headgroup + PO4 - (ROOH) ]\n
  Fragment for headgroup neutral loss, excluding phosphate and loss of free-fatty acid\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MH_PO4_s_HG_FAkx(lipid, adduct, intensity, MH_PO4_s_HG_FAk, tail)

class MH_PO4_s_HG_FAkx(MH_PO4_s_HG):
  '''Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - (self.tail.mass-Masses['H2O'])
  def Formula(self):
    formula = super().Formula()
    formula.update({'H':2 ,'O':1})
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

# ~ #

def HG_FAk_NL_B(lipid, adduct, intensity):
  '''[ M(+/-)H - Headgroup + PO4 - (R=O) ]\n
  Fragment for headgroup neutral loss, excluding phosphate and loss of a fatty acid ketene\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield HG_FAk_NL_Bx(lipid, adduct, intensity, HG_FAk_NL_B, tail)

class HG_FAk_NL_Bx(MH_PO4_s_HG):
  '''Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - (self.tail.mass-Masses['H2O'])
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

class HG_sn1k_NL_B(MH_PO4_s_HG):
  '''[ M(+/-)H - Headgroup + PO4 - (R=O) ] (sn1)\n
  Fragment for headgroup neutral loss, excluding phosphate and loss of a fatty acid ketene'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[0].mass-Masses['H2O'])
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

class HG_sn2k_NL_B(MH_PO4_s_HG):
  '''[ M(+/-)H - Headgroup + PO4 - (R=O) ] (sn2)\n
  Fragment for headgroup neutral loss, excluding phosphate and loss of a fatty acid ketene'''
  def MZ(self):
    return super().MZ() - (self.lipid.tails[1].mass-Masses['H2O'])
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment

# ~ # 

def HG_FA_NL_H2O_B(lipid, adduct, intensity):
  '''[ M(+/-)H - Headgroup + PO4 - H2O - (ROOH) ]\n
  Fragment for headgroup neutral loss, excluding phosphate and loss of free-fatty acid and water\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield HG_FA_NL_H2O_Bx(lipid, adduct, intensity, HG_FA_NL_H2O_B, tail)

class HG_FA_NL_H2O_Bx(MH_PO4_s_HG_H2O):
  '''Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - self.tail.mass
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

class HG_sn1_NL_H2O_B(MH_PO4_s_HG_H2O):
  '''[ M(+/-)H - Headgroup + PO4 - (ROOH) ] (sn1)\n
  Fragment for a 'clean' headgroup neutral loss and loss of free-fatty acid and water\n
  i.e. loss of headgroup, including phospate and bridging -OH from adducted molecular ion'''
  def MZ(self):
    return super().MZ() - self.lipid.tails[0].mass
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[0].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

class HG_sn2_NL_H2O_B(MH_PO4_s_HG_H2O):
  '''[ M(+/-)H - Headgroup + PO4 - (ROOH) ] (sn2)\n
  Fragment for a 'clean' headgroup neutral loss and loss of free-fatty acid and water\n
  i.e. loss of headgroup, including phospate and bridging -OH from adducted molecular ion'''
  def MZ(self):
    return super().MZ() - self.lipid.tails[1].mass
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.lipid.tails[1].formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment

# ~ #  Sometimes the headgroup leaves with the adduct, leaving M+/-H

class MH_s_HG(MH):
  def __init__(self, lipid, adduct, intensity, fragmentType=None):
      for sn in lipid.tails: # ie, sn1, sn2, or sn3
        if sn.type == 'Headgroup':
          self.headgroup = sn
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
      return super().MZ() - (self.headgroup.mass - Masses['H2O'])
  def Formula(self):
      formula = super().Formula()
      formula.subtract(self.headgroup.formula)
      formula.update({'H':2, 'O':1})
      return formula
  def Comment(self):
    comment = super().Comment()
    headgroup = ''.join(''.join((key, str(val))) for (key, val) in self.headgroup.formula.items())
    comment = comment.replace('M', 'M-'+headgroup)
    return comment

class MH_s_HG_H2O(MH_s_HG):
  def MZ(self):
      return super().MZ() - Masses['H2O']
  def Formula(self):
      formula = super().Formula()
      formula.subtract({'H':2, 'O':1})
      return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-H2O]')
    return comment  

class MH_s_HG_2H2O(MH_s_HG):
  def MZ(self):
      return super().MZ() - 2*Masses['H2O']
  def Formula(self):
      formula = super().Formula()
      formula.subtract({'H':4, 'O':2})
      return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-2H2O]')
    return comment  

# ~ # 

def MH_s_HG_FA(lipid, adduct, intensity):
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MH_s_HG_FAx(lipid, adduct, intensity, MH_s_HG_FA, tail)

class MH_s_HG_FAx(MH_s_HG):
  '''Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - self.tail.mass
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

# ~ #

def MH_s_HG_FA_H2O(lipid, adduct, intensity):
  '''[ MA - Headgroup - H2O - (ROOH) ]\n
  Fragment for a 'clean' headgroup neutral loss and loss of free-fatty acid and water\n
  i.e. loss of headgroup, including phospate and bridging -OH from adducted molecular ion\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MH_s_HG_FA_H2Ox(lipid, adduct, intensity, MH_s_HG_FA_H2O, tail)

class MH_s_HG_FA_H2Ox(MH_s_HG_H2O):
  '''Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - (self.tail.mass/abs(Masses[self.adduct][2]))
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(ROOH)')
    comment += ' ('+self.tail.name+')'
    return comment  

# ~ # 

def MH_s_HG_FAk(lipid, adduct, intensity):
  '''[ MH - Headgroup - (R=O) ]\n
  Fragment for a 'clean' headgroup neutral loss and loss of a fatty acid ketene\n
  For nonspecific sn position\n
  Method used to generate multiple objects'''
  for tail in lipid.tails:
    if tail.type in ['Acyl']:
      yield MH_s_HG_FAkx(lipid, adduct, intensity, MH_s_HG_FAk, tail)

class MH_s_HG_FAkx(MH_s_HG):
  '''Do not use this class, intended for use in loop'''
  def __init__(self, lipid, adduct, intensity, fragmentType, tail):
      self.tail = tail
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return super().MZ() - (self.tail.mass-Masses['H2O'])/abs(Masses[self.adduct][2])
  def Formula(self):
    formula = super().Formula()
    formula.subtract(self.tail.formula)
    formula.update({'H':2 ,'O':1})
    return formula
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace('M', 'M-(R=O)')
    comment += ' ('+self.tail.name+')'
    return comment  

# ~ #  Sometimes the headgroup flies off the the adduct!

class HGA(Fragment):  # Headgroup + Adduct
  def __init__(self, lipid, adduct, intensity, fragmentType=None):
      for sn in lipid.tails: # ie, sn1, sn2, or sn3
        if sn.type == 'Headgroup':
          self.headgroup = sn
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    return (self.headgroup.mass + Masses[self.adduct][0])/abs(Masses[self.adduct][2])
  def Formula(self):
    formula = Counter(self.headgroup.formula)
    formula.update(Masses[self.adduct][3])
    return formula
  def Charge(self):
    return Masses[self.adduct][2]
  def Comment(self):
    comment = self.adduct
    headgroup = ''.join(''.join((key, str(val))) for (key, val) in self.headgroup.formula.items())
    comment = comment.replace('M', headgroup)
    return comment

class HGA_s_H2O(HGA):  # Headgroup + Adduct - H2O
  def MZ(self):
    return super().MZ() - (Masses['H2O'])/abs(Masses[self.adduct][2])
  def Formula(self):
    formula = Counter(self.headgroup.formula)
    formula.subtract({'H':2, 'O':1})
    return formula
  def Charge(self):
    return Masses[self.adduct][2]
  def Comment(self):
    comment = super().Comment()
    comment = comment.replace(']', '-H2O]')
    return comment

class HGH(Fragment):  # Headgroup + Adduct
  def __init__(self, lipid, adduct, intensity, fragmentType=None):
      for sn in lipid.tails: # ie, sn1, sn2, or sn3
        if sn.type == 'Headgroup':
          self.headgroup = sn
      super().__init__(lipid, adduct, intensity, fragmentType)

  def MZ(self):
    if Masses[self.adduct][1] == 'Positive':
      return (self.headgroup.mass + Masses['H+'])
    else:
      return (self.headgroup.mass - Masses['H+'])
  def Formula(self):
    formula = Counter(self.headgroup.formula)
    if Masses[self.adduct][1] == 'Positive':
      formula.update({'H':1})
    else:
      formula.subtract({'H':1})
    return formula
  def Charge(self):
    if Masses[self.adduct][1] == 'Positive':
      return 1
    else:
      return -1
  def Comment(self):
    if Masses[self.adduct][1] == 'Positive':
      comment = '[M+H]+'
    else:
      comment = '[M-H]-'
    headgroup = ''.join(''.join((key, str(val))) for (key, val) in self.headgroup.formula.items())
    comment = comment.replace('M', headgroup)
    return comment

class HGH_s_H2O(HGH):  # Headgroup + Adduct - H2O
  def MZ(self):
    return super().MZ() - (Masses['H2O'])
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'H':2, 'O':1})
    return formula
  def Charge(self):
    return Masses[self.adduct][2]
  def Comment(self):
    comment = super().Comment()
    headgroup = ''.join(''.join((key, str(val))) for (key, val) in self.headgroup.formula.items())
    comment = comment.replace('M', headgroup+'-H2O')
    return comment

class HGH_s_PO4(HGH):  # Headgroup + Adduct - H2O
  def MZ(self):
    return super().MZ() - (Masses['PO4H3'])
  def Formula(self):
    formula = super().Formula()
    formula.subtract({'H':3, 'P':1, 'O':4})
    return formula
  def Charge(self):
    return super().Charge()
  def Comment(self):
    if Masses[self.adduct][1] == 'Positive':
      comment = '[M+H]+'
    else:
      comment = '[M-H]-'
    chnops = Counter(self.headgroup.formula)
    chnops.subtract({'H':3, 'P':1, 'O':4})
    new_chnops = {k: v for k, v in chnops.items() if v != 0}
    headgroup = ''.join(''.join((key, str(val))) for (key, val) in new_chnops.items())
    comment = comment.replace('M', headgroup)
    return comment

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# ~ # ~ # ~ # Characteristic Fragments

class C10H16N5O10P2(Fragment):
  '''X+H Fragment common to AcylCoA under Positive ESI\n
  MZ: 428.03669'''
  def MZ(self):
      return 428.03669
  def Formula(self):
    formula = Counter({'C':10, 'H':16, 'N':5, 'O':10, 'P':2,})
    return formula
  def Charge(self):
      return 1
  def Comment(self):
    return '[C10H16N5O10P2]+' 

class C15H27O13(Fragment):
  '''X-H Fragment common to DGDG under negative ESI\n
  MZ: 415.145715'''
  def MZ(self):
      return 415.145715
  def Formula(self):
    formula = Counter({'C':15, 'H':27, 'O':13})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C15H27O13]-' 

class C15H25O12(Fragment):
  '''X-H Fragment common to DGDG under negative ESI\n
  MZ: 397.13515'''
  def MZ(self):
      return 397.13515
  def Formula(self):
    formula = Counter({'C':15, 'H':25, 'O':12})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C15H25O12]-' 

class C15H23O11(Fragment):
  '''X-H Fragment common to DGDG under negative ESI\n
  MZ: 379.124585'''
  def MZ(self):
      return 379.124585
  def Formula(self):
    formula = Counter({'C':15, 'H':23, 'O':11})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C15H23O11]-' 

class C6H10NaO11P2(Fragment):
  '''X+Na-2H Fragment for sodiated PIP\n
  MZ: 342.960152'''
  def MZ(self):
    return 342.960152
  def Formula(self):
    formula = Counter({'C':6, 'H':10, 'Na':1, 'O':11, 'P':2})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C6H10O11P2Na]-'

class C6H8NaO10P2(Fragment):
  '''X+Na-2H Fragment for sodiated PIP\n
  MZ: 324.949587'''
  def MZ(self):
    return 324.949587
  def Formula(self):
    formula = Counter({'C':6, 'H':8, 'Na':1, 'O':10, 'P':2})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C6H8O10P2Na]-'

class C6H11O11P2(Fragment):
  '''X-H Fragment common for PIP\n
  MZ: 320.978207'''
  def MZ(self):
    return 320.978207
  def Formula(self):
    formula = Counter({'C':6, 'H':11, 'O':11, 'P':2})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C6H11O11P2]-'

class C9H16O10P(Fragment):
  '''X-H Fragment common for PI\n
  MZ: 315.048656'''
  def MZ(self):
    return 315.048656
  def Formula(self):
    formula = Counter({'C':9, 'H':16, 'O':10, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C9H16O10P]-'

class C6H9O10P2(Fragment):
  '''X-H Fragment common for PIP\n
  MZ: 302.967642'''
  def MZ(self):
    return 302.967642
  def Formula(self):
    formula = Counter({'C':6, 'H':9, 'O':10, 'P':2})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C6H9O10P2]-'

class C9H14O9P(Fragment):
  '''X-H Fragment common for PI\n
  MZ: 297.038092'''
  def MZ(self):
    return 297.038092
  def Formula(self):
    formula = Counter({'C':9, 'H':14, 'O':9, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C9H14O9P]-'

class C6H12O9P(Fragment):
  '''X-H Fragment common for PI\n
  MZ: 259.022442'''
  def MZ(self):
    return 259.022442
  def Formula(self):
    formula = Counter({'C':6, 'H':12, 'O':9, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C6H12O9P]-'

class C6H10O8P(Fragment):
  '''X-H Fragment common for PI\n
  MZ: 241.011876845'''
  def MZ(self):
    return 241.011876845
  def Formula(self):
    formula = Counter({'C':6, 'H':10, 'O':8, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C6H10O8P]-'

class C6H9O8S(Fragment):
  '''X-H Fragment common for PI\n
  MZ: 241.002360812'''
  def MZ(self):
    return 241.002360812
  def Formula(self):
    formula = Counter({'C':6, 'H':9, 'O':8, 'S':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C6H9O8S]-'

class C9H15O7(Fragment):
  '''X-H Fragment common to DGDG under negative ESI\n
  MZ: 235.082326'''
  def MZ(self):
      return 235.082326
  def Formula(self):
    formula = Counter({'C':9, 'H':15, 'O':7})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C9H15O7]-' 

class C6H12O7P(Fragment):
  '''X-H Headgroup fragment for PG\n
  MZ: 227.032612'''
  def MZ(self):
      return 227.032612
  def Formula(self):
    formula = Counter({'C':6, 'H':12, 'O':7, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C6H12O7P]-'

class C6H9O7S(Fragment):
  '''X-H Fragment common to Sulphoquinovosyl diacylglycerol under negative ESI\n
  MZ: 225.007446'''
  def MZ(self):
      return 225.007446
  def Formula(self):
    formula = Counter({'C':6, 'H':9, 'O':7, 'S':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C6H9O7S]-'

class C6H8O7P(Fragment):
  '''X-H Fragment common for PI\n
  MZ: 223.001312'''
  def MZ(self):
    return 223.001312
  def Formula(self):
    formula = Counter({'C':6, 'H':8, 'O':7, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C6H8O7P]-'

class C5H11NO5P(Fragment):
  '''X-H Headgroup fragment for PE\n
  MZ: 196.038032'''
  def MZ(self):
      return 196.038032
  def Formula(self):
    formula = Counter({'C':5, 'H':11, 'N':1, 'O':5, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C5H11NO5P]-'

class C3H7NaO6P(Fragment):
  '''X+Na-2H Fragment common to sodiated Glycerophospholipids under negative ESI\n
  MZ: 192.988343'''
  def MZ(self):
    return 192.988343
  def Formula(self):
    formula = Counter({'C':3, 'H':7, 'Na':1, 'O':6, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C3H7O6PNa]-'

class C4H10O6P(Fragment):
  '''X-H Fragment common to MPA under negative ESI\n
  MZ: 185.022048'''
  def MZ(self):
      return 185.022048
  def Formula(self):
    formula = Counter({'C':4, 'H':10, 'O':6, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C4H10O6P]-' 

class C5H15NO4P(Fragment):
  '''X+H Headgroup fragment for PC\n
  MZ: 184.0733204'''
  def MZ(self):
    return 184.0733204
  def Formula(self):
    formula = Counter({'C':5, 'H':15, 'N':1, 'O':4, 'P':1})
    return formula
  def Charge(self):
      return 1
  def Comment(self):
    return '[C5H15NO4P]+'

class C3H8O6P(Fragment):
  '''X+H Headgroup fragment for PG or PA\n
  MZ: 171.006397541'''
  def MZ(self):
    return 171.006397541
  def Formula(self):
    formula = Counter({'C':3, 'H':8, 'O':6, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C3H8O6P]-'

class C4H8O5P(Fragment):
  '''X-H Fragment common to MPA under negative ESI\n
  MZ: 167.011483'''
  def MZ(self):
      return 167.011483
  def Formula(self):
    formula = Counter({'C':4, 'H':8, 'O':5, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C4H8O5P]-' 

class HO6P2(Fragment):
  '''X-H Headgroup fragment for PPA\n
  MZ: 158.925383317'''
  def MZ(self):
      return 158.925383317
  def Formula(self):
    formula = Counter({'H':1, 'O':6, 'P':2})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[HP2O6]-'

class C3H6O5P(Fragment):
  '''X-H Fragment common to Glycerophospholipids under negative ESI\n
  MZ: 152.995832857'''
  def MZ(self):
    return 152.995832857
  def Formula(self):
    formula = Counter({'C':3, 'H':6, 'O':5, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C3H6O5P]-'

class C2H5NaO4P(Fragment):
  '''X+Na Fragment common to PC under positive ESI\n
  MZ: 146.981766'''
  def MZ(self):
    return 146.981766
  def Formula(self):
    formula = Counter({'C':2, 'H':5, 'Na':1, 'O':4, 'P':1})
    return formula
  def Charge(self):
      return 1
  def Comment(self):
    return '[C2H5O4PNa]+' 

class C7H14N1O2(Fragment):
  '''X+H Fragment common to N-trimethylhomoserine diacylglycerol under positive ESI\n
  MZ: 144.101905128'''
  def MZ(self):
    return 144.101905128
  def Formula(self):
    formula = Counter({'C':7, 'H':14, 'N':1, 'O':2})
    return formula
  def Charge(self):
      return 1
  def Comment(self):
    return '[C7H14NO2]+' 

class C2H7NO4P(Fragment):
  '''X-H Headgroup fragment for PE\n
  MZ: 140.011817274'''
  def MZ(self):
      return 140.011817274
  def Formula(self):
    formula = Counter({'C':2, 'H':7, 'N':1, 'O':4, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C2H7NO4P]-' 

class C6H14NO2(Fragment):
  '''X+H Fragment common to DGCC under positive ESI\n
  MZ: 132.101905'''
  def MZ(self):
    return 132.101905
  def Formula(self):
    formula = Counter({'C':6, 'H':14, 'N':1, 'O':2})
    return formula
  def Charge(self):
      return 1
  def Comment(self):
    return '[C6H14NO2]+' 

class C6H5O3(Fragment):
  '''X-H Fragment common to DGDG under negative ESI\n
  MZ: 125.024418'''
  def MZ(self):
      return 125.024418
  def Formula(self):
    formula = Counter({'C':6, 'H':5, 'O':3})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[C6H5O3]-' 

class CH4O4P(Fragment):
  '''X-H Fragment common to MPA under negative ESI\n
  MZ: 110.985268'''
  def MZ(self):
      return 110.985268
  def Formula(self):
    formula = Counter({'C':1, 'H':4, 'O':4, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[CH4O4P]-' 

class C5H14NO(Fragment):
  '''X+H Fragment common to DGCC under positive ESI\n
  MZ: 104.10699'''
  def MZ(self):
    return 104.10699
  def Formula(self):
    formula = Counter({'C':5, 'H':14, 'N':1, 'O':1})
    return formula
  def Charge(self):
      return 1
  def Comment(self):
    return '[C5H14NO]+' 

class H2O4P(Fragment):
  '''X-H Fragment common to phospholipids under negative ESI\n
  MZ: 96.969618109'''
  def MZ(self):
      return 96.969618109
  def Formula(self):
    formula = Counter({'H':2, 'O':4, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[H2PO4]-' 

class O4SH(Fragment):
  '''X-H Fragment common to sulfatide under negative ESI\n
  MZ: 96.960102077'''
  def MZ(self):
      return 96.960102077
  def Formula(self):
    formula = Counter({'H':1, 'O':4, 'S':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[HSO4]-' 

class O3SH(Fragment):
  '''X-H Fragment common to Sulphoquinovosyl diacylglycerol under negative ESI\n
  MZ: 80.965187457'''
  def MZ(self):
      return 80.965187457
  def Formula(self):
    formula = Counter({'H':1, 'O':3, 'S':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[HSO3]-' 

class O3P(Fragment):
  '''X-H Fragment common to phospholipids under negative ESI\n
  MZ: 78.959053425'''
  def MZ(self):
      return 78.959053425
  def Formula(self):
    formula = Counter({'O':3, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[PO3]-' 

class H2O2P(Fragment):
  '''X-H Fragment common to CerP under negative ESI\n
  According to LipiDex\n
  MZ: 64.979789'''
  def MZ(self):
      return 64.979789
  def Formula(self):
    formula = Counter({'H':2, 'O':2, 'P':1})
    return formula
  def Charge(self):
      return -1
  def Comment(self):
    return '[H2O2P]-' 

class C3H10N(Fragment):
  '''X-H Fragment common to phospholipids under negative ESI\n
  MZ: 60.080776'''
  def MZ(self):
      return 60.080776
  def Formula(self):
    formula = Counter({'C':3, 'H':10, 'N':1})
    return formula
  def Charge(self):
      return 1
  def Comment(self):
    return '[C3H10N]+' 

# ~ # ~ # Cholesteryl ring fragments:

class C13H19(Fragment):
  '''X+H Fragment common to cholesteryl under positive ESI\n
  MZ: 175.148127043'''
  def MZ(self):
      return 175.148127043
  def Formula(self):
    formula = Counter({'C':13, 'H':19})
    return formula
  def Charge(self):
      return 1
  def Comment(self):
    return '[C13H19]+'  

class C12H17(Fragment):
  '''X+H Fragment common to cholesteryl under positive ESI\n
  MZ: 161.132476979'''
  def MZ(self):
      return 161.132476979
  def Formula(self):
    formula = Counter({'C':12, 'H':17})
    return formula
  def Charge(self):
      return 1
  def Comment(self):
    return '[C12H17]+'  

class C11H15(Fragment):
  '''X+H Fragment common to cholesteryl under positive ESI\n
  MZ: 147.116826915'''
  def MZ(self):
      return 147.116826915
  def Formula(self):
    formula = Counter({'C':11, 'H':15})
    return formula
  def Charge(self):
      return 1
  def Comment(self):
    return '[C11H15]+'

class C10H15(Fragment):
  '''X+H Fragment common to cholesteryl under positive ESI\n
  MZ: 135.116826915'''
  def MZ(self):
      return 135.116826915
  def Formula(self):
    formula = Counter({'C':10, 'H':15})
    return formula
  def Charge(self):
      return 1
  def Comment(self):
    return '[C10H15]+'

# ~ # Goodness gracious that's a lot of fragments...


