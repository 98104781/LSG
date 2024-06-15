import sys
import csv

import Lipids.GenerateLipids as GL
from PySide6.QtCore import QObject, Signal
from itertools import combinations_with_replacement as cwr, product
from collections import Counter

class Generator(QObject):
    
    finished = Signal()
    progress = Signal(GL.Lipid)
    fileError = Signal()

    progress_bar_increment = Signal()

    def __init__(self, file_name, filter, classes_to_generate, tails_to_generate, bases_to_generate, 
                 isomerism, lipidSpecifics, lipidList, tailSpecifics, tailList, specificOrganisation):
        super().__init__()

        self.file_name = file_name
        self.filter = filter

        self.classes_to_generate = classes_to_generate
        self.tails_to_generate = tails_to_generate
        self.bases_to_generate = bases_to_generate

        self.isomerism = isomerism
        self.lipidSpecifics = lipidSpecifics
        self.lipidList = lipidList
        self.tailSpecifics = tailSpecifics
        self.tailList = tailList
        self.specificOrganisation = specificOrganisation

        self.count = 0


    def run(self):
        try:
            #import pydevd;pydevd.settrace(suspend=False)
            with open(self.file_name, 'x', newline='') as self.save_file:
                if self.lipidSpecifics: self.lipid_data = self.generate_specific()
                else: self.lipid_data = self.generate_range()
                if self.filter == "MSP (*.msp)": self.as_msp()
                elif self.filter == "Orbitrap Inclusion (*.csv)": self.as_orb()
                elif self.filter == "Skyline Transition (*.csv)": self.as_sky()
                else: self.fileError.emit()
        except: 
            (type, value, traceback) = sys.exc_info()
            sys.excepthook(type, value, traceback)
            self.fileError.emit()

    # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

    # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

    def addTailNames(self, lipid):
        c =  sum(snx.c  for snx in lipid.tails if snx.type != 'Headgroup')
        d =  sum(snx.d  for snx in lipid.tails if snx.type != 'Headgroup')
        me = sum(snx.me for snx in lipid.tails if snx.type != 'Headgroup')
        oh = sum(snx.oh for snx in lipid.tails if snx.type != 'Headgroup')
        dt = sum(snx.dt for snx in lipid.tails if snx.type != 'Headgroup')

        if 'Ether' in [snx.type for snx in lipid.tails if snx.type != 'Headgroup']:
            name = f"O-{c}:{d}"
        elif 'Vinyl' in [snx.type for snx in lipid.tails if snx.type != 'Headgroup']:
            name = f"P-{c}:{d}"
        else: name = f"{c}:{d}"
        
        if me > 0: # Methyl branching of fatty acid
            name += f";{me}-M" 
        if oh > 0: # Hydroxy functionalisation of fatty acid
            name += f";O{oh}"
        if dt > 0: # deuterium labelled fatty acids
            name += f"(D{dt})" # Deuterium doesn't update smiles currently.

        return name

    def checklipidAmbiguity(self, lipid, adduct):
        if not next((k.__name__ for k in lipid.adducts[adduct].keys() if 'FA' in k.__name__ or 'Cer' in k.__name__), False):
            return False
        else: return True

    ambiguousLipids = []
    def redifineAmbiguousLipid(self, lipid, adduct):
        lipid.ambiguousName = f"{lipid.lipid_class} {self.addTailNames(lipid)}" 
        if [lipid.ambiguousName, adduct] in self.ambiguousLipids:
            return False
        else: self.ambiguousLipids.append([lipid.ambiguousName, adduct])
        lipid.ambiguoussmiles = ' '
        return True

    def generate_range(self):

        self.acyls = []
        self.ethers = []
        self.vinyls = []
        self.bases = {}
        defaultBases = GL.generate_base_tails(self.bases_to_generate)

        if self.tailSpecifics:

            for tail in self.tailList:
                if tail.type == 'Acyl':
                    self.acyls.append(tail)
                elif tail.type == 'Ether':
                    self.ethers.append(tail)
                elif tail.type == 'Vinyl':
                    self.vinyls.append(tail)
                elif tail.type in GL.baseTypes:
                    if tail.type not in self.bases.keys(): self.bases[tail.type] = []
                    self.bases[tail.type].append(tail)

        else:

            self.acyls = GL.generate_tails(self.tails_to_generate, 'Acyl')
            self.ethers = GL.generate_tails(self.tails_to_generate, 'Ether')
            self.vinyls = GL.generate_tails(self.tails_to_generate, 'Vinyl')
            self.bases = defaultBases


        for cls in self.classes_to_generate:    
            cls.ambiguousSpectra = []             
            for adduct in cls.adducts: # Remove all ions in spectra with an intensity of 0
                cls.adducts[adduct] = {k: v for k, v in cls.adducts[adduct].items() if v != 0}
                cls.ambiguousSpectra.append(adduct if not self.checklipidAmbiguity(cls, adduct) else None)

            if self.specificOrganisation == True:
                try:constituents = cls.specificTailOrganisation
                except: constituents = cls.tailOrganisation
            else: constituents = cls.tailOrganisation # List of tails and organisation: ['B', 'AA', 'A'], indicates a sphingoid Base,
            constituentList = [] # a combination of two tails, and another tail independent of the previous combination is needed.

            for x in constituents: # ie, for 'B', 'AA', 'A' in ['B', 'AA', 'A']
                group = dict(Counter(x)) # ie, {'B':1}, {'A':2}, {'A':1}
                for key in group:
                    if key == 'B': # In the case of 'B', it indicates a base is needed. Get base list!
                        constituentList.append(self.flatten([self.bases[basetype] if basetype in self.bases.keys() else defaultBases[basetype] for basetype in cls.base_types ]))
                    elif key == 'A': # In the case of 'A', it indicates an acyl tail is needed. generate tail combination!
                        constituentList.append(cwr(self.acyls, r=group[key]))
                    elif key == 'O': # In the case of 'O', it indicates an ether tail is needed. generate tail combination!
                        constituentList.append(cwr(self.ethers, r=group[key]))
                    elif key == 'P': # In the case of 'P', it indicates a vinyl tail is needed. generate tail combination!
                        constituentList.append(cwr(self.vinyls, r=group[key]))
                    else: pass

            combs = list(product(*constituentList))
            for combination in combs:
                combination = self.flatten(combination)
                try:yield cls(*combination)
                except:pass
            self.progress.emit(cls)

    def flatten(self, data):
        if isinstance(data, tuple) or isinstance(data, list):
            for x in data: yield from self.flatten(x)
        else: yield data

    def generate_specific(self):
        lipidList = self.lipidList
        for lipid, selected_adduct in lipidList:
            adducts = lipid.adducts
            lipid.ambiguousSpectra = []
            lipid.adducts = {k: v for k, v in adducts.items() if k is selected_adduct}  # Remove all ions in spectra with an intensity of 0
            lipid.adducts[selected_adduct] = {k: v for k, v in adducts[selected_adduct].items() if v != 0}
            lipid.ambiguousSpectra.append(selected_adduct if not self.checklipidAmbiguity(lipid, selected_adduct) else None)
        return [x[0] for x in lipidList]

    # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

    # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

    def as_msp(self):
        '''
        Defines how to export data when saved as .MSP.
        Contains lipid fragmentation informaiton.
        '''
        
        self.noun = 'spectra' # Noun is used in Page 3 console when generation is completed
        string = ''
        for lipid in self.lipid_data:
            
            if self.specificOrganisation == True: 
                try: lipid.name = lipid.specificname
                except:pass

            for adduct in lipid.adducts:

                lipid.ambiguousName = False
                lipid.ambiguoussmiles = False
                if adduct in lipid.ambiguousSpectra:
                    generate = self.redifineAmbiguousLipid(lipid, adduct)
                else: generate = True

                if generate:
                    try: # Seems to be ever so slightly faster to batch print them, if uses a bit more memory.
                        lipid.resolve_spectra(adduct, lipid.adducts[adduct])
                        spectrum = lipid.spectra[adduct]
                        string += (f"NAME: {lipid.ambiguousName if lipid.ambiguousName else lipid.name} {adduct}\n"
                                f"IONMODE: {GL.adducts[adduct][1]}\n"
                                f"MW: {lipid.mass}\n"
                                f"PRECURSORMZ: {GL.MA(lipid, adduct, 0).mass}\n"
                                f"COMPOUNDCLASS: {lipid.lipid_class}\n"
                                f"FORMULA: {''.join(''.join((key, str(val))) for (key, val) in lipid.formula.items())}\n"
                                f"SMILES: {lipid.ambiguoussmiles if lipid.ambiguoussmiles else lipid.smiles}\n"
                                f"COMMENT: LSG in-silico\n" 
                                f"RETENTIONTIME: 0.00\n" # Pointless
                                f"PRECURSORTYPE: {adduct}\n"
                                f"Num Peaks: {len(spectrum)}\n")
                        for peak in spectrum:
                            string += f'{peak.mass} {peak.intensity} "{peak.Comment()}" \n'
                        string += '\n'
                        self.count += 1
                    except: pass
            del lipid

            if (self.count % 500 == 0): # Every 500, batch print to file
                self.save_file.write(string)
                string = ''
        self.save_file.write(string) # Batch print remaining to file
        string = ''
        self.finished.emit()

    # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

    # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

    def as_orb(self):
        '''
        Defines how to export data when saved as .CSV.
        Specifically for use in QE 'Orbitrap' inclusion list.
        '''

        self.noun = 'precursors' # Noun is used in Page 3 console when generation is completed
        writer = csv.writer(self.save_file)
        writer.writerow(['Mass [m/z]','Formula [M]','Formula type',
                         'Species','CS [z]','Polarity','Start [min]',
                         'End [min]','(N)CE','(N)CE type','MSX ID','Comment'])

        unique_mass = []
        for lipid in self.lipid_data:

            if self.specificOrganisation == True: 
                try: lipid.name = lipid.specificname
                except:pass

            for adduct in lipid.adducts:
                prec = GL.MA(lipid, adduct, 0)
                if prec.mass not in unique_mass: # This can take some lot of time
                    unique_mass.append(prec.mass) # Removes all the duplicate precursor masses

                    writer.writerow([prec.mass,'','',
                                     type(lipid).__name__ ,GL.adducts[adduct][2],GL.adducts[adduct][1],'',
                                     '','','','',adduct])
                    self.count +=1
                else: continue
            del lipid

        self.finished.emit()

    # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

    # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

    def as_sky(self):
        '''
        Defines how to export data when saved as .CSV.
        Specifically for use in Skyline Transition list.
        '''

        self.noun = 'transitions'  # Noun is used in Page 3 console when generation is completed
        writer = csv.writer(self.save_file)
        writer.writerow(['Molecule List Name', 'Precursor Name', 'Precursor Formula',
                         'Precursor Adduct', 'Precursor m/z', 'Precursor Charge', 'Product Formula',
                         'Product m/z', 'Product Charge', 'Explicit Retention Time', 'Explicit Collision Energy'])

        for lipid in self.lipid_data:

            if self.specificOrganisation == True: 
                try: lipid.name = lipid.specificname
                except:pass

            for adduct in lipid.adducts:

                lipid.ambiguousName = False
                lipid.ambiguoussmiles = False
                if adduct in lipid.ambiguousSpectra:
                    generate = self.redifineAmbiguousLipid(lipid, adduct)
                else: generate = True

                if generate:

                    lipid.resolve_spectra(adduct, lipid.adducts[adduct])

                    prec_mz = GL.MA(lipid, adduct, 0).mass
                    prec_formula = ''.join(''.join((key, str(val))) for (key, val) in lipid.formula.items())       

                    written_masses = []
                    for prod in lipid.spectra[adduct]:
                        prod_formula = ''.join(''.join((key, str(val))) for (key, val) in prod.Formula().items())
                        if prod.mass not in written_masses and prod.intensity > 0 and prod.mass != prec_mz:

                            writer.writerow([lipid.lipid_class, (lipid.ambiguousName if lipid.ambiguousName else lipid.name), prec_formula,
                                            adduct, prec_mz, GL.adducts[adduct][2], prod_formula,
                                            prod.mass, prod.Charge(), '', ''])

                            written_masses.append(prod.mass)
                            self.count +=1            
            del lipid

        self.finished.emit()