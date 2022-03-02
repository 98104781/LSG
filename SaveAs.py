import csv
import GenerateLipids as GL
from PySide6.QtCore import QObject, Signal
from itertools import combinations_with_replacement as cwr, product

class Generator(QObject):
    
    finished = Signal()
    progress = Signal()
    fileError = Signal()

    progress_bar_increment = Signal()

    def __init__(self, file_name, filter, classes_to_generate, tails_to_generate, bases_to_generate, 
                 isomerism, specifics, lipidList):
        super().__init__()

        self.file_name = file_name
        self.filter = filter
        self.classes_to_generate = classes_to_generate
        self.tails_to_generate = tails_to_generate
        self.bases_to_generate = bases_to_generate
        self.isomerism = isomerism
        self.specifics = specifics
        self.lipidList = lipidList
        self.count = 0


    def run(self):
        self.save_file = open(self.file_name, 'x', newline='')
        if self.specifics:
            self.lipid_data = self.generate_specific()
        else: self.lipid_data = self.generate_range()
        if self.filter == "MSP (*.msp)": self.as_msp()
        elif self.filter == "Orbitrap Inclusion (*.csv)": self.as_orb()
        elif self.filter =="Skyline Transition (*.csv)": self.as_sky()
        else: self.fileError.emit()
        self.save_file.close()

    # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

    def generate_range(self):
        self.tails = GL.generate_acyl_tails(self.tails_to_generate)
        self.bases = GL.generate_base_tails(self.bases_to_generate)
        for cls in self.classes_to_generate:
            for adduct in cls.adducts: cls.adducts[adduct] = {k: v for k, v in cls.adducts[adduct].items() if v != 0}
            if not issubclass(cls, GL.Sphingolipid):
                if self.isomerism == True:
                    for comb in product(self.tails, repeat = cls.No_Tails):
                        yield cls(*comb)
                else: 
                    for comb in cwr(self.tails, cls.No_Tails):
                        yield cls(*comb)
            elif issubclass(cls, GL.Sphingolipid):
                for type in cls.base_types:
                    for comb in product(self.bases[type], self.tails):
                        yield cls(*comb)
                else: pass


    def generate_specific(self):
        lipidList = self.lipidList
        for lipid, selected_adduct in lipidList:
            adducts = lipid.adducts
            lipid.adducts = {k: v for k, v in adducts.items() if k is selected_adduct}
            lipid.adducts[selected_adduct] = {k: v for k, v in adducts[selected_adduct].items() if v != 0}
        return [x[0] for x in lipidList]

    # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #
    
    def as_msp(self):
        '''
        Defines how to export data when saved as .MSP.
        Contains lipid fragmentation informaiton.
        '''

        self.noun = 'spectra'
        for lipid in self.lipid_data:
            for adduct in lipid.adducts:
                
                lipid.resolve_spectra(adduct, lipid.adducts[adduct])

                self.save_file.write(f"NAME: {lipid.name} {adduct}\n"
                                     f"IONMODE: {GL.Masses[adduct][1]}\n"
                                     f"PRECURSORMZ: {GL.MA(lipid, adduct, 0).mass}\n"
                                     f"COMPOUNDCLASS: {lipid.lipid_class}\n"
                                     f"FORMULA: {''.join(''.join((key, str(val))) for (key, val) in lipid.formula.items())} \n"
                                     f"RETENTIONTIME: 0.00\n" # Pointless
                                     f"PRECURSORTYPE: {adduct}\n")
                spectrum = lipid.spectra[adduct]

                self.save_file.write(f"Num Peaks: {len(spectrum)}\n")
                self.save_file.writelines(f"{peak.mass} {peak.intensity}\n" for peak in spectrum)
                self.save_file.write("\n")
                self.count += 1


            del lipid

        self.finished.emit()

    # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #
    
    def as_orb(self):
        '''
        Defines how to export data when saved as .CSV.
        Specifically for use in QE 'Orbitrap' inclusion list.
        '''

        self.noun = 'precursors'

        writer = csv.writer(self.save_file)
        writer.writerow(['Mass [m/z]','Formula [M]','Formula type','Species','CS [z]','Polarity','Start [min]','End [min]','(N)CE','(N)CE type','MSX ID','Comment'])

        unique_mass = []
        for lipid in self.lipid_data:
            for adduct in lipid.adducts:
                prec = GL.MA(lipid, adduct, 0)
                if prec.mass not in unique_mass: # This takes a lot of time!
                    unique_mass.append(prec.mass) # Removes all the duplicates, orbitrap inclusion looks for the mass, doesn't care about anything else!
                    writer.writerow([prec.mass,'','',type(lipid).__name__ ,GL.Masses[adduct][2],GL.Masses[adduct][1],'','','','','',adduct])
                    self.count +=1
                else: continue
            del lipid

        self.finished.emit()

    # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #
    
    def as_sky(self):
        '''
        Defines how to export data when saved as .CSV.
        Specifically for use in Skyline Transition list.
        '''

        self.noun = 'transitions'

        writer = csv.writer(self.save_file)
        writer.writerow(['Molecule List Name', 'Precursor Name', 'Precursor Formula', 'Precursor Adduct', 'Precursor m/z', 'Precursor Charge', 'Product Formula', 'Product m/z', 'Product Charge', 'Explicit Retention Time', 'Explicit Collision Energy'])

        for lipid in self.lipid_data:

            for adduct in lipid.adducts:
                lipid.resolve_spectra(adduct, lipid.adducts[adduct])

                prec_mz = GL.MA(lipid, adduct, 0).mass
                prec_formula = ''.join(''.join((key, str(val))) for (key, val) in lipid.formula.items())        
            
                written_masses = []
                for prod in lipid.spectra[adduct]:

                    prod_formula = ''.join(''.join((key, str(val))) for (key, val) in prod.Formula().items())
                    if prod.mass not in written_masses and prod.intensity > 0 and prod.mass != prec_mz:
                        writer.writerow([lipid.lipid_class, lipid.name, prec_formula, adduct, prec_mz, GL.Masses[adduct][2], prod_formula, prod.mass, prod.Charge(), '', ''])
                        written_masses.append(prod.mass)
                        self.count +=1            
            del lipid

        self.finished.emit()