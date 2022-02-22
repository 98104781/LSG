import csv
import GenerateLipids as GL
from memory_profiler import profile

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #
#@profile
def as_msp(self, save_file, lipid_data):
    '''
    Defines how to export data when saved as .MSP.
    Contains lipid fragmentation informaiton.
    '''
    count = 0
    for lipid in lipid_data:

        for adduct in lipid.adducts:
            lipid.resolve_spectra(adduct, lipid.adducts[adduct])
    
            save_file.write(f"NAME: {lipid.name} {adduct}\n")
            save_file.write(f"IONMODE: {GL.Masses[adduct][1]}\n")
            save_file.write(f"PRECURSORMZ: {GL.MA(lipid, adduct, 0).mass}\n")
            save_file.write(f"COMPOUNDCLASS: {lipid.lipid_class}\n")
            save_file.write(f"FORMULA: {''.join(''.join((key, str(val))) for (key, val) in lipid.formula.items())} \n")
            save_file.write(f"RETENTIONTIME: 0.00\n") # Pointless
            save_file.write(f"PRECURSORTYPE: {adduct}\n")

            spectrum = lipid.spectra[adduct]
            save_file.write(f"Num Peaks: {len(spectrum)}\n")
            save_file.writelines(f"{peak.mass} {peak.intensity}\n" for peak in spectrum)

            save_file.write("\n")
            count += 1

        del lipid

    self.progress_bar.setValue(self.progress_bar.value()+1)
    return count

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

def as_orb(self, save_file, lipid_data):
    '''
    Defines how to export data when saved as .CSV.
    Specifically for use in QE 'Orbitrap' inclusion list.
    '''
    count = 0
    writer = csv.writer(save_file)

    writer.writerow(['Mass [m/z]','Formula [M]','Formula type','Species','CS [z]','Polarity','Start [min]','End [min]','(N)CE','(N)CE type','MSX ID','Comment'])
    unique_mass = []

    for lipid in lipid_data:
        for adduct in lipid.adducts:
            prec = GL.MA(lipid, adduct, 0)
            if prec.mass not in unique_mass: # This takes a lot of time!
                unique_mass.append(prec.mass) # Removes all the duplicates
                writer.writerow([prec.mass,'','',type(lipid).__name__ ,GL.Masses[adduct][2],GL.Masses[adduct][1],'','','','','',adduct])
                count +=1
            else: continue
        del lipid

    self.progress_bar.setValue(self.progress_bar.value()+1)
    return count

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

def as_sky(self, save_file, lipid_data):
    '''
    Defines how to export data when saved as .CSV.
    Specifically for use in Skyline Transition list.
    '''
    count = 0
    writer = csv.writer(save_file)

    writer.writerow(['Molecule List Name', 'Precursor Name', 'Precursor Formula', 'Precursor Adduct', 'Precursor m/z', 'Precursor Charge', 'Product Formula', 'Product m/z', 'Product Charge', 'Explicit Retention Time', 'Explicit Collision Energy'])

    for lipid in lipid_data:

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
                    count +=1            
        del lipid

    self.progress_bar.setValue(self.progress_bar.value()+1)
    return count