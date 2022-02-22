import csv
import GenerateLipids as GL

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

def as_msp(self, save_file, lipid_data):
    '''
    Defines how to export data when saved as .MSP.
    Contains lipid fragmentation informaiton.
    '''
    count = 0
    for lipid in lipid_data:

        for adduct in lipid.adduct_dictionary:
            lipid.resolve_spectra(adduct, lipid.adduct_dictionary[adduct])
    
            written_masses = []
            for frag in lipid.spectra[adduct]:
                mz = round(frag.MZ(), 6)
                intensity = frag.intensity
                if [mz, intensity] not in written_masses and intensity > 0:
                    written_masses.append([mz, intensity])

            save_file.write(f"NAME: {lipid.name} {adduct}\n")
            save_file.write(f"IONMODE: {GL.Masses[adduct][1]}\n")
            save_file.write(f"PRECURSORMZ: {round(GL.MA(lipid, adduct, 0).MZ(), 6)}\n")
            save_file.write(f"COMPOUNDCLASS: {lipid.lipid_class}\n")
            save_file.write(f"FORMULA: {''.join(''.join((key, str(val))) for (key, val) in lipid.formula.items())} \n")
            save_file.write(f"RETENTIONTIME: 0.00\n") # Pointless
            save_file.write(f"PRECURSORTYPE: {adduct}\n")

            save_file.write(f"Num Peaks: {len(written_masses)}\n")
            save_file.writelines(f"{x[0]} {x[1]}\n" for x in written_masses)
            
            save_file.write("\n")
            count += 1

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
            ma = round(GL.MA(lipid, adduct, 0).MZ(), 6) # lipid.resolve_frag(adduct, {GL.MA:1})[0][0]
            if ma not in unique_mass: # This takes a lot of time!
                unique_mass.append(ma) # Removes all the duplicates
                writer.writerow([ma,'','',type(lipid).__name__ ,GL.Masses[adduct][2],GL.Masses[adduct][1],'','','','','',adduct])
                count +=1
            else: continue

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

        for adduct in lipid.adduct_dictionary:
            lipid.resolve_spectra(adduct, lipid.adduct_dictionary[adduct])

            prec_mz = round(GL.MA(lipid, adduct, 0).MZ(), 6)
            prec_formula = ''.join(''.join((key, str(val))) for (key, val) in lipid.formula.items())        
        
            written_masses = []
            for prod in lipid.spectra[adduct]:

                prod_mz = round(prod.MZ(), 6)
                prod_formula = ''.join(''.join((key, str(val))) for (key, val) in prod.Formula().items())

                if prod_mz not in written_masses and prod.intensity > 0 and prod_mz != prec_mz:
                    writer.writerow([lipid.lipid_class, lipid.name, prec_formula, adduct, prec_mz, GL.Masses[adduct][2], prod_formula, prod_mz, prod.Charge(), '', ''])
                    written_masses.append(prod_mz)
                    count +=1

    self.progress_bar.setValue(self.progress_bar.value()+1)
    return count

