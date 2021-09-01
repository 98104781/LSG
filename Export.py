import GlycerolipidClasses as GC
import GenerateLipids as GL

from PySide6.QtWidgets import QFileDialog
import time
import csv
import os

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

def as_msp(save_file, lipid_data):

    count = 0
    for lipid in lipid_data:
        lipid.generate_spectra()
    
        for spectrum in lipid.spectra:

            save_file.write(f"NAME: {lipid.name} {spectrum}\n")
            save_file.write(f"IONMODE: {GL.Masses[spectrum][1]}\n")
            save_file.write(f"PRECURSORMZ: {round(GL.MA(lipid, GL.Masses[spectrum]), 6)}\n")
            save_file.write(f"PRECURSORTYPE: {spectrum}\n")
            save_file.write(f"Num Peaks: {len(lipid.spectra[spectrum])}\n")
            save_file.writelines(f"{peak[0]} {peak[1]}\n" for peak in lipid.spectra[spectrum])
            save_file.write("\n")
            count += 1

    return count

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

def as_csv(save_file, lipid_data): # For orbitrap inclusion list!

    count = 0
    writer = csv.writer(save_file)

    writer.writerow(['Mass [m/z]','Formula [M]','Formula type','Species','CS [z]','Polarity','Start [min]','End [min]','(N)CE','(N)CE type','MSX ID','Comment'])
    unique_mass = []
    for lipid in lipid_data:
        for adduct in lipid.adducts:
            MA = round(GL.Adduct_Spectra(lipid, adduct,{GL.MA:1}).spectrum[0][0], 6)
            if MA not in unique_mass: # This takes a lot of time!
                unique_mass.append(MA) 
                writer.writerow([MA,'','',type(lipid).__name__ ,GL.Masses[adduct][2],GL.Masses[adduct][1],'','','','','',adduct])
                count += 1

    return count

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

def save_as():

    # Generate lipids
    lipid_list = GC.Generate_Lipids()
    # Create save location
    file_name, _ = QFileDialog.getSaveFileName(filter="MSP (*.msp);;CSV (*.csv)")

    if file_name:

        t0 = time.time()

        # If save location exists, override
        if os.path.exists(file_name):
            os.remove(file_name) # Removes if exists
            print('Overwriting {}'.format(file_name))
        else: print('Creating {}'.format(file_name))

        save_file = open(file_name, 'x')

        # Based on file extension, export differently
        if '.msp' in file_name:
            count = as_msp(save_file, lipid_list)
        elif '.csv' in file_name:
            count = as_csv(save_file, lipid_list)
        else: print('Unsupported file type')

        t1 = time.time()        
        print(f"Generated {count} spectra in {t1-t0:.4f} seconds!")

        save_file.close()

    else: pass

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #


