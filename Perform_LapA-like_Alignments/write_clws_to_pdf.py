from fpdf import FPDF
import os
import sys

def main(inDir, outFile):
	all_species = []
	out = ""
	files = os.listdir(inDir)
	for inFile in files:

		LapAorMapA = inFile[-8:-4]
		number = inFile[-10]
		species = inFile[:-10].replace('_', ' ')
		if species not in all_species:
			all_species.append(species)

		with open(str(inDir + inFile), 'r') as f:
			seq = ""
			id_count = 0
			sim_count = 0
			entry =[]
			for line in f.readlines():
				entry.append(line)
				if "_" in line:
					seq += line.split()[1]
				elif "MUSCLE" in line:
					continue
				elif "MapA" in line or 'LapA' in line:
					continue
				elif "*" in line:
					id_count += line.count("*")
					sim_count += line.count("*")
					sim_count += line.count(".")
					sim_count += line.count(":")
				elif "." in line:
					sim_count += line.count(".")
					sim_count += line.count(":")
				elif ":" in line:
					sim_count += line.count(":")
			
			ID = str(round(100 * id_count / len(seq), 1))
			Sim = str(round(100 * sim_count / len(seq), 1))
			header = str("%s adhesin %s aligned with Pf0-1 %s: %s %% identity; %s %% similarity" %(species, number, LapAorMapA, ID, Sim))
			entry ="".join(entry)
			entry = entry.replace("MUSCLE (3.8) multiple sequence alignment", header)
			out = out + '\n\n' + entry
	contents = "\n".join(all_species)
	out = "Species present in this file:\n\n" + contents + '\n\n\nAlignments:\n\n' + out

	pdf = FPDF()
	pdf.add_page()
	pdf.set_xy(0, 0)
	pdf.set_font('courier', 'B', 9.5)
	pdf.multi_cell(h=5.0, w=0, txt=out)
	pdf.output(outFile, 'F')

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])