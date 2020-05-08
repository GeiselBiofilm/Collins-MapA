from fpdf import FPDF
import os
import sys

def main(inDir, outFile):
	"""
	Reads in clustal format alignments produced by run_muscle_on_dir.py. Creates .txt and .pdf files with 
	pairwise alignments of all LapA-like proteins with LapA and MapA. Also writes a contents table with 
	all species present in file at top of file. Each alignment has a title of the format:

	Pseudomonas alcaligenes adhesin 1 aligned with Pf0-1 LapA: 21.0 % identity; 57.4 % similarity

	inDir: Directory containing clustal format pairwise alignments
	outFile: Filename for output files. Don't add an extension as this script will output a .pdf and a .txt
	with whatever filename you provide.

	Example usage: python3 write_clws_to_file.py Alignments Output/LapA-like_protein_alignments
	"""

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

	txtoutFile = outFile + ".txt"
	pdfoutFile = outFile + ".pdf"

	with open(txtoutFile, "w+") as outf:
		outf.write(out)
	outf.close()

	pdf = FPDF()
	pdf.add_page()
	pdf.set_xy(0, 0)
	pdf.set_font('courier', 'B', 9.5)
	pdf.multi_cell(h=5.0, w=0, txt=out)
	pdf.output(pdfoutFile, 'F')

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])