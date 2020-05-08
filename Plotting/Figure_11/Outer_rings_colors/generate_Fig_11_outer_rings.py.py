import os
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt

class Lapmap():
	"""object to store information about the LapA and MapA alignments of each LapA-like protein analysed."""
	def __init__(self):
		self.species = ""
		self.LapAorMapA = ""
		self.number = 0
		self.LapApercent_id = 0
		self.LapApercent_sim = 0
		self.MapApercent_id = 0
		self.MapApercent_sim = 0



def main(inDir, outDir):
	"""
	Given a directory containing clustalw format alignments generated by run_muscle_on_dir.py from the "Perform_LapA-like_alignments" dir
	output files into the given outDir containing color maps corresponding to either the similarity or identity of each LapA-like 
	protein with LapA or MapA. Outputs a number of files equal to two times the largest number of LapA-like proteins encoded by any 
	organism in the dataset. The files correspond to one file representing the identity and one file representing the similarity with 
	LapA or MapA. The first ID or similarity file corresponds to the largest LapA-like protein found in each organism. The second 
	file corresponds to the second largest LapA-like protein and so on. Each file only contains color mapping for species with that 
	many LapA-like proteins. i.e. ID3.txt conains identity score color mapping for the third largest LapA-like proteins of all 
	species found to encode at least 3 LapA-like proteins. The resulting similarity .txt files were used as color strip datasets for 
	the 3 outer rings in Fig 11 in the manuscript. In order to fork with a diverging color map, MapA identity and similarity scores
	are converted to negative numbers and the color map is scaled from -100 to 100, with -100 being high MapA identity or similarity
	and 100 being high LapA identity or similarity. For each LapA-like protein, only the highest score of the LapA and MapA alignments
	is kept. e.g. if a given protein is 20% identical to LapA and 40% identical to MapA, the MapA score will be kept as -40 and then 
	converted into an RGBA color using matplotlib cmap.

	inDir: Directory containing ClustalW format alignments generated by run_muscle_on_dir.py from the "Perform_LapA-like_alignments" dir.
	outDir: Directory into which you would like the color map score files to be placed.

	Example usage: python3 generate_Fig_11_outer_rings.py Alignments Score_files
	"""
	results = {}
	files = os.listdir(inDir)
	protein_max = 0
	for file in files:
		Homolog = Lapmap()
		Homolog.LapAorMapA = file[-8:-4]
		Homolog.number = file[-10]
		Homolog.species = file[:-10]
		if int(Homolog.number) > protein_max:
			protein_max = int(Homolog.number)

		seq = ""
		id_count = 0
		sim_count = 0
		with open("%s%s" %(inDir, file), 'r') as f:
			for line in f.readlines():
				if "Pseudomonas" in line:
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

		f.close()

		if Homolog.LapAorMapA == 'LapA':
			Homolog.LapApercent_id = 100 * id_count / len(seq)
			Homolog.LapApercent_sim = 100 * sim_count / len(seq)
		else:
			Homolog.MapApercent_id = 100 * id_count / len(seq)
			Homolog.MapApercent_sim = 100 * sim_count / len(seq)

		if Homolog.species in results.keys():
			if len(results[Homolog.species]) == int(Homolog.number):
				current = results[Homolog.species][int(Homolog.number)-1]
				if current.number == Homolog.number:
					if current.LapApercent_id == 0:
						results[Homolog.species][int(Homolog.number)-1].LapApercent_id = Homolog.LapApercent_id
						results[Homolog.species][int(Homolog.number)-1].LapApercent_sim = Homolog.LapApercent_sim
					elif current.MapApercent_id == 0:
						results[Homolog.species][int(Homolog.number)-1].MapApercent_id = Homolog.MapApercent_id
						results[Homolog.species][int(Homolog.number)-1].MapApercent_sim = Homolog.MapApercent_sim
				else:
					results[Homolog.species].append(Homolog)
			else:
				results[Homolog.species].append(Homolog)
		else:
			results[Homolog.species] = [Homolog]
	
	norm = mpl.colors.Normalize(vmin=-100, vmax=100)
	cmap=plt.cm.PiYG
	
	for i in range(1, protein_max+1):
		with open("%sID%i.txt" % (outDir, i), 'w+') as outfile:
				for species in results.keys():
					for entry in results[species]:
						if int(entry.number) == i:
							if entry.LapApercent_id > entry.MapApercent_id:
								ID = entry.LapApercent_id
							else:
								ID = entry.MapApercent_id * -1
							outfile.write("%s\trgba%s\t%i\n" %(species, cmap(norm(ID), bytes=True), ID))

	for i in range(1, protein_max+1):
		with open("%ssim%i.txt" % (outDir, i), 'w+') as outfile:
				for species in results.keys():
					for entry in results[species]:
						if int(entry.number) == i:
							if entry.LapApercent_sim > entry.MapApercent_sim:
								sim = entry.LapApercent_sim
							else:
								sim = entry.MapApercent_sim * -1
							outfile.write("%s\trgba%s\t%i\n" %(species, cmap(norm(sim), bytes=True), sim))


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])