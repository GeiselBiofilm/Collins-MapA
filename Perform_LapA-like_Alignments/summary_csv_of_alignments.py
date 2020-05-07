import os
import sys

class Lapmap():
	"""dobject to store """
	def __init__(self):
		self.species = ""
		self.LapAorMapA = ""
		self.number = 0
		self.LapApercent_id = 0
		self.LapApercent_sim = 0
		self.MapApercent_id = 0
		self.MapApercent_sim = 0



def main(inDir, outFile):
	results = {}
	files = os.listdir(inDir)
	for file in files:
		Homolog = Lapmap()
		Homolog.LapAorMapA = file[-8:-4]
		Homolog.number = file[-10]
		Homolog.species = file[:-10]

		seq = ""
		id_count = 0
		sim_count = 0
		with open("%s%s" %(inDir, file), 'r') as f:
			for line in f.readlines():
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
	with open(outFile, 'w+') as outcsv:
		outcsv.write("Species,ID with LapA (%), Similarity with LapA (%), ID with MapA (%), Similarity with MapA (%)")
		for species in results.keys():
			for i in results[species]:
				outcsv.write("\n%s Protein %i,%i,%i,%i,%i" %(species.replace('_',' '), int(i.number), i.LapApercent_id, i.LapApercent_sim, i.MapApercent_id, i.MapApercent_sim))
	outcsv.close()

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])