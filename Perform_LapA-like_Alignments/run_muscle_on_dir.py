import sys
import os
import subprocess



def main(inDir, outDir):

	"""
	N.B. Requires a local installation of the MUSCLE alignment software available here (http://drive5.com/muscle/downloads.htm)

	Perform MUSLCE alignments of the paired FASTA files produced by the create_paired_fastas.py script. Outputs clustalW format alignments to
	specified outdir. Each file is named for example: "Pseudomonas_alcaligenes1_LapA.clw" where "1" indicates that this is the largest protein 
	identified for that species and "LapA" indicates that the sequence of the protein was aligned with that of LapA.
	inDir: Directory containing the paired FASTA files output by the create_paired_fastas.py script.
	outDir: desired destination folder for alignments.

	Example usage: python3 run_muscle_on_dir.py Paired_FASTAs Alignments
	"""

	for file in os.listdir(inDir):
		outfilename = file[:-4] + ".clw"

		bashCommand = "muscle -in %s/%s -out %s/%s -clw" %(inDir, file, outDir, outfilename)
		process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
		output, error = process.communicate()




if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])