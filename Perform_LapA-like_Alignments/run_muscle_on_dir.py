import sys
import os
import subprocess



def main(inDir, outDir):
	for file in os.listdir(inDir):
		outfilename = file[:-4] + ".clw"

		bashCommand = "muscle -in %s/%s -out %s/%s -clw" %(inDir, file, outDir, outfilename)
		process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
		output, error = process.communicate()




if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])