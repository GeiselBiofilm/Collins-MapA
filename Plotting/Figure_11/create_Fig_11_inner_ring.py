import re
import sys

def main(TreeFile, No_LapDGFile, LapACountFile, OutFile):

	"""
	Creates color strip dataset for iToL given a treefile to annotate and files containing information about the presence/absence of LapD and 
	LapG and the number of LapA-like proteins encoded. The input files can contain taxa not included in the tree.
	TreeFile: a treefile in Newick format
	No_LapDGFile: .csv or .txt with one species per row that does not encode LapD and LapG
	LapACountFile: .csv with two columns: Species, Number of LapA-like proteins encoded
	OutFile: name of file to output .txt format dataset for iToL annotation
	Example usage:	python3 create_itol_color_dataset Data/my_tree.tre Data/Pseudomonas_no_lapDG.txt Data/Pseudomonas_LapA_counts.csv Output/tree_lap_status.txt
	"""
	
	with open(TreeFile, "r") as treefile:
		
		leaf_dict = dict.fromkeys(re.findall(r'([^(),;\s]+)', treefile.read())) 


	treefile.close()

	with open(No_LapDGFile, 'r') as laplessfile:
		nolapDG_pseudos = []
		for line in laplessfile.readlines():
			if "Candidatus" in line:
				nolapDG_pseudos.append("_".join(line.split()[1:3]).strip())
			else:
				nolapDG_pseudos.append("_".join(line.split()[:2]).strip().replace("[", "").replace("]",""))

	laplessfile.close()

	with open(LapACountFile, 'r') as LapAcountsfile:

		LapAorgs = []
		for line in LapAcountsfile.readlines():
			if not "organism" in line:
				elements = line.split(",")
				Species = "_".join(elements[0].split()[:2])
				if "sp." not in Species:
					LapAorgs.append((Species, elements[1].strip()))
	LapAcountsfile.close()

	for i in nolapDG_pseudos:
		if i in leaf_dict.keys():
			leaf_dict[i] = "No_LapD/G_found"

	for i in LapAorgs:
		Species = i[0].strip('"')
		if Species in leaf_dict.keys():
			if leaf_dict[Species] == None:
				leaf_dict[Species] = i[1]
			if leaf_dict[Species] == "No_LapD/G_found":
				leaf_dict[Species] = i[1]
			elif leaf_dict[Species] == "Multiple":
				continue
			elif leaf_dict[Species] == "One" and i[1] == "Multiple":
				leaf_dict[Species] = i[1]
			elif leaf_dict[Species] == "Zero" and i[1] in ["Single", "Multiple"]:
				leaf_dict[Species] = i[1]
			else:
				continue
		else:
			print("Taxon %s was not found in the tree file." % Species)
		leaf_dict["Cellvibrio_japonicus"] = "Zero"


	colour_dict = {"No_LapD/G_found": "rgb(0,0,0)", "Zero_adhesins": "rgb(255,255,255)", "One_adhesin": "rgb(254,178,76)", "Multiple_adhesin": "rgb(240,59,32)"}
	with open(OutFile, 'w+') as outfile:
		outfile.write("Species\tColour\tLap_status\n")
		for k, v in leaf_dict.items():
			outfile.write("%s\t%s\t%s\n" %(k, colour_dict[v.strip('"')], v))
	outfile.close()


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])