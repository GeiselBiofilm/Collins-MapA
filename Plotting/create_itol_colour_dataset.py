import re

def main():
	
	with open("new_mlst_ML_tree.tre", "r") as treefile:
		
		leaf_dict = dict.fromkeys(re.findall(r'([^(),;\s]+)', treefile.read())) 
		# find all instances of stuff not containing (due to '^') (),; or whitespace (\s)


	treefile.close()

	with open("Pseudomonas_no_lapDG.csv", 'r') as laplessfile:
		nolapDG_pseudos = []
		for line in laplessfile.readlines():
			if "Candidatus" in line:
				nolapDG_pseudos.append("_".join(line.split()[1:3]).strip())
			else:
				nolapDG_pseudos.append("_".join(line.split()[:2]).strip().replace("[", "").replace("]",""))

	laplessfile.close()

	with open("Pseudomonas_LapA_counts.csv", 'r') as LapAcountsfile:

		LapAorgs = []
		for line in LapAcountsfile.readlines():
			if not "organism" in line:
				elements = line.split(",")
				Species = "_".join(elements[1].split()[:2])
				if "sp." not in Species:
					LapAorgs.append((Species, elements[2].strip()))
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
				#print("%s is the query and %s is the dict entry" %(i, leaf_dict[Species]))
		else:
			print(Species)
	leaf_dict["Azotobacter_vinelandii"] = "No_LapD/G_found"


	colour_dict = {"No_LapD/G_found": "rgb(0,0,0)", "Zero": "rgb(255,255,255)", "One": "rgb(254,178,76)", "Multiple": "rgb(240,59,32)"}
	with open("tree_lap_status.txt", 'w+') as outfile:
		outfile.write("Species\tColour\tLap_status\n")
		for k, v in leaf_dict.items():
			outfile.write("%s\t%s\t%s\n" %(k, colour_dict[v.strip('"')], v))
	outfile.close()


if __name__ == '__main__':
	main()