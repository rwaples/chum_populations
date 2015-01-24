#script to convert haploid mstmap file to lepmap
#11/3/2014
#wes larson
#wlarson1@uw.edu
#script name: haploid_mst_to_lepmap.py


#this script takes a haploid mstmap file and converts it to lepmap format

#name of cross and parent for file naming later
cross_name="Kok_hap"
female_parent="SISSAQ12X03_0001F"
male_parent="Kok_male"

#read in MSTmap, make each line an arrary entry
haploid_mst_file=open("kokanee_haps_mstmap.txt","r")
haploid_mst_array=haploid_mst_file.readlines()
haploid_mst_file.close()

def recode_gens(genotype):
    new_gen=""
    if genotype == "a":
        new_gen="1 1"
    if genotype == "b":
        new_gen="1 2"
    if genotype == "-":
        new_gen="0 0"  
    return(new_gen)        


#get ind id line
mst_header_line=haploid_mst_array[13].rstrip().split("\t")
list_of_inds=mst_header_line[1:len(mst_header_line)]
num_inds=len(list_of_inds)


#trim off header of mst_map file and transpose matrix
mst_array_no_header=haploid_mst_array[13:len(haploid_mst_array)]

#read through genotypes lines, make a dictionary with keys as SNP names
#and values as genotypes for those SNPs
good_locus_dict={}
for i in mst_array_no_header:
    current_locus=i.rstrip().split("\t")
    current_locus_name=current_locus[0]
    #skip first row
    if "locus" in current_locus[0]:
        continue
    #if not first row iterate through genotypes for each locus and convert them to lepmap format
    new_gens_by_locus=[]
    for j in range(1,num_inds+1):
       new_gens_by_locus.append(recode_gens(current_locus[j]))
       
    good_locus_dict[current_locus_name]=new_gens_by_locus
 
        
##############build lepmap output file
lepmap_filename=cross_name+"_lepmap.linkage"
lepmap_out_file=open(lepmap_filename,"w")

#this is sorted by not necessarily in proper numeric order
list_of_snps_in_order=sorted(good_locus_dict.keys())
#print(list_of_snps_in_order)

#write header line
lepmap_out_file.write("#Family\tSample\tSire\tDam\tSex\tSomething\t")
for i in list_of_snps_in_order:
    lepmap_out_file.write(i+"\t")
lepmap_out_file.write("\n")

#write full file
family=cross_name
female=female_parent
male=male_parent

#write female line
lepmap_out_file.write(family+"\t"+female+"\t"+"0"+"\t"+"0"+"\t"+"2"+"\t"+"0"+"\t")
for i in range(1,len(list_of_snps_in_order)+1):
    lepmap_out_file.write("0 0"+"\t")
lepmap_out_file.write("\n")


#write male line
lepmap_out_file.write(family+"\t"+male+"\t"+"0"+"\t"+"0"+"\t"+"1"+"\t"+"0"+"\t")
for i in range(1,len(list_of_snps_in_order)+1):
    lepmap_out_file.write("1 2"+"\t")
lepmap_out_file.write("\n")

#write offspring genotypes
for i in range(0,num_inds):
    ind_name=list_of_inds[i]
    #print(ind_name)
    lepmap_out_file.write(family+"\t"+ind_name+"\t"+male+"\t"+female+"\t"+"0"+"\t"+"0"+"\t")
    for j in list_of_snps_in_order:
        current_gens=good_locus_dict[j]
        #print(current_gens)
        lepmap_out_file.write(current_gens[i]+"\t")
        #print(current_gens)
    lepmap_out_file.write("\n")
        
lepmap_out_file.close()
