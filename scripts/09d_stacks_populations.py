# genotypes


#
my_markers_file = "/media/Shared/Data/chum/populations/pstacks/batch_1.ALL_markers.tsv"
my_whitelist_file = "/media/Shared/Data/chum/populations/pstacks/batch_1.whitelist.txt"
def generate_whitelist(markers_file, whitelist_file, n_min):
    with open(markers_file) as INFILE:
        with open(whitelist_file, 'w') as OUTFILE:
            header = next(INFILE).strip().split("\t")
            line_count = 0
            for line in INFILE:
                line_count += 1
                n_genotyped = int(line.strip().split("\t")[4])
                catID = line.strip().split("\t")[2]
                if line_count % 1000 == 0:
                    print(line_count,catID, n_genotyped)
                if n_genotyped > n_min:
                    OUTFILE.write(catID + "\n")
                

generate_whitelist(my_markers_file, my_whitelist_file,150)

"""
populations -b 1 -P /media/Shared/Data/chum/populations/pstacks -s -t 6 -r .5 -p 6 -a .05 --fstats --plink  --genepop --plink --vcf -M /media/Shared/Data/chum/populations/stacks/pop_map/pop_map.txt -W /media/Shared/Data/chum/populations/pstacks/batch_1.whitelist.txt
"""


"""
populations -b 2 -P /media/Shared/Data/chum/populations/pstacks -s -t 6 -r .5 -p 6 -a .05 --fstats --plink --genepop --vcf -M /media/Shared/Data/chum/populations/stacks/pop_map/pop_map.txt

populations -b 2 -P /media/Shared/Data/chum/populations/pstacks -s -t 6 -r .5 -p 6 -a .05 --fstats --plink --genepop --vcf -M /media/Shared/Data/chum/populations/stacks/pop_map/pop_map.txt -W /home/ipseg/Desktop/waples/chum_populations/data/batch_2/on_map.txt

"""



-W /media/Shared/Data/chum/populations/stacks_output/populations/t
