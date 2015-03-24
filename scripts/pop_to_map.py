import pandas as pd

map_data = pd.read_csv('/home/ipseg/Desktop/waples/chum_populations/linkage_map/LEPmap/with_paralogs/collapsed.linkagemap', sep = "\t")
pop_data = pd.read_csv('/home/ipseg/Desktop/waples/chum_populations/data/batch_2/batch_2.final.fst', sep = "\t")

pop_start = pd.read_csv('/home/ipseg/Desktop/waples/chum_populations/data/batch_2/batch_2.plink.map', sep = "\t", skiprows = 1, header = None)




map_data = map_data.iloc[:,[0,1,2]]

map_data.head()
pop_data.head()
pop_start.head()


map_data['catID'] = map_data['Marker'].apply(lambda x: (x.split('_', 1)[0]))
pop_data['catID'] = pop_data['SNP'].apply(lambda x: (x.split('_', 1)[0]))
pop_data['catID']  = pop_data['catID'].apply(lambda x: (x[1:]))
pop_start['catID']  = pop_start[0].apply(lambda x: (x[1:]))


len(set(map_data['catID']).intersection(set(pop_data['catID'])))

pd.merge(pop_data, map_data, on ='catID', how = 'inner')
