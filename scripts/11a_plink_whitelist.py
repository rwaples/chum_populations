import pandas as pd
from __future__ import division

bim = pd.read_csv('/home/ipseg/Desktop/waples/chum_populations/data/batch_2/batch_2.final.bim', sep = "\t", header = None, index_col=False)

bim.columns = ['chr', 'varID', 'cM', 'bp', 'allele_1', 'allele_2']

bim['catID'],bim['pos']  = zip(*bim['varID'].apply(lambda x: x.split('_', 1)))
bim = bim.convert_objects(convert_numeric=True)
bim = bim.sort(columns=(['catID', 'pos']))

whitelist = bim.loc[bim['pos'] != 84]

whitelist.drop_duplicates(subset= 'catID', inplace= True)

whitelist['varID'].to_csv('/home/ipseg/Desktop/waples/chum_populations/data/batch_2/whitelist', index = False)

###########
from scipy.stats import nbinom
def test_minor_homo(aa, ab, bb):
    n = aa+ab+bb
    a = 2*aa + ab
    b = 2*bb + ab
    p = a/(2*n)
    q = b/(2*n)
    p_aa = p**2
    p_not_aa = 2*p*q + q**2
    print(n,(a,b), (p,q), (p_aa,p_not_aa ))
    rv = nbinom(n = n, p = p_not_aa)
    plot(-rv.logpmf(np.arange(10)))
    return(rv.pmf(aa), rv.cdf(aa), -rv.logpmf(aa))
    
test_minor_homo(1,62,100)

def 



rv = nbinom(n = 12, p = .4)
plot(rv.pmf(np.arange(30)))
ax.vlines(x, 0, rv.pmf(x), colors='k', linestyles='-', lw=1,label='frozen pmf')


np.arange(20)
test_minor_homo(0,40,120)




bim.head()
bim.dtypes




bim['pos'].hist(bins = max(bim['pos'])-min(bim['pos']))

bim['catID'].value_counts().hist()



# whitelist criteria
not in pos 84
