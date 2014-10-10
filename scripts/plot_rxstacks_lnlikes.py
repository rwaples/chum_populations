#plot rxstacks

import pandas



lnls = pandas.DataFrame.from_csv(path = "/media/Shared/Data/chum/populations/stacks/rxstacks_b10/batch_10.rxstacks_lnls.tsv", sep = "\t")

my_means = array(lnls.Mean)
my_medians = array(lnls.Median)

figure()
hist(my_means, log = 'y', bins = 100)
title("mean ln_like")


figure()
hist(my_medians, log = 'y', bins = 100)
title("median ln_like")
ylim(.01, 10000000)


figure()
hist(my_means[my_means > -100], log = 'y', bins = 100)
title("mean ln_like")
