import pickler
import result_sorting
import result_plots
import params
from itertools import product

time = list(range(240))

print('reading...')
total_results = pickler.read('res_pickle')

print('searching...')
runs = []
settings = list(product(params.plasmid_copies_BamA, params.init_mRNA_BamA))
for each in settings:
    runs.append(result_sorting.find_runs(each[0], 0, each[1], 15, total_results))

runs_results = []
for each in runs:
    runs_results.append([total_results[i] for i in each])

print('sorting...')
avgs = []
for each in runs_results:
    avgs.append(result_sorting.avg_proteins(each)[0])

print('plotting...')
result_plots.plt_runs_of_interest(avgs, time)
