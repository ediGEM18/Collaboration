import systems
import pickler
import result_sorting
import result_plots
from tqdm import tqdm
from params import settings
from scipy.integrate import odeint

time = list(range((240))) #per minute for 24hrs

total_results = []

#screen
print('simulating...')
for param_set in tqdm(settings):

  model = systems.Systems(param_set[:15])

  #run
  result = odeint(model.sys, param_set[15:], time)

  #save
  total_results.append((param_set, result))

#serialize
print('serializing results for later')
pickler.write('res_pickle', total_results)

#result sorting
print('sorting results')
average_mRNAs = result_sorting.avg_mRNAs(total_results)
average_proteins = result_sorting.avg_proteins(total_results)

min_mRNAs = result_sorting.min_mRNAs(total_results)
min_proteins = result_sorting.min_proteins(total_results)

max_mRNAs = result_sorting.max_mRNAs(total_results)
max_proteins = result_sorting.max_proteins(total_results)

first = result_sorting.get_run(0, total_results)

mRNAs = result_sorting.get_mRNAs(total_results)

proteins = result_sorting.get_proteins(total_results)

#result plotting
print('plotting results')
result_plots.plt_avg_min_max(average_mRNAs, average_proteins, min_mRNAs, min_proteins, max_mRNAs, max_proteins, time)

result_plots.plt_run(first, time)

result_plots.end_ratios(average_proteins, min_proteins, max_proteins)
