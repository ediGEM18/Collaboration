import pickler
import result_sorting
import result_plots

time = list(range(240))

print('reading...')
total_results = pickler.read('res_pickle')

#result sorting
print('sorting...')
average_mRNAs = result_sorting.avg_mRNAs(total_results)
average_proteins = result_sorting.avg_proteins(total_results)

min_mRNAs = result_sorting.min_mRNAs(total_results)
min_proteins = result_sorting.min_proteins(total_results)

max_mRNAs = result_sorting.max_mRNAs(total_results)
max_proteins = result_sorting.max_proteins(total_results)

first = result_sorting.get_run(0, total_results)

every_20 = [(i * 20) - 1 for i in range(1, 13)]
mRNAs_per_20 = [result_sorting.get_mRNAs_at_time(total_results, i) for i in every_20]
proteins_per_20 = [result_sorting.get_proteins_at_time(total_results, i) for i in every_20]

#result plotting
print('plotting max, min, avg')
result_plots.plt_avg_min_max(average_mRNAs, average_proteins, min_mRNAs, min_proteins, max_mRNAs, max_proteins, time)
print('plotting run')
result_plots.plt_run(first, time)
print('finding ratios')
result_plots.end_ratios(average_proteins, min_proteins, max_proteins)
print('dist mRNA per 20')
result_plots.plt_BamA_mRNA_dists(mRNAs_per_20)
result_plots.plt_OmpA_mRNA_dists(mRNAs_per_20)
result_plots.plt_lgA_mRNA_dists(mRNAs_per_20)
print('dist proteins per 20')
result_plots.plt_BamA_protein_dists(proteins_per_20)
result_plots.plt_OmpA_protein_dists(proteins_per_20)
result_plots.plt_lgA_protein_dists(proteins_per_20)
