from matplotlib import pyplot as plt
from matplotlib import ticker

# headers of aggregate.out:
# env	py	comm	weather		jobid		nodes	ranks	r/gpu	evts

ranks_per_gpu_100000_evts = {1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[]}
# data stored in each list: (nnodes, scalable time)
colorlist = ["bright lilac", "pure blue", "greenish blue", "dark mint green", \
             "green yellow", "sunflower yellow", "tangerine", "cerise"]
colorlist = ["xkcd:" + color for color in colorlist]

with open("../work/aggregate.out", "r") as aggfile:
  for line in aggfile.readlines():
    if line.startswith("env"): continue
    parts = line.split(",")
    if len(parts) < 9: continue
    if parts[8].strip() != "100000": continue
    rpg = int(parts[7])
    nodes = int(parts[5])
    envtime = float(parts[0])
    ranks_per_gpu_100000_evts[rpg].append((nodes, envtime))

fig, ax = plt.subplots()
for rpg in range(1,9):
  data = ranks_per_gpu_100000_evts[rpg]
  nnodes, envtime = zip(*data)
  nnodes_sorted = sorted(nnodes)
  order = [nnodes.index(i) for i in nnodes_sorted]
  envtime_sorted = [envtime[i] for i in order]
  ax.plot(nnodes_sorted, envtime_sorted, label="{} ranks/gpu".format(rpg), color=colorlist.pop(0), marker="s")
ax.plot([10,20,40,80], [1650,825,412.5,206.25], label="ideal strong scaling", color="xkcd:goldenrod", linestyle="--")
ax.axvline(x=77, linewidth=1, color='k', linestyle=':')#, color='darkgray')
ax.set_xscale('log',base=10)
ax.set_yscale('log',base=10)
ax.set_xlabel("Nodes requested")
ax.set_ylabel("Scalable wall time (s)")
plt.text(77,2000,"5% of Perlmutter ",horizontalalignment='right')#,color='darkgray',fontweight='bold')
plt.title("Strong scaling on t > 0 execution time")
plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter())
plt.gca().yaxis.set_minor_formatter(ticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:n}"))
plt.gca().xaxis.set_minor_formatter(ticker.ScalarFormatter())
ax.legend(loc='lower left')
plt.show()
