from matplotlib import pyplot as plt
from matplotlib import ticker

# headers of aggregate.out:
# env	py	comm	weather		jobid		nodes	ranks	r/gpu	evts

ranks_per_gpu_100000_evts_kokkos = {1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[]}
ranks_per_gpu_100000_evts_cuda   = {1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[]}
ranks_per_gpu_100000_evts = {"adse13-241":ranks_per_gpu_100000_evts_kokkos, "adse13-249":ranks_per_gpu_100000_evts_cuda}
# data stored in each list: (nnodes, scalable time)
#colorlist = ["bright lilac", "pure blue", "greenish blue", "dark mint green", "green yellow", "sunflower yellow", "tangerine", "cerise"]
colorlist = ["greenish blue", "dark mint green", "green yellow", "sunflower yellow", "tangerine", "cerise"]
colorlist = colorlist + colorlist
colorlist = ["xkcd:" + color for color in colorlist]
linetypes = {"adse13-241":"-","adse13-249":"--"}

for target in ("adse13-241","adse13-249"):
  with open("../../{}/work/aggregate.out".format(target), "r") as aggfile:
    for line in aggfile.readlines():
      if line.startswith("env"): continue
      parts = line.split(",")
      if len(parts) < 9: continue
      if parts[8].strip() != "100000": continue
      rpg = int(parts[7])
      nodes = int(parts[5])
      envtime = float(parts[0])
      ranks_per_gpu_100000_evts[target][rpg].append((nodes, envtime))

fig, ax = plt.subplots()
for target in ("adse13-241","adse13-249"):
  for rpg in range(3,9):
    code = "Kokkos" if target == "adse13-241" else "CUDA"
    data = ranks_per_gpu_100000_evts[target][rpg]
    nnodes, envtime = zip(*data)
    nnodes_sorted = sorted(nnodes)
    order = [nnodes.index(i) for i in nnodes_sorted]
    envtime_sorted = [envtime[i] for i in order]
    ax.plot(nnodes_sorted, envtime_sorted, label="{} ranks/gpu, {}".format(rpg, code), color=colorlist.pop(0), marker="s", linestyle=linetypes[target])
  if target == "adse13-241":
    ax.plot([10,20,40,80], [1600,800,400,200], label="ideal strong scaling", color="k", linestyle="-")
ax.axvline(x=77, linewidth=1, color='k', linestyle=':')#, color='darkgray')
ax.set_xscale('log',base=10)
ax.set_yscale('log',base=10)
ax.set_ylim(ymin=100)
ax.set_xlabel("Nodes requested")
ax.set_ylabel("Scalable wall time (s)")
plt.text(77,2000,"5% of Perlmutter ",horizontalalignment='right')#,color='darkgray',fontweight='bold')
plt.title("Strong scaling on t > 0 execution time")
plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter())
plt.gca().yaxis.set_minor_formatter(ticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:n}"))
plt.gca().xaxis.set_minor_formatter(ticker.ScalarFormatter())
ax.legend(loc='lower left',ncol=2)
plt.show()
