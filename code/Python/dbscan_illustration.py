import matplotlib as mpl
mpl.rcParams['font.family'] = "Liberation Sans"

import matplotlib.pyplot as plt

from matplotlib.patches import FancyArrowPatch

from collections import defaultdict
#core
fig = plt.gcf()
#fig = plt.figure(figsize=(4, 4), dpi=100)
#fig.set_size_inches(4, 4)
ax = fig.gca()

ax.set_xlim((0,10))
ax.set_ylim((0,6))

ax.set_xticks([])
ax.set_yticks([])

ax.set_aspect(1)

points = []
points.append(['r','Core',     [(3,2.1),(4,1.9),(3.8,3.2),(5,3),(5.5,1.5)]])
points.append(['y','Reachable',[(6.8,1.2),(6,4),(1.7,1.8)]])
points.append(['b','Noise',    [(1,4.2),(8.5,3)]])

L = defaultdict(list)
for col, label, pts in points:
    for pt in pts:
        cir=plt.Circle(pt,.2,facecolor=col, edgecolor="k")
        ax.add_artist(cir)
        L[label].append(cir)
        
        cir=plt.Circle(pt,1.7,color=col, fill=False)
        ax.add_artist(cir)

#opts = {'head_width':0.3, 'head_length':0.2, 'width':0.1, 'length_includes_head':True, 'color':'k'}

#reachables
#plt.arrow(3.8, 3, -0.9, 0, **opts)
#plt.arrow(6.15, 4.15, 0.7, 0.7, **opts)
#plt.arrow(6.7, 2.5, 0.6, 0, **opts)

#core
#plt.arrow(4.3, 3, 0.5, 0, **opts)
#plt.arrow(4.7, 3, -0.5, 0, **opts)

#plt.arrow(4, 3.2, 0.5, 1.5, **opts)
#plt.arrow(4.5, 4.5, -0.5, -1.5, **opts)

#ax.add_patch(FancyArrowPatch(posA=(9,1), posB=(1,9), arrowstyle='<|-|>', mutation_scale=30., linewidth=0.9))

plt.legend([L["Core"][0], L["Noise"][0], L["Reachable"][0]],["Core", "Noise", "Reachable"])


#circle1=plt.Circle((0,0),.2,color='r')
#circle2=plt.Circle((.5,.5),.2,color='b')
#circle3=plt.Circle((1,1),.2,color='g',clip_on=False)

#fig = plt.gcf()

#fig.gca().add_artist(circle1)
#fig.gca().add_artist(circle2)
#fig.gca().add_artist(circle3)

#plt.legend([circle1],["red"])
#plt.show()
plt.savefig("test.pdf",bbox_inches='tight', pad_inches=0)

#fig.show()
