#!/usr/bin/python
# -*- coding: utf8 -*-
import os.path
from scipy import arange
import sys

from interface import *

ger_small_file="/home/seydanator/Desktop/SVN/egtdata/tweets.germany.gt"
#  2MB
# 57.554
ger_small= "a56c48c3df72452177dce28efd790ddc"
ger_small_d = "a56c48c3df72452177dce28efd790ddcdisplace"


ger_file = "/home/seydanator/Desktop/SVN/egtdata/complete_2008-2011.germany.egt"
#140MB
# 1.254.298
ger = "90bda4cc98289c9d3c231127f8253189"

usa_small_file="/home/seydanator/Desktop/SVN/egtdata/complete_2008-2011.usa.egt"
#719MB
# 5.976.723
usa_small= "3f87f7a7aa3d8ba083407ff956f1ada1"
usa_small_1M= "3f87f7a7aa3d8ba083407ff956f1ada1_1M"
#1.086.678

usa_file = "/home/seydanator/Desktop/SVN/egtdata/tweets_110716_110815_text.us.egt"
#1200MB
# 14.531.622
usa = "32f974a4edb4eb8f1e2aa366fa1259b1"
# 279.889.256
#  16.777.216

dir_col = {
    ger_small : "GER_SMALL",
    ger_small_d : "GER_SMALL_D",
    ger       : "GER",
    usa_small : "USA_SMALL",
    usa_small_1M : "USA_SMALL_1M",
    usa       : "USA",
    "blobs"    : "",
    "noisy_circles":"",
    "noisy_moons":"",
    "random":""
}

#
# Max 1.333.000 points because of 16MB file limit in MongoDB
#

#insert_file(ger_small_file)
#insert_file(ger_file)
#insert_file(usa_small_file)
#insert_file(usa_file)

#new_collection(ger_small, ger_small+"displace", displace=True)

#triang_combiner(ger_small)
#triang_combiner(ger_small+"displace")
#triang_combiner(ger)
#triang_combiner(usa_small)
#triang_combiner(usa)

#sys.exit()

#new_collection(usa_small, usa_small+"_1M", bbox=[0,0,0,0], sampling=(2,9), limit=0, count=False, combine=True, displace=False)
#triang_combiner(usa_small+"_1M")
#sys.exit()


#same_point("ger_test200k")
#   .207
#200.000

#same_point("ger_test100k")
#   .609
#114.028

#same_point("ger_test2")
#  1.457
#418.100

#ger
#uniques:    3.229
#points: 1.254.298

#usa_small
#13757
#5976723

# 
#frequency_map(ger_small, borders_by_collection(ger_small), legend=True)
#frequency_map(ger, borders_by_collection(ger), legend=True)
#frequency_map(usa_small_1M, borders_by_collection(usa_small_1M), legend=True)
#edge_pic(ger_small, cm="winter")
#edge_pic(ger, cm="winter")
#edge_pic(usa_small_1M, cm="winter") 

#sys.exit()



#cluster_histogram([home+dir_col[col]+"/"+efile])

#from datetime import datetime
#startTime = datetime.now()
#print eval_clustering(col, home+dir_col[col]+"/"+efile, 950, 8000)
#print eval_clustering(col, home+dir_col[col]+"/"+efile, 10, 1600, True)

#print eval_clustering(col, home+dir_col[col]+"/"+efile, 10, 1600, False)
#eval_clustering(col, home+dir_col[col]+"/"+efile, 980, 8000)
#eval_clustering(col, home+dir_col[col]+"/"+efile, 990, 8000)

#delta = datetime.now()-startTime
#if (delta.seconds > 120):
#    print delta.seconds/float(60), "minutes"
#else:
#    print delta.seconds, "seconds"


#D
#DB
#C
#SW

#Geo
#Jac

home = "/home/seydanator/Desktop/SVN/code/"

#different clusters
import matplotlib as mpl
mpl.rcParams['font.family'] = "Liberation Serif"

import matplotlib.pyplot as plt
import numpy as np

def get_scores(make, col, efile):
    stats = {'Geo': {'SW': [], 'C': [], 'DB': [], 'D': []}, 'JBi': {'SW': [], 'C': [], 'DB': [], 'D': []}}
        
    for m in make:
        skip = 0
        count_limit = 1000000000
        count = count_limit
        count_ = 0
        equal = False
        
        #while(count >= count_limit):
        #    skip += 500
        #    count_ = int(eval_clustering(col, home+dir_col[col]+"/"+efile, m, skip, True)["count"]["count"])
            
        #    if count == count_:
        #        equal = True
        #        break
        #    count = count_
        
        count = int(eval_clustering(col, home+dir_col[col]+"/"+efile, m, skip, True)["count"]["count"])
        
        if count > count_limit:
            equal = True
        
        if not equal:
            stat = eval_clustering(col, home+dir_col[col]+"/"+efile, m, skip, False)
            for typ, val in stat.iteritems():
                for qm, score in val.iteritems():
                    stats[typ][qm].append(float(score))
    return stats


def eval_comp(col, files, make = [10, 20, 30, 40, 50]):
    #f, (geo_ax, jac_ax) = plt.subplots(2, sharex=True, sharey=True)
    ##f, geo_ax= plt.subplots()
    
    #plt.xlabel('Quality Measure')
    #geo_ax.set_ylabel('Location Scores')
    #jac_ax.set_ylabel('Jaccard Scores')
    #
    #geo_ax.grid()
    #jac_ax.grid()

    #n_groups = 4
    #index = np.arange(n_groups)
    #bar_width = 1/float(len(files)+1)
    
    #opacity = 0.4
    #error_config = {'ecolor': '0.3'}
    ##color=["b","g","r","y"]
    #cm = plt.cm.get_cmap("rainbow")#cubehelix")
    #
    j_means = []
    j_std = []
    
    g_means = []
    g_std = []
    
    for i, efile in enumerate(files):
        j_means.append([])
        j_std.append([])
         
        g_means.append([])
        g_std.append([])
        
        #break
        
        stats = get_scores(make, col, efile)
        #print stats
        
        
        for qm, scores in stats["Geo"].iteritems():
            mean = np.ma.masked_invalid(scores).mean()
            if not np.ma.count_masked(mean):
                g_means[i].append(mean)
                g_std[i].append(np.ma.masked_invalid(scores).std())
            else:
                g_means[i].append(0.0)
                g_std[i].append(0.0)
        
        for qm, scores in stats["JBi"].iteritems():
            mean = np.ma.masked_invalid(scores).mean()
            if not np.ma.count_masked(mean):
                j_means[i].append(mean)
                j_std[i].append(np.ma.masked_invalid(scores).std())
            else:
                j_means[i].append(0.0)
                j_std[i].append(0.0)
    
    
    print j_means
    print j_std
    
    print g_means
    print g_std
    
    print "location"
    for i, efile in enumerate(files):
        print "%s & %.2e & %.2e & %.2e & %.2e\\\\" % (efile, g_means[i][0], g_means[i][1], g_means[i][2], g_means[i][3])
    print "text"
    for i, efile in enumerate(files):
        print "%s & %.2e & %.2e & %.2e & %.2e\\\\" % (efile, j_means[i][0], j_means[i][1], j_means[i][2], j_means[i][3])
    
    #DB C SW D
    #efile & 1.35e-02 & 3.18e-09 & 9.87e-01& 3.13e+00
    
    
    #j_means=[[0.0, 0.62174666666666667, -0.053039099999999999, 0.0], [0.0, 0.62304533333333334, -0.057713966666666672, 0.0], [0.0, 0.728271, -0.118867, 0.0], [0.0, 0.50574433333333335, -0.090671500000000002, 0.0]]
    
    #j_std=[[0.0, 0.017993807552105858, 0.011821196722272524, 0.0], [0.0, 0.019107415773870481, 0.0096527306023851209, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.10828102195162771, 0.011682913534730969, 0.0]]
    
    #g_means=[[0.13513900000000001, 0.00911818, 0.32162033333333334, 0.00073883399999999995], [0.148448, 0.0085566400000000008, 0.30092066666666667, 0.00068666933333333326], [37.1616, 0.00557852, 0.43602200000000002, 3.2653099999999999e-12], [3.7070700000000003, 0.004697743333333333, 0.59381033333333333, 6.491913333333332e-12]]
    
    #g_std=[[0.001865104286628497, 0.00086073188268279328, 0.07201798690759291, 0.00019194617108449967], [0.0069335815179939024, 0.00067720244196251987, 0.073248796760690155, 1.6565789835950764e-05], [0.0, 0.0, 0.0, 0.0], [1.4463635771824457, 0.0050081985185915128, 0.018190218916280868, 2.6178531075461219e-12]]
    
    # copy lists in order to normalize to 0...1 and show also the original values
    #import copy
    #j_plot = copy.deepcopy(j_means)
    #g_plot = copy.deepcopy(g_means)
    #print g_means, g_plot
    
    #normalizing DB values
    #g_db = 0.0
    #g_d = 0.0
    #for line in g_means:
    #    g_db = max(g_db, line[0])
    #    g_d = max(g_d, line[3])
    #
    #j_db = 0.0
    #j_d = 0.0
    #for line in j_means:
    #    j_db = max(j_db, line[0])
    #    j_d = max(j_d, line[3])
    #
    #if g_db != 0.0:
    #    for line in g_plot:
    #        line[0]/=g_db
    #if j_db != 0.0:
    #    for line in j_plot:
    #        line[0]/=j_db
    #
    #if g_d != 0.0:
    #    for line in g_plot:
    #        line[3]/=g_d
    #if j_d != 0.0:
    #    for line in j_plot:
    #        line[3]/=j_d
    
    
    #print g_means, g_plot
    
    #for i, efile in enumerate(files):
    #    rects1 = geo_ax.bar(index+i*bar_width, g_plot[i], bar_width,
    #                     alpha=opacity,
    #                     color=cm(float(i)/len(files)),
    #                     #yerr=g_std[i],
    #                     error_kw=error_config,
    #                     label=efile)
    #    for j, rect in enumerate(rects1):
    #        height = rect.get_height()
    #        geo_ax.text(rect.get_x()+rect.get_width()/2., 0.1, '%.2e'%g_means[i][j], ha='center', va='bottom', rotation='vertical')#1.05*height
    #                
    #    rects2 = jac_ax.bar(index+i*bar_width, j_plot[i], bar_width,
    #                     alpha=opacity,
    #                     color=cm(float(i)/len(files)),
    #                     #yerr=j_std[i],
    #                     error_kw=error_config,
    #                     label=efile)
    #    for j, rect in enumerate(rects2):
    #        height = rect.get_height()
    #        jac_ax.text(rect.get_x()+rect.get_width()/2., 0.1, '%.2e'%j_means[i][j], ha='center', va='bottom', rotation='vertical')#1.05*height
    #
    #plt.xticks(index + (len(files)*bar_width) / 2.0, [u'0 ≤ DB ≤ ∞\nminimize', u'0 ≤ C ≤ 1\nminimize', u'-1 ≤ SW ≤ 1\nmaximize', u'0 ≤ D ≤ ∞\nmaximize'])
    #plt.legend(ncol=2, loc=8, mode="expand", bbox_to_anchor=(0., -0.75, 1., .02))
    
    #plt.tight_layout()
    ##plt.show()
    #plt.savefig("compare.pdf", bbox_inches='tight', pad_inches=0)
    #plt.close()


# GER Details

#inttt = dir_col[ger]+"/RW_0.986000_4_2.000000_0.100000_35000"
#outtt = inttt + "_testttt.pdf"
#bor = [52.3,52.6, 12.9,13.8]
#plot(ger, clusters_in_borders(ger, bor, inttt), show=False, output=outtt, cluster_count=27, borders=None, scatter=False, colored=True, shape_size=1000, legend=False)

#inttt = dir_col[ger]+"/COM_4_1.0_1000_0.1_b"
#outtt = inttt + "_testttt.pdf"
#bor = [52.13,52.683043, 12.6,14.0]
##min_lat, max_lat, min_lon, max_lon = borders
#plot(ger, parse_clustering(inttt), show=False, output=outtt, cluster_count=27, borders=bor, scatter=True, colored=True, shape_size=2000, legend=False)

# GER SMALL Details
col = ger_small
todo = ["/RW_0.969400_4_1.000000_0.100000_40000", "/COM_4_1.0_300_0.1_b"]
#COM_4_1.0_300_0.1_b
#clusters 888
#clustered 52418
#unclustered 5136
#
#RW_0.969400_4_1.000000_0.100000_40000
#clusters 2082
#clustered 43175
#unclustered 14379

for ttt in todo:
    inttt = dir_col[col]+ttt
    
    for key in ["acakfilm"]:#["berlin", "brandenburg", "hamburg", "potsdam", "hessen", "stuttgart", "frankfurt"]:
        print col, key
        
        frequency_map(col, borders_by_collection(col), clustering=parse_clustering(inttt), keywords=[key], legend=True, out=ttt[1:]+"_"+key, bin_scale=4)
sys.exit()

# USA Details
col = usa_small_1M
todo = ["/RW_0.951000_4_0.100000_0.100000_300000", "/COM_4_1.0_7000_0.5_b"]
for ttt in todo:
    inttt = dir_col[col]+ttt
    
    for key in ["nature", "desert", "coast"]:
        print col, key
        
        frequency_map(col, borders_by_collection(col), clustering=parse_clustering(inttt), keywords=[key], legend=True, out=ttt[1:]+"_"+key)
        continue
        
        outtt = dir_col[col]+"_"+key+"_"+ttt[:3]
    #bor = [36.8, 43.8, -79.0,-69.5]#general, shapesize = 20000
    #bor = [38.706946, 39.049052,-77.311707,-76.784363]#washington, shapesize = 500
    #bor = [40.229218, 40.984045,-74.48822,-73.675232]#new york, ss=500
    
        plot(col, parse_clustering(inttt), show=False, output=outtt, cluster_count=None, borders=borders_by_collection(col), scatter=False, colored=True, shape_size=50000, legend=False, keywords = [key])

    #plot(col, clusters_in_borders(col, bor, inttt, 5000), show=False, output=outtt, cluster_count=27, borders=bor, scatter=False, colored=True, shape_size=10000, legend=False)




# GER Details
col = ger
todo = ["/RW_0.972000_4_1.000000_0.100000_35000", "/COM_4_1.0_1000_0.1_b"]
#COM_4_1.0_1000_0.1_b
#clusters 200
#clustered 1253613
#unclustered 685
#
#RW_0.972000_4_1.000000_0.100000_35000
#clusters 194
#clustered 955138
#unclustered 299160

for ttt in todo:
    inttt = dir_col[col]+ttt
    
    for key in ["schluchsee", "sanssouci", "oster", "berlin", "hamburg", "autobahn"]:
        print col, key
        
        frequency_map(col, borders_by_collection(col), clustering=parse_clustering(inttt), keywords=[key], legend=True, out=ttt[1:]+"_"+key, bin_scale=4)
        continue



""""""
ger_files = [\
#"COM_4_1.0_1000_0.9_b", \
#"COM_4_1.0_1000_0.1_b", \
#"COM_4_1.0_500_0.5_b"]
#"RW_0.951000_4_0.100000_0.100000_35000",\
#"RW_0.955000_4_0.100000_0.100000_35000",\
#"RW_0.965000_4_0.100000_0.100000_35000",\
#"RW_0.972000_4_1.000000_0.100000_35000",\
#"RW_0.976000_4_1.000000_0.100000_35000",\
#"RW_0.981000_4_1.000000_0.100000_35000"
#"RW_0.983000_4_2.000000_0.100000_35000",\
#"RW_0.986000_4_2.000000_0.100000_35000",\
#"RW_0.992000_4_2.000000_0.100000_35000",\
#"random", "random2"
]

"""
location                                    min         min     max         max
COM_4_1.0_1000_0.9_b                    & 2.88e-02 & 3.76e-06 & 9.82e-01 & 7.19e-02\\
COM_4_1.0_1000_0.1_b                    & 1.20e-01 & 2.52e-05 & 9.05e-01 & 2.48e-03\\
COM_4_1.0_500_0.5_b                     & 1.71e-02 & 7.98e-06 & 9.76e-01 & 4.44e-02\\
RW_0.951000_4_0.100000_0.100000_35000   & 3.35e-01 & 3.27e-05 & 8.77e-01 & 5.65e-08\\
RW_0.955000_4_0.100000_0.100000_35000   & 6.09e-01 & 3.21e-05 & 8.89e-01 & 2.05e-09\\
RW_0.965000_4_0.100000_0.100000_35000   & 9.78e-01 & 8.62e-05 & 7.51e-01 & 3.29e-10\\
RW_0.972000_4_1.000000_0.100000_35000   & 2.20e-01 & 3.52e-05 & 8.93e-01 & 5.45e-08\\
RW_0.976000_4_1.000000_0.100000_35000   & 1.97e+03 & 6.84e-05 & 8.39e-01 & 5.45e-08\\
RW_0.981000_4_1.000000_0.100000_35000   & 7.22e+01 & 2.06e-04 & 7.59e-01 & 2.95e-09\\
RW_0.983000_4_2.000000_0.100000_35000   & 2.08e+00 & 1.26e-04 & 7.49e-01 & 7.51e-09\\
RW_0.986000_4_2.000000_0.100000_35000   & 2.43e+01 & 3.54e-04 & 7.05e-01 & 1.61e-09\\
RW_0.992000_4_2.000000_0.100000_35000   & 1.59e+04 & 1.07e-03 & 5.73e-01 & 4.35e-10\\
random                                  & 0.00e+00 & 0.00e+00 & 0.00e+00 & 0.00e+00\\
random2                                 & 1.07e+08 & 2.64e-01 & -1.01e-01 & 0.00e+00\\

text
COM_4_1.0_1000_0.9_b                    & 1.29e+00 & 1.95e-01 & 5.34e-01 & 0.00e+00\\
COM_4_1.0_1000_0.1_b                    & 1.17e+00 & 1.90e-01 & 5.60e-01 & 0.00e+00\\
COM_4_1.0_500_0.5_b                     & 1.22e+00 & 1.94e-01 & 5.47e-01 & 0.00e+00\\
RW_0.951000_4_0.100000_0.100000_35000   & 1.18e+00 & 2.84e-01 & 3.19e-01 & 0.00e+00\\
RW_0.955000_4_0.100000_0.100000_35000   & 1.30e+00 & 2.77e-01 & 4.04e-01 & 0.00e+00\\
RW_0.965000_4_0.100000_0.100000_35000   & 1.59e+00 & 2.91e-01 & 3.86e-01 & 0.00e+00\\
RW_0.972000_4_1.000000_0.100000_35000   & 1.13e+00 & 2.78e-01 & 4.23e-01 & 0.00e+00\\
RW_0.976000_4_1.000000_0.100000_35000   & 1.25e+00 & 2.83e-01 & 3.41e-01 & 0.00e+00\\
RW_0.981000_4_1.000000_0.100000_35000   & 1.62e+00 & 3.02e-01 & 3.54e-01 & 0.00e+00\\
RW_0.983000_4_2.000000_0.100000_35000   & 1.47e+00 & 3.37e-01 & 2.93e-01 & 0.00e+00\\
RW_0.986000_4_2.000000_0.100000_35000   & 1.56e+00 & 3.52e-01 & 2.54e-01 & 0.00e+00\\
RW_0.992000_4_2.000000_0.100000_35000   & 1.42e+00 & 3.37e-01 & 3.10e-01 & 0.00e+00\\
random                                  & 0.00e+00 & 0.00e+00 & 0.00e+00 & 0.00e+00\\
random2                                 & 1.57e+00 & 8.68e-01 & -3.34e-02 & 0.00e+00\\
"""



ger_small_files = [\
#"COM_4_1.0_300_0.5_b", \
#"COM_4_1.0_300_0.1_b", \
#"COM_4_1.0_200_0.5_b", \
#"COM_4_1.0_200_0.1_b"]
#"RW_0.960000_4_0.100000_0.100000_40000",\
#"RW_0.964000_4_0.100000_0.100000_40000",\
#"RW_0.965000_4_0.100000_0.100000_40000",\
#"RW_0.970000_4_1.000000_0.100000_40000",\
#"RW_0.976000_4_2.000000_0.100000_40000"
#"RW_0.969400_4_1.000000_0.100000_40000",
#"RW_0.969600_4_1.000000_0.100000_40000",
#"RW_0.969900_4_1.000000_0.100000_40000",
#"random", "random2"
]

"""
#DB C SW D
location                                    min         min     max         max
COM_4_1.0_300_0.5_b                     & 1.58e-01 & 2.72e-04 & 7.02e-01 & 3.54e-04\\
COM_4_1.0_300_0.1_b                     & 2.25e-01 & 1.17e-04 & 7.14e-01 & 2.77e-04\\
COM_4_1.0_200_0.5_b                     & 1.37e-01 & 6.21e-05 & 7.51e-01 & 1.20e-03\\
COM_4_1.0_200_0.1_b                     & 1.63e-01 & 5.07e-05 & 7.57e-01 & 6.95e-04\\
RW_0.960000_4_0.100000_0.100000_40000   & 1.14e+00 & 1.10e-04 & 7.64e-01 & 7.59e-10\\
RW_0.964000_4_0.100000_0.100000_40000   & 2.29e+00 & 3.42e-04 & 5.09e-01 & 7.99e-12\\
RW_0.965000_4_0.100000_0.100000_40000   & 2.48e+01 & 6.57e-04 & 4.20e-01 & 1.15e-10\\
RW_0.970000_4_1.000000_0.100000_40000   & 1.20e+00 & 1.78e-04 & 6.86e-01 & 3.60e-09\\
RW_0.976000_4_2.000000_0.100000_40000   & 2.60e+00 & 6.34e-04 & 4.72e-01 & 2.21e-10\\
RW_0.969400_4_1.000000_0.100000_40000   & 4.19e+00 & 1.66e-04 & 7.32e-01 & 9.30e-08\\
RW_0.969600_4_1.000000_0.100000_40000   & 6.92e+00 & 1.72e-04 & 7.18e-01 & 4.24e-09\\
RW_0.969900_4_1.000000_0.100000_40000   & 2.74e+00 & 1.79e-04 & 6.93e-01 & 3.90e-09\\
random                                  & 3.57e+08 & 2.77e-01 & -2.70e-01 & 0.00e+00\\
random2                                 & 1.63e+08 & 2.79e-01 & -2.34e-01 & 0.00e+00\\

text
COM_4_1.0_300_0.5_b                     & 2.75e+00 & 5.48e-01 & -4.83e-03 & 0.00e+00\\
COM_4_1.0_300_0.1_b                     & 2.83e+00 & 5.29e-01 & 1.87e-03 & 0.00e+00\\
COM_4_1.0_200_0.5_b                     & 2.82e+00 & 5.54e-01 & -9.90e-03 & 0.00e+00\\
COM_4_1.0_200_0.1_b                     & 3.09e+00 & 5.35e-01 & -2.79e-03 & 0.00e+00\\
RW_0.960000_4_0.100000_0.100000_40000   & 3.54e+00 & 4.40e-01 & -7.30e-03 & 0.00e+00\\
RW_0.964000_4_0.100000_0.100000_40000   & 3.80e+00 & 5.92e-01 & -5.64e-02 & 0.00e+00\\
RW_0.965000_4_0.100000_0.100000_40000   & 3.63e+00 & 6.19e-01 & -5.05e-02 & 0.00e+00\\
RW_0.970000_4_1.000000_0.100000_40000   & 3.57e+00 & 4.00e-01 & 4.55e-03 & 0.00e+00\\
RW_0.976000_4_2.000000_0.100000_40000   & 3.52e+00 & 4.63e-01 & -2.01e-02 & 0.00e+00\\
RW_0.969400_4_1.000000_0.100000_40000   & 3.46e+00 & 3.71e-01 & 1.89e-02 & 0.00e+00\\
RW_0.969600_4_1.000000_0.100000_40000   & 3.38e+00 & 3.81e-01 & 1.28e-02 & 0.00e+00\\
RW_0.969900_4_1.000000_0.100000_40000   & 3.45e+00 & 3.95e-01 & 8.18e-03 & 0.00e+00\\
random                                  & 6.93e+00 & 9.65e-01 & -2.53e-02 & 0.00e+00\\
random2                                 & 7.28e+00 & 9.63e-01 & -2.20e-02 & 0.00e+00\\
"""


usa_small_files = [\
#"COM_4_1.0_3000_0.9_b",\
#"COM_4_1.0_3000_0.5_b",\
#"COM_4_1.0_3000_0.1_b",\
#"RW_0.951000_4_0.100000_0.100000_300000",\
#"RW_0.955000_4_0.100000_0.100000_300000",\
#"RW_0.965000_4_0.100000_0.100000_300000",\
#"RW_0.981000_4_1.000000_0.100000_300000",\
#"RW_0.985000_4_1.000000_0.100000_300000",\
#"RW_0.986000_4_2.000000_0.100000_300000",\
#"RW_0.987000_4_2.000000_0.100000_300000",\
#"random", "random2"
]

"""
location                                    min         min     max         max
COM_4_1.0_3000_0.9_b                   & 4.40e-02 & 4.68e-06 & 9.16e-01 & 2.18e-02\\
COM_4_1.0_3000_0.5_b                   & 3.19e-02 & -3.10e-07 & 9.51e-01 & 9.11e-03\\
COM_4_1.0_3000_0.1_b                   & 6.96e-01 & -5.73e-07 & 8.78e-01 & 6.15e-04\\
RW_0.951000_4_0.100000_0.100000_300000 & 7.09e+00 & 7.01e-05 & 6.79e-01 & 2.17e-12\\
RW_0.955000_4_0.100000_0.100000_300000 & 5.52e+03 & 1.35e-04 & 6.51e-01 & 1.94e-12\\
RW_0.965000_4_0.100000_0.100000_300000 & 1.18e+03 & 2.00e-04 & 5.94e-01 & 2.68e-12\\
RW_0.981000_4_1.000000_0.100000_300000 & 2.51e+04 & 1.12e-03 & 4.65e-01 & 4.04e-13\\
RW_0.985000_4_1.000000_0.100000_300000 & 1.93e+04 & 3.14e-03 & 3.15e-01 & 1.29e-12\\
RW_0.986000_4_2.000000_0.100000_300000 & 1.31e+04 & 1.81e-03 & 3.28e-01 & 3.70e-13\\
RW_0.987000_4_2.000000_0.100000_300000 & 6.09e+04 & 3.58e-03 & 2.89e-01 & 2.78e-13\\
random                                 & 0.00e+00 & 0.00e+00 & 0.00e+00 & 0.00e+00\\
random2                                & 1.29e+03 & 2.35e-01 & -1.41e-01 & 0.00e+00\\

text
COM_4_1.0_3000_0.9_b                   & 1.35e+00 & 1.96e-01 & 4.82e-01 & 0.00e+00\\
COM_4_1.0_3000_0.5_b                   & 1.37e+00 & 1.74e-01 & 5.20e-01 & 0.00e+00\\
COM_4_1.0_3000_0.1_b                   & 1.37e+00 & 1.42e-01 & 6.69e-01 & 0.00e+00\\
RW_0.951000_4_0.100000_0.100000_300000 & 1.28e+00 & 2.70e-01 & 4.13e-01 & 0.00e+00\\
RW_0.955000_4_0.100000_0.100000_300000 & 1.41e+00 & 2.69e-01 & 4.08e-01 & 0.00e+00\\
RW_0.965000_4_0.100000_0.100000_300000 & 1.44e+00 & 2.31e-01 & 4.04e-01 & 0.00e+00\\
RW_0.981000_4_1.000000_0.100000_300000 & 1.54e+00 & 3.92e-01 & 3.25e-01 & 0.00e+00\\
RW_0.985000_4_1.000000_0.100000_300000 & 1.52e+00 & 3.22e-01 & 2.60e-01 & 0.00e+00\\
RW_0.986000_4_2.000000_0.100000_300000 & 1.48e+00 & 4.06e-01 & 2.47e-01 & 0.00e+00\\
RW_0.987000_4_2.000000_0.100000_300000 & 1.48e+00 & 4.69e-01 & 2.20e-01 & 0.00e+00\\
random                                 & 0.00e+00 & 0.00e+00 & 0.00e+00 & 0.00e+00\\
random2                                & 3.90e+00 & 8.69e-01 & -2.05e-02 & 0.00e+00\\

"""

#fff = ["blobs","noisy_circles","noisy_moons","random"]
#eval_comp(col, fff)

#t_files = ["COM_4_1.0_200_0.5_b", "COM_4_1.0_200_0.1_b", "COM_4_1.0_300_0.1_b", "RW_0.960000_4_0.100000_0.100000_40000"]

#eval_comp(ger_small, ger_small_files, [6,7,8,9])
#eval_comp(ger, ger_files, [6,7,8,9])
#eval_comp(usa_small_1M, usa_small_files, [6,7,8,9])
sys.exit()

def random_clustering(nodes, clusters, noise=0.1):
    import random
    for nid in xrange(nodes):
        if random.random() > noise:
            print nid, 1+nid%(clusters-1)

# 57.554 ger_small

# 1.254.298 ger

#1.086.678 usa_small_1M

#random_clustering(57554, 150)
#random_clustering(1254298, 150)
#random_clustering(1086678, 150)

sys.exit()

#col_ger_small = "a56c48c3df72452177dce28efd790ddc"
#col_ger = "90bda4cc98289c9d3c231127f8253189"
#clustering_str = "/home/seydanator/Desktop/SVN/code/usa"#LocF_10000_8.txt"

#lf = LocF(10000, 8)
#jac = Jaccard(0.5, 8, "w")
#graph = RWalk(0.5, 8, 1, 1)
#com = Comb(lf, jac)

#print lf
#print jac
#print com

#args = ['./dbscan', 'scan', "", ""]
#args.extend(str(com).split())
#print args
    
#sys.exit()
#print lf
#print jac
#print com
#usa_jac = "USA_"+collection_usa+"_jac"
#usa_lf =  "USA_"+collection_usa+"_lf"
#usa_com = "USA_"+collection_usa+"_com"

#gsg = "GER_SMALL"
#gs = "GER"
#borders = borders_by_collection(ger)
#print borders
#plot(ger_small, parse_clustering(gsg + "/LF_2000_6"), borders=borders, output=gsg+"/LF_2000_6.png", cluster_count=16)

from itertools import product

def cluster_RWalk(col, points, text_w, jump, epses, dists, text_t=100, ss=15000):
    borders = borders_by_collection(col)
    for pts, w, c, eps, dist in product(points, text_w, jump, epses, dists):
        #cl_st = dir_col[col]+"/"+"RW_"+str(eps)+"_"+str(pts)+"_"+str(w)+"_"+str(c)+"_"+str(dist)
        cl_st = "{}/RW_{:f}_{:d}_{:f}_{:f}_{:d}".format(dir_col[col], eps, pts, w, c, dist)
        if not os.path.isfile(cl_st):
            print cl_st
            stat = cluster(col, RWalk(eps, pts, w, c, dist, text_t), cl_st)
            
            if 3*stat["unclustered"] <= stat["clustered"]:
                clusters = parse_clustering(cl_st)
                
                for count in [25]:
                    print "plotting", count, ss
                    plot(col, clusters, show=False, output=cl_st+"_"+str(count), cluster_count=count, borders=borders, legend=False, shape_size=ss)
        
def cluster_LocF(col, points, epses, ss=15000):
    borders = borders_by_collection(col)
    for eps, pts in product(epses, points):
        cl_st = dir_col[col]+"/"+"LF_"+str(eps)+"_"+str(pts)
        if not os.path.isfile(cl_st):
            stat = cluster(col, LocF(eps, pts), cl_st)
            
            if 3*stat["unclustered"] <= stat["clustered"]:
                clusters = parse_clustering(cl_st)
                for count in [25]:
                    plot(col, clusters, show=False, output=cl_st+"_"+str(count), cluster_count=count, borders=borders, shape_size=ss)

def cluster_Jaccard(col, points, epses, kind, ss=15000):
    borders = borders_by_collection(col)
    for eps, pts in product(epses, points):
        #cl_st = dir_col[col]+"/"+"JW_"+str(eps)+"_"+str(pts)+"_"+kind
        cl_st = "{}/JW_{:f}_{:d}_{}".format(dir_col[col], eps, pts, kind)
        if not os.path.isfile(cl_st):
            stat = cluster(col, Jaccard(eps, pts, kind), cl_st)
            
            if 3*stat["unclustered"] <= stat["clustered"]:
                clusters = parse_clustering(cl_st)
                for count in [25]:
                    plot(col, clusters, show=False, output=cl_st+"_"+str(count), cluster_count=count, borders=borders, shape_size=ss)

def cluster_Combined(col, points, epses, loc_epses, jac_epses, kind, ss=15000):
    borders = borders_by_collection(col)
    for pts, eps, loc_eps, jac_eps in product(points, epses, loc_epses, jac_epses):
        cl_st = dir_col[col]+"/"+"COM_"+str(pts)+"_"+str(eps)+"_"+str(loc_eps)+"_"+str(jac_eps)+"_"+kind
        
        if not os.path.isfile(cl_st):
            l = LocF(loc_eps, pts)
            j = Jaccard(jac_eps, pts, kind)
            stat = cluster(col, Comb(eps, pts, l, j), cl_st)
            
            if 3*stat["unclustered"] <= stat["clustered"]:
                clusters = parse_clustering(cl_st)
                for count in [25]:
                    plot(col, clusters, show=False, output=cl_st+"_"+str(count), cluster_count=count, borders=borders, shape_size=ss)


#print borders_by_collection(ger)
#[47.17342, 55.05963, 5.89202, 15.02338]

#print borders_by_collection(usa_small)
#[25.3716, 50.06222, -124.83338, -52.62152]

#print borders_by_collection(usa)
#[23.5, 52.99994, -129.72311, -51.17668]

#shape_size usa = 30000
#ger = 15000

#cluster_Combined(usa_small_1M, [4], [1.0], arange(2000, 100000 , 2000), [0.5], "b")
#cluster_Combined(ger_small, [4], [1.0],[100,200, 300, 500, 1000], [0.5], "b")


#cluster_LocF(ger_small, [8], arange(100, 501, 100))
#sys.exit(0)
#plot(usa_small_1M, parse_clustering("USA_SMALL_1M/LF_15000_8"), show=False, cluster_count=30, borders=borders_by_collection(usa_small_1M), shape_size=30000, legend=True)




#cluster_RWalk(usa_small_1M, [4], [0.1], [0.1], arange(0.95, .97, 0.001), [300000], 5, ss=50000)
#cluster_RWalk(usa_small_1M, [4], [1.0], [0.1], arange(0.98, .99, 0.001), [300000], 5, ss=50000)
#cluster_RWalk(usa_small_1M, [4], [2.0], [0.1], arange(0.98, .99, 0.001), [300000], 5, ss=50000)

#cluster_Combined(usa_small_1M, [4], [1.0],[2000, 3000, 5000, 7000], [0.1, 0.5, 0.9], "b", ss=50000)


#cluster_RWalk(ger, [4], [0.1], [0.1], arange(0.95, .97, 0.001), [35000], 5)
#cluster_RWalk(ger, [4], [1.0], [0.1], arange(0.97, .99, 0.001), [35000], 5)
#cluster_RWalk(ger, [4], [2.0], [0.1], arange(0.98, 1.0, 0.001), [35000], 5)


#cluster_Combined(ger, [4], [1.0],[100,200, 300, 500, 1000], [0.1, 0.5, 0.9], "b")

#cluster_RWalk(ger_small, [4], [0.1], [0.1], arange(0.96, .98, 0.001), [40000], 5)

#cluster_RWalk(ger_small, [4], [1.0], [0.1], arange(0.969, .97, 0.0001), [40000], 5)
#cluster_RWalk(ger_small, [4], [2.0], [0.1], arange(0.97, .98, 0.001), [40000], 5)

#cluster_Combined(ger_small, [4], [1.0],[100,200, 300, 500, 1000], [0.1, 0.5, 0.9], "b")


#try other scopes
#cluster_RWalk(usa_small_1M, [4], [2.0, 3.0], [0.01], arange(0.998, .999, 0.0001), [300000], 5, ss=50000)
#cluster_Combined(usa_small_1M, [4], [0.5, 0.7, 1.0, 1.2, 1.5, 2.0],[100000], [0.1], "b", ss=50000)

#cluster_Jaccard(usa_small_1M, [8], [0.4, 0.5, 0,6 ,0.7], "b")#arange(0.9, 1.01, 0.01)
sys.exit()
#cluster_histogram(["/home/seydanator/Desktop/SVN/code/USA_SMALL_1M/LF_60000_8","/home/seydanator/Desktop/SVN/code/USA_SMALL_1M/LF_10000_8","/home/seydanator/Desktop/SVN/code/USA_SMALL_1M/LF_40000_8","/home/seydanator/Desktop/SVN/code/USA_SMALL_1M/LF_20000_8"])

#cluster_RWalk(ger_small_d, [4], [1.0], [0.1], arange(0.965, 0.97, 0.0005))

#triang_combiner(ger)
#



#from matplotlib import cm
#maps=[m for m in cm.datad if not m.endswith("_r")]
#maps.sort()
#l=len(maps)+1
#for i, m in enumerate(maps):
#    edge_pic(usa_small_1M, cm=m)

#edge_pic(ger_small)
#edge_pic(ger)

sys.exit()


#from datetime import datetime
#for pair in [(10, 90),(100, 0), (90, 10), (900, 100), (50, 50), (500, 500), (75, 25), (25, 75)]:
#    startTime = datetime.now()
#    ev = eval_clustering(ger_small, "/home/seydanator/Desktop/SVN/code/GER_SMALL/LF_2000_4", pair[0], pair[1])
#    print pair, datetime.now()-startTime
#    print ev
    #print ev['Geo']
    #print ev['JBi']
 



#locf

#jac

def add_collection_points(arr, name):
    def add_id(gen):
        _id = 0
        for en in gen:
            en["_id"] = _id
            yield en
            _id += 1
    
    tweet_list = []
    tweet_counter = 0
    write_buffer = 1000

    for tweet in add_id(arr):
        tweet_counter += 1
        tweet_list.append(tweet)

        if tweet_counter == write_buffer:
            tweets[name].insert(tweet_list)
            tweet_list = []
            tweet_counter = 0

    if tweet_counter > 0:
        tweets[name].insert(tweet_list)




#new_collection(ger, "ger_test")# sampling=(2,1))#count=30000)#
#new_collection(ger, "ger_test90k", sampling=(1,12), count=False)
#frequency_map("ger_test100k", borders_by_collection("ger_test100k"))
#cluster_RWalk(ger_small, [4], [1.1], [0.1], arange(0.95, 0.99, 0.01))
cluster_RWalk("ger_test90k", [4], [1.0], [0.1], arange(0.95, 0.99, 0.01))

sys.exit()

#cluster(collection_ger_small, graph, ger_small_graph)
#plot(collection_ger_small, ger_small_graph, show=False, output=ger_small_graph, cluster_count=30)


#cluster(collection_usa, jac, usa_jac)
#plot(collection_usa, usa_jac, show=False, output=usa_jac, cluster_count=8)

#cluster(collection_usa, lf, usa_lf)
#plot(collection_usa, usa_lf, show=False, output=usa_lf, cluster_count=8)

#cluster(collection_usa, com, usa_com)
#plot(collection_usa, usa_com, show=False, output=usa_com, cluster_count=8)

#evaluation_lf  = eval_clustering(collection_usa, usa_lf)
#print evaluation_lf

#evaluation_jac = eval_clustering(collection_usa, usa_jac)
#print evaluation_jac

#evaluation_com = eval_clustering(collection_usa, usa_com)
#print evaluation_com

#print evaluation

sys.exit()

from sklearn import datasets
from sklearn.metrics import euclidean_distances
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler

# Generate datasets. We choose the size big enough to see the scalability
# of the algorithms, but not too big to avoid too long running times
n_samples = 1500

noisy_circles = datasets.make_circles(n_samples=n_samples, factor=.5, noise=.05)
noisy_moons = datasets.make_moons(n_samples=n_samples, noise=.05)
blobs = datasets.make_blobs(n_samples=n_samples, random_state=8)

np.random.seed(0)
no_structure = np.random.rand(n_samples, 2), None

#sys.exit()

import pymongo
connection = pymongo.MongoClient()
tweets = connection['tweets']

def add_points(arr, name):
    def add_id(gen):
        _id = 0
        for pair in gen:
            elem = {}
            elem["_id"] = _id
            x,y = pair
            elem["loc"] = [x,y]
            elem["tag"] = []
            
            yield elem
            _id += 1
    
    tweet_list = []
    tweet_counter = 0
    write_buffer = 1000

    for tweet in add_id(arr):
        tweet_counter += 1
        tweet_list.append(tweet)

        if tweet_counter == write_buffer:
            tweets[name].insert(tweet_list)
            tweet_list = []
            tweet_counter = 0

    if tweet_counter > 0:
        tweets[name].insert(tweet_list)


add_points(StandardScaler().fit_transform(noisy_circles[0]), "noisy_circles")
add_points(StandardScaler().fit_transform(noisy_moons[0]), "noisy_moons")
add_points(StandardScaler().fit_transform(blobs[0]), "blobs")
add_points(StandardScaler().fit_transform(no_structure[0]), "random")

with open("random", "w") as ran:
    for i in xrange(n_samples):
        ran.write("%s %s\n" % (i,(i%2)))

for name in ["noisy_circles", "noisy_moons", "blobs"]:
    cluster(name, LocF(5000.0, 2), name)
for name in ["noisy_circles", "noisy_moons", "blobs", "random"]:
    ev = eval_clustering(name, name, 1, 0)
    print ev
    plot(name, parse_clustering(name), show=False, scatter=False, colored=True, output="test/"+name+"_color_2000",shape_size=2000)

sys.exit()
