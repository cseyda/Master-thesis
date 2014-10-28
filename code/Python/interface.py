#!/usr/bin/python
# -*- coding: utf-8 -*-
import pymongo
import bson
import codecs
import file_input

import sys
from collections import defaultdict, Counter
from itertools import islice, izip, combinations, cycle

from subprocess import Popen, PIPE

from math import log

import matplotlib as mpl
#mpl.use('Cairo')

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from matplotlib.collections import PolyCollection, PatchCollection

from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import ScalarMappable

#import matplotlib.font_manager
#print matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')


#import matplotlib.font_manager as font_manager

#path = '/usr/share/fonts/TTF/LiberationSerif-Regular.ttf'
#prop = font_manager.FontProperties(fname=path)
#print prop.get_name()
#mpl.rcParams['font.family'] = prop.get_name()

mpl.rcParams['font.family'] = "Liberation Serif"

from shapely.ops import cascaded_union
from shapely.geometry import Point
from numpy import asarray

import itertools

import numpy as np
import networkx as nx

# geo: {loc: [lon, lat]}
connection = pymongo.MongoClient()
tweets = connection['tweets']
#stats  = connection['stats']

# Database: 1
# Collection: files
# Documents: tweets

def borders_by_collection(collection):
    "returns [min_lat, max_lat, min_lon, max_lon]"
    max_lat = tweets[collection].find_one(sort=[("loc.1", -1)])["loc"][1]
    min_lat = tweets[collection].find_one(sort=[("loc.1", 1)])["loc"][1]
    
    max_lon = tweets[collection].find_one(sort=[("loc.0", -1)])["loc"][0]
    min_lon = tweets[collection].find_one(sort=[("loc.0", 1)])["loc"][0]
    
    return [min_lat, max_lat, min_lon, max_lon]

def clusters_in_borders(collection, borders, clustering_file, minPts=30):
    
    min_lat, max_lat, min_lon, max_lon = borders
    out_file = clustering_file + "{:f}.{:f}.{:f}.{:f}".format(*borders)
    
    import os
    if not os.path.isfile(out_file):
        print "not"
        clusters = parse_clustering(clustering_file)
        
        cid = 1
        for cluster in clusters:
            if len(cluster) < 5:
                continue
            
            lats = []
            lons = []
            lon_lat = defaultdict(set)
            tags = Counter()
            
            for entry in tweets[collection].find({"_id":{"$in":cluster}}):
                lon, lat = entry["loc"]
                tag = entry["tag"]
    
                lon_lat[lon].add(lat)
                
                tags.update(tag)
            
            for lon, lat_set in lon_lat.iteritems():
                for lat in lat_set:
                    lats.append(lat)
                    lons.append(lon)
            
            #print min_lat, max(lats), min(lats), max_lat
            
            if not any((min_lat <= lat <= max_lat) for lat in lats):
                #print "break 1"
                continue
            
            if not any((min_lon <= lon <= max_lon) for lon in lons):
                #print "break 2"
                continue
            
            with codecs.open(out_file, 'a') as new_clustering:
                for nid in cluster:
                    new_clustering.write(u"{:d} {:d}\n".format(nid, cid))
            cid+=1
    else:
        print "in"
    return parse_clustering(out_file)
    
def plot(collection, clusters, show=False, output="", cluster_count=8, borders=[], scatter=False, colored=True, shape_size=15000, legend=False, keywords=[]):
    if borders:
        min_lat, max_lat, min_lon, max_lon = borders
        
    else:
        max_lat = []
        min_lat = []
        
        max_lon = []
        min_lon = []
        
        for cluster in clusters:
            for entry in tweets[collection].find({"_id":{"$in":cluster}}):
                lon, lat = entry["loc"]
                
                max_lat.append(lat)
                max_lon.append(lon)
                
                min_lat.append(lat)
                min_lon.append(lon)
        
        max_lat = max(max_lat)
        max_lon = max(max_lon)
        
        min_lat = min(min_lat)
        min_lon = min(min_lon)
        
    lat_border = (max_lat - min_lat) / 20
    lon_border = (max_lon - min_lon) / 20
    
    m = Basemap(projection='merc',\
        llcrnrlat=min_lat-lat_border, urcrnrlat=max_lat+lat_border,\
        llcrnrlon=min_lon-lon_border, urcrnrlon=max_lon+lon_border,\
        lat_ts=20,resolution='l')
    
    #m.drawcoastlines()
    #m.drawcountries()
    m.drawmapboundary()
    
    min_lon -= lon_border
    min_lat -= lat_border
    max_lon += lon_border
    max_lat += lat_border
    min_x, min_y = m(min_lon, min_lat)
    max_x, max_y = m(max_lon, max_lat)
    
    if colored:
        colors = cycle(['#FF0000','#00FF00','#0000FF','#FF00FF','#00FFFF','#FFFF00','#70DB93','#B5A642','#5F9F9F','#B87333','#2F4F2F','#9932CD','#871F78','#855E42','#545454','#8E2323','#F5CCB0','#238E23','#CD7F32','#DBDB70','#C0C0C0','#527F76','#9F9F5F','#8E236B','#2F2F4F','#EBC79E','#CFB53B','#FF7F00','#DB70DB','#D9D9F3','#5959AB','#8C1717','#238E68','#6B4226','#8E6B23','#007FFF','#00FF7F','#236B8E','#38B0DE','#DB9370','#ADEAEA','#5C4033','#4F2F4F','#CC3299','#99CC32','#FFFFFF','#000000'])
    else:
        colors =   cycle(['k'])
    #markers =  cycle(['o', 'D', 'h', 'H', '*', 's', 'v', 'p'])
    #patterns = cycle([''])#,'--', '++', 'xx', '\\\\', '**', 'oo', 'OO', '..'])
    #patterns = cycle(['-', '+', 'x', '\\', '*', 'o', 'O', '.'])
    
    color_name = {'#FF0000':'Red','#00FF00':'Green','#0000FF':'Blue','#FF00FF':'Magenta','#00FFFF':'Cyan','#FFFF00':'Yellow','#70DB93':'Aquamarine','#B5A642':'Brass','#5F9F9F':'Cadet Blue','#B87333':'Copper','#2F4F2F':'Dark Green','#9932CD':'Dark Orchid','#871F78':'Dark Purple','#855E42':'Dark Wood','#545454':'Dim Grey','#8E2323':'Firebrick','#F5CCB0':'Flesh','#238E23':'Forest Green','#CD7F32':'Gold','#DBDB70':'Goldenrod','#C0C0C0':'Grey','#527F76':'Green Copper','#9F9F5F':'Khaki','#8E236B':'Maroon','#2F2F4F':'Midnight Blue','#EBC79E':'New Tan','#CFB53B':'Old Gold','#FF7F00':'Orange','#DB70DB':'Orchid','#D9D9F3':'Quartz','#5959AB':'Rich Blue','#8C1717':'Scarlet','#238E68':'Sea Green','#6B4226':'Semi-Sweet Chocolate','#8E6B23':'Sienna','#007FFF':'Slate Blue','#00FF7F':'Spring Green','#236B8E':'Steel Blue','#38B0DE':'Summer Sky','#DB9370':'Tan','#ADEAEA':'Turquoise','#5C4033':'Very Dark Brown','#4F2F4F':'Violet','#CC3299':'Violet Red','#99CC32':'Yellow Green','#FFFFFF':'White','#000000':'Black'}
    
    fname = ''
    if not show:
        if output != '':
            fname = output
        else:
            fname = 'noname' + '_' + collection
    
    handles=[]
    legends=[]
    
    fig = plt.gcf()
    ax  = plt.gca()
    
    stats = tweets["stats"].find_one({"_id":collection})
    word_hist = dict(izip(stats["words"], stats["words_histogram"]))
    
    node_count = tweets[collection].count()
    
    x_ = []
    y_ = []
    z_ = []
    c_ = ['w']
    c_count = 1
    for cluster in clusters[:cluster_count]:
        lats = []
        lons = []
        lon_lat = defaultdict(set)
        tags = Counter()
        
        for entry in tweets[collection].find({"_id":{"$in":cluster}}):
            lon, lat = entry["loc"]
            tag = entry["tag"]

            lon_lat[lon].add(lat)
            
            tags.update(tag)
        
        for lon, lat_set in lon_lat.iteritems():
            for lat in lat_set:
                lats.append(lat)
                lons.append(lon)
        
        #print min_lat, max(lats), min(lats), max_lat
        
        if not any((min_lat <= lat <= max_lat) for lat in lats):
            continue
        
        if not any((min_lon <= lon <= max_lon) for lon in lons):
            continue
        
        if keywords:
            for key, tag in itertools.product(keywords, tags.keys()):
                if key in tag:
                    break
            else:
                continue
        
        x, y = m(lons,lats)
        
        x_.extend(x)
        y_.extend(y)
        z_.extend(len(x) * [c_count])
        c_count+=1
        

        c = next(colors)
        c_.append(c)
        #mark = next(markers)
        #pat = next(patterns)
        
        
        if show:
            print color_name[c], len(cluster), len(cluster)/float(node_count)
            print tags.most_common(20)
        else:
            with codecs.open(fname+'txt', 'a') as not_show:
                #not_show.write(color_name[c] + u"\n")
                #for key, val in tags.most_common(20):
                #    not_show.write(key.encode('utf8'))
                #    not_show.write(u" ")
                #    not_show.write(str(val))
                #    not_show.write(u"\n")
                #not_show.write(u"\n")
                not_show.write(u"{} {:d} {:.2f}\n".format(color_name[c], len(cluster), len(cluster)/float(node_count)))
                
                kv = u" {:d} {:.2f}\n"
                for key, val in tags.most_common(20):
                    not_show.write(key.encode('utf8'))
                    not_show.write(kv.format(val, val/float(word_hist[key])))
                not_show.write(u"\n")
        
        if scatter:
            plt.scatter(x,y, c=c, zorder=10)
        else:
            shapes = []
            for lon,lat in izip(x,y):
                shapes.append(Point(lon, lat).buffer(shape_size, resolution=10))
            merged = cascaded_union(shapes)
            
            patches = []
            
            try:
                for part in merged:
                    a = asarray(part.exterior)
                    poly = Polygon(a)#, facecolor=c)
                    patches.append(poly)
            except:
                a = asarray(merged.exterior)
                poly = Polygon(a)#, facecolor=c)
                patches.append(poly)
            
            ax.add_collection(PatchCollection(patches, facecolor=c, zorder=4))
            
            #ProxyArtist
            handles.append(Rectangle((0, 0), 1, 1, fc=c))
            legends.append(color_name[c])
    
    #http://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
    #from mpl_toolkits.axes_grid1 import make_axes_locatable
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    
    if legend:
        plt.legend(handles, legends, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)#, cax=cax)
    
    """
    x_w = np.linspace(min_x, max_x,5)
    y_w = np.linspace(min_y, max_y,5)
    x_w, y_w = np.meshgrid(x_w, y_w)
    
    x_ = np.r_[x_,[item for sublist in x_w for item in sublist]]
    y_ = np.r_[y_,[item for sublist in y_w for item in sublist]]
    z_ = np.r_[z_,[0]*25]

    xi = np.linspace(min_x, max_x,len(x_))
    yi = np.linspace(min_y, max_y,len(y_))
    xi, yi = np.meshgrid(xi, yi)
    
    cc = mpl.colors.ColorConverter()
    for i, col in enumerate(c_[1:], 1):
        print i
        zi = mpl.mlab.griddata(x_,y_,np.array([1.0 if elem==i else 0.0 for elem in z_]),xi,yi)
        cmap = mpl.colors.LinearSegmentedColormap.from_list("c",[cc.to_rgba(col, alpha=0),cc.to_rgba(col, alpha=0.2), cc.to_rgb(col)])
        
        plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
        plt.pcolormesh(xi,yi,zi,cmap=cmap)#,norm=norm,)
    """

    """    
>>> att = set()
>>> for d in m.natural_info:
...     att.add(d["type"])
>>> att
>>> set(['water', 'park', 'riverbank', 'forest'])

roads
set(['proposed', 'primary', 'pedestrian', 'elevator', 'bridleway', 'secondary_link', 'tertiary', 'primary_link', 'escalator', 'service', 'residential', 'motorway_link', 'cycleway', 'platform', 'corridor', 'secondary', 'living_street', 'track', 'access_ramp', 'motorway', 'construction', 'tertiary_link', 'trunk', 'path', 'layby', 'trunk_link', 'none', 'rest_area', 'abandoned', 'footway', 'wheelchair', 'unclassified', 'bus_stop', 'steps', 'crossing', 'footway; track', 'road'])
    
    """
    import shapefile
    from shapely.geometry import Polygon as Polygon2
    #from shapely.geometry import MultiPolygon
    import shapely
    #import math
    
    def bb(bbox):
        min_lonbb, min_latbb, max_lonbb, max_latbb = bbox
        if any((min_lat <= lat <= max_lat) for lat in [min_latbb, max_latbb]):
            if any((min_lon <= lon <= max_lon) for lon in [min_lonbb, max_lonbb]):
                return True
        if any((min_latbb <= lat <= max_latbb) for lat in [min_lat, max_lat]):
            if any((min_lonbb <= lon <= max_lonbb) for lon in [min_lon, max_lon]):
                return True
        
        return False
    
    def proc_points(shape_pts):
        lons,lats = zip(*shape_pts)
        data = np.array(m(lons, lats)).T
        
        try:
            poly=Polygon2(data)
            data = list(poly.simplify(0.1).exterior.coords)
        except:
            pass
        
        #except:
        #    pass
        return Polygon(data, False)
    
    if True:
        #place = "ne_10m_ocean"
        #print place
        #sf = shapefile.Reader('shapefiles/'+place+'/'+place)
        #records = sf.records()
        #patches = []
        #for i, info in enumerate(records):
        #    sh = sf.shape(i)
        #    if not bb(sh.bbox):
        #        #continue
        #        pass
        #    
        #    parts = sh.parts
        #    parts.append(len(sh.points))
        #    
        #    for j in xrange(len(parts)-1):
        #        patch = proc_points(sh.points[parts[j]:parts[j+1]])
        #        
        #        patch.set_facecolor('blue')
        #        patch.set_edgecolor('blue')
        #        patch.set_linewidth(10.4)
        #        patch.set_alpha(0.3)
        #    
        #        patches.append(patch)
        #if patches:
        #    ax.add_collection(PatchCollection(patches, match_original=True))
        #del sf
        #del records
        ocean = [[min_lon, min_lat], [min_lon, max_lat], [max_lon, max_lat], [max_lon, min_lat], [min_lon, min_lat]]
        ocean = proc_points(ocean)
        ocean.set_facecolor('#0066CC')
        ocean.set_zorder(-1)
        plt.gca().add_patch(ocean)
        
        place = "ne_10m_admin_1_states_provinces"
        print place
        sf = shapefile.Reader('shapefiles/'+place+'/'+place)
        records = sf.records()
        patches = []
        for i, info in enumerate(records):
            sh = sf.shape(i)
            if not bb(sh.bbox):
                continue
            
            parts = sh.parts
            parts.append(len(sh.points))
            
            for j in xrange(len(parts)-1):
                patch = proc_points(sh.points[parts[j]:parts[j+1]])
                
                patch.set_facecolor('green')
                patch.set_edgecolor('k')
                patch.set_linewidth(1.)
                patch.set_alpha(.7)
                
                patches.append(patch)
        if patches:
            ax.add_collection(PatchCollection(patches, match_original=True))
        del sf
        del records
        
        place = "ne_10m_urban_areas"
        print place
        sf = shapefile.Reader('shapefiles/'+place+'/'+place)
        records = sf.records()
        patches = []
        for i, info in enumerate(records):
            sh = sf.shape(i)
            if not bb(sh.bbox):
                continue
            
            parts = sh.parts
            parts.append(len(sh.points))
            
            for j in xrange(len(parts)-1):
                patch = proc_points(sh.points[parts[j]:parts[j+1]])
                
                patch.set_facecolor('brown')
                patch.set_edgecolor('k')
                patch.set_linewidth(.8)
                patch.set_alpha(.5)
            
                patches.append(patch)
        if patches:
            ax.add_collection(PatchCollection(patches, match_original=True))
        del sf
        del records
        """
        places = "ne_10m_roads", "ne_10m_roads_north_america"
        for place in places:
            print place
            sf = shapefile.Reader('shapefiles/'+place+'/'+place)
            records = sf.records()
            patches = []
            for i, info in enumerate(records):
                sh = sf.shape(i)
                if not bb(sh.bbox):
                    continue
                
                parts = sh.parts
                parts.append(len(sh.points))
                
                for j in xrange(len(parts)-1):
                    patch = proc_points(sh.points[parts[j]:parts[j+1]])
                    
                    patch.set_facecolor('none')
                    patch.set_edgecolor('grey')
                    patch.set_linewidth(0.6)
                    patch.set_alpha(0.8)
                
                    patches.append(patch)
            if patches:
                ax.add_collection(PatchCollection(patches, match_original=True))
            del sf
            del records
        """
        
        """
        places = "ne_10m_rivers_lake_centerlines", "ne_10m_rivers_north_america", "ne_10m_rivers_europe"
        for place in places:
            print place
            sf = shapefile.Reader('shapefiles/'+place+'/'+place)
            records = sf.records()
            patches = []
            for i, info in enumerate(records):
                sh = sf.shape(i)
                if not bb(sh.bbox):
                    continue
                
                parts = sh.parts
                parts.append(len(sh.points))
                
                for j in xrange(len(parts)-1):
                    patch = proc_points(sh.points[parts[j]:parts[j+1]])
                    
                    #patch.set_facecolor('blue')
                    #patch.set_edgecolor('none')
                    #patch.set_linewidth(0.4)
                    #patch.set_alpha(0.7)
                
                    patches.append(patch)
            if patches:
                patches = PatchCollection(patches, facecolors='none', edgecolors='blue')
                patches.set_alpha(.7)
                ax.add_collection(patches)
            del sf
            del records
        """
        
        """
        places = "ne_10m_lakes", "ne_10m_lakes_north_america", "ne_10m_lakes_europe"
        for place in places:
            print place
            sf = shapefile.Reader('shapefiles/'+place+'/'+place)
            records = sf.records()
            patches = []
            for i, info in enumerate(records):
                sh = sf.shape(i)
                if not bb(sh.bbox):
                    continue
                
                parts = sh.parts
                parts.append(len(sh.points))
                
                for j in xrange(len(parts)-1):
                    patch = proc_points(sh.points[parts[j]:parts[j+1]])
                    #patch.set_facecolor('none')
                    #patch.set_edgecolor('blue')
                    #patch.set_linewidth(0.6)
                    #patch.set_alpha(0.7)
                    patches.append(patch)
            if patches:
                patches = PatchCollection(patches, facecolors='blue', edgecolors='none', linewidths=.6)
                patches.set_alpha(.7)
                ax.add_collection(patches)
            del sf
            del records
        """
        
        """
        place = "ne_10m_populated_places"#[['osm_id'], ['name'], ['type'], ['population']
        print place
        import matplotlib.patheffects as PathEffects

        sf = shapefile.Reader('shapefiles/'+place+'/'+place)
        records = sf.records()
        patches = []
        for i, info in enumerate(records):
            sh = sf.shape(i)
            if not bb([sh.points[0][0], sh.points[0][1], sh.points[0][0], sh.points[0][1]]):
                continue
            
            #if type(info[3]) == int:
            if info[0] <= 5:#scalerank
                x,y = m(*sf.shape(i).points[0])
                ax.text(x,y,info[4].decode("utf8"), fontsize=10, weight='bold', path_effects=[PathEffects.withStroke(linewidth=.5, foreground="w")], zorder=100, alpha=1.)
        del sf
        del records
        """
        
        
        
    else:
        places =  []#"brandenburg",#"berlin",
        
        for place in places:
            f = "natural"#[4317997, 'Schlachtensee', 'water']
            sf = shapefile.Reader('shapefiles/'+place+'/'+f)
            records = sf.records()
            
            print place, f
            
            patches = []
            for i, info in enumerate(records):
                sh = sf.shape(i)
                if not bb(sh.bbox):
                    continue
                
                lines = proc_points(sh)
                
                if info[2] in ['water']:
                    lines.set_facecolor('blue')
                    #lines.set_edgecolor('k')
                    #lines.set_linewidth(0.3)
                    lines.set_alpha(0.4)
                    
                elif info[2] in ['park']:
                    #lines = proc_points(sh)
                    lines.set_facecolor('green')
                    #lines.set_edgecolor('k')
                    #lines.set_linewidth(0.3)
                    lines.set_alpha(0.2)
                    
                elif info[2] in ['forest']:
                    #lines = proc_points(sh)
                    lines.set_facecolor('green')
                    #lines.set_edgecolor('k')
                    #lines.set_linewidth(0.3)
                    lines.set_alpha(0.5)
                    
                else:
                    continue
                patches.append(lines)
            ax.add_collection(PatchCollection(patches, match_original=True))
            
            del sf
            del records
            
            f = "waterways"#[id, name, type, width]
            sf = shapefile.Reader('shapefiles/'+place+'/'+f)
            records = sf.records()
            
            print place, f
            
            patches = []
            for i, info in enumerate(records):
                sh = sf.shape(i)
                if not bb(sh.bbox):
                    continue
                
                lines = proc_points(sh)
                lines.set_edgecolor('blue')
                lines.set_facecolor('none')
                
                if info[2] in ['stream', 'canal']:#,'river', ]:
                    lines.set_linewidth(0.6)
                    lines.set_alpha(0.6)
                elif info[2] in ['drain', 'ditch']:
                    lines.set_linewidth(0.4)
                    lines.set_alpha(0.4)
                else:
                    continue
                patches.append(lines)
            ax.add_collection(PatchCollection(patches, match_original=True))
            
            del sf
            del records
            
            f = "roads"#[3996955, '', 'A 115', 'motorway', 1, 0, 0, 120]
            
            print place, f
            sf = shapefile.Reader('shapefiles/'+place+'/'+f)
            records = sf.records()
            
            patches = []
            for i, info in enumerate(records):
                sh = sf.shape(i)
                if not bb(sh.bbox):
                    continue
                lines = proc_points(sh)
                lines.set_facecolor('none')
                
                if info[3] in ['motorway', 'motorway_link']:
                    #lines.set_facecolor((0,0,0,0))
                    lines.set_edgecolor('k')
                    lines.set_linewidth(1.2)
                    lines.set_alpha(0.7)
                elif info[3] in ['primary','primary_link']:
                    #lines.set_facecolor((0,0,0,0))
                    lines.set_edgecolor('k')
                    lines.set_linewidth(1.)
                    lines.set_alpha(0.6)
                elif info[3] in ['secondary','secondary_link']:
                    lines.set_edgecolor('k')
                    lines.set_linewidth(.8)
                    lines.set_alpha(0.5)
                elif info[3] in ['residential','track']:
                    lines.set_edgecolor('gray')
                    lines.set_linewidth(.6)
                    lines.set_alpha(0.5)
                else:
                    continue
                patches.append(lines)
            ax.add_collection(PatchCollection(patches, match_original=True))
            del sf
            del records
            
            #f = "railways"#[['osm_id', 'N', 11, 0], ['name', 'C', 48, 0], ['type', 'C', 16, 0]]
            
            #print place, f
            #sf = shapefile.Reader('shapefiles/'+place+'/'+f)
            #records = sf.records()
            #patches=[]
            #for i, info in enumerate(records):
            #    sh = sf.shape(i)
            #    if not bb(sh.bbox):
            #        continue
            #    lines = proc_points(sh)
            #    
            #    if info[2] in ['rail']:
            #        lines.set_facecolor('none')
            #        lines.set_edgecolor('red')
            #        lines.set_linewidth(.8)
            #        lines.set_alpha(0.7)
            #    else:
            #        continue
            #    patches.append(lines)
            #ax.add_collection(PatchCollection(patches, match_original=True))
            #del sf
            #del records
            
            #f = "buildings"#[['osm_id', 'N', 11, 0], ['name', 'C', 48, 0], ['type', 'C', 16, 0]]
            
            #print place, f
            #sf = shapefile.Reader('shapefiles/'+place+'/'+f)
            #records = sf.records()
            
            #patches=[]
            #for i, info in enumerate(records):
            #    sh = sf.shape(i)
            #    if not bb(sh.bbox):
            #        continue
            #    lines = proc_points(sh)
            #    
            #    lines.set_facecolor('blue')
            #    lines.set_edgecolor('yellow')
            #    lines.set_linewidth(.3)
            #    
            #    patches.append(lines)
            #ax.add_collection(PatchCollection(patches, match_original=True))
            #del sf
            #del records
            
            
            
            ### TODO all in left down?
            
            f = "places"#[['osm_id'], ['name'], ['type'], ['population']
            
            print place, f
            
            sf = shapefile.Reader('shapefiles/'+place+'/'+f)
            records = sf.records()
            
            import matplotlib.patheffects as PathEffects
            
            for i, info in enumerate(records):
                sh = sf.shape(i)
                if not bb([sh.points[0][0], sh.points[0][1], sh.points[0][0], sh.points[0][1]]):
                    continue
                
                if type(info[3]) == int:
                    if info[3] > 20000:
                        x,y = m(*sf.shape(i).points[0])
                        ax.text(x,y,info[1].decode("utf8"), fontsize=8, path_effects=[PathEffects.withStroke(linewidth=1, foreground="w")], zorder=100, alpha=0.9)
            del sf
            del records
            
            print "rdy"
    
        
    if show:
        plt.show()
    else:
        ftype = 'pdf'
        plt.savefig(fname + "." + ftype, format=ftype, bbox_inches='tight', pad_inches=0, dpi=600)
    plt.close()

def parse_clustering(input_str):
    #input_str = "/home/seydanator/Downloads/graphlabapi/release/toolkits/graph_analytics/"
    #files = 4
    
    #name = "out"
    
    clustering = defaultdict(list)
    
    with open(input_str, "r") as source:
        for line in source:
            vid, cid = line.split(" ")#vertex, component ID
            clustering[int(cid)].append(int(vid))
            #print vid, cid
            
    cluster_list = []
    for k,v in clustering.iteritems():
        cluster_list.append((k,len(v)))
    cluster_list.sort(key = lambda x: x[1], reverse=True)
    
    return [clustering[cid[0]] for cid in cluster_list]

def cluster_histogram(c_list):
    #plt.hist([len(cluster) for cluster in c1], 100)
    #plt.hist([len(cluster) for cluster in c2], 100)
    #plt.show()
    #plt.close()
    
    #plt.plot([len(cluster) for cluster in c1])
    #plt.plot([len(cluster) for cluster in c2])
    #plt.show()
    #plt.close()
    clusterings = []
    for name in c_list:
        clusterings.append([parse_clustering(name), name])
    
    for clustering in clusterings:
        cl = [len(cluster) for cluster in clustering[0]]
        print cl, len(cl)
        bins=[pow(2,x) for x in xrange(20) if pow(2,x-1)<max(cl)]
        #print bins
        plt.hist(cl, bins=bins, cumulative=True)
        plt.hist(cl, bins=bins)
        plt.gca().set_xscale('log',basex=2)
        plt.gca().xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%d'))
        plt.show()
    
    l=[]
    legends=[]
    for clustering in clusterings:
        c=0
        p=[]
        for cluster in reversed(clustering[0]):
            c+=len(cluster)
            p.append(c)
        plt.plot(p, label=clustering[1])
    
    plt.gca().legend()
    plt.show()
    plt.close()
    
def edge_pic(collection, cm="Reds"):
    import pymongo
    connection = pymongo.MongoClient()
    graphs = connection['graphs']
    
    g_stat = graphs["graph"].find({"_id":collection})[0]
    
    orig2triang = g_stat["orig2triang"]
    triang2orig = []
    for i in range(g_stat["combined_points"]):
        triang2orig.append([])
    
    for i, tr in enumerate(orig2triang):
        triang2orig[tr].append(i)
    
    e_gen = ((triang2orig[entry["edge"][0]][0], triang2orig[entry["edge"][1]][0], entry["dist"]) for entry in graphs[collection].find())
    
    triang(collection, borders_by_collection(collection), legend=False, edges=e_gen, name="tri_"+collection+"_"+cm+".pdf", cm=cm)











    """
    stats = tweets["stats"].find_one({"_id":collection})
    word_hist = dict(izip(stats["words"], stats["words_histogram"]))
    
    node_count = tweets[collection].count()
    
    for cluster in clusters:
        lats = []
        lons = []
        tags = Counter()
        
        for entry in (e for e in tweets[collection].find({"_id":{"$in":cluster}}) for cluster in clusters):
            lon, lat = entry["loc"]
            tag = entry["tag"]

            
            tags.update(tag)
        
        
            lats.append(lat)
            lons.append(lon)
        
        if keywords:
            for key, tag in itertools.product(keywords, tags.keys()):
                if key in tag:
                    break
            else:
                continue
        
        x, y = m(lons,lats)
        
        c_count+=1
        
            with codecs.open(fname+'txt', 'a') as not_show:
                #not_show.write(color_name[c] + u"\n")
                #for key, val in tags.most_common(20):
                #    not_show.write(key.encode('utf8'))
                #    not_show.write(u" ")
                #    not_show.write(str(val))
                #    not_show.write(u"\n")
                #not_show.write(u"\n")
                not_show.write(u"{} {:d} {:.2f}\n".format(color_name[c], len(cluster), len(cluster)/float(node_count)))
                
                kv = u" {:d} {:.2f}\n"
                for key, val in tags.most_common(20):
                    not_show.write(key.encode('utf8'))
                    not_show.write(kv.format(val, val/float(word_hist[key])))
                not_show.write(u"\n")
       """
        
def frequency_map(collection, borders, legend=False, clustering=None, keywords=[], out="", bin_scale=1):
    min_lat, max_lat, min_lon, max_lon = borders
    
    lat_border = (max_lat - min_lat) / 20
    lon_border = (max_lon - min_lon) / 20
    
    m = Basemap(projection='merc',\
        llcrnrlat=min_lat-lat_border, urcrnrlat=max_lat+lat_border,\
        llcrnrlon=min_lon-lon_border, urcrnrlon=max_lon+lon_border,\
        lat_ts=20,resolution='l')
    
    m.drawcoastlines()
    m.drawcountries()
    
    lons = []
    lats = []
    
    file_out = "freq_"+out+collection
    
    if clustering:
        stats = tweets["stats"].find_one({"_id":collection})
        word_hist = dict(izip(stats["words"], stats["words_histogram"]))
        
        #node_count = tweets[collection].count()
        
        D_tags = Counter()
        D_tags_in_d = Counter()
        
        d_tags_list = []
        
        for cluster in clustering:
            l_lons = []
            l_lats = []
            l_tags = Counter()
            
            for entry in tweets[collection].find({"_id":{"$in":cluster}}):
                lon, lat = entry["loc"]
                tag = entry["tag"]
                
                l_tags.update(tag)
                
                l_lats.append(lat)
                l_lons.append(lon)
            
            #print keywords
            #p_key = 0.01
            if keywords:
                #c = 0
                #c_all = float(len(cluster))
                for key, tag in itertools.product(keywords, l_tags.keys()):
                    if key in tag:
                        #c+=l_tags[tag]
                        #if c/c_all > p_key:
                        break
                else:
                    continue
            
            D_tags.update(l_tags)
            D_tags_in_d.update(l_tags.keys())
            
            d_tags_list.append(l_tags)
            
            lons.extend(l_lons)
            lats.extend(l_lats)
        
        """
        # averaged tf idf
        tfidf = []
        
        def tf(term, document):
            return document[term]
        
        def idf(term, frequency):
            return log(len(d_tags_list) / float(frequency))
        
        for term, frequency in D_tags.iteritems():
            idf_ = idf(term, frequency)
            s = 0
            for document in d_tags_list:
                s += tf(term, document) * idf_
            
            s /= float(len(d_tags_list))
            tfidf.append([term, s])
        
        tfidf.sort(key=lambda a: a[1], reverse=True)
        #print tfidf
        """
        #for key, val in D_tags.iteritems():
        #    print key, val, D_tags_in_d[val], float(len(d_tags_list)), val*D_tags_in_d[key] / float(len(d_tags_list))
        
        tfidf = [(key, val*D_tags_in_d[key]/float(len(d_tags_list))*val/float(word_hist[key])) for key, val in D_tags.iteritems()]
        tfidf.sort(key=lambda a: a[1], reverse=True)
        
        with codecs.open(file_out+'txt', 'a') as not_show:
            not_show.write(str(len(clustering)) + " " + str(len(d_tags_list)) + "\n")
            kv = u" {:f}\n"
            for key, val in tfidf:#tags.most_common(20):
                not_show.write(key.encode('utf8'))
                not_show.write(kv.format(val))
            not_show.write(u"\n")
        
        #print D_tags.items()
        #print sorted(D_tags.items(), key=lambda a: a[1], reverse=True)
        
        with codecs.open(file_out+'txtfreq', 'a') as not_show:
            for key, val in sorted(D_tags.items(), key=lambda a: a[1], reverse=True):
                not_show.write(key.encode('utf8'))
                not_show.write(" {:d} {:.2f}\n".format(val, float(word_hist[key])/val))
            not_show.write(u"\n")
        
    else:
        for entry in tweets[collection].find():
            lon, lat = entry["loc"]
            
            lats.append(lat)
            lons.append(lon)
    
    if len(lons) == 0:
        print "nothing found for ", keywords
        return
    
    x, y = m(lons,lats)
    
    x_max, y_max = m(max_lon, max_lat)
    x_min, y_min = m(min_lon, min_lat)
    
    
    lat_bin = 15.0
    lat_dif = y_max-y_min#max_lat-min_lat
    lon_dif = x_max-x_min#max_lon-min_lon
    
    lon_bin =  lon_dif / lat_dif * lat_bin
    
    bin_size_lat = lat_dif / lat_bin
    bin_size_lon = lon_dif / lon_bin
    
    lon_bin = max(1, lon_bin * (max(x)-min(x)) / lon_dif)
    lat_bin = max(1, lat_bin * (max(y)-min(y)) / lat_dif)
    
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=[int(lon_bin), int(lat_bin)])#, normed=True)
    
    extent = [xedges.min(), xedges.max(), yedges.min(), yedges.max()]
    if (extent[1]-extent[0]) < bin_size_lon:
        extent[1] = extent[0] + bin_size_lon
    if (extent[3]-extent[2]) < bin_size_lat:
        extent[3] = extent[2] + bin_size_lat
    
    cm = plt.cm.get_cmap("winter_r")#Reds")
    cm.set_under("w", 0.0)
    
    bounds = [heatmap.max()]#-min_d]
    steps = 8
    divider = pow(heatmap.max() / 1.0, 1.0/steps)
    
    for i in xrange(steps-1):
        bounds.append(bounds[i]/divider)#float(3.2))
    bounds.append(1.0)
    bounds = list(reversed(bounds))
    
    norm = mpl.colors.BoundaryNorm(bounds, cm.N)
    
    
    im = plt.gca().imshow(heatmap.T, cmap=cm, norm=norm, interpolation='none', extent=extent, origin="lower", zorder=4, alpha=0.9)#,  #http://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
    
    #im = plt.pcolormesh(xedges, yedges, heatmap.T, cmap=cm, norm=norm, zorder=4)
    #d_y = max(lats)-min(lats)
    #print d_y
    #im = plt.hexbin(x,y,cmap=cm, norm=norm, zorder=4, gridsize=int(d_y))
    
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    if legend:
        cbar = plt.gcf().colorbar(im, cax=cax, ticks=bounds, boundaries=bounds, format='%.1e')
        cbar.set_label('Number of points')

    plt.savefig(file_out, format="pdf", bbox_inches='tight', pad_inches=0)
    plt.close()

def triang(collection, borders, legend=False, edges=False, name=False, cm="Reds"):
    print "triang"
    from pyproj import Geod
    g = Geod(ellps='WGS84')
    
    min_lat, max_lat, min_lon, max_lon = borders
    
    lat_border = (max_lat - min_lat) / 20
    lon_border = (max_lon - min_lon) / 20
    
    m = Basemap(projection='merc',\
        llcrnrlat=min_lat-lat_border, urcrnrlat=max_lat+lat_border,\
        llcrnrlon=min_lon-lon_border, urcrnrlon=max_lon+lon_border,\
        lat_ts=20,resolution='l')
    m.drawcoastlines()
    m.drawcountries()
    
    
    
    #scal = ScalarMappable(cmap=plt.get_cmap('gray'),norm=Normalize(vmax=100000))
    scal = 0
    
    dists = []
    edge_points = []
    mst = []
    
    if not edges:
        lons = []
        lats = []
        pts = defaultdict(lambda : defaultdict(list))
        
        for entry in tweets[collection].find():
            lon, lat = entry["loc"]
            pts[lon][lat].append(entry["_id"])
    
        for lon, lat_pts in pts.iteritems():
            for lat in lat_pts.iterkeys():
                lats.append(lat)
                lons.append(lon)
    
        x, y = m(lons,lats)
        
        print "no edges, generating"
        import matplotlib.tri as tri
    
        triang = tri.Triangulation(x, y)
        dists.append('b')
        for id1, id2 in triang.edges:
            az12,az21,dist = g.inv(lons[id1], lats[id1],  lons[id2], lats[id2])
            if dist < 40000:
                edge_points.append([[x[id1],y[id1]], [x[id2],y[id2]]])
        scal = ScalarMappable(cmap=plt.get_cmap('Reds_r'),norm=Normalize(vmax=100000, vmin=0))
    
    if edges:
        print "edges"
        
        for edge in edges:
            #print list(tweets[collection].find({"_id":{"$in":edge[0:2]}}))
            id1, id2 = list(tweets[collection].find({"_id":{"$in":edge[0:2]}}))
            x1, y1 = m(*id1["loc"])
            x2, y2 = m(*id2["loc"])
            
            if edge[2] == 0:
                mst.append([[x1,y1], [x2,y2]])
            else:
                edge_points.append([[x1,y1], [x2,y2]])
                dists.append(edge[2])
        
        print min(dists), np.median(dists), max(dists)
        
        #scal = ScalarMappable(cmap=plt.get_cmap('Reds_r'),norm=LogNorm(vmax=max(colors), vmin=np.median(colors)))
    
    fig = plt.gcf()
    ax  = plt.gca()
    
    #http://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    #col_map = plt.cm.Accent#is ok
    #col_map = plt.cm.Set3#is not really ok
    #col_map = plt.cm.Reds_r#
    col_map = plt.cm.get_cmap(cm)
    
    
    
    #bounds = np.linspace(0,max(dists),16)
    #min_d = 13000#min(dists)
    bounds = [max(dists)]#-min_d]
    for i in xrange(11):
        bounds.append(bounds[i]/float(1.6))
    
    #print bounds
    bounds = list(reversed(bounds))
    #bounds = [i+min_d for i in reversed(bounds)]
    #print bounds
    
    #from numpy import histogram
    #print histogram(dists, bins=bounds)
    
    norm = mpl.colors.BoundaryNorm(bounds, col_map.N)
    scal = ScalarMappable(cmap=col_map,norm=norm)
    
    lc = LineCollection(edge_points, colors=scal.to_rgba(dists))
    
    ax.add_collection(lc)
    ax.add_collection(LineCollection(mst, colors=scal.to_rgba([0])))
    
    
    scal._A = []
    cbar = fig.colorbar(scal, cax=cax, ticks=bounds, boundaries=bounds, format='%.1e')#'%1i')
    # spacing='proportional'
    cbar.set_label('Length of edges in km')
    plt.tight_layout()
    
    if not name:
        plt.show()
    else:
        plt.savefig(name, bbox_inches='tight', pad_inches=0)
    plt.close()
    
def eval_clustering(collection, clustering_file, sample_make, sample_skip, count = False):
    i_count = 0
    if count:
        i_count = 1
    
    args = ['./dbscan', 'eval', collection, clustering_file, str(sample_make), str(sample_skip), str(i_count)]
    print args
    process = Popen(args, stdout=PIPE)
    stdout, stderr = process.communicate()
    
    score = defaultdict(dict)
    for line in stdout.splitlines():
        #print line
        line = line.split()
        if line:
            print line
            if not count:
                if line[0] == "Eval":
                    metric, evaluation, rating = line[1:]
                    score[metric][evaluation] = rating
            else:
                if line[0] == "Counter:":
                    score["count"]["count"] = line[1]
    
    return score


class LocF:
    def __init__(self, threshold, minPts):
        self.threshold = threshold
        self.minPts = minPts
        self.metric = 1
    def __str__(self):
        return str(self.metric) + " " + str(self.threshold) + " " + str(self.minPts)

class Jaccard:
    def __init__(self, threshold, minPts, kind):
        self.threshold = threshold
        self.minPts = minPts
        self.kind = kind
        self.metric = 2
    def __str__(self):
        return str(self.metric) + " " + str(self.threshold) + " " + str(self.minPts) + " " + self.kind

class Comb:
    def __init__(self, threshold, minPts, metric1, metric2):
        self.metric1 = metric1
        self.metric2 = metric2
        self.metric = 4
        self.minPts = minPts
        self.threshold = threshold
    
    def __str__(self):
        return str(self.metric) + " " + str(self.threshold) + " " + str(self.minPts) + " " + str(self.metric1) + " " + str(self.metric2)

class RWalk:
    def __init__(self, threshold, minPts, text_weight, c, dist, text_t):
        self.threshold = threshold
        self.minPts = minPts
        self.metric = 3
        self.text_weight = text_weight
        self.c = c
        self.dist = dist
        self. text_t = text_t
        
    def __str__(self):
        return str(self.metric) + " " + str(self.threshold) + " " + str(self.minPts) + " " + str(self.text_weight) + " " + str(self.c) + " " + str(self.dist) + " " + str(self.text_t)

from datetime import datetime

def cluster(collection, metric, output_str):
    #1 LocF: threshold, minPts
    #2 J:    threshold, minPts, type (bigram, word)
        
    #3 Random Walk:
    #4 Comb:
    args = ['./dbscan', 'scan', collection, output_str]
    args.extend(str(metric).split())
    print args
    stats = dict()
    startTime = datetime.now()
    process = Popen(args, stdout=PIPE)
    
    #for line in process.stdout.readlines():
    with open(output_str+".out", "w") as output:
        for line in iter(process.stdout.readline, b''):
        
            if "%" not in line:
                output.write(line)
        
            #process.stdout.flush()
            print line
            line = line.split()
            
            if line:
                if line[0] == "clusters":
                    stats["clusters"] = int(line[1])
                
                elif line[0] == "unclustered":
                    stats["unclustered"] = int(line[1])
                
                elif line[0] == "clustered":
                    stats["clustered"] = int(line[1])
        
    
        delta = datetime.now()-startTime
        if (delta.seconds > 120):
            print delta.seconds/float(60), "minutes"
            output.write(str(delta.seconds/float(60))+" minutes\n")
        else:
            print delta.seconds, "seconds"
            output.write(str(delta.seconds)+" seconds\n")
    return stats


def insert_file(file_path):
    gen_list = []
    hex_str,gen = file_input.open_file(file_path)
    gen_list.append(gen)
    
    hex_str,gen2 = file_input.open_file(file_path)
    gen_list.append(gen2)
    
    hex_str,gen3 = file_input.open_file(file_path)
    gen_list.append(gen3)
    
    add_collection_to_mongo(hex_str, gen_list)
        
def add_collection_to_mongo(collection_name, gen_list):
    def add_id(gen, bi_map, wo_map):
        for _id, elem in enumerate(gen):
            elem["_id"] = _id
            
            wo_set = set()
            #wo_set = list()
            bi_set = set()
            #bi_set = list()
                    
            for tag in elem["tag"]:
                wo_set.add(wo_map[tag])
                #wo_set.append(wo_map[tag])
                
                tokens = u" " + tag + u" "
                n_tokens = len(tokens)

                for i in xrange(n_tokens-1):
                    bi_set.add(bi_map[tokens[i:i+2]])
                    #bi_set.append(bi_map[tokens[i:i+2]])
                    
            elem["bigram"] = sorted(list(bi_set))
            elem["words"]  = sorted(list(wo_set))
            yield elem
    
    print "open"
    #hex_str,gen = file_input.open_file(file_path)
    hex_str = collection_name
    gen = gen_list[0]
    
    elements=0
    for i in gen:
        elements+=1
    print elements
    
    #hex_str,gen = file_input.open_file(file_path)
    gen = gen_list[1]
    
    print hex_str
    
    if hex_str not in tweets.collection_names():
        bigrams = Counter()
        words = Counter()
        
        for elem in gen:
            words.update(elem["tag"])
            
            for tag in elem["tag"]:
                tokens = u" " + tag + u" "
                n_tokens = len(tokens)
                
                bigrams.update((tokens[i:i+2] for i in xrange(n_tokens-1)))

        bi_keys, bi_values = zip(*(item for item in bigrams.most_common()))
        wo_keys, wo_values = zip(*(item for item in words.most_common()))
        
        bi_map = dict()
        
        for count, key in enumerate(bi_keys, 1):
            bi_map[key] = count
        
        wo_map = dict()
        
        for count, key in enumerate(wo_keys, 1):
            wo_map[key] = count
        
        #hex_str,gen = file_input.open_file(file_path)
        gen = gen_list[2]
    
        tweet_list = []
        tweet_counter = 0
        write_buffer = 100
        percent = (write_buffer * 100) / float(elements)
        print percent
        wrote = 0
        print_count = 5
        for tweet in add_id(gen, bi_map, wo_map):
            tweet_counter += 1
            tweet_list.append(tweet)
        
            if tweet_counter == write_buffer:
                tweets[hex_str].insert(tweet_list)
                tweet_list = []
                tweet_counter = 0
                wrote+=1
                if (wrote * percent) > print_count:
                    print wrote * percent, "%"
                    print_count+=5
        if tweet_counter > 0:
            tweets[hex_str].insert(tweet_list)
        
        tweets[hex_str].ensure_index([("loc", pymongo.GEO2D)])
        
        print "insert stats"
        tweets["stats"].insert({"_id": hex_str, "bigrams":bi_keys, "bigrams_histogram":bi_values, "bigrams_count":len(bi_keys), "words": wo_keys, "words_histogram": wo_values, "words_count":len(wo_keys)})
    else:
        print "%s already inserted" % hex_str
    return hex_str

def new_collection(base, target, bbox=[0,0,0,0], sampling=(0,0), limit=0, count=False, combine=False, displace=False):
    """
    sampling: parameter 0: take, parameter 1: skip
    limit:    maximum points created in new collection
    bbox:     points have to be in bbox [min_lat, max_lat, min_lon, max_lon]
    count:    return only the amount of tweets
    """
    import pymongo
    connection = pymongo.MongoClient()
    tweets = connection['tweets']
    from math import ceil

    b_count = tweets[base].count()
    if count > 0:
        if b_count <= count:
            print "new count is more than old collection has"
            return False
    
    if target in tweets.collection_names():
        print target,"is already in tweets"
        return False
    
    
    def combiner(source):
        #if combine:
            #lons = []
            #lats = []
        pts = defaultdict(lambda : defaultdict(list))
        tags = defaultdict(list)
        
        for entry in source:#tweets[collection].find():
            lon, lat = entry["loc"]
            pts[lon][lat].append(entry["_id"])
            tags[lon, lat].extend(entry["tag"])
            
        for lon, lat_pts in pts.iteritems():
            for lat in lat_pts.iterkeys():
                #lats.append(lat)
                #lons.append(lon)
                entry = {}
                entry["loc"] = [lon, lat]
                entry["tag"] = tags[lon, lat]
                entry["o_ids"] = pts[lon][lat]
                yield entry
            
            
        #else:
        #    for entry in tweets[base].find():
        #        yield entry
    
    def displacer(source):
        from random import randint
        
        coords = defaultdict(set)
        
        for entry in source:
            x,y = entry["loc"]
            
            case = randint(0,7) #Inclusive
            while True:
                if y not in coords[x]:
                    coords[x].add(y)
                    entry["loc"] = [x,y]
                    #unique = True
                    break
                    
                else:
                    if case == 0:
                        x += 0.00001
                    if case == 1:
                        y += 0.00001
                    if case == 2:
                        x -= 0.00001
                    if case == 3:
                        y -= 0.00001
                    if case == 4:
                        y += 0.00001
                        x += 0.00001
                    if case == 5:
                        y -= 0.00001
                        x -= 0.00001
                    if case == 6:
                        y += 0.00001
                        x -= 0.00001
                    if case == 7:
                        y -= 0.00001
                        x += 0.00001
            yield entry
    
    def bboxer(source):
        for entry in source:
            lon, lat = entry["loc"]#[lon,lat]
            if bbox[0] <= lat <= bbox[1]:
                if bbox[2] <= lon <= bbox[3]:
                    yield entry
    
    def sampler(gen):
        while True:
            for i in xrange(sampling[0]):
                yield gen.next()
            for i in xrange(sampling[1]):
                gen.next()
    
    def counter(gen, limit):
        for i in xrange(limit):
            yield gen.next()
    
    gen_list = []
    
    source1 = tweets[base].find()
    source2 = tweets[base].find()
    source3 = tweets[base].find()
    
    bbox_filter = True
    if bbox[0] == bbox[1]:
        if bbox[2] == bbox[3]:
            bbox_filter = False
        
    if bbox_filter:
        source1 = bboxer(source1)
        source2 = bboxer(source2)
        source3 = bboxer(source3)
    
    if sampling[1]>0:
        source1 = sampler(source1)
        source2 = sampler(source2)
        source3 = sampler(source3)
    
    #if combine:
     #   source1 = combiner(source1)
      #  source2 = combiner(source2)
       # source3 = combiner(source3)
    
    #else:
    if displace:
        source1 = displacer(source1)
        source2 = displacer(source2)
        source3 = displacer(source3)

    if limit>0:
        source1 = counter(source1, limit)
        source2 = counter(source2, limit)
        source3 = counter(source3, limit)
    
    gen_list.append(source1)
    gen_list.append(source2)
    gen_list.append(source3)
    
    if count:
        print "base_count:", b_count
        print "after filtering:", sum(1 for i in gen_list[0])
    
    else:
        #add_collection_points(arr, target)
        add_collection_to_mongo(target, gen_list)

def triang_combiner(collection):
    from pyproj import Geod
    g = Geod(ellps='WGS84')
    
    lons = []
    lats = []
    pts = defaultdict(lambda : defaultdict(list))
    
    #points on the same position
    for entry in tweets[collection].find():
        lon, lat = entry["loc"]
        pts[lon][lat].append(int(entry["_id"]))
    
    triang2orig = list()
    
    for lon, lat_pts in pts.iteritems():
        for lat in lat_pts.iterkeys():
            lats.append(lat)
            lons.append(lon)
            
            triang2orig.append(pts[lon][lat])

    collection_count = tweets[collection].count()
    orig2triang = [None] * collection_count
    
    for i, elem_list in enumerate(triang2orig):
        for elem in elem_list:
            orig2triang[elem] = i
    
    points_count = len(lats)
    
    
    
    
    ##points = np.array([[0, 0], [0, 1.1], [1, 0], [1, 1]])
    points = np.array(zip(lons, lats))
    from scipy.spatial import Delaunay
    tri = Delaunay(points, qhull_options="QJ")
    tri.simplices
    
    sci_edges = []
    for i in xrange(tri.nsimplex):
        if i > tri.neighbors[i,2]:
            sci_edges.append((tri.vertices[i,0], tri.vertices[i,1]))
        if i > tri.neighbors[i,0]:
            sci_edges.append((tri.vertices[i,1], tri.vertices[i,2]))
        if i > tri.neighbors[i,1]:
            sci_edges.append((tri.vertices[i,2], tri.vertices[i,0]))



    #import matplotlib.tri as tri
    #triangulation = tri.Triangulation(lons, lats)
    
    G = nx.Graph()
    
    edges = dict()
    #for id1, id2 in triangulation.edges:
    for id1, id2 in sci_edges:
        if id1 < id2:
            id1, id2 = id2, id1
        az12,az21,dist = g.inv(lons[id1], lats[id1],  lons[id2], lats[id2])
        #edge_points.append([int(id1), int(id2), int(dist)])
        edges[int(id1), int(id2)] = dist
        G.add_edge(int(id1), int(id2),weight=dist)
    
    T=nx.minimum_spanning_tree(G)
    #print(sorted(T.edges(data=True)))
    
    for id1, id2 in T.edges():
        if id1 < id2:
            id1, id2 = id2, id1
        edges[id1, id2] = 0
    
    triang(collection, borders_by_collection(collection), edges=((triang2orig[edge[0]][0], triang2orig[edge[1]][0], dist) for edge, dist in edges.iteritems()), name="triang_"+str(collection)+".pdf")
    
    return
    
    graphs = connection['graphs']
    
    graphs["graph"].insert({"_id": collection,\
                            "combined_points" : points_count,\
                            #"orig2triang":orig2triang,\
                            })
    
    update_buffer = 70000
    update_list = []
    for i, nid in enumerate(orig2triang):
        update_list.append(nid)
        
        if i%update_buffer==0:
            graphs["graph"].update({"_id": collection},\
            {'$push': { "orig2triang": { "$each": update_list } }},True)
            update_list=[]
    if len(update_list)>0:
        graphs["graph"].update({"_id": collection},\
        {'$push': { "orig2triang": { "$each": update_list } }},True)
    
    e_list = []
    e_counter = 0
    buffer_limit = 1000
    
    for edge, dist in edges.iteritems():
        #print edge, dist
        e_counter += 1
        e_list.append({"edge": edge, "dist":int(dist)})
    
        if e_counter == buffer_limit:
            graphs[collection].insert(e_list)
            e_list = []
            e_counter = 0
            
    if e_counter > 0:
        graphs[collection].insert(e_list)
    
