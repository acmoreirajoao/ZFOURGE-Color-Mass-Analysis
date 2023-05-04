# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

"""

#libraries used:
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from itertools import count
from textwrap import wrap


# data lists

data = ['CDFS1.tsv', 'COSMOS1.tsv', 'UDS1.tsv', 
        'CDFS2.tsv', 'COSMOS2.tsv', 'UDS2.tsv',
        'CDFS3.tsv', 'COSMOS3.tsv', 'UDS3.tsv',
        'CDFS4.tsv', 'COSMOS4.tsv', 'UDS4.tsv']

datacdfs = ['CDFS1.tsv', 'CDFS2.tsv', 'CDFS3.tsv', 'CDFS4.tsv']
datacosmos = ['UDS1.tsv', 'UDS2.tsv', 'UDS3.tsv' , 'UDS4.tsv']
datauds = ['COSMOS1.tsv', 'COSMOS2.tsv', 'COSMOS3.tsv', 'COSMOS4.tsv']

data1 = ['CDFS1.tsv', 'COSMOS1.tsv', 'UDS1.tsv']
data2 = ['CDFS2.tsv', 'COSMOS2.tsv', 'UDS2.tsv']
data3 = ['CDFS3.tsv', 'COSMOS3.tsv', 'UDS3.tsv']
data4 = ['CDFS4.tsv', 'COSMOS4.tsv', 'UDS4.tsv']

#some house-keeping, to keep track of some relevant values

quantity_per_sample = {}
max_min_colors = {}

# =============================================================================
# Green valley as defined in Schawinski et al 2014
# =============================================================================

xval = np.linspace(0,15)

def y_max(x):
    return -0.24 + 0.25*x

def y_min(x):
    return -0.75 + 0.25*x


# =============================================================================
# Analysis 1: data in individual samples (all 12 seperate ones)
# =============================================================================

for i in data:
    sample = pd.read_csv(i,  sep='\t', on_bad_lines='skip', dtype = 'str')
    sample = sample.dropna(axis='rows')
    sample = sample.drop_duplicates(keep='first')
    sample = sample.apply(pd.to_numeric, errors='coerce')
    quantity_per_sample[i] = len(sample.index)
    #calculate color
    Flux_u = sample['Fu'] * 1e-6*3500
    Flux_r = sample['Fr'] *1e-6* 6750
    sample['color'] = -2.5 * np.log10(Flux_u/Flux_r)
    max_min_colors[max(sample['color'])] =  min(sample['color'])
    
    #plot the data
    color = sample['color'].values.tolist()
    mass = sample['lmass'].values.tolist()
    plt.scatter(mass, color, marker = '.', s = 3, color = '#322c2c')
    plt.ylim(-5,7)
    plt.xlim(0,14)
    plt.title(i)
    plt.suptitle("Color-Mass Plot for Sample")
    plt.xlabel(u"Stellar Mass log $M_*$($M_\u2609$)")
    plt.ylabel("u-r color")
    plt.plot(xval,y_min(xval), 'k--')
    plt.plot(xval,y_max(xval), 'k--')
    plt.savefig(i+ ".png",dpi=300)
    plt.show()
    
# =============================================================================
# Analysis 2: grouping by region
# =============================================================================

#list for CDFS data
lmasscdfs = []
colorcdfs = []

#list for COSMOS data
lmasscosmos = []
colorcosmos = []
  
#list for UDS data
lmassuds = []
coloruds = []


for q in datacdfs:
    dataframecdfs = pd.read_csv(q,  sep='\t', on_bad_lines='skip', dtype = 'str')
    dataframecdfs = dataframecdfs.dropna(axis='rows')
    dataframecdfs = dataframecdfs.drop_duplicates(keep='first')
    dataframecdfs = dataframecdfs.apply(pd.to_numeric, errors='coerce')
    Flux_ucdfs = dataframecdfs['Fu'] * 1e-6*3500
    Flux_rcdfs = dataframecdfs['Fr'] * 1e-6*6750
    dataframecdfs['color'] = -2.5 * np.log10(Flux_ucdfs/Flux_rcdfs)
    max_min_colors[max(dataframecdfs['color'])] =  min(dataframecdfs['color'])
    lmass_add = dataframecdfs['lmass'].values.tolist()
    color_add = dataframecdfs['color'].values.tolist()
    lmasscdfs.extend(lmass_add)
    colorcdfs.extend(color_add)    


for w in datacosmos:
    dataframecosmos = pd.read_csv(w,  sep='\t', on_bad_lines='skip', dtype = 'str')
    dataframecosmos = dataframecosmos.dropna(axis='rows')
    dataframecosmos = dataframecosmos.drop_duplicates(keep='first')
    dataframecosmos = dataframecosmos.apply(pd.to_numeric, errors='coerce')
    Flux_ucosmos = dataframecosmos['Fu'] * 1e-6*3500
    Flux_rcosmos = dataframecosmos['Fr'] * 1e-6*6750
    dataframecosmos['color'] = -2.5 * np.log10(Flux_ucosmos/Flux_rcosmos)
    max_min_colors[max(dataframecosmos['color'])] =  min(dataframecosmos['color'])
    lmass_add = dataframecosmos['lmass'].values.tolist()
    color_add = dataframecosmos['color'].values.tolist()
    lmasscosmos.extend(lmass_add)
    colorcosmos.extend(color_add)

for e in datauds:
    dataframeuds = pd.read_csv(e,  sep='\t', on_bad_lines='skip', dtype = 'str')
    dataframeuds = dataframeuds.dropna(axis='rows')
    dataframeuds = dataframeuds.drop_duplicates(keep='first')
    dataframeuds = dataframeuds.apply(pd.to_numeric, errors='coerce')
    Flux_uuds = dataframeuds['Fu'] * 1e-6*3500
    Flux_ruds = dataframeuds['Fr'] *1e-6*6750
    dataframeuds['color'] = -2.5 * np.log10(Flux_uuds/Flux_ruds)
    max_min_colors[max(dataframecosmos['color'])] =  min(dataframecosmos['color'])
    lmass_add = dataframeuds['lmass'].values.tolist()
    color_add = dataframeuds['color'].values.tolist()
    lmassuds.extend(lmass_add)
    coloruds.extend(color_add)
    
#plot data
plt.scatter(lmasscdfs, colorcdfs, marker = '.', s = 3, color = '#322c2c')
plt.title('Color-Mass Plot for CDFS data')
plt.xlabel(u"Stellar Mass log $M_*$($M_\u2609$)")
plt.ylabel("u-r color")
plt.plot(xval,y_min(xval), 'k--')
plt.plot(xval,y_max(xval), 'k--')
plt.ylim(-5,7)
plt.xlim(0,14)
plt.figure()

plt.scatter(lmasscosmos, colorcosmos, marker = '.', s = 3, color = '#322c2c')
plt.title('Color-Mass Plot for COSMOS data')
plt.xlabel(u"Stellar Mass log $M_*$($M_\u2609$)")
plt.ylabel("u-r color")
plt.plot(xval,y_min(xval), 'k--')
plt.plot(xval,y_max(xval), 'k--')
plt.ylim(-5,7)
plt.xlim(0,14)
plt.figure()

plt.scatter(lmassuds, coloruds, marker = '.', s = 3, color = '#322c2c')
plt.title('Color-Mass Plot for UDS data')
plt.xlabel(u"Stellar Mass log $M_*$($M_\u2609$)")
plt.ylabel("u-r color")
plt.plot(xval,y_min(xval), 'k--')
plt.plot(xval,y_max(xval), 'k--')
plt.ylim(-5,7)
plt.xlim(0,14)
plt.figure()

# =============================================================================
# Analysis 3: data in redshift categories    
# =============================================================================

#list for z<1
lmass1 = []
color1 = []
redshift1 = []

#list for 1<z<2
lmass2 = []
color2 = []
redshift2 = []
  
#list for 2<z<3
lmass3 = []
color3 = []
redshift3 = []
  
#list for z>3
lmass4 = []
color4 = []
redshift4 = []

#now to fill the lists

for a in data1:
    dataframe1 = pd.read_csv(a,  sep='\t', on_bad_lines='skip', dtype = 'str')
    dataframe1 = dataframe1.dropna(axis='rows')
    dataframe1 = dataframe1.drop_duplicates(keep='first')
    dataframe1 = dataframe1.apply(pd.to_numeric, errors='coerce')
    Flux_u1 = dataframe1['Fu'] * 1e-6*6750
    Flux_r1 = dataframe1['Fr'] * 1e-6*6750
    dataframe1['color'] = -2.5 * np.log10(Flux_u1/Flux_r1)
    max_min_colors[max(dataframe1['color'])] =  min(dataframe1['color'])
    lmass_add = dataframe1['lmass'].values.tolist()
    color_add = dataframe1['color'].values.tolist()
    redshift_add = dataframe1['zp'].values.tolist()
    lmass1.extend(lmass_add)
    color1.extend(color_add)
    redshift1.extend(redshift_add)

for b in data2:
    dataframe2 = pd.read_csv(b,  sep='\t', on_bad_lines='skip', dtype = 'str')
    dataframe2 = dataframe2.dropna(axis='rows')
    dataframe2 = dataframe2.drop_duplicates(keep='first')
    dataframe2 = dataframe2.apply(pd.to_numeric, errors='coerce')
    Flux_u2 = dataframe2['Fu'] * 1e-6*3500
    Flux_r2 = dataframe2['Fr'] * 1e-6*6750
    dataframe2['color'] = -2.5 * np.log10(Flux_u2/Flux_r2)
    max_min_colors[max(dataframe2['color'])] =  min(dataframe2['color'])     
    lmass_add = dataframe2['lmass'].values.tolist()
    color_add = dataframe2['color'].values.tolist()
    redshift_add = dataframe2['zp'].values.tolist()
    lmass2.extend(lmass_add)
    color2.extend(color_add)
    redshift2.extend(redshift_add)
    
for c in data3:
    dataframe3 = pd.read_csv(c,  sep='\t', on_bad_lines='skip', dtype = 'str')
    dataframe3 = dataframe3.dropna(axis='rows')
    dataframe3 = dataframe3.drop_duplicates(keep='first')
    dataframe3 = dataframe3.apply(pd.to_numeric, errors='coerce')
    Flux_u3 = dataframe3['Fu'] *1e-6* 3500
    Flux_r3 = dataframe3['Fr'] * 1e-6*6750
    dataframe3['color'] = -2.5 * np.log10(Flux_u3/Flux_r3)
    max_min_colors[max(dataframe3['color'])] =  min(dataframe3['color'])
    lmass_add = dataframe3['lmass'].values.tolist()
    color_add = dataframe3['color'].values.tolist()
    redshift_add = dataframe3['zp'].values.tolist()
    lmass3.extend(lmass_add)
    color3.extend(color_add)
    redshift3.extend(redshift_add)

for d in data4:
    dataframe4 = pd.read_csv(d,  sep='\t', on_bad_lines='skip', dtype = 'str')
    dataframe4 = dataframe4.dropna(axis='rows')
    dataframe4 = dataframe4.drop_duplicates(keep='first')
    dataframe4 = dataframe4.apply(pd.to_numeric, errors='coerce')
    Flux_u4 = dataframe4['Fu'] * 1e-6*3500
    Flux_r4 = dataframe4['Fr'] * 1e-6*6750
    dataframe4['color'] = -2.5 * np.log10(Flux_u4/Flux_r4)
    max_min_colors[max(dataframe4['color'])] =  min(dataframe4['color'])
    lmass_add = dataframe4['lmass'].values.tolist()
    color_add = dataframe4['color'].values.tolist()
    redshift_add = dataframe4['zp'].values.tolist()
    lmass4.extend(lmass_add)
    color4.extend(color_add)
    redshift4.extend(redshift_add)
    
    
        
#plot the data

iterator_graphtitle =(count(start = 1, step = 1))
#above iterator is just used for titles, to keep track of graphs

plt.scatter(lmass1, color1, marker = '.', s = 3, color = '#bb1279')
plt.suptitle("Color-Mass Plot for Redshift Bin")
plt.title(next(iterator_graphtitle))
plt.xlabel(u"Stellar Mass log $M_*$($M_\u2609$)")
plt.ylabel("u-r color")
plt.plot(xval,y_min(xval), 'k--')
plt.plot(xval,y_max(xval), 'k--')
plt.ylim(-5,7)
plt.xlim(0,14)
# plt.savefig("bin1.png",dpi=300)
plt.figure()

plt.scatter(lmass2, color2, marker = '.', s = 3, color = '#1221bb')
plt.suptitle("Color-Mass Plot for Redshift Bin")
plt.title(next(iterator_graphtitle))
plt.xlabel(u"Stellar Mass log $M_*$($M_\u2609$)")
plt.ylabel("u-r color")
plt.plot(xval,y_min(xval), 'k--')
plt.plot(xval,y_max(xval), 'k--')
plt.ylim(-5,7)
plt.xlim(0,14)
# plt.savefig("bin2.png",dpi=300)
plt.figure()

plt.scatter(lmass3, color3, marker = '.', s = 3, color = '#4eb313')
plt.suptitle("Color-Mass Plot for Redshift Bin")
plt.title(next(iterator_graphtitle))
plt.xlabel(u"Stellar Mass log $M_*$($M_\u2609$)")
plt.ylabel("u-r color")
plt.plot(xval,y_min(xval), 'k--')
plt.plot(xval,y_max(xval), 'k--')
plt.ylim(-5,7)
plt.xlim(0,14)
# plt.savefig("bin3.png",dpi=300)
plt.figure()

plt.scatter(lmass4, color4, marker = '.', s = 3, color = '#fb853a')
plt.suptitle("Color-Mass Plot for Redshift Bin")
plt.title(next(iterator_graphtitle))
plt.xlabel(u"Stellar Mass log $M_*$($M_\u2609$)")
plt.ylabel("u-r color")
plt.plot(xval,y_min(xval), 'k--')
plt.plot(xval,y_max(xval), 'k--')
plt.ylim(-5,7)
plt.xlim(0,14)
# plt.savefig("bin4.png",dpi=300)
plt.figure()


# =============================================================================
# Analysis 4: visualize all data categorized by z ranges
# =============================================================================

z1 = plt.scatter(lmass1, color1, marker = '1', s = 3, color = '#bb1279', alpha = 0.8)
z2 = plt.scatter(lmass2, color2, marker = '2', s = 3, color = '#1221bb', alpha = 0.8)
z3 = plt.scatter(lmass3, color3, marker = '3', s = 3, color = '#4eb313', alpha = 0.8)
z4 = plt.scatter(lmass4, color4, marker = '4', s = 3, color = '#fb853a', alpha = 0.8)
plt.title('Color-Mass Diagram for All Data')
plt.xlabel(u"Stellar Mass log $M_*$($M_\u2609$)")
plt.ylabel("u-r color")
plt.legend((z1,z2,z3,z4),
           ('z ≤ 1', '1 < z < 2', '2 < z <3', 'z ≥ 4'),
           scatterpoints=1,
           markerscale = 3,
           loc='upper left',
           ncol=1,
           fontsize=10)
# plt.savefig("All Data wout green valley.png",dpi=300)
plt.figure()


# =============================================================================
# Further analysis based on all data
# =============================================================================

#adding y axis limits, acknowledging possible outliers
colormass1 = plt.scatter(lmass1, color1, marker = '1', s = 3, color = '#bb1279', alpha = 0.8)
colormass2 = plt.scatter(lmass2, color2, marker = '2', s = 3, color = '#1221bb', alpha = 0.8)
colormass3 = plt.scatter(lmass3, color3, marker = '3', s = 3, color = '#4eb313', alpha = 0.8)
colormass4 = plt.scatter(lmass4, color4, marker = '4', s = 3, color = '#fb853a', alpha = 0.4)
plt.plot(xval,y_min(xval), 'k--')
plt.plot(xval,y_max(xval), 'k--')
plt.ylim(-1,5)
plt.title('Color-Mass Diagram for All Data')
plt.xlabel(u"Stellar Mass log $M_*$($M_\u2609$)")
plt.ylabel("u-r color")
plt.legend((colormass1, colormass2, colormass3, colormass4),
           ('z ≤ 1', '1 < z < 2', '2 < z <3', 'z ≥ 4'),
           scatterpoints=1,
           markerscale = 3,
           loc='upper left',
           ncol=1,
           fontsize=10)
#plt.savefig("CColor-Mass Diagram for All Data angstrom.png",dpi=300)
plt.figure()


#removing z bin 1 data, as it seems analogous compared to the others
plt.scatter(lmass2, color2, marker = '2', s = 3, color = '#1221bb', alpha = 0.8)
plt.scatter(lmass3, color3, marker = '3', s = 3, color = '#4eb313', alpha = 0.8)
plt.scatter(lmass4, color4, marker = '4', s = 3, color = '#fb853a', alpha = 0.8)
plt.plot(xval,y_min(xval), 'k--')
plt.plot(xval,y_max(xval), 'k--')
plt.ylim(-1,5)
plt.title('Color-Mass Diagram for All Data')
plt.xlabel(u"Stellar Mass log $M_*$($M_\u2609$)")
plt.ylabel("u-r color")
plt.legend((colormass2, colormass3, colormass4),
           ('1 < z < 2', '2 < z <3', 'z ≥ 4'),
           scatterpoints=1,
           markerscale = 3,
           loc='upper left',
           ncol=1,
           fontsize=10)
plt.figure()

#zooming into the previous graph
plt.scatter(lmass2, color2, marker = '2', s = 3, color = '#1221bb', alpha = 0.9)
plt.scatter(lmass3, color3, marker = '3', s = 3, color = '#4eb313', alpha = 0.8)
plt.scatter(lmass4, color4, marker = '4', s = 3, color = '#fb853a', alpha = 0.4)
plt.plot(xval,y_min(xval), 'k--')
plt.plot(xval,y_max(xval), 'k--')
plt.ylim(-1,3.5)
plt.xlim(6,14)
plt.title('Color-Mass Diagram for All Data')
plt.xlabel(u"Stellar Mass log $M_*$($M_\u2609$)")
plt.ylabel("u-r color")
plt.legend((colormass2, colormass3, colormass4),
           ('1 < z < 2', '2 < z <3', 'z ≥ 4'),
           scatterpoints=1,
           markerscale = 3,
           loc='upper left',
           ncol=1,
           fontsize=10)
#plt.savefig("zoomzoom.png",dpi=300)
plt.figure()

# =============================================================================
# 3d plot
# =============================================================================

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# ax.scatter(lmass1, color1, redshift1, marker = '1', s = 3, color = '#bb1279', alpha = 0.8)
# ax.scatter(lmass2, color2, redshift2, marker = '2', s = 3, color = '#1221bb', alpha = 0.9)
# ax.scatter(lmass3, color3, redshift3, marker = '3', s = 3, color = '#4eb313', alpha = 0.8)
# ax.scatter(lmass4, color4, redshift4, marker = '4', s = 3, color = '#fb853a', alpha = 0.4)
# ax.xlabel(u"Stellar Mass log $M_*$($M_\u2609$)")
# ax.ylabel("u-r color")
# ax.zlabel = ('z')
# plt.show()


# =============================================================================
# Calculating how many are in the green valley, above and below
# =============================================================================

#redshift bin1
red1 = []
green1 =[]
blue1 = []

for o,p in zip(lmass1,color1):
    if p >= y_max(o):
        red1.append(p)
    elif p <= y_min(o):
        blue1.append(p)
    else:
        green1.append(p)

print('bin 1:')
print('red',len(red1))
print('green',len(green1))
print('blue', len(blue1))
print('total', len(red1)+len(green1)+len(blue1))
print('-----------------------------------')


#redshift bin2
red2 = []
green2 =[]
blue2 = []

for k,l in zip(lmass2,color2):
    if l >= y_max(k):
        red2.append(l)
    elif l <= y_min(k):
        blue2.append(l)
    else:
        green2.append(l)
    
print('bin 2:')
print('red',len(red2))
print('green',len(green2))
print('blue', len(blue2))
print('total', len(red2)+len(green2)+len(blue2))
print('-----------------------------------')


#redshift bin3
red3 = []
green3 =[]
blue3 = []

for n,m in zip(lmass3,color3):
    if m >= y_max(n):
        red3.append(m)
    elif m <= y_min(n):
        blue3.append(m)
    else:
        green3.append(m)
    
print('bin 3:')
print('red',len(red3))
print('green',len(green3))
print('blue', len(blue3))
print('total', len(red3)+len(green3)+len(blue3))
print('-----------------------------------')


#redshift bin4
red4 = []
green4 =[]
blue4 = []

for g,h in zip(lmass4,color4):
    if h >= y_max(g):
        red4.append(h)
    elif h <= y_min(g):
        blue4.append(h)
    else:
        green4.append(h)
        
print('bin 4:')
print('red',len(red4))
print('green',len(green4))
print('blue', len(blue4))
print('total', len(red4)+len(green4)+len(blue4))

#graphing that data
zbins = ['z ≤ 1', '1 < z < 2', '2 < z <3', 'z ≥ 4']

prop_red1 = len(red1)/(len(red1)+len(green1)+len(blue1))
prop_blue1 = len(blue1)/(len(red1)+len(green1)+len(blue1))
prop_green1 = len(green1)/(len(red1)+len(green1)+len(blue1))

prop_red2 = len(red2)/(len(red2)+len(green2)+len(blue2))
prop_blue2 = len(blue2)/(len(red2)+len(green2)+len(blue2))
prop_green2 = len(green2)/(len(red2)+len(green2)+len(blue2))

prop_red3 = len(red3)/(len(red3)+len(green3)+len(blue3))
prop_blue3 = len(blue3)/(len(red3)+len(green3)+len(blue3))
prop_green3 = len(green3)/(len(red3)+len(green3)+len(blue3))

prop_red4 = len(red4)/(len(red4)+len(green4)+len(blue4))
prop_blue4 = len(blue4)/(len(red4)+len(green4)+len(blue4))
prop_green4 = len(green4)/(len(red4)+len(green4)+len(blue4))

prop_red = [prop_red1, prop_red2, prop_red3, prop_red4]
prop_blue = [prop_blue1, prop_blue2, prop_blue3, prop_blue4]
prop_green = [prop_green1, prop_green2, prop_green3, prop_green4]

p1 = plt.bar(zbins, prop_red, width = 0.7, color = '#d96161')
p2 = plt.bar(zbins, prop_green, width = 0.7,bottom = prop_red, color= '#37c278')
bottom_blue = [prop_red1+prop_green1, prop_red2+prop_green2, 
               prop_red3+prop_green3, prop_red4+prop_green4]
p3 = plt.bar(zbins, prop_blue, width = 0.7, bottom =bottom_blue, color = '#2f74c7')
plt.legend(['Red Sequence','Green Valley','Blue Cloud'])
plt.xlabel('Redshift Bins')
plt.ylabel('Percentage of sample')
plt.title("\n".join(wrap("Proportions of Red Sequence, Green Valley, and Blue Cloud Galaxies over Different z bins", 60)))
plt.savefig("percentages.png",dpi=300)
plt.show()


# =============================================================================
# Below: drafts, brainstorming mess
# =============================================================================

    
# CDFS1 = pd.read_csv('CDFS1.tsv',  sep='\t', on_bad_lines='skip', dtype = 'str')
# CDFS1.dropna(axis='rows')
# CDFS1 = CDFS1.drop_duplicates(keep='first')
# CDFS1 = CDFS1.apply(pd.to_numeric, errors='coerce')
# # print(CDFS1.dtypes)


# CDFS2 = pd.read_csv('CDFS2.tsv',  sep='\t', on_bad_lines='skip', dtype = 'str')
# CDFS2.dropna(axis='rows')
# CDFS2 = CDFS2.apply(pd.to_numeric, errors='coerce')
                  

# CDFS3 = pd.read_csv('CDFS3.tsv',  sep='\t', on_bad_lines='skip', dtype = 'str')
# CDFS3.dropna(axis='rows')
# CDFS3 = CDFS3.apply(pd.to_numeric, errors='coerce')


# CDFS4 = pd.read_csv('CDFS4.tsv',  sep='\t', on_bad_lines='skip', dtype = 'str')
# CDFS4 = CDFS4.apply(pd.to_numeric, errors='coerce')
# CDFS4.dropna(axis='rows')



# calculate U-R color
# CDFS = [CDFS1, CDFS2, CDFS3, CDFS4]
# for n in CDFS:
#     n['color'] = -2.5*np.log10(n['Fu']/n['Fr'])
#     # print(max(n['color']), min(n['color']))
    
#     #graph data
#     color = n['color'].values.tolist()
#     mass = n['lmass'].values.tolist()

#     plt.scatter(mass, color, marker = '.', s = 10)
#     # plt.ylim(1.4816,1.4818)
#     plt.show()





