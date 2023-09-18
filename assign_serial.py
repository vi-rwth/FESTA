import numpy as np
import os
import sys
import time 
import itertools as it
import shapely.geometry
import copy
import tqdm
import MDAnalysis as mda
import matplotlib.pyplot as plt
from matplotlib import ticker
import warnings
warnings.filterwarnings('ignore')
import operator as op
 
try:
  topo = sys.argv[4]
except IndexError:
  topo = None
T, thresh, traj = sys.argv[1:]

def index_search(values, l):
    indices = []
    for i, x in enumerate(l):
        if values == x:
            indices.append(i)
    return indices

def have_common_elem(l1, l2):
    for elem in l2:
        if op.countOf(l1, elem) > 0:
            return True
            break
    return False

def group_numbers(numbers, max_diff):
    separate_groups, subgroup = [], []
    tmplist = copy.deepcopy(numbers)
    seed_elem = tmplist[0]
    while any(tmplist):  
        min_distance = max_diff
        found = False
        tmplist.remove(seed_elem)
        for compare_elem in tmplist: 
            if ((seed_elem[0]-compare_elem[0])**2 + (seed_elem[1]-compare_elem[1])**2)**0.5< min_distance: 
                found = True
                min_distance = ((seed_elem[0]-compare_elem[0])**2 + (seed_elem[1]-compare_elem[1])**2)**0.5
                min_elem = compare_elem
        if found == True and any(subgroup):
            subgroup.append(seed_elem)
            seed_elem = min_elem                 
        else:
            if any(subgroup):  
                separate_groups.append(subgroup)
            subgroup = []
            sec_run = False
            min_distance = max_diff
            for group in separate_groups:
                for elem in group:
                    if ((seed_elem[0]-elem[0])**2 + (seed_elem[1]-elem[1])**2)**0.5< min_distance:
                        group.append(seed_elem)
                        sec_run = True
                        break
                if sec_run == True:
                    break
            if sec_run == False:
                subgroup.append(seed_elem)
            if any(tmplist):
                seed_elem = tmplist[0]
    return separate_groups

print('working on directory: ' + T)
outline, bins, coords = [], [], []
count1, count2 = it.count(0), it.count(0)
thresh_val, tot_min_frames, dimX, dimY = 0, 0, 0, 1
os.chdir(T)

data_fes = np.genfromtxt('fes.dat')
data_fes_T = np.transpose(data_fes)
a_fes, b_fes, ener = data_fes_T[0].copy(), data_fes_T[1].copy(), data_fes_T[2].copy()

if thresh == 'auto':
    thresh_val = abs(min(ener)) - (abs(min(ener))-abs(max(ener)))/3.5
    if min(ener) < 0:
        thresh_val = thresh_val * (-1)
    print('automatically determined', end =' ') 
else:
    thresh_val = int(thresh)
print('threshold value: ' + str(thresh_val))
    
b_count = b_fes[0]
for elem in b_fes:
    if not elem == b_count:
        dimY += 1
        b_count = elem    
while b_fes[next(count1)] == b_fes[0]:
    dimX += 1
    
low_max_a, high_max_a, low_max_b, high_max_b = a_fes[0], a_fes[dimX-1], b_fes[0], b_fes[-1]
outline_show_a, outline_show_b, edge = [], [], []

for i,elem in enumerate(ener):
    try:
        if elem<thresh_val and (a_fes[i] == low_max_a or a_fes[i] == high_max_a or b_fes[i] == low_max_b or b_fes[i] == high_max_b or ener[i-1]>thresh_val or ener[i+1]>thresh_val or ener[i-dimX]>thresh_val or ener[i+dimX]>thresh_val):
            if a_fes[i] == low_max_a or a_fes[i] == high_max_a or b_fes[i] == low_max_b or b_fes[i] == high_max_b:
                edge.append([a_fes[i],b_fes[i]])
            outline.append([a_fes[i],b_fes[i]])
            outline_show_a.append((a_fes[i]+abs(low_max_a))*(dimX/(abs(low_max_a)+abs(high_max_a))))
            outline_show_b.append(abs(b_fes[i]-abs(high_max_b))*(dimY/(abs(low_max_b)+abs(high_max_b))))
    except IndexError:
        pass

tolerance = abs(a_fes[0]-a_fes[1])/2    
data_colvar = np.genfromtxt('COLVAR')
data_colvar_T = np.transpose(data_colvar)
step, a, b = data_colvar_T[0].copy(), data_colvar_T[1].copy(), data_colvar_T[2].copy()

start1 = time.perf_counter()
all_points, sorted_coords = [], []
for i, elem in enumerate(a):
    all_points.append(shapely.geometry.Point(elem,b[i]))
grouped_points = group_numbers(outline, 20*np.sqrt(8)*tolerance)
#if edge:
#    print('periodicity detected: boundaries will be considered periodic')
#    pbc = []
#    grouped_edges = group_numbers(edge, 20*np.sqrt(8)*tolerance)
#    for i in range(len(grouped_edges)-1):
#        if (abs(grouped_edges[i][0][0])-abs(grouped_edges[i+1][0][0]))**2 + (abs(grouped_edges[i][0][0])-abs(grouped_edges[i+1][1][0]))**2 < 0.01:
#            pbc_group = []
#            for j, groups in enumerate(grouped_points):
#                if have_common_elem(grouped_edges[i],groups) or have_common_elem(grouped_edges[i+1],groups):
#                    pbc_group.append(j)
#                if len(pbc_group) == 2:
#                    pbc.append(pbc_group)
#                if len(pbc)*2 == len(grouped_edges):
#                    break
#    print(pbc)
#quit()
polygons = []
for groups in grouped_points:
    polygons.append(shapely.geometry.Polygon(groups))
    
print(str(len(polygons)) + ' minima identified')
for polygon in tqdm.tqdm(polygons, desc='determining minima frames', leave=False):
    ab, coords = [], []
    for point in all_points:
        if polygon.distance(point) <= tolerance:
            ab.append(True)
        else:
            ab.append(False)
    for i,elem2 in enumerate(ab):
        if elem2 == True:
            coords.append([a[i],b[i]])
    tot_min_frames += len(coords)
    sorted_coords.append(coords)
print('processed ' + str(len(a)) + ' frames')
print('found ' + str(tot_min_frames) + ' minima frames')
print('time needed for minima frames identification step: ' + str(round(time.perf_counter() - start1,3)) + ' s')

if not os.path.isdir('minima'):
    os.mkdir('minima')

start3 = time.perf_counter()

if topo == None:
    if traj.endswith('.pdb'):
        f = open(traj, 'r+')
        lines = f.readlines()
        f.close()
        atom_count, head = 0, 0
        for line in lines:
            if line.startswith('ATOM'):
                atom_count += 1
            elif line.startswith('END'):
                break
            elif line.startswith('AUTHOR') or line.startswith('TITLE'):
                head += 1
    else:
        u = mda.Universe(traj, in_memory=True)
else:
    u = mda.Universe(topo, traj, in_memory=True)
    
os.chdir('minima')
with open('min_overview.txt', 'w') as overviewfile:
    pass

for i in range(dimY):
    bins.append(np.zeros(dimX))
for i in range(dimY):
    for l in range(dimX):
        bins[-1-i][l] = ener[next(count2)]

plt.figure(figsize=(8,6), dpi=100)
plt.imshow(bins, interpolation='gaussian', cmap='nipy_spectral')
plt.xticks(np.linspace(-0.5,dimX-0.5,5),np.round(np.linspace(low_max_a,high_max_a, num=5),2))
plt.yticks(np.linspace(-0.5,dimY-0.5,5),np.round(np.linspace(high_max_b,low_max_b, num=5),2))
plt.xlabel('CV1 [a.U.]')
plt.ylabel('CV2 [a.U.]')
plt.axis('tight')
plt.title('threshold: ' + str(round(thresh_val,3)) + ' a.U.')
plt.plot(outline_show_a, outline_show_b, '.', color='white')
cb = plt.colorbar(label='free energy [a.U.]', format="{x:.0f}")
tick_locator = ticker.MaxNLocator(nbins=8)
cb.locator = tick_locator
cb.update_ticks()    
plt.savefig('fes_visual.png',bbox_inches='tight')

for i,elem in enumerate(sorted_coords):
    with open('min_overview.txt', 'a') as overviewfile:
        try:
            overviewfile.writelines('min_' + str(i) + ': <CV1> = ' + str(round(np.mean(sorted_coords[i], axis=0)[0],4)) + ', <CV2> = ' + str(round(np.mean(sorted_coords[i], axis=0)[1],4)) + '\n')
        except IndexError:
            overviewfile.writelines('min_' + str(i) + ': No minima-frames found for this minimum\n')
    if traj.endswith('.pdb'):
        tempfile = open('min_' + str(i) + '.pdb', 'w')
        ref_point = [0,0,0]
    else:
        indx_list = []
    with tqdm.tqdm(total=len(elem), desc='min ' + str(i) + ': ' + str(len(elem)) + ' frames', leave=False) as progress_bar:
        for o in range(len(elem)):
            indx_a = index_search(sorted_coords[i][o][0],a)
            indx_b = index_search(sorted_coords[i][o][1],b)
            both = set(indx_a).intersection(indx_b)
            indx = both.pop()
            if traj.endswith('.pdb'):
                min_coords, shift_vect, print_out = [], [0,0,0], []
                count = it.count(0)
                for k in range(head+(indx*(atom_count+3)), head+((indx+1)*(atom_count+3))):
                    if lines[k].startswith('ATOM'):
                        tmp_count = next(count)
                        print_lines = []
                        str_part_1 = lines[k].split('                 ')[0]
                        str_part_2 = lines[k].split('0.00  0.00')[1]
                        if not len(" ".join(lines[k].split()).split(' ')) < 9:
                            atom_coords =[float(" ".join(lines[k].split()).split(' ')[3]),float(" ".join(lines[k].split()).split(' ')[4]),float(" ".join(lines[k].split()).split(' ')[5])]
                        else:
                            if len(" ".join(lines[k].split()).split(' ')[3]) > 8:
                                tmp_str = " ".join(lines[k].split()).split(' ')[3].split('-')
                                if len(tmp_str) == 2:
                                    atom_coords = [float(tmp_str[0]),float(tmp_str[1])*(-1),float(" ".join(lines[k].split()).split(' ')[4])]
                                elif len(tmp_str) == 4:
                                    atom_coords = [float(tmp_str[1])*(-1),float(tmp_str[2])*(-1),float(tmp_str[3])*(-1)]
                                else:
                                    if len(tmp_str[0]) == 0:
                                        atom_coords = [float(tmp_str[1])*(-1),float(tmp_str[2])*(-1),float(" ".join(lines[k].split()).split(' ')[4])]
                                    else:
                                        atom_coords = [float(tmp_str[0]),float(tmp_str[1])*(-1),float(tmp_str[2])*(-1)]
                            elif len(" ".join(lines[k].split()).split(' ')[4]) > 8:
                                tmp_str = " ".join(lines[k].split()).split(' ')[4].split('-')
                                atom_coords = [float(" ".join(lines[k].split()).split(' ')[3]),float(tmp_str[0]),float(tmp_str[1])*(-1)]
                        if o == 0 and tmp_count == 0:
                            ref_point = atom_coords
                        if not o == 0 and tmp_count == 0:
                            shift_vect = (np.array(ref_point) - np.array(atom_coords))
                        min_coords.append((np.array(atom_coords) + np.array(shift_vect)).tolist())
                        print_lines = ''.join((str_part_1,' '*(25-len(str(round(min_coords[tmp_count][0],3)))),str(round(min_coords[tmp_count][0],3)),' '*(8-len(str(round(min_coords[tmp_count][1],3)))),str(round(min_coords[tmp_count][1],3)),' '*(8-len(str(round(min_coords[tmp_count][2],3)))),str(round(min_coords[tmp_count][2],3)),'  0.00  0.00',str_part_2))  
                        print_out.append(print_lines)
                    elif lines[k].startswith('REMARK'):
                        print_out.append(lines[k])
                    elif lines[k].startswith('CRYST1'):
                        print_out.append(lines[k])                
                print_out.append('END\n')
                tempfile.writelines(print_out)
            else:
                indx_list.append(indx)
            progress_bar.update(1)
    if traj.endswith('.pdb'):        
        tempfile.close()
    else:
        ag =u.select_atoms('all')
        try:
            ag.write('min_' + str(i) + '.' + traj.split('.')[1], frames=u.trajectory[indx_list])
        except IndexError:
            raise Exception('Multiple frames are not supported with this trajectory-format.')
        except Exception:
            print('MDAnalysis does not support writing of ' + traj.split('.')[1] + ' files, writing in xyz')
            ag.write('min_' + str(i) + '.xyz', frames=u.trajectory[indx_list])
            
print('time needed for postprocessing step: ' + str(round(time.perf_counter() - start3,3)) + ' s')
