import numpy as np
import os
import sys
import time 
import itertools as it
import shapely.geometry
import copy
import tqdm

default_args=['','','']
args=sys.argv[3:]
args += default_args[len(args) :]

T = args[0]
thresh = args[1]
frozen_atom_number = args[2]

def index_search(values, l):
    indices = []
    for i, x in enumerate(l):
        if values == x:
            indices.append(i)
    return indices

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
outline, coords = [], []
count1 = it.count(0)
thresh_val, max_size, tot_min_frames, dimX = 0, 0, 0, 0
cur_dir = os.listdir('.')
for elem in cur_dir:
    if elem.startswith('fes_' + T):
        filename = elem

data_fes = np.genfromtxt(filename)
data_fes_T = np.transpose(data_fes)
a_fes, b_fes, ener = data_fes_T[0].copy(), data_fes_T[1].copy(), data_fes_T[2].copy()

if thresh == 'auto':
    thresh_val = abs(min(ener)) - (abs(min(ener))-abs(max(ener)))/3.5
    if min(ener) < 0:
        thresh_val = thresh_val * (-1)
    print('automaticly determined threshold value: ' + str(thresh_val))
else:
    thresh_val = int(thresh)
    
while b_fes[next(count1)] == b_fes[0]:
    dimX += 1
low_max_a, high_max_a, low_max_b, high_max_b = a_fes[0], a_fes[dimX-1], b_fes[0], b_fes[-1]

for i,elem in enumerate(ener):
    try:
        if elem<thresh_val and (a_fes[i] == low_max_a or a_fes[i] == high_max_a or b_fes[i] == low_max_b or b_fes[i] == high_max_b or ener[i-1]>thresh_val or ener[i+1]>thresh_val or ener[i-dimX]>thresh_val or ener[i+dimX]>thresh_val):
            outline.append([a_fes[i],b_fes[i]])
    except:
        pass
    
tolerance = abs(a_fes[0]-a_fes[1])/2
os.chdir(T)     
data_colvar = np.genfromtxt('COLVAR')
data_colvar_T = np.transpose(data_colvar)
step, a, b = data_colvar_T[0].copy(), data_colvar_T[1].copy(), data_colvar_T[2].copy()

start1 = time.perf_counter()
all_points, sorted_coords = [], []
for i, elem in enumerate(a):
    all_points.append(shapely.geometry.Point(elem,b[i]))
grouped_points = group_numbers(outline, 20*np.sqrt(8)*tolerance)
polygons = []
for groups in grouped_points:
    polygons.append(shapely.geometry.Polygon(groups))
    
print(str(len(polygons)) + ' minima identified')
print('determining minima frames:')
for polygon in tqdm.tqdm(polygons):
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

in_dir = os.listdir('.')

for elem2 in in_dir:
	if elem2.endswith('.pdb'):
		size = os.stat(elem2).st_size
		if size > max_size:
			max_size = size
			max_file = elem2
            
f = open(max_file, 'r+')
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
 
if not os.path.isdir('minima'):
    os.mkdir('minima')
os.chdir('minima')

start3 = time.perf_counter()
for i,elem in enumerate(sorted_coords):
    tempfile = open('min_' + str(i) + '.pdb', 'w')
    tempfile.writelines('TITLE     <a> = ' + str(round(np.mean(elem, axis=0)[0],4)) + ', <b> = ' + str(round(np.mean(elem, axis=0)[1],4)) + '\n')
    tot_lines, ref_point = [], [0,0,0]
    with tqdm.tqdm(total=len(elem), desc='min ' + str(i) + ': ' + str(len(elem)) + ' frames', leave=True) as progress_bar:
        for o in range(len(elem)):
            indx_a = index_search(sorted_coords[i][o][0],a)
            indx_b = index_search(sorted_coords[i][o][1],b)
            both = set(indx_a).intersection(indx_b)
            indx = both.pop()
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
                    if o == 0 and tmp_count == frozen_atom_number:
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
            tot_lines.append(print_out)
            progress_bar.update(1)
    for q in range(len(tot_lines)):
        tempfile.writelines(tot_lines[q])
    tempfile.close()
print('time needed for postprocessing step: ' + str(round(time.perf_counter() - start3,3)) + ' s')
    
