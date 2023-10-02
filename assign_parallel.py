import numpy as np
import os
import sys
import time 
import itertools as it
import shapely.geometry
import copy
import tqdm
import multiprocessing as mp
import MDAnalysis as mda
import matplotlib.pyplot as plt
from matplotlib import ticker
import operator as op
import shutil
import warnings
warnings.filterwarnings('ignore')

try:
  topo = sys.argv[4]
except IndexError:
  topo = None
T, thresh, traj = sys.argv[1:]

def init_worker_pdb(sorted_coords,a,b,lines,atom_count,head,desc,barargs):
    tqdm.tqdm.set_lock(barargs)
    global Gsorted_coords
    global Ga
    global Gb
    global Gatom_count
    global Ghead
    global Glines
    global Gdesc
    Gatom_count, Ghead, Glines, Gsorted_coords, Ga, Gb, Gdesc = atom_count, head, lines, sorted_coords, a, b, desc

def init_worker(sorted_coords,a,b,u,ag,desc,barargs):
    tqdm.tqdm.set_lock(barargs)
    global Gsorted_coords
    global Ga
    global Gb 
    global Gu
    global Gag
    global Gdesc
    Gsorted_coords, Ga, Gb, Gu, Gag, Gdesc = sorted_coords, a, b, u, ag, desc

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
        try:
            tmplist.remove(seed_elem)
            new_group_found = False
        except ValueError:
            pass
        for compare_elem in tmplist: 
            if ((seed_elem[0]-compare_elem[0])**2 + (seed_elem[1]-compare_elem[1])**2)**0.5 < min_distance: 
                found = True
                min_distance = ((seed_elem[0]-compare_elem[0])**2 + (seed_elem[1]-compare_elem[1])**2)**0.5
                min_elem = compare_elem
        if found == True and any(subgroup):
            if new_group_found == False:
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
                    if  ((seed_elem[0]-elem[0])**2 + (seed_elem[1]-elem[1])**2)**0.5 < min_distance:
                        group.append(seed_elem)
                        sec_run = True
                        break
                if sec_run == True:
                    break
            if sec_run == False:
                subgroup.append(seed_elem)
                new_group_found = True
            elif any(tmplist):
                seed_elem = tmplist[0]
    return separate_groups

def sort(i):
    indx_list = []
    with open('min_overview.txt', 'a') as overviewfile:
        if Gdesc:
            overviewfile.writelines('min_' + str(i) + ': ' + Gdesc[i] + '\n')
        else:
            overviewfile.writelines('min_' + str(i) + ': CV1: ' + str(round(np.mean(grouped_points[i], axis=0)[0],4)) + ' CV2: ' + str(round(np.mean(grouped_points[i], axis=0)[1],4)) + '\n')
    pos = mp.current_process()._identity[0]-1
    with tqdm.tqdm(total=len(Gsorted_coords[i]), desc='min ' + str(i) + ': ' + str(len(Gsorted_coords[i])) + ' frames', position=pos, leave=False) as progress_bar:
        for o in range(len(Gsorted_coords[i])):
            indx_a = index_search(Gsorted_coords[i][o][0],Ga)
            indx_b = index_search(Gsorted_coords[i][o][1],Gb)
            both = set(indx_a).intersection(indx_b)
            indx = both.pop()
            indx_list.append(indx)
            progress_bar.update(1)
    try:
        Gag.write('min_' + str(i) + '.' + traj.split('.')[1], frames=Gu.trajectory[indx_list])
    except (IndexError, NameError):
        if traj.endswith('.pdb'):
            tempfile = open('min_' + str(i) + '.pdb', 'w')
            ref_point = [0,0,0]
            for o,elem_inner in enumerate(indx_list):
                min_coords, shift_vect, print_out = [], [0,0,0], []
                count = it.count(0)
                for k in range(Ghead+(elem_inner*(Gatom_count+3)), Ghead+((elem_inner+1)*(Gatom_count+3))):
                    if Glines[k].startswith('ATOM'):
                        tmp_count = next(count)
                        print_lines = []
                        str_part_1 = Glines[k].split('                 ')[0]
                        str_part_2 = Glines[k].split('0.00  0.00')[1]
                        if not len(" ".join(Glines[k].split()).split(' ')) < 9:
                            atom_coords =[float(" ".join(Glines[k].split()).split(' ')[3]),float(" ".join(Glines[k].split()).split(' ')[4]),float(" ".join(Glines[k].split()).split(' ')[5])]
                        else:
                            if len(" ".join(Glines[k].split()).split(' ')[3]) > 8:
                                tmp_str = " ".join(Glines[k].split()).split(' ')[3].split('-')
                                if len(tmp_str) == 2:
                                    atom_coords = [float(tmp_str[0]),float(tmp_str[1])*(-1),float(" ".join(Glines[k].split()).split(' ')[4])]
                                elif len(tmp_str) == 4:
                                    atom_coords = [float(tmp_str[1])*(-1),float(tmp_str[2])*(-1),float(tmp_str[3])*(-1)]
                                else:
                                    if len(tmp_str[0]) == 0:
                                        atom_coords = [float(tmp_str[1])*(-1),float(tmp_str[2])*(-1),float(" ".join(Glines[k].split()).split(' ')[4])]
                                    else:
                                        atom_coords = [float(tmp_str[0]),float(tmp_str[1])*(-1),float(tmp_str[2])*(-1)]
                            elif len(" ".join(Glines[k].split()).split(' ')[4]) > 8:
                                tmp_str = " ".join(Glines[k].split()).split(' ')[4].split('-')
                                atom_coords = [float(" ".join(Glines[k].split()).split(' ')[3]),float(tmp_str[0]),float(tmp_str[1])*(-1)]
                        if o == 0 and tmp_count == 0:
                            ref_point = atom_coords
                        if not o == 0 and tmp_count == 0:
                            shift_vect = (np.array(ref_point) - np.array(atom_coords))
                        min_coords.append((np.array(atom_coords) + np.array(shift_vect)).tolist())
                        print_lines = ''.join((str_part_1,' '*(25-len(str(round(min_coords[tmp_count][0],3)))),str(round(min_coords[tmp_count][0],3)),' '*(8-len(str(round(min_coords[tmp_count][1],3)))),str(round(min_coords[tmp_count][1],3)),' '*(8-len(str(round(min_coords[tmp_count][2],3)))),str(round(min_coords[tmp_count][2],3)),'  0.00  0.00',str_part_2))  
                        print_out.append(print_lines)
                    elif Glines[k].startswith('REMARK'):
                        print_out.append(Glines[k])
                    elif Glines[k].startswith('CRYST1'):
                        print_out.append(Glines[k])                
                print_out.append('END\n')           
                tempfile.writelines(print_out) 
            tempfile.close()
        else:
            raise Exception('Multiple frames are not supported with this trajectory-format.')
    except (TypeError, ValueError):
        print('MDAnalysis does not support writing in ' + traj.split('.')[1] + '-format, writing in xyz-format instead')
        Gag.write('min_' + str(i) + '.xyz', frames=Gu.trajectory[indx_list])

def sort_pdb_cp2k_prework():
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
    return head, atom_count, lines

if __name__ == '__main__':
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
        thresh_val = float(thresh)
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
    polygons = []
    for groups in grouped_points:
        polygons.append(shapely.geometry.Polygon(groups))
        
    periodicity = False
    if edge:
        edge_points, pbc = [], []
        grouped_edges = group_numbers(edge, 10*np.sqrt(8)*tolerance)
        for i,group1 in enumerate(grouped_edges):
            if sum(list(map(len, pbc))) >= len(grouped_edges):
                break
            expect_group, tmp_lst = [], []
            for elem in group1:
                tmp_pt = copy.deepcopy(elem)
                if elem[0] == high_max_a:
                    tmp_pt[0] = low_max_a
                elif elem[0] == low_max_a:
                    tmp_pt[0] = high_max_a
                if elem[1] == high_max_b:
                    tmp_pt[1] = low_max_b
                elif elem[1] == low_max_b:
                    tmp_pt[1] = high_max_b
                expect_group.append(tmp_pt)
            found_periodic = False
            for j,group2 in enumerate(grouped_points):
                if have_common_elem(group2, expect_group) or have_common_elem(group2, group1):
                    periodicity = True
                    found_periodic = True
                    tmp_lst.append(j)
            if found_periodic == True:
                if len(tmp_lst) == 1:
                    break
                elif i == 0:
                    print('periodicity detected: boundaries will be considered periodic')
                pbc.append(tmp_lst)
    print(str(len(polygons)), end = ' ')
    if periodicity == True:
        print('distinctive areas identified')
    else:
        print('minima identified')

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
    
    desc = []
    if periodicity == True:
        sorted_coords_period, tot_pbc  = [], []
        for elem in pbc:
            desc.append(' + '.join(('CV1: ' + str(round(np.mean(grouped_points[j], axis=0)[0],4)) + ' CV2: ' + str(round(np.mean(grouped_points[j], axis=0)[1],4))) for j in elem))
            help_list = []
            for i in elem:
                tot_pbc.append(i)
                help_list += sorted_coords[i]
            sorted_coords_period.append(help_list)
        for i,elem in enumerate(sorted_coords):
            if not i in tot_pbc:
                desc.append('CV1: ' + str(round(np.mean(grouped_points[i], axis=0)[0],4)) + ' CV2: ' + str(round(np.mean(grouped_points[i], axis=0)[1],4)))
                sorted_coords_period.append(elem)
        sorted_coords = sorted_coords_period
        print(str(len(sorted_coords)) + ' minima identified')

    try:
        os.mkdir('minima')
    except FileExistsError:
        shutil.rmtree('minima')
        os.mkdir('minima')

    try:
        if topo == None:
            u = mda.Universe(traj, in_memory=True)
        else:
            u = mda.Universe(topo, traj, in_memory=True)
        ag =u.select_atoms('all')
    except IndexError:
        if traj.endswith('.pdb'):
            head, atom_count, lines = sort_pdb_cp2k_prework()
        else:
            raise Exception('MDAnalysis does not support the topology- or trajectory-file.')
    except FileNotFoundError:
        raise
    except Exception:
        raise Exception('MDAnalysis does not support the topology- or trajectory-file.')
        
    for i in range(dimY):
        bins.append(np.zeros(dimX))
    for i in range(dimY):
        for l in range(dimX):
            bins[-1-i][l] = ener[next(count2)]

    plt.figure(figsize=(8,6), dpi=100)
    plt.imshow(bins, interpolation='gaussian', cmap='nipy_spectral')
    plt.xticks(np.linspace(low_max_a,dimX-np.round(high_max_a,1),5),np.round(np.linspace(low_max_a,high_max_a, num=5),2))
    plt.yticks(np.linspace(low_max_b,dimY-np.round(high_max_b,1),5),np.round(np.linspace(high_max_b,low_max_b, num=5),2))
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
    
    if len(sorted_coords) > os.cpu_count()-1:
        usable_cpu = os.cpu_count()-1
    else:
        usable_cpu = len(sorted_coords)
    
    start3 = time.perf_counter()
    try:
        pool = mp.Pool(processes = usable_cpu, initializer=init_worker, initargs=(sorted_coords,a,b,u,ag,desc,mp.RLock(),)) 
    except NameError:
        pool = mp.Pool(processes = usable_cpu, initializer=init_worker_pdb, initargs=(sorted_coords,a,b,lines,atom_count,head,desc,mp.RLock(),))
    out = pool.map(sort, range(len(sorted_coords)))
    pool.close()
    pool.join()

    print('time needed for postprocessing step: ' + str(round(time.perf_counter() - start3,3)) + ' s')
