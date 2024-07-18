import numpy as np
import os
import time 
import itertools as it
import shapely.geometry
import copy
import tqdm
import multiprocessing as mp
import matplotlib.pyplot as plt
from matplotlib import ticker
import operator as op
import shutil
import warnings
import argparse
import MDAnalysis as mda

warnings.filterwarnings('ignore', category=UserWarning)

parser = argparse.ArgumentParser()

parser.add_argument('-traj', dest='traj', required=True, help='!REQUIRED! MD trajectory-file name in the MD-output-directory. Format is also used for output-files.', type=str)
parser.add_argument('-md', dest='md_dir', default=os.getcwd(), help='MD-output-directory path. DEFAULT: Current directory path.')
parser.add_argument('-thresh', dest='thresh', default=None, help='Specifies threshold for assigning. Input value has to correspond with values in FES-file. DEFAULT: Lowest 1/12 of the energy span.', type=float)
parser.add_argument('-topo', dest='topo', default=None, help='MD topology-file name in the MD-output-directory, if trajectory-file does not specify topology. DEFAULT: None.', type=str)
parser.add_argument('-fes', dest='fes', default='fes.dat', help='FES-file name in the MD-output-directory. DEFAULT: "fes.dat".', type=str)
parser.add_argument('-colv', dest='colvar', default='COLVAR', help='COLVAR-file in the MD-output-directory. DEFAULT: "COLVAR".')
parser.add_argument('-png', dest='fes_png', default=True, help='Specifies whether a PNG-visualization of the FES should be created. Expects True/False. DEFAULT: True.', type=bool)
parser.add_argument('-nopbc', dest='nopbc', default=False, help='Suppresses the automatic periodicity search (triggered when minima touch the edges). Expects True/False. DEFAULT: False.', type=bool)
parser.add_argument('-mindist', dest='mindist', default=10, help='Smallest allowed distance at which areas are considered separate minima (unit: bins of FES-histogram). Must be larger than 1. DEFAULT: 10.', type=float)
parser.add_argument('-stride', dest='stride', default=1, help='Reads only every n-th frame of trajectory. DEFAULT: 1.', type=int)

args = parser.parse_args()
    
def init_polygon(Ga,Gb,Gtolerance,Gmp_sorted_coords,Ggrouped_points, barargs):
    tqdm.tqdm.set_lock(barargs)
    global all_points
    global tolerance
    global managed_list
    global grouped_points
    global a
    global b
    a,b, tolerance, managed_list, grouped_points = Ga,Gb, Gtolerance, Gmp_sorted_coords, Ggrouped_points
    
def init_worker_pdb(Gsorted_indx, barargs):
    tqdm.tqdm.set_lock(barargs)
    global sorted_indx
    sorted_indx = Gsorted_indx
    
def det_min_frames(j):
    all_points = [shapely.geometry.Point(a[i],b[i]) for i in range(len(a))]
    polygon = shapely.geometry.Polygon(grouped_points[j])
    convex_hull = shapely.geometry.Polygon(polygon.convex_hull)
    if abs(1-(polygon.area/convex_hull.area)) > 0.4:
        polygon = convex_hull
        print('polygon ' + str(j) + ' did not initialize properly, using convex-hull')
    indxes = []
    with tqdm.tqdm(total=len(all_points), desc='min ' + str(j), position=tqdm.tqdm._get_free_pos()+j, leave=False) as pbar:
        for i,point in enumerate(all_points):
            if polygon.distance(point) <= tolerance:
                indxes.append(i)
            pbar.update(1)
        managed_list[j] = indxes

def have_common_elem(l1, l2):
    for elem in l2:
        if op.countOf(l1, elem) > 0:
            return True
            break
    return False

def ex3(numbers, max_diff):
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
                dists = np.empty(len(group))
                for i,elem in enumerate(group):
                    dists[i] = ((seed_elem[0]-elem[0])**2 + (seed_elem[1]-elem[1])**2)**0.5
                    if dists[i] < min_distance:
                        sec_run = True
                if sec_run == True:
                    try:
                        nih = np.empty(len(dists)-1)
                        for j in range(len(dists)-1):
                                nih[j] = dists[j]+dists[j+1]
                        group.insert(np.argmin(nih)+1, seed_elem)
                    except ValueError:
                        group.append(seed_elem)
                    break
            if sec_run == False:
                subgroup.append(seed_elem)
                new_group_found = True
            elif any(tmplist):
                seed_elem = tmplist[0]
    
    connect_groups = []
    min_distance = max_diff          
    for g1 in range(0,len(separate_groups)):
        for g2 in range(g1+1, len(separate_groups)):
            dists = np.empty(len(separate_groups[g1])*len(separate_groups[g2]))
            indx = it.count(0)
            for e1 in separate_groups[g1]:
                for e2 in separate_groups[g2]:
                    dists[next(indx)] = (((e1[0]-e2[0])**2+(e1[1]-e2[1])**2)**0.5)
            if np.min(dists) <= min_distance:
                connect_groups.append([g1,g2])
    grouped_connected_groups = []
    while len(connect_groups)>0:
        first, *rest = connect_groups
        first = set(first)
        lf = -1
        while len(first)>lf:
            lf = len(first)
            rest2 = []
            for r in rest:
                if len(first.intersection(set(r)))>0:
                    first |= set(r)
                else:
                    rest2.append(r)     
            rest = rest2
        grouped_connected_groups.append(first)
        connect_groups = rest
        
    fin_sep_groups, tot  = [], []
    for elem in grouped_connected_groups:
        help_list = []
        for i in elem:
            tot.append(i)
            help_list += separate_groups[i]
        fin_sep_groups.append(help_list)
    for i,elem in enumerate(separate_groups):
        if not i in tot:
            fin_sep_groups.append(elem)
        
    return fin_sep_groups
    
def printout(i):
    if args.stride > 1:
        sorted_indx[i] = [args.stride*q for q in sorted_indx[i]]
    try:
        ag.write('minima/min_' + str(i) + '.' + args.traj.split('.')[-1], frames=u.trajectory[sorted_indx[i]])
    except (TypeError, ValueError):
        print('MDAnalysis does not support writing in ' + args.traj.split('.')[-1] + '-format, writing in xyz-format instead')
        ag.write('minima/min_' + str(i) + '.xyz', frames=u.trajectory[sorted_indx[i]])
        
def printout_pdb_cp2k(i):
    if args.stride > 1:
        sorted_indx[i] = [args.stride*q for q in sorted_indx[i]]
    linecount, printcount = 0, 0
    with open('minima/min_' + str(i) + '.pdb', 'w') as minfile:
        with open(args.traj, 'r') as ftraj:
            with tqdm.tqdm(total=len(sorted_indx[i]), desc='writing min ' + str(i), position=tqdm.tqdm._get_free_pos()+i, leave=False) as pbar:
                for line in ftraj:
                    try:
                        if sorted_indx[i][printcount] == linecount:
                            minfile.write(line)
                            if line.startswith('END'):
                                printcount += 1
                                pbar.update(1)
                    except IndexError:
                        break
                    if line.startswith('END'):
                        linecount += 1

def printout_pdb_cp2k_prework():
    with open(args.traj, 'rb') as bfile:
        line_count = sum(1 for _ in bfile)
    with open(args.traj, 'r') as tfile:
        atom_count, head = 0, 0
        for line in tfile:
            if line.startswith('ATOM'):
                atom_count += 1
            elif line.startswith('END'):
                break
            elif line.startswith('AUTHOR') or line.startswith('TITLE'):
                head += 1
    return int((line_count-head)/(atom_count+3))

if __name__ == '__main__':
    title = '.: Free Energy Surface Trajectory Analysis - FESTA :.'
    termin = '.: terminated successfully :.'
    try:
        terminal_size = os.get_terminal_size()[0]
    except OSError:
        terminal_size = len(title)
    print(' '*terminal_size)
    print(' '*int((terminal_size-len(title))/2) + title + ' '*int((terminal_size-len(title))/2))
    print(' '*terminal_size)
    print('working on directory: ' + args.md_dir)
    if args.mindist <= 1:
        print('Invalid minimal separation distance specified. Using default.')
        args.mindist = 10
    outline = []
    count1, count2 = it.count(0), it.count(0)
    thresh_val, tot_min_frames, dimX, dimY = 0, 0, 0, 1
    os.chdir(args.md_dir)
    
    with open(args.colvar, 'r') as colvar_file:
        col_var = colvar_file.readline().split(' ')
        col_var[-1] = col_var[-1][:-1]
    with open(args.fes, 'r') as fes_file:
        fes_var = fes_file.readline().split(' ')
        fes_var[-1] = fes_var[-1][:-1]

    pos_ener = fes_var.index('file.free')-2
    com = list(set(col_var).intersection(fes_var))
    pos_cvs_fes, pos_cvs_col = [], []
    for elem in com:
        if not elem == '#!' and not elem == 'FIELDS':
            pos_cvs_fes.append(fes_var.index(elem)-2)
            pos_cvs_col.append(col_var.index(elem)-2)            
    pos_cvs_fes.sort()
    pos_cvs_col.sort()
    
    if not len(pos_cvs_fes) == 2 or not len(pos_cvs_col) == 2:
        raise Exception('Only MD-runs with 2 CVs supported')
    
    a_fes, b_fes, ener = np.loadtxt(args.fes, unpack=True, usecols=(pos_cvs_fes[0], pos_cvs_fes[1], pos_ener))
    
    for i in range(len(ener)):
        if not np.isfinite(ener[i]):
            raise Exception('Non-finite value (NaN or inf) discovered in FES-file')
    
    if args.thresh == None:
        thresh_val = max(ener) - abs(max(ener)-min(ener))*(1-1/12)
        print('automatically determined', end =' ') 
    else:
        thresh_val = args.thresh
    print('threshold value: ' + str(round(thresh_val,3)) + ' a.U.')
        
    b_central = True
    if (b_fes[0] == b_fes[1]):
        while b_fes[next(count1)] == b_fes[0]:
            dimX += 1
        b_count = b_fes[0]
        for elem in b_fes:
            if not elem == b_count:
                dimY += 1
                b_count = elem
        high_max_a, high_max_b  = a_fes[dimX-1], b_fes[-1]
        tolX = abs(a_fes[0]-a_fes[1])/2
        tolY = abs(b_fes[0]-b_fes[dimY])/2
    else:
        b_central = False
        while a_fes[next(count1)] == a_fes[0]:
            dimX += 1
        a_count = a_fes[0]
        for elem in a_fes:
            if not elem == a_count:
                dimY += 1
                a_count = elem
        high_max_a, high_max_b  = a_fes[-1], b_fes[dimX-1]
        tolX = abs(a_fes[0]-a_fes[dimY])/2
        tolY = abs(b_fes[0]-b_fes[1])/2
    low_max_a, low_max_b = a_fes[0], b_fes[0]
    outline_show_a, outline_show_b, edge = [], [], []
    
    for i in tqdm.tqdm(range(len(ener)),desc='collecting outline', leave=False):
        try:
            if ener[i]<thresh_val and (a_fes[i] == low_max_a or a_fes[i] == high_max_a or b_fes[i] == low_max_b or b_fes[i] == high_max_b or ener[i-1]>thresh_val or ener[i+1]>thresh_val or ener[i-dimY]>thresh_val or ener[i+dimY]>thresh_val or ener[i+1+dimY]>thresh_val or ener[i-1+dimY]>thresh_val or ener[i+1-dimY]>thresh_val or ener[i-1-dimY]>thresh_val):
                if a_fes[i] == low_max_a or a_fes[i] == high_max_a or b_fes[i] == low_max_b or b_fes[i] == high_max_b:
                    edge.append([a_fes[i],b_fes[i]])
                outline.append([a_fes[i],b_fes[i]])
                outline_show_a.append(abs((a_fes[i]-low_max_a)/((high_max_a-low_max_a)/dimX)))
                outline_show_b.append(dimY-abs((b_fes[i]-low_max_b)/((high_max_b-low_max_b)/dimY)))
        except IndexError:
            pass
    
    a, b = np.loadtxt(args.colvar, unpack=True, usecols=(pos_cvs_col[0], pos_cvs_col[1]))
    if args.stride > 1:
        a, b = a[0::args.stride], b[0::args.stride]
    
    print('reading trajectory in ... ' , end='', flush=True)
    cp2k_pdb = False
    try:
        if args.topo == None:
            if args.traj.split('.')[-1] == 'lammpstrj':
                u = mda.Universe(args.traj, topology_format='LAMMPSDUMP')
            else:
                u = mda.Universe(args.traj)
        else:
            u = mda.Universe(args.topo, args.traj)
        ag = u.select_atoms('all')
        if not int((len(u.trajectory)-1)/args.stride+1) == len(a):
            raise Exception('COLVAR-file and trajectory-file must have similar step length, here: ' + str(len(a)) + ' vs ' + str(int((len(u.trajectory)-1)/args.stride+1)))
    except IndexError:
        if args.traj.endswith('.pdb'):
            frame_count = printout_pdb_cp2k_prework()
            cp2k_pdb = True
            if not int((frame_count-1)/args.stride+1) == len(a):
                raise Exception('COLVAR-file and trajectory-file must have similar step length, here: ' + str(len(a)) + ' vs ' + str(int((frame_count-1)/args.stride+1)))
        else:
            raise Exception('MDAnalysis does not support this topology- or trajectory-file')
    except FileNotFoundError:
        raise
    print('done')
    
    start0 = time.perf_counter() 
    grouped_points = ex3(outline, args.mindist*2*np.sqrt(tolX**2+tolY**2))
    grouped_points = [groups for groups in grouped_points if len(groups)>3]
    print('time needed for CCL step: ' + str(round(time.perf_counter() - start0,3)) + ' s')

    start1 = time.perf_counter()
    periodicity = False
    if edge and args.nopbc == False:
        edge_points, pbc = [], []
        grouped_edges = ex3(edge, 10*2*np.sqrt(tolX**2+tolY**2))
        for i in range(len(grouped_edges)):
            if sum(list(map(len, pbc))) >= len(grouped_edges):
                break
            expect_group, tmp_lst = [], []
            for elem in grouped_edges[i]:
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
                if have_common_elem(group2, expect_group) or have_common_elem(group2, grouped_edges[i]):
                    periodicity = True
                    found_periodic = True
                    tmp_lst.append(j)
            if found_periodic == True:
                if len(tmp_lst) == 1:
                    break
                elif i == 0:
                    print('periodicity detected: boundaries will be considered periodic')
                pbc.append(tmp_lst)
    print(str(len(grouped_points)), end = ' ')
    if periodicity == True:
        print('distinctive areas identified')
    else:
        if len(grouped_points) == 1:
            print('minimum identified')
        else:
            print('minima identified')
    
    tot_min_frames =0
    if len(grouped_points) > os.cpu_count()-1:
        usable_cpu = os.cpu_count()-1
    else:
        usable_cpu = len(grouped_points)
    mp_sorted_coords = mp.Manager().list([[] for _ in range(len(grouped_points))])
    with mp.Pool(processes=usable_cpu, initializer=init_polygon, initargs=(a,b,np.sqrt(tolY**2+tolX**2), mp_sorted_coords,grouped_points,mp.RLock(),)) as pool:
        pool.map(det_min_frames, range(len(grouped_points)))
    sorted_indx = list(mp_sorted_coords)
    for lists in sorted_indx:
        tot_min_frames += len(lists)
    print('processed ' + str(len(a)) + ' frames')
    print('found ' + str(tot_min_frames) + ' minima frames')
    print('time needed for minima frames identification step: ' + str(round(time.perf_counter() - start1,3)) + ' s')
    
    desc = []
    if periodicity == True:
        sorted_coords_period, tot_pbc  = [], []
        for elem in pbc:
            desc.append(' + '.join((fes_var[pos_cvs_fes[0]+2] + ': ' + str(round(np.mean(grouped_points[j], axis=0)[0],4)) + ' '+fes_var[pos_cvs_fes[1]+2]+': ' + str(round(np.mean(grouped_points[j], axis=0)[1],4))) for j in elem))
            help_list = []
            for i in elem:
                tot_pbc.append(i)
                help_list += sorted_indx[i]
            if cp2k_pdb == True:
                help_list.sort()
            sorted_coords_period.append(help_list)
        for i,elem in enumerate(sorted_indx):
            if not i in tot_pbc:
                desc.append(fes_var[pos_cvs_fes[0]+2]+': ' + str(round(np.mean(grouped_points[i], axis=0)[0],4)) + ' '+fes_var[pos_cvs_fes[1]+2]+': ' + str(round(np.mean(grouped_points[i], axis=0)[1],4)))
                sorted_coords_period.append(elem)
        sorted_indx = sorted_coords_period
        print(str(len(sorted_indx)), end=' ')
        if len(sorted_indx) == 1:
            print('minimum identified')
        else:
            print('minima identified')
    try:
        os.mkdir('minima')
    except FileExistsError:
        shutil.rmtree('minima')
        os.mkdir('minima')

    if args.fes_png == True:
        bins = np.empty((dimY,dimX))
        if b_central == True:
            for i in range(dimY):
                for j in range(dimX):
                    bins[-1-i,j] = ener[next(count2)]
        else:
            for i in range(dimX):
                for j in range(dimY):
                    bins[-1-j,i] = ener[next(count2)]
    
        plt.figure(figsize=(8,6), dpi=300)
        plt.imshow(bins, interpolation='gaussian', cmap='nipy_spectral')
        plt.xticks(np.linspace(-0.5,dimX-0.5,5),np.round(np.linspace(low_max_a,high_max_a, num=5),3))
        plt.yticks(np.linspace(-0.5,dimY-0.5,5),np.round(np.linspace(high_max_b,low_max_b, num=5),3))
        plt.xlabel(fes_var[pos_cvs_fes[0]+2] + ' [a.U.]')
        plt.ylabel(fes_var[pos_cvs_fes[1]+2] + ' [a.U.]')
        plt.axis('tight')
        plt.title('threshold: ' + str(round(thresh_val,3)) + ' a.U.')
        plt.plot(outline_show_a, outline_show_b, '.', color='white', markersize=2)
        cb = plt.colorbar(label='free energy [a.U.]', format="{x:.0f}")
        tick_locator = ticker.MaxNLocator(nbins=8)
        cb.locator = tick_locator
        cb.update_ticks()    
        plt.savefig('minima/fes_visual.png',bbox_inches='tight')
    
    start3 = time.perf_counter()
    with open('minima/min_overview.txt', 'w') as overviewfile:
        for i in range(len(sorted_indx)):
            if desc:
                overviewfile.writelines('min_' + str(i) + ': ' + desc[i] + '\n')
            else:
                overviewfile.writelines('min_' + str(i) + ': '+ fes_var[pos_cvs_fes[0]+2] +': ' + str(round(np.mean(grouped_points[i], axis=0)[0],4)) + ' '+fes_var[pos_cvs_fes[1]+2]+': ' + str(round(np.mean(grouped_points[i], axis=0)[1],4)) + '\n')
    if cp2k_pdb == False:
        for i in tqdm.tqdm(range(len(sorted_indx)), desc='printing to file', leave=False):
            printout(i)
    else:
        with mp.Pool(processes = usable_cpu, initializer=init_worker_pdb, initargs=(sorted_indx,mp.RLock(),)) as pool:
            pool.map(printout_pdb_cp2k, range(len(sorted_indx)))

    print('time needed for postprocessing step: ' + str(round(time.perf_counter() - start3,3)) + ' s')
    print(' '*terminal_size)
    print(' '*int((terminal_size-len(termin))/2) + termin + ' '*int((terminal_size-len(termin))/2))
