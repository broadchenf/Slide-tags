"""
SlidePipe

matching bead barcode 
"""
import os
from collections import Counter
import editdistance
from multiprocessing import Pool
import logging
import argparse
import numpy as np

bead_barcode_dir = "/path/to/Barcodes"

def get_barcode_position(bead_bc_file,bead_pos_file):
    def low_cpx_bc(bc,max_homo):
        return Counter(list(bc)).most_common(1)[0][1]>max_homo
    tmp = open(bead_pos_file).readlines()
    pos_x = [float(it) for it in tmp[0].strip().split(",")]
    pos_y =  [float(it) for it in tmp[1].strip().split(",")]
    bc_list = [it.strip().replace(",","") for it in open(bead_bc_file).readlines()]
    max_homo = int(0.8*len(list(bc_list[0])))  # assume all have the same length
    bc_pos_dict = {}
    for bc,coordx,coordy in zip(bc_list,pos_x,pos_y):
        if not low_cpx_bc(bc,max_homo):
            bc_pos_dict[bc] = (coordx,coordy)
    return bc_pos_dict, bc_list


def build_6mer_dist(bc_list):
    start_km = {}
    mid_km = {}
    end_km = {}
    for bc in bc_list:
        start_km.setdefault(bc[:6] , []).append(bc)
        mid_km.setdefault(bc[4:10], []).append(bc)
        end_km.setdefault(bc[-6:] , []).append(bc)
    return start_km,mid_km,end_km


def barcode_matching(bc_pos_dict, spatial_bc_list, max_dist=1):
    bc_matching_dict = {}
    bc_ref_list = list(bc_pos_dict.keys())
    start_km, mid_km, end_km = build_6mer_dist(bc_ref_list)
    exact_match = 0
    fuzzy_match = 0

    bc_list_array = np.array(bc_ref_list)
    spatial_bc_array = np.array(spatial_bc_list)

    start_km_array = np.array(list(start_km.keys()))
    mid_km_array = np.array(list(mid_km.keys()))
    end_km_array = np.array(list(end_km.keys()))

    for bc in spatial_bc_array:
        if bc in bc_pos_dict:
            exact_match += 1
            bc_matching_dict[bc] = bc
        else:
            sel_bc = get_sel_bc(bc, start_km_array, mid_km_array, end_km_array)
            if len(sel_bc) > 0:
                fz = np.array([(it, editdistance.eval(it, bc)) for it in sel_bc])
                fz = fz[fz[:, 1] <= max_dist]
                fz = fz[fz[:, 1].argsort(kind='mergesort')]
                if len(fz) == 0:
                    continue
                if len(fz) > 1 and fz[0][1] == fz[1][1]:
                    if editdistance.eval(fz[0][0][:-1], bc[-1]) > editdistance.eval(fz[1][0][:-1], bc[-1]):
                        fuzzy_match += 1
                        bc_matching_dict[bc] = fz[1][0]
                    elif editdistance.eval(fz[0][0][:-1], bc[-1]) < editdistance.eval(fz[1][0][:-1], bc[-1]):
                        fuzzy_match += 1
                        bc_matching_dict[bc] = fz[0][0]
                else:
                    fuzzy_match += 1
                    bc_matching_dict[bc] = fz[0][0]

    return bc_matching_dict, exact_match, fuzzy_match


def get_sel_bc(bc, start_km_array, mid_km_array, end_km_array):
    sel_bc = []
    if bc[:6] in start_km_array:
        sel_bc += list(start_km[bc[:6]])
    if bc[-6:] in end_km_array:
        sel_bc += list(end_km[bc[-6:]])
    if bc[4:10] in mid_km_array:
        sel_bc += list(mid_km[bc[4:10]])
    return set(sel_bc)

def split_list(li, trunk_size=3000):
    new_li = []
    new_st = 0
    if len(li)<trunk_size:
        return [li]
    while True:
        if (new_st+trunk_size)<len(li):
            new_li.append(li[new_st:(new_st+trunk_size)])
            new_st += trunk_size
        else:
            new_li.append(li[new_st:] )
            break
    return new_li


def write_barcode_match(bc_matching_dict,bc_pos_dict,output_file):
    """_summary_

    Args:
        bc_matching_dict (dict): key:Illumina barcode, value: matched barcode
        bc_pos_dict (dict): key:puck barcode, value: (x, y)
        output_file (string): absolute path of output file
    """
    with open(output_file,"w") as fo:
        fo.write("Illumina_barcode,matched_beadbarcode,xcoord,ycoord\n")
        for bc in bc_matching_dict:
            fo.write("{},{},{},{}\n".format(bc,
                                          bc_matching_dict[bc],
                                          bc_pos_dict[bc_matching_dict[bc]][0],
                                          bc_pos_dict[bc_matching_dict[bc]][1]))


def main(puck_id,barcode_file,output_file):
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    bead_bc_file = os.path.join(bead_barcode_dir,puck_id,"BeadBarcodes.txt")
    bead_pos_file = os.path.join(bead_barcode_dir,puck_id,"BeadLocations.txt")
    bc_pos_dict,_ = get_barcode_position(bead_bc_file,bead_pos_file)
    logging.info("Read {} bead barcodes from puckid folder.".format(len(bc_pos_dict)) )
    spatial_bc_list = open(barcode_file).readline().strip().split(",")
    logging.info("Number of barcode from Illumina reads: {}".format(len(spatial_bc_list)))
    logging.info("First 5 barcode: {}".format(spatial_bc_list[:5]))
    logging.info("Start matching...")
    spatial_bc_list_split = split_list(spatial_bc_list)
    with Pool(8) as p:
        res_list = p.starmap(barcode_matching, [(bc_pos_dict, it) for it in spatial_bc_list_split ] )  # match bead barcode from Illumina reads to bead barcode from puck sequencing
        bc_matching_dict = {}
        for d,_,_ in res_list:
            for k in d:
                bc_matching_dict[k] = d[k]
        exact_match = sum(it[1] for it in res_list)
        fuzzy_match = sum(it[2] for it in res_list)
    logging.info("Matching finished...")
    logging.info("Number of exact match: {}".format(exact_match))
    logging.info("Number of fuzzy match: {}".format(fuzzy_match))
    logging.info("Writing output to file: {}".format(output_file))
    write_barcode_match(bc_matching_dict,bc_pos_dict,output_file)
    logging.info("Done!")

def get_args():
    parser = argparse.ArgumentParser(description='Match bead barcode acquired from Illumina reads to the referecne list from puck sequencing.')
    parser.add_argument(
        "-b", "--barcode_file",
        help="the file that contain barcodes from Illumina reads, the first line is barcode separated by comma.",
        type=str,
        required=True
        )
    parser.add_argument(
        "-o", "--outcsv",
        help="output csv file .",
        type=str,
        required=True
        )
    parser.add_argument(
        "-i", "--puckid",
        help="the puck_id such as Puck_211004_34. will be looking for folder with the same name in /broad/macosko/data/Slideseq/Barcode.",
        type=str,
        required=True
        )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    main(args.puckid,args.barcode_file,args.outcsv)
    # python bead_matching.py -b rc_dge.csv -i Puck_220201_09 -o matching_result.csv


