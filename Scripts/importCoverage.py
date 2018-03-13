#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 20:59:56 2017

@author: login
"""

import os,  sys
import cPickle as pickle
import pandas as pd 
import numpy as np

def loadFullMetadata(verbose=True):
    code_dir = 'MetaData/countedCodes'
    code_fs = sorted(os.listdir(code_dir))
    code_ps = [os.path.join(os.getcwd(), code_dir, i ) for i in code_fs]
    
               
               
    source_files = ['100714_B/100702_MystikLake_amaterna_Alm_L1_1_sequence.txt',
                    '121114Alm_dir/121114Alm_D12-4491_1_sequence.fastq',
                    '130719Alm_dir/130719Alm_D13-3437_1_sequence.fastq',
                    '130823Alm_dir/130823Alm_D13-4258_1_sequence.fastq',
                    '131001Alm_dir/131001Alm_D13-4961_phiX_best_1.fastq',
                    '131011Alm_dir/131011Alm_D13-5267_phiX_best_1.fastq',
                    '131114Alm_dir/131114Alm_D13-6069_1_sequence.fastq',
                    '131126Alm_dir/131126Alm_D13-6337_1_sequence.fastq',
                    'BGI_092012_dir/newsplit.Sarah_ML.1']
                    
    source_dirs = [os.path.dirname(i) for i in source_files]
    #============================================================================#
    if verbose:
        print "Loading Barcode Counts & Verifying Integrity"               
    
    code_dicts = []
    for ps in code_ps:
        with open(ps, 'r') as f:
            code_dicts.append(pickle.load(f))
    
    sev_check = code_dicts[-2]
    sev_check_2 = code_dicts[-3]
    for k in sev_check.keys():
        assert sev_check[k] == sev_check_2[k]
    del code_dicts[-2]
    
    code_dict2 = {i:j for i, j in zip(source_dirs, code_dicts)}
                  
    #===============================================================================#
    if verbose:
        print "Loading Metadata Table and creating new placeholder columns"
    
    metadata_d = 'MetaData/Master_mapping_file_remastered_4'
    metadata_f = 'Master_remastered_mapping.txt-Table 1.csv'
    metadata_p = os.path.join(os.getcwd(), metadata_d, metadata_f)
    metadata_df = pd.read_csv(metadata_p)
    metadata_df['Coverage'] = np.full(metadata_df.shape[0], np.nan)
    metadata_df['TotalSeqs'] = np.full(metadata_df.shape[0], np.nan)
    metadata_df['BadBarcodes'] = np.full(metadata_df.shape[0], np.nan)
    metadata_df['GoodBarcodes'] = np.full(metadata_df.shape[0], np.nan)
    
    #===============================================================================#
    if verbose:
        print "Processing Group 1"
    
    run_1_key, run_1_codes = '100714_B', code_dict2['100714_B']
    sub_df = metadata_df[metadata_df['Sequencing Date'] == run_1_key]
    if verbose:
        print "{} samples detected for the given date".format(sub_df.shape[0])
    
    all_matched = []
    run_1_idx = list(sub_df.index)
    for grp, idx in zip(list(sub_df.BarcodeSequence), run_1_idx):
        code_choices = run_1_codes.keys()
        matches = [hit for hit in code_choices if hit[-7:] == grp[-7:]]
        all_matched += matches
        tally = np.array([run_1_codes[match] for match in matches])
        metadata_df.ix[idx, 'Coverage'] = tally.sum()
        
    all_tally = np.array([run_1_codes[match] for match in all_matched])
    all_codes = np.array(run_1_codes.values())
    total_cd = all_codes.sum()
    bad_codes = total_cd - all_tally.sum()
    good_codes = all_tally.sum()
    gd_pc = good_codes / float(total_cd)
    bad_pc = bad_codes / float(total_cd)
    metadata_df.ix[run_1_idx, 'TotalSeqs'] = total_cd
    metadata_df.ix[run_1_idx, 'GoodBarcodes'] = gd_pc
    metadata_df.ix[run_1_idx, 'BadBarcodes'] = bad_pc
    if verbose:
        print "{} good, {} bad, {} total".format(gd_pc, bad_pc, total_cd)
    
    #===============================================================================#
    if verbose:
        print "Processing Group 2"
    
    run_2_key, run_2_codes = '130719Alm', code_dict2['130719Alm_dir']
    sub_df = metadata_df[metadata_df['Sequencing Date'] == run_2_key]
    if verbose:
        print "{} samples detected for the given date".format(sub_df.shape[0])
    
    all_matched = []
    run_2_idx = list(sub_df.index)
    for grp, idx in zip(list(sub_df.BarcodeSequence), run_2_idx):
        if 'N' in grp:
            sys.exit("N detected")
        code_choices = run_2_codes.keys()
        matches = [hit for hit in code_choices if hit == grp]
        all_matched += matches
        tally = np.array([run_2_codes[match] for match in matches])
        metadata_df.ix[idx, 'Coverage'] = tally.sum()
        
    all_tally = np.array([run_2_codes[match] for match in all_matched])
    all_codes = np.array(run_2_codes.values())
    total_cd = all_codes.sum()
    bad_codes = total_cd - all_tally.sum()
    good_codes = all_tally.sum()
    gd_pc = good_codes / float(total_cd)
    bad_pc = bad_codes / float(total_cd)
    
    metadata_df.ix[run_2_idx, 'TotalSeqs'] = total_cd
    metadata_df.ix[run_2_idx, 'GoodBarcodes'] = gd_pc
    metadata_df.ix[run_2_idx, 'BadBarcodes'] = bad_pc
    if verbose:
        print "{} good, {} bad, {} total".format(gd_pc, bad_pc, total_cd)
    
    #===============================================================================#
    if verbose:
        print "Processing Group 3"
    
    run_3_key, run_3_codes = '131126Alm', code_dict2['131126Alm_dir']
    sub_df = metadata_df[metadata_df['Sequencing Date'] == run_3_key]
    if verbose:
        print "{} samples detected for the given date".format(sub_df.shape[0])
    
    all_matched = []
    run_3_idx = list(sub_df.index)
    for grp, idx in zip(list(sub_df.BarcodeSequence), run_3_idx):
        if 'N' in grp:
            sys.exit("N detected")
        code_choices = run_3_codes.keys()
        matches = [hit for hit in code_choices if hit == grp]
        all_matched += matches
        tally = np.array([run_3_codes[match] for match in matches])
        metadata_df.ix[idx, 'Coverage'] = tally.sum()
        
    all_tally = np.array([run_3_codes[match] for match in all_matched])
    all_codes = np.array(run_3_codes.values())
    total_cd = all_codes.sum()
    bad_codes = total_cd - all_tally.sum()
    good_codes = all_tally.sum()
    gd_pc = good_codes / float(total_cd)
    bad_pc = bad_codes / float(total_cd)
    
    metadata_df.ix[run_3_idx, 'TotalSeqs'] = total_cd
    metadata_df.ix[run_3_idx, 'GoodBarcodes'] = gd_pc
    metadata_df.ix[run_3_idx, 'BadBarcodes'] = bad_pc
    if verbose:
        print "{} good, {} bad, {} total".format(gd_pc, bad_pc, total_cd)
    
    #===============================================================================# 
    if verbose:
        print "Processing Group 4"
    
    run_4_key, run_4_codes = '131114Alm', code_dict2['131114Alm_dir']
    sub_df = metadata_df[metadata_df['Sequencing Date'] == run_4_key]
    if verbose:
        print "{} samples detected for the given date".format(sub_df.shape[0])
    
    all_matched = []
    run_4_idx = list(sub_df.index)
    for grp, idx in zip(list(sub_df.BarcodeSequence), run_4_idx):
        if 'N' in grp:
            sys.exit("N detected")
        code_choices = run_4_codes.keys()
        matches = [hit for hit in code_choices if hit == grp]
        all_matched += matches
        tally = np.array([run_4_codes[match] for match in matches])
        metadata_df.ix[idx, 'Coverage'] = tally.sum()
        
    all_tally = np.array([run_4_codes[match] for match in all_matched])
    all_codes = np.array(run_4_codes.values())
    total_cd = all_codes.sum()
    bad_codes = total_cd - all_tally.sum()
    good_codes = all_tally.sum()
    gd_pc = good_codes / float(total_cd)
    bad_pc = bad_codes / float(total_cd)
    
    metadata_df.ix[run_4_idx, 'TotalSeqs'] = total_cd
    metadata_df.ix[run_4_idx, 'GoodBarcodes'] = gd_pc
    metadata_df.ix[run_4_idx, 'BadBarcodes'] = bad_pc
    if verbose:
        print "{} good, {} bad, {} total".format(gd_pc, bad_pc, total_cd)
    
    #===============================================================================# 
    if verbose:
        print "Processing Group 5"
    
    run_5_key, run_5_codes = '121114Alm', code_dict2['121114Alm_dir']
    sub_df = metadata_df[metadata_df['Sequencing Date'] == run_5_key]
    if verbose:
        print "{} samples detected for the given date".format(sub_df.shape[0])
    
    all_matched = []
    run_5_idx = list(sub_df.index)
    for grp, idx in zip(list(sub_df.BarcodeSequence), run_5_idx):
        if 'N' in grp:
            sys.exit("N detected")
        code_choices = run_5_codes.keys()
        matches = [hit for hit in code_choices if hit == grp]
        all_matched += matches
        tally = np.array([run_5_codes[match] for match in matches])
        metadata_df.ix[idx, 'Coverage'] = tally.sum()
        
    all_tally = np.array([run_5_codes[match] for match in all_matched])
    all_codes = np.array(run_5_codes.values())
    total_cd = all_codes.sum()
    bad_codes = total_cd - all_tally.sum()
    good_codes = all_tally.sum()
    gd_pc = good_codes / float(total_cd)
    bad_pc = bad_codes / float(total_cd)
    
    metadata_df.ix[run_5_idx, 'TotalSeqs'] = total_cd
    metadata_df.ix[run_5_idx, 'GoodBarcodes'] = gd_pc
    metadata_df.ix[run_5_idx, 'BadBarcodes'] = bad_pc
    if verbose:
        print "{} good, {} bad, {} total".format(gd_pc, bad_pc, total_cd)
    
    #===============================================================================# 
    if verbose:
        print "Processing Group 5"
    
    run_5_key, run_5_codes = '121114Alm', code_dict2['121114Alm_dir']
    sub_df = metadata_df[metadata_df['Sequencing Date'] == run_5_key]
    if verbose:
        print "{} samples detected for the given date".format(sub_df.shape[0])
    
    all_matched = []
    run_5_idx = list(sub_df.index)
    for grp, idx in zip(list(sub_df.BarcodeSequence), run_5_idx):
        if 'N' in grp:
            sys.exit("N detected")
        code_choices = run_5_codes.keys()
        matches = [hit for hit in code_choices if hit == grp]
        all_matched += matches
        tally = np.array([run_5_codes[match] for match in matches])
        metadata_df.ix[idx, 'Coverage'] = tally.sum()
        
    all_tally = np.array([run_5_codes[match] for match in all_matched])
    all_codes = np.array(run_5_codes.values())
    total_cd = all_codes.sum()
    bad_codes = total_cd - all_tally.sum()
    good_codes = all_tally.sum()
    gd_pc = good_codes / float(total_cd)
    bad_pc = bad_codes / float(total_cd)
    
    metadata_df.ix[run_5_idx, 'TotalSeqs'] = total_cd
    metadata_df.ix[run_5_idx, 'GoodBarcodes'] = gd_pc
    metadata_df.ix[run_5_idx, 'BadBarcodes'] = bad_pc
    if verbose:
        print "{} good, {} bad, {} total".format(gd_pc, bad_pc, total_cd)
    
    #===============================================================================# 
    if verbose:
        print "Processing Group 6"
    
    run_6_key, run_6_codes = 'BGI_092012', code_dict2['BGI_092012_dir']
    sub_df = metadata_df[metadata_df['Sequencing Date'] == run_6_key]
    if verbose:
        print "{} samples detected for the given date".format(sub_df.shape[0])
    
    all_matched = []
    run_6_idx = list(sub_df.index)
    for grp, idx in zip(list(sub_df.BarcodeSequence), run_6_idx):
        if 'N' in grp:
            sys.exit("N detected")
        code_choices = run_6_codes.keys()
        matches = [hit for hit in code_choices if hit == grp]
        all_matched += matches
        tally = np.array([run_6_codes[match] for match in matches])
        metadata_df.ix[idx, 'Coverage'] = tally.sum()
        
    all_tally = np.array([run_6_codes[match] for match in all_matched])
    all_codes = np.array(run_6_codes.values())
    total_cd = all_codes.sum()
    bad_codes = total_cd - all_tally.sum()
    good_codes = all_tally.sum()
    gd_pc = good_codes / float(total_cd)
    bad_pc = bad_codes / float(total_cd)
    
    metadata_df.ix[run_6_idx, 'TotalSeqs'] = total_cd
    metadata_df.ix[run_6_idx, 'GoodBarcodes'] = gd_pc
    metadata_df.ix[run_6_idx, 'BadBarcodes'] = bad_pc
    if verbose:
        print "{} good, {} bad, {} total".format(gd_pc, bad_pc, total_cd)
    
    
    #===============================================================================# 
    if verbose:
        print "Processing Group 7"
    
    run_7_key, run_7_codes = '131011Alm', code_dict2['131011Alm_dir']
    sub_df = metadata_df[metadata_df['Sequencing Date'] == run_7_key]
    if verbose:
        print "{} samples detected for the given date".format(sub_df.shape[0])
    
    all_matched = []
    run_7_idx = list(sub_df.index)
    for grp, idx in zip(list(sub_df.BarcodeSequence), run_7_idx):
        if 'N' in grp:
            sys.exit("N detected")
        code_choices = run_7_codes.keys()
        matches = [hit for hit in code_choices if hit == grp]
        all_matched += matches
        tally = np.array([run_7_codes[match] for match in matches])
        metadata_df.ix[idx, 'Coverage'] = tally.sum()
        
    all_tally = np.array([run_7_codes[match] for match in all_matched])
    all_codes = np.array(run_7_codes.values())
    total_cd = all_codes.sum()
    bad_codes = total_cd - all_tally.sum()
    good_codes = all_tally.sum()
    gd_pc = good_codes / float(total_cd)
    bad_pc = bad_codes / float(total_cd)
    
    metadata_df.ix[run_7_idx, 'TotalSeqs'] = total_cd
    metadata_df.ix[run_7_idx, 'GoodBarcodes'] = gd_pc
    metadata_df.ix[run_7_idx, 'BadBarcodes'] = bad_pc
    if verbose:
        print "{} good, {} bad, {} total".format(gd_pc, bad_pc, total_cd)
    
    #===============================================================================# 
    if verbose:
        print "Processing Group 8"
    
    run_8_key, run_8_codes = '131001Alm', code_dict2['131001Alm_dir']
    sub_df = metadata_df[metadata_df['Sequencing Date'] == run_8_key]
    if verbose:
        print "{} samples detected for the given date".format(sub_df.shape[0])
    
    all_matched = []
    run_8_idx = list(sub_df.index)
    for grp, idx in zip(list(sub_df.BarcodeSequence), run_8_idx):
        if 'N' in grp:
            sys.exit("N detected")
        code_choices = run_8_codes.keys()
        matches = [hit for hit in code_choices if hit == grp]
        all_matched += matches
        tally = np.array([run_8_codes[match] for match in matches])
        metadata_df.ix[idx, 'Coverage'] = tally.sum()
        
    all_tally = np.array([run_8_codes[match] for match in all_matched])
    all_codes = np.array(run_8_codes.values())
    total_cd = all_codes.sum()
    bad_codes = total_cd - all_tally.sum()
    good_codes = all_tally.sum()
    gd_pc = good_codes / float(total_cd)
    bad_pc = bad_codes / float(total_cd)
    
    metadata_df.ix[run_8_idx, 'TotalSeqs'] = total_cd
    metadata_df.ix[run_8_idx, 'GoodBarcodes'] = gd_pc
    metadata_df.ix[run_8_idx, 'BadBarcodes'] = bad_pc
    if verbose:
        print "{} good, {} bad, {} total".format(gd_pc, bad_pc, total_cd)
    
    
    #===============================================================================# 
    if verbose:
        print "Processing Group 9"
    
    run_9_key, run_9_codes = '130823Alm', code_dict2['130823Alm_dir']
    sub_df = metadata_df[metadata_df['Sequencing Date'] == run_9_key]
    if verbose:
        print "{} samples detected for the given date".format(sub_df.shape[0])
    
    all_matched = []
    run_9_idx = list(sub_df.index)
    for grp, idx in zip(list(sub_df.BarcodeSequence), run_9_idx):
        if 'N' in grp:
            sys.exit("N detected")
        code_choices = run_9_codes.keys()
        matches = [hit for hit in code_choices if hit == grp]
        all_matched += matches
        tally = np.array([run_9_codes[match] for match in matches])
        metadata_df.ix[idx, 'Coverage'] = tally.sum()
        
    all_tally = np.array([run_9_codes[match] for match in all_matched])
    all_codes = np.array(run_9_codes.values())
    total_cd = all_codes.sum()
    bad_codes = total_cd - all_tally.sum()
    good_codes = all_tally.sum()
    gd_pc = good_codes / float(total_cd)
    bad_pc = bad_codes / float(total_cd)
    
    metadata_df.ix[run_9_idx, 'TotalSeqs'] = total_cd
    metadata_df.ix[run_9_idx, 'GoodBarcodes'] = gd_pc
    metadata_df.ix[run_9_idx, 'BadBarcodes'] = bad_pc
    if verbose:
        print "{} good, {} bad, {} total".format(gd_pc, bad_pc, total_cd)
    
    
    to_drop = ['BarcodeSequence', 'LinkerPrimerSequence', 'BioSampleID', 'Replicate', 
               'DateMMDDYY', 'Date', 'DepthName', 'Depth', 'Sample Type', 
               '16S region', 'Forward  16S primer', 'Reverse 16S primer', 
               'Treatment', 'Sequence notes', 'Reverse read length', 'Description', 
               'OldSampleID', 'Notes']
               
               
    metadata_df.drop(to_drop, axis=1, inplace=True)
    
    return metadata_df
