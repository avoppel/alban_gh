#!/usr/bin/env python
# encoding: utf-8
"""
analyze_7T_S1.py

Created by Tomas HJ Knapen on 2009-11-26.
Copyright (c) 2009 TK. All rights reserved.

To run this script headlessly and hopefully avoiding X server errors in the process, execute this script using 'xvfb-run python subjects.py'
"""

import os, sys, datetime
import subprocess, logging

import scipy as sp
import scipy.stats as stats
import numpy as np
import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pylab as pl

from IPython import embed as shell

this_project_folder = '/home/shared/Weibull_PRF/'
this_raw_folder = '/home/raw_data/Weibull_PRF/'

sys.path.append( os.environ['ANALYSIS_HOME'] )
	
sys.path.append( os.path.join('/home', 'voppel', 'PRF_2_analysis') )


from PopulationReceptiveFieldMappingSession import *


from Tools.Sessions import *
from Tools.Subjects.Subject import *
from Tools.Run import *
from Tools.Projects.Project import *
# -----------------
# Comments / to do
# alle physlog kunnen weg, gaan we niets mee doen. (nog steeds doen? volgens mij hebben alle runs deze data wel, tot nu toe)
# raw behavior zat niet in alle runs; wel in de goede. Toch weg doen?(idem, zie boven?)
# check dat ook de .txt met logfiles meegenomen worden uit de raw mappen.
# -----------------

# subject_initials = ['DvE', 'JWdG', 'MK', 'NM']
subject_initials = ['DvE']
# subjects = ['DvE', 'JWdG', 'MK', 'NM']
subjects = ['DvE']
run_arrays = []
projects = []
session_dates = []

for which_subject in subject_initials:
	

	def runWholeSession( rA, session ):
		for ri,r in enumerate(rA):
			thisRun = Run( **r )
			session.addRun(thisRun)
		session.parcelateConditions()
		session.parallelize = True

		# ----------------------
		# Functions to execute -
		# ----------------------

		#session.setupFiles(rawBase = presentSubject.initials, process_eyelink_file = False)
		# session.motionCorrectFunctionals(use_ref_file=False)
		# ^done for DvE
		session.create_moco_check_gifs()

	# ----------------------
	# Initialise session   -
	# ----------------------
	
	if which_subject == 'DvE':
		
		# ----------------------
		# Subject information  -
		# ----------------------
		
		initials = 'DvE'
		firstName = 'DvE'
		standardFSID = 'DE_110412'						# look up
		birthdate = datetime.date(1989,07,14)					# look up, yyyymmdd
		labelFolderOfPreference = 'visual areas'
		presentSubject = Subject(initials, firstName, birthdate, 
								 standardFSID, labelFolderOfPreference)
		presentProject = Project('Weibull_PRF', subject = presentSubject,
								 base_dir = os.path.join(this_project_folder, 'data'))
		sessionID = 'Weibull_PRF' + presentSubject.initials
		sessionDate = datetime.date(2014, 07, 24)
		sj_session = 'DvE_240714'
		
		subject_session = PopulationReceptiveFieldMappingSession(sessionID, sessionDate, 
													   presentProject, presentSubject, this_project_folder)
		
		try:
			os.mkdir(os.path.join(this_project_folder, 'data', initials))
		except OSError:
			subject_session.logger.debug('output folders already exist')
		
		# ----------------------
		# Run array			   -
		# ----------------------
		
		subject_run_array = [
			# 10 runs met de PRF taak, 1 run T2, 1 run T1.
			{
				'ID': 1, 'scanType': 'epi_bold', 'condition': '06', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF01_SENSE_4_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_01_cat_06_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_01_cat_06_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF01.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE01.edf')				
			},
			{
				'ID': 2, 'scanType': 'epi_bold', 'condition': '03', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF02_SENSE_5_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_02_cat_03_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag', 												
												 'DvE_240714_run_02_cat_03_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF02.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE02.edf')				
			},
			{
				'ID': 3, 'scanType': 'epi_bold', 'condition': '07', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF03_SENSE_6_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_03_cat_07_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_03_cat_07_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF03.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE03.edf')				
			},
			{
				'ID': 4, 'scanType': 'epi_bold', 'condition': '08', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF04_SENSE_7_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_04_cat_08_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_04_cat_08_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF04.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE04.edf')		
			},
			{
				'ID': 5, 'scanType': 'epi_bold', 'condition': '05', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF05_SENSE_9_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_05_cat_05_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_05_cat_05_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF05.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE05.edf')				
			},
			{
				'ID': 6, 'scanType': 'epi_bold', 'condition': '01', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF06_SENSE_10_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_06_cat_01_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_06_cat_01_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF06.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE06.edf')				
			},
			{
				'ID': 7, 'scanType': 'epi_bold', 'condition': '02', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF07_SENSE_11_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_07_cat_02_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_07_cat_02_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF07.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE07.edf')				
			},
			{
				'ID': 8, 'scanType': 'epi_bold', 'condition': '04', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF08_SENSE_13_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_08_cat_04_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_08_cat_04_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF08.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE08.edf')				
			},
			{
				'ID': 9, 'scanType': 'epi_bold', 'condition': '09', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF09_SENSE_14_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_09_cat_09_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_09_cat_09_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF09.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE09.edf')				
			},
			{
				'ID': 10, 'scanType': 'epi_bold', 'condition': '10', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF10_SENSE_15_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_10_cat_10_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_10_cat_10_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF10.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE10.edf')				
			},
			{
				'ID': 11, 'scanType': 'T1', 'condition': 'mapper', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_sT13DTFE_P25_S2_3m_SENSE_12_1.nii.gz'),
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_T1.phy')
			},
			{
				'ID': 12, 'scanType': 'T2', 'condition': 'mapper', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_T2W_RetMap_1.25_CLEAR_8_1.nii.gz'),
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_T2.phy')
			}
		]
		
		runWholeSession(subject_run_array, subject_session)
	
	# ----------------------
	# Initialise session   -
	# ----------------------
	
	elif which_subject == 'JWdG':
		
		# ----------------------
		# Subject information  -
		# ----------------------

		
		initials = 'JWdG'
		firstName = 'JWdG'
		standardFSID = 'JW_310312'						# look up
		birthdate = datetime.date(1900,01,01)					# look up, yyyymmdd
		labelFolderOfPreference = 'visual areas'
		presentSubject = Subject(initials, firstName, birthdate, 
								 standardFSID, labelFolderOfPreference)
		presentProject = Project('Weibull_PRF', subject = presentSubject,
								 base_dir = os.path.join(this_project_folder, 'data'))
		sessionID = 'Weibull_PRF' + presentSubject.initials
		sessionDate = datetime.date(2014, 07, 14)
		# sessionDate = datetime.date(2014, 08, 05)
		sj_session = 'JWdG_140714'
		
		subject_session = PopulationReceptiveFieldMappingSession(sessionID, sessionDate, 
													   presentProject, presentSubject, this_project_folder)
		
		try:
			os.mkdir(os.path.join(this_project_folder, 'data', initials))
		except OSError:
			subject_session.logger.debug('output folders already exist')
		
		# ----------------------
		# Run array			   -
		# ----------------------

	
		subject_run_array = [
			# 10 runs met de PRF taak, 1 run T2, 1 run T1.
			{
				'ID': 1, 'scanType': 'epi_bold', 'condition': '04', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'prf_jwg_WIP_pRF01_SENSE_5_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_01_cat_04_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_01_cat_04_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'WBpRF_jwg_01.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'jwg1.edf')
			},
			{
				'ID': 2, 'scanType': 'epi_bold', 'condition': '01', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'prf_jwg_WIP_pRF02_SENSE_6_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_02_cat_01_stim.txt'), 
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag', 
												 'JWdG_140714_run_02_cat_01_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'WBpRF_jwg_02.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'jwg2.edf')
			},
			{
				'ID': 3, 'scanType': 'epi_bold', 'condition': '02', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'prf_jwg_WIP_pRF03_SENSE_7_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_03_cat_02_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_03_cat_02_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'WBpRF_jwg_03.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'jwg3.edf')
			},
			{
				'ID': 4, 'scanType': 'epi_bold', 'condition': '09', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'prf_jwg_WIP_pRF04_SENSE_9_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_04_cat_09_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_04_cat_09_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'WBpRF_jwg_04.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'jwg4.edf')
			},
			{
				'ID': 5, 'scanType': 'epi_bold', 'condition': '03', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'prf_jwg_WIP_pRF05_SENSE_11_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_05_cat_03_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_05_cat_03_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'WBpRF_jwg_05.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'jwg5.edf')
			},
			{
				'ID': 6, 'scanType': 'epi_bold', 'condition': '08', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'prf_jwg_WIP_pRF06_SENSE_12_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_06_cat_08_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_06_cat_08_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'WBpRF_jwg_06.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'jwg6.edf')
			},
			{
				'ID': 7, 'scanType': 'epi_bold', 'condition': '07', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'prf_jwg_WIP_pRF07_SENSE_13_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_07_cat_07_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_07_cat_07_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'WBpRF_jwg_07.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'jwg7.edf')
			},
			{
				'ID': 8, 'scanType': 'epi_bold', 'condition': '10', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'prf_jwg_WIP_pRF08_SENSE_15_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_08_cat_10_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_08_cat_10_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'WBpRF_jwg_08.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'jwg8.edf')
			},
			{
				'ID': 9, 'scanType': 'epi_bold', 'condition': '05', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'prf_jwg_WIP_pRF09_SENSE_17_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_09_cat_05_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_09_cat_05_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'WBpRF_jwg_09.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'jwg9.edf')
			},
			{
				'ID': 10, 'scanType': 'epi_bold', 'condition': '06', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'prf_jwg_WIP_pRF10_SENSE_18_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_10_cat_06_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_140714_run_10_cat_06_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'WBpRF_jwg_10.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'jwg10.edf')
			},
			{
				'ID': 11, 'scanType': 'T1', 'condition': 'mapper', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'prf_jwg_WIP_sT13DTFE_P25_S2_3m_SENSE_14_1.nii.gz')
			},
			{
				'ID': 12, 'scanType': 'T2', 'condition': 'mapper', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'prf_jwg_WIP_T2W_RetMap_1.25_CLEAR_10_1.nii.gz')
			},
		]
		runWholeSession(subject_run_array, subject_session)


	
	elif which_subject == 'MK':
		
		# ----------------------
		# Subject information  -
		# ----------------------
		
		initials = 'MK'
		firstName = 'MK'
		standardFSID = 'MK_170714_4'
		birthdate = datetime.date(1990, 9, 24)
		labelFolderOfPreference = 'visual areas'
		presentSubject = Subject(initials, firstName, birthdate, 
								 standardFSID, labelFolderOfPreference)
		presentProject = Project('Weibull_PRF', subject = presentSubject,
								 base_dir = os.path.join(this_project_folder, 'data'))
		sessionID = 'Weibull_PRF' + presentSubject.initials
		sessionDate = datetime.date(2014, 07, 17)
		# sessionDate = datetime.date(2014, 08, 05)
		sj_session = 'MK_050814'
		
		subject_session = PopulationReceptiveFieldMappingSession(sessionID, sessionDate, 
													   presentProject, presentSubject, this_project_folder)
		
		try:
			os.mkdir(os.path.join(this_project_folder, 'data', initials))
		except OSError:
			subject_session.logger.debug('output folders already exist')
		
		# ----------------------
		# Run array			   -
		# ----------------------
		
		
		subject_run_array = [
			# 10 runs met de PRF taak, 1 run T2, 1 run T1.
			# Twee extra T1 runs.
			{
				'ID': 1, 'scanType': 'epi_bold', 'condition': '07', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_pRF01_SENSE_5_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'MK_050814_run_01_cat_07_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'MK_050814_run_01_cat_07_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'MK2_pRF01.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'MK2_1.edf')
			},
			{
				'ID': 2, 'scanType': 'epi_bold', 'condition': '03', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_pRF02_SENSE_6_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'MK_050814_run_02_cat_03_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'MK_050814_run_02_cat_03_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'MK2_pRF02.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'MK2_2.edf')
			},
			{
				'ID': 3, 'scanType': 'epi_bold', 'condition': '10', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_pRF03_SENSE_7_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'MK_050814_run_03_cat_10_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'MK_050814_run_03_cat_10_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'MK2_pRF03.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'MK2_3.edf')
			},
			{
				'ID': 4, 'scanType': 'epi_bold', 'condition': '02', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_pRF04_SENSE_8_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'MK_050814_run_04_cat_02_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'MK_050814_run_04_cat_02_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'MK2_pRF04.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'MK2_4.edf')
			},
			{
				'ID': 5, 'scanType': 'epi_bold', 'condition': '01', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_pRF05_SENSE_10_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'MK_050814_run_05_cat_01_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'MK_050814_run_05_cat_01_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'MK2_pRF05.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'MK2_5.edf')
			},
			{
				'ID': 6, 'scanType': 'epi_bold', 'condition': '04', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_pRF06_SENSE_5_11.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'MK_050814_run_06_cat_04_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'MK_050814_run_06_cat_04_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'MK2_pRF06.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'MK2_6.edf')
			},
			{
				'ID': 7, 'scanType': 'epi_bold', 'condition': '09', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_pRF07_SENSE_12_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'MK_050814_run_07_cat_09_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'MK_050814_run_07_cat_09_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'MK2_pRF07.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'MK2_7.edf')
			},
			{
				'ID': 8, 'scanType': 'epi_bold', 'condition': '05', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_pRF08_SENSE_14_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'MK_050814_run_08_cat_05_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'MK_050814_run_08_cat_05_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'MK2_pRF08.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'MK2_8.edf')
			},
			{
				'ID': 9, 'scanType': 'epi_bold', 'condition': '08', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_pRF09_SENSE_15_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'MK_050814_run_09_cat_08_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'MK_050814_run_09_cat_08_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'MK2_pRF09.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'MK2_9.edf')
			},
			{
				'ID': 10, 'scanType': 'epi_bold', 'condition': '06', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_pRF10_SENSE_16_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'MK_050814_run_10_cat_06_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'MK_050814_run_10_cat_06_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'MK2_pRF10.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'MK2_10.edf')
			},
			{
				'ID': 11, 'scanType': 'T1', 'condition': 'mapper', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_sT13DTFE_P25_S2_3m_SENSE_13_1.nii.gz')
			},
			{
				'ID': 12, 'scanType': 'T1', 'condition': 'mapper', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_sT13DTFE_P25_S2_3m_SENSE_17_1.nii.gz')
			},
			{
				'ID': 13, 'scanType': 'T1', 'condition': 'mapper', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_sT13DTFE_P25_S2_3m_SENSE_18_1.nii.gz')
			},
			# weird stuff is going on here. Why are there 3 T1 scans?
			# There are three files called _13,, _17 en 18 in the folder.
			# All 3 look kinda ok to my (very untrained) eye.
			# just pick one, delete the rest?
			{
				'ID': 14, 'scanType': 'T2', 'condition': 'mapper', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_T2W_RetMap_125_CLEAR_9_1.nii.gz')
			}
		]
		runWholeSession(subject_run_array, subject_session)

	
	elif which_subject == 'NM':
		
		# ----------------------
		# Subject information  -
		# ----------------------

		initials = 'NM'
		firstName = 'NM'
		standardFSID = 'NM_290314'
		birthdate = datetime.date(1990, 01, 29)
		labelFolderOfPreference = 'visual areas'
		presentSubject = Subject(initials, firstName, birthdate, 
								 standardFSID, labelFolderOfPreference)
		presentProject = Project('Weibull_PRF', subject = presentSubject,
								 base_dir = os.path.join(this_project_folder, 'data'))
		sessionID = 'Weibull_PRF' + presentSubject.initials
		sessionDate = datetime.date(2014, 07, 06)
		# sessionDate = datetime.date(2014, 07, 22)
		sj_session = 'NM_220714'
		
		subject_session = PopulationReceptiveFieldMappingSession(sessionID, sessionDate, 
													   presentProject, presentSubject, this_project_folder)
		
		try:
			os.mkdir(os.path.join(this_project_folder, 'data', initials))
		except OSError:
			subject_session.logger.debug('output folders already exist')
		
		# ----------------------
		# Run array			   -
		# ----------------------
		subject_run_array = [
			# 10 runs met de PRF taak, 1 run T2, 1 run T1.
			{
				'ID': 1, 'scanType': 'epi_bold', 'condition': '03', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'NM_WIP_pRF01_SENSE_5_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'NM_220714_run_01_cat_03_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'NM_220714_run_01_cat_03_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'NM_pRF01.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'NM01.edf')
			},
			{
				'ID': 2, 'scanType': 'epi_bold', 'condition': '09', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'NM_WIP_pRF02_SENSE_6_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'NM_220714_run_02_cat_09_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'NM_220714_run_02_cat_09_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'NM_pRF02.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'NM02.edf')
			},
			{
				'ID': 3, 'scanType': 'epi_bold', 'condition': '10', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'NM_WIP_pRF03_SENSE_7_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'NM_220714_run_03_cat_10_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'NM_220714_run_03_cat_10_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'NM_pRF03.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'NM03.edf')
			},
			{
				'ID': 4, 'scanType': 'epi_bold', 'condition': '02', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'NM_WIP_pRF04_SENSE_8_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'NM_220714_run_04_cat_02_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'NM_220714_run_04_cat_02_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'NM_pRF04.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'NM04.edf')
			},
			{
				'ID': 5, 'scanType': 'epi_bold', 'condition': '05', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'NM_WIP_pRF05_SENSE_10_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'NM_220714_run_05_cat_05_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'NM_220714_run_05_cat_05_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'NM_pRF05.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'NM05.edf')
			},
			{
				'ID': 6, 'scanType': 'epi_bold', 'condition': '08', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'NM_WIP_pRF06_SENSE_11_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'NM_220714_run_06_cat_08_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'NM_220714_run_06_cat_08_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'NM_pRF06.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'NM06.edf')
			},
			{
				'ID': 7, 'scanType': 'epi_bold', 'condition': '07', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'NM_WIP_pRF07_SENSE_12_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'NM_220714_run_07_cat_07_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'NM_220714_run_07_cat_07_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'NM_pRF07.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'NM07.edf')
			},
			{
				'ID': 8, 'scanType': 'epi_bold', 'condition': '01', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'NM_WIP_pRF08_SENSE_14_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'NM_220714_run_08_cat_01_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'NM_220714_run_08_cat_01_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'NM_pRF08.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'NM08.edf')
			},
			{
				'ID': 9, 'scanType': 'epi_bold', 'condition': '04', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'NM_WIP_pRF08_SENSE_15_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'NM_220714_run_09_cat_04_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'NM_220714_run_09_cat_04_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'NM_pRF09.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'NM09.edf')
			},
			{
				'ID': 10, 'scanType': 'epi_bold', 'condition': '06', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'NM_WIP_pRF10_SENSE_16_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'NM_220714_run_10_cat_06_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'NM_220714_run_10_cat_06_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'NM_pRF10.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'NM10.edf')
			},
			{
				'ID': 11, 'scanType': 'T1', 'condition': 'mapper', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'NM_WIP_sT13DTFE_P25_S2_3m_SENSE_13_1.nii.gz')
			},
			{
				'ID': 12, 'scanType': 'T2', 'condition': 'mapper', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'NM_WIP_T2W_RetMap_1.25_CLEAR_9_1.nii.gz')
			}
		]
		runWholeSession(subject_run_array, subject_session)

		
	# subjects.append(presentSubject)
	# run_arrays.append(runWBPRFArray)
	# projects.append(presentProject)
	# session_dates.append(sessionDate)