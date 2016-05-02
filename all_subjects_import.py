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
from Tools.Sessions.AlbanPRF import *
from Tools.Sessions.AlbanPRF import WeibullPopulationReceptiveFieldMappingSession

# -----------------
# Comments / to do
# alle physlog kunnen weg, gaan we niets mee doen. (nog steeds doen? volgens mij hebben alle runs deze data wel, tot nu toe)
# raw behavior zat niet in alle runs; wel in de goede. Toch weg doen?(idem, zie boven?)
# check dat ook de .txt met logfiles meegenomen worden uit de raw mappen. *!!!!!!!!!!!!!!!!!!!!!!!!!!!* 
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

		## SETUP FILES: 
		#session.setupFiles(rawBase = presentSubject.initials, process_eyelink_file = False)
		#session.weibull_prf_trialtimes(task_conditions=['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10'])

		## WE'LL FIRST MOTION CORRECT THE EPIS TO THEMSELVES
		# session.motionCorrectFunctionals(use_ref_file=False)
		# session.create_moco_check_gifs()

		# ## NOW, LET'S FLIRT ALL MEAN MOTION CORRECTED VOLUMES TO THE SESSIONS T2 
		# # SO THAT WE CAN VISUALLY CHECK WHICH ONE HAS THE LEAST B0 DISTORTION
		# session.flirt_mean_moco_to_session_T2()
		# session.create_B0_distortion_check_gifs()

		# ## WITH A TARGET EPI VOLUME SELECTED AND MARKED IN THE SUBJECT DEFINITION
		# # WE CAN NOW REGISTER TO THE FREESURFER T1
		# session.registerSession(input_type='target_meanvol_moco_epi', MNI = False)
		
		"""
		oppassen; dit is een poging zonder MNI. de originele staat hieronder.
		MNI is uitgezet, omdat ik geen rechten heb om te schrijven voor deze folder.
		session.registerSession(input_type='target_meanvol_moco_epi')

		"""

		# ## WITH A TARGET EPI VOLUME SELECTED AND MARKED IN THE SUBJECT DEFINITION
		# # WE CAN NOW FLIRT ALL MEAN MOTION CORRECTED EPIS TO THAT EPI
		# # AND CREATE VISUAL SANITY CHECKS
		
		# session.flirt_mean_moco_to_mean_target_EPI()
		# session.check_EPI_alignment(postFix=['mcf','meanvol','NB','flirted2targetEPI'])

		# # ## FOR THE FINAL TOUCH, WE'LL NOW FNIRT THE MEAN MOTION CORRECTED AND FLIRTED
		# # EPI TO THE TARGET MEAN MOTION CORRECTED EPI

		# session.fnirt_mean_moco_to_mean_target_EPI()
		# session.check_EPI_alignment(postFix=['mcf','meanvol','fnirted2targetEPI'])

		# # NOW COMBINE MOCO / FLIRT / FNIRT AND APPLY TO ALL DATA
		# session.applywarp_to_moco_data()
		# session.create_mean_vol(postFix=['mcf','fnirted'])
		# session.check_EPI_alignment(postFix=['mcf','fnirted','meanvol'])

		## MASKS
		# session.dilate_and_move_func_bet_mask()
		# session.createMasksFromFreeSurferLabels(annot = False, annotFile = 'aparc.a2009s', labelFolders = ['retmap_PRF'], cortex = False)


		## SGTF
	 	# for condition in ['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10']:
			# session.rescaleFunctionals(condition=condition,operations = ['sgtf'],filterFreqs={'highpass':120}, funcPostFix = ['mcf','fnirted'], mask_file = os.path.join(session.stageFolder('processed/mri/masks/anat'), 'bet_mask_dilated.nii.gz'))
			# session.rescaleFunctionals(condition=condition,operations = ['percentsignalchange'], funcPostFix = ['mcf','fnirted','sgtf'])


		## Design Matrices

		# session.design_matrices_for_concatenated_data(n_pixel_elements_raw = 101,n_pixel_elements_convolved=31,
		# 					task_conditions=['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10'])
		# session.design_matrices_for_concatenated_data(n_pixel_elements_raw = 101,n_pixel_elements_convolved=31,
		# 					task_conditions=['PRF_01'])

		# session.design_matrices_for_averaged_data()


		
		## SETUP FIT PARAMETERS:
		# for condition in ['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10']:
		# moet dit wel in een for loop?
		# 	#task_conditions = ['PRF']
		# 	task_conditions = condition

		"""
		time to start fitting.
		de retmap / PRF folder in labels bestaat nog niet voor NM_290314 en MK_170714_4
		Deze moeten dus gemaakt worden. time to get this code running!

		n_jobs = 20 max.
		"""


		n_jobs = 20
		mask = 'lh.V1' #or any other mask here. visible in FSL
		postFix = ['mcf','fnirted','sgtf','psc']
		model = 'OG'# OG or DoG
		hrf_type = 'canonical' #'median'
		#slice_no was eerst 0; dan kreeg ik not iterable.
		slice_no = 12

		session.setup_fit_PRF_on_concatenated_data(
			mask_file_name = mask, 
			n_jobs = n_jobs, 
			postFix = postFix, 
			plotbool = True,
			model = model,
			hrf_type = hrf_type,
			fit_on_all_data = True,
			slice_no = slice_no
			)

		# session.combine_seperate_slice_niftis(mask,postFix,model,task_conditions=['All'],hrf_type=hrf_type)
		# session.convert_to_surf(mask_file = mask,postFix=postFix,model=model,hrf_type=hrf_type,depth_min=-1.0,depth_max=2.0,depth_step=0.25,task_conditions=['Fix'],sms=[0])
		# session.combine_surfaces(mask_file = mask,postFix=postFix,model=model,hrf_type=hrf_type,depth_min=-1.0,depth_max=2.0,depth_step=0.25,task_conditions=['Fix'],sms=[0])

		# session.design_matrices_for_concatenated_data(n_pixel_elements_raw = 101,n_pixel_elements_convolved=31,
		# 	change_type=change_type,run_num=run_num,task_conditions=['Fix','Color','Speed'])
		# r_squared_threshold = 0.0005 # threshold for which voxels to fit in the non-ALL condition
		# session.setup_fit_PRF_on_concatenated_data(
		# 	mask_file_name = mask, 
		# 	n_jobs = n_jobs, 
		# 	postFix = postFix, 
		# 	plotbool = True,
		# 	model = model,
		# 	hrf_type = 'median',
		# 	fit_on_all_data = False,
		# 	r_squared_threshold = r_squared_threshold,
		# 	slice_no = slice_no,
		# 	change_type = change_type,
		# 	run_num = run_num,
		# 	)	

		# session.combine_seperate_slice_niftis(mask,postFix,model,task_conditions = ['Fix','Color','Speed'],hrf_type=hrf_type)

		# CV:
		# session.combine_seperate_slice_niftis(mask,postFix,model,task_conditions = ['Fix','Color','Speed'],hrf_type=hrf_type,change_type='leave_one_out')
		# session.cross_predictions_concatenated_data(mask,postFix,model,hrf_type,n_jobs,run_num)
		# session.combine_predictions_concatenated_data(hrf_type,model,postFix,mask)


		# session.compare_fits(mask)
		
		# task_conditions = ['Stim']#['All','Fix','Color','Speed']
		# session.mask_stats_to_hdf(mask_file = mask , postFix = postFix, task_conditions = task_conditions,model=model,hrf_type=hrf_type)
		# session.fit_diagnostics(task_conditions = task_conditions, maskfile=mask, ecc_thresh=[0.0,7.0],model=model, 
		# 						hists = False, 
		# 						eccen_surf = True, 
		# 						r_squared_threshold = 0.2,
		# 						condition_comparison_plots = False, 
		# 						correlation_of_ecc_fwhm_diff_plot = False,
		# 						polar_plot = False,
		# 						simple_size_plot = False,
		# 						polar_imshow=False,
		# 						fba_attent_corr_plots=False,
		# 						postFix=postFix)

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
		#subject_session = PopulationReceptiveFieldMappingSession(sessionID, sessionDate, 
		#											   presentProject, presentSubject, this_project_folder, targetEPIID=3)
		subject_session = WeibullPopulationReceptiveFieldMappingSession(sessionID, sessionDate, 
		 											   presentProject, presentSubject, this_project_folder, targetEPIID=3)

		
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
				'ID': 1, 'scanType': 'epi_bold', 'condition': 'PRF_06', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF01_SENSE_4_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_01_cat_06_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_01_cat_06_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF01.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE01.edf'), 
				'thisSessionT2ID':12,				
			},
			{
				'ID': 2, 'scanType': 'epi_bold', 'condition': 'PRF_03', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF02_SENSE_5_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_02_cat_03_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag', 												
												 'DvE_240714_run_02_cat_03_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF02.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE02.edf'), 
				'thisSessionT2ID':12,				
			},
			{
				'ID': 3, 'scanType': 'epi_bold', 'condition': 'PRF_07', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF03_SENSE_6_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_03_cat_07_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_03_cat_07_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF03.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE03.edf'), 
				'thisSessionT2ID':12,				
			},
			{
				'ID': 4, 'scanType': 'epi_bold', 'condition': 'PRF_08', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF04_SENSE_7_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_04_cat_08_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_04_cat_08_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF04.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE04.edf'), 
				'thisSessionT2ID':12,		
			},
			{
				'ID': 5, 'scanType': 'epi_bold', 'condition': 'PRF_05', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF05_SENSE_9_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_05_cat_05_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_05_cat_05_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF05.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE05.edf'), 
				'thisSessionT2ID':12,				
			},
			{
				'ID': 6, 'scanType': 'epi_bold', 'condition': 'PRF_01', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF06_SENSE_10_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_06_cat_01_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_06_cat_01_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF06.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE06.edf'), 
				'thisSessionT2ID':12,				
			},
			{
				'ID': 7, 'scanType': 'epi_bold', 'condition': 'PRF_02', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF07_SENSE_11_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_07_cat_02_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_07_cat_02_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF07.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE07.edf'), 
				'thisSessionT2ID':12,				
			},
			{
				'ID': 8, 'scanType': 'epi_bold', 'condition': 'PRF_04', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF08_SENSE_13_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_08_cat_04_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_08_cat_04_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF08.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE08.edf'), 
				'thisSessionT2ID':12,				
			},
			{
				'ID': 9, 'scanType': 'epi_bold', 'condition': 'PRF_09', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF09_SENSE_14_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_09_cat_09_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_09_cat_09_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF09.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE09.edf'), 
				'thisSessionT2ID':12,				
			},
			{
				'ID': 10, 'scanType': 'epi_bold', 'condition': 'PRF_10', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_pRF10_SENSE_15_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_10_cat_10_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'DvE_240714_run_10_cat_10_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_pRF10.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'DvE10.edf'), 
				'thisSessionT2ID':12,				
			},
			{
				'ID': 11, 'scanType': 'T1', 'condition': 'T1', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'DvE_1_WIP_sT13DTFE_P25_S2_3m_SENSE_12_1.nii.gz'),
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'DvE_1_T1.phy')
			},
			{
				'ID': 12, 'scanType': 'inplane_anat', 'condition': 'T2_anat', 'session': 1,
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
		
		subject_session = WeibullPopulationReceptiveFieldMappingSession(sessionID, sessionDate, 
													   presentProject, presentSubject, this_project_folder, targetEPIID=3)
		
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
				'ID': 1, 'scanType': 'epi_bold', 'condition': 'PRF_04', 'session': 1,
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
				'ID': 2, 'scanType': 'epi_bold', 'condition': 'PRF_01', 'session': 1,
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
				'ID': 3, 'scanType': 'epi_bold', 'condition': 'PRF_02', 'session': 1,
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
				'ID': 4, 'scanType': 'epi_bold', 'condition': 'PRF_09', 'session': 1,
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
				'ID': 5, 'scanType': 'epi_bold', 'condition': 'PRF_03', 'session': 1,
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
				'ID': 6, 'scanType': 'epi_bold', 'condition': 'PRF_08', 'session': 1,
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
				'ID': 7, 'scanType': 'epi_bold', 'condition': 'PRF_07', 'session': 1,
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
				'ID': 8, 'scanType': 'epi_bold', 'condition': 'PRF_10', 'session': 1,
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
				'ID': 9, 'scanType': 'epi_bold', 'condition': 'PRF_05', 'session': 1,
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
				'ID': 10, 'scanType': 'epi_bold', 'condition': 'PRF_06', 'session': 1,
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
				'ID': 11, 'scanType': 'T1', 'condition': 'T1', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'prf_jwg_WIP_sT13DTFE_P25_S2_3m_SENSE_14_1.nii.gz')
			},
			{
				'ID': 12, 'scanType': 'inplane_anat', 'condition': 'T2_anat', 'session': 1,
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
		
		subject_session = WeibullPopulationReceptiveFieldMappingSession(sessionID, sessionDate, 
													   presentProject, presentSubject, this_project_folder, targetEPIID=3)
		
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
				'ID': 1, 'scanType': 'epi_bold', 'condition': 'PRF_07', 'session': 1,
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
				'ID': 2, 'scanType': 'epi_bold', 'condition': 'PRF_03', 'session': 1,
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
				'ID': 3, 'scanType': 'epi_bold', 'condition': 'PRF_10', 'session': 1,
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
				'ID': 4, 'scanType': 'epi_bold', 'condition': 'PRF_02', 'session': 1,
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
				'ID': 5, 'scanType': 'epi_bold', 'condition': 'PRF_01', 'session': 1,
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
				'ID': 6, 'scanType': 'epi_bold', 'condition': 'PRF_04', 'session': 1,
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
				'ID': 7, 'scanType': 'epi_bold', 'condition': 'PRF_09', 'session': 1,
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
				'ID': 8, 'scanType': 'epi_bold', 'condition': 'PRF_05', 'session': 1,
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
				'ID': 9, 'scanType': 'epi_bold', 'condition': 'PRF_08', 'session': 1,
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
				'ID': 10, 'scanType': 'epi_bold', 'condition': 'PRF_06', 'session': 1,
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
				'ID': 11, 'scanType': 'T1', 'condition': 'T1', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_sT13DTFE_P25_S2_3m_SENSE_13_1.nii.gz')
			},
			{
				'ID': 12, 'scanType': 'T1', 'condition': 'T1', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_sT13DTFE_P25_S2_3m_SENSE_17_1.nii.gz')
			},
			{
				'ID': 13, 'scanType': 'T1', 'condition': 'T1', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'MK_2_WIP_sT13DTFE_P25_S2_3m_SENSE_18_1.nii.gz')
			},
			# weird stuff is going on here. Why are there 3 T1 scans?
			# There are three files called _13,, _17 en 18 in the folder.
			# All 3 look kinda ok to my (very untrained) eye.
			# just pick one, delete the rest?
			{
				'ID': 14, 'scanType': 'inplane_anat', 'condition': 'T2_anat', 'session': 1,
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
		
		subject_session = WeibullPopulationReceptiveFieldMappingSession(sessionID, sessionDate, 
													   presentProject, presentSubject, this_project_folder, targetEPIID=3)
		
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
				'ID': 1, 'scanType': 'epi_bold', 'condition': 'PRF_03', 'session': 1,
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
				'ID': 2, 'scanType': 'epi_bold', 'condition': 'PRF_09', 'session': 1,
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
				'ID': 3, 'scanType': 'epi_bold', 'condition': 'PRF_10', 'session': 1,
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
				'ID': 4, 'scanType': 'epi_bold', 'condition': 'PRF_02', 'session': 1,
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
				'ID': 5, 'scanType': 'epi_bold', 'condition': 'PRF_05', 'session': 1,
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
				'ID': 6, 'scanType': 'epi_bold', 'condition': 'PRF_08', 'session': 1,
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
				'ID': 7, 'scanType': 'epi_bold', 'condition': 'PRF_07', 'session': 1,
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
				'ID': 8, 'scanType': 'epi_bold', 'condition': 'PRF_01', 'session': 1,
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
				'ID': 9, 'scanType': 'epi_bold', 'condition': 'PRF_04', 'session': 1,
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
				'ID': 10, 'scanType': 'epi_bold', 'condition': 'PRF_06', 'session': 1,
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
				'ID': 11, 'scanType': 'T1', 'condition': 'T1', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'NM_WIP_sT13DTFE_P25_S2_3m_SENSE_13_1.nii.gz')
			},
			{
				'ID': 12, 'scanType': 'inplane_anat', 'condition': 'T2_anat', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'NM_WIP_T2W_RetMap_1.25_CLEAR_9_1.nii.gz')
			}
		]
		runWholeSession(subject_run_array, subject_session)

		
	# subjects.append(presentSubject)
	# run_arrays.append(runWBPRFArray)
	# projects.append(presentProject)
	# session_dates.append(sessionDate)