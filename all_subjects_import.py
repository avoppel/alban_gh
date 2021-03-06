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

sys.path.append( os.path.join('/home', 'voppel', 'anaconda2', 'lib', 'python2.7', 'site-packages') )

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
#from Tools.Sessions.AlbanPRF import *
#from Tools.Sessions.AlbanPRF import WeibullPopulationReceptiveFieldMappingSession
from AlbanPRF import *
from AlbanPRF import WeibullPopulationReceptiveFieldMappingSession
from ModelExperiment import *
import time


subject_initials = ['DvE', 'JWdG', 'MK', 'NM']
# subject_initials = ['NM']
subjects = ['DvE', 'JWdG', 'MK', 'NM']
# subjects = ['NM']
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

		# ## SETUP FILES: 
		# session.setupFiles(rawBase = presentSubject.initials, process_eyelink_file = False)
		# session.weibull_prf_trialtimes(task_conditions=['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10'])

		# # WE'LL FIRST MOTION CORRECT THE EPIS TO THEMSELVES
		# session.motionCorrectFunctionals(use_ref_file=False)
		# session.create_moco_check_gifs()

		# NOW, LET'S FLIRT ALL MEAN MOTION CORRECTED VOLUMES TO THE SESSIONS T2 
		# SO THAT WE CAN VISUALLY CHECK WHICH ONE HAS THE LEAST B0 DISTORTION
		# session.flirt_mean_moco_to_session_T2()
		# session.create_B0_distortion_check_gifs()

		# ## WITH A TARGET EPI VOLUME SELECTED AND MARKED IN THE SUBJECT DEFINITION
		# # WE CAN NOW REGISTER TO THE FREESURFER T1
		# session.registerSession(input_type='target_meanvol_moco_epi', MNI = False)
		# original version: #	session.registerSession(input_type='target_meanvol_moco_epi')

		# ## WITH A TARGET EPI VOLUME SELECTED AND MARKED IN THE SUBJECT DEFINITION
		# # WE CAN NOW FLIRT ALL MEAN MOTION CORRECTED EPIS TO THAT EPI
		# # # AND CREATE VISUAL SANITY CHECKS
		
		# session.flirt_mean_moco_to_mean_target_EPI()
		# session.check_EPI_alignment(postFix=['mcf','meanvol','NB','flirted2targetEPI'])

		# # # # ## FOR THE FINAL TOUCH, WE'LL NOW FNIRT THE MEAN MOTION CORRECTED AND FLIRTED
		# # # # EPI TO THE TARGET MEAN MOTION CORRECTED EPI

		# session.fnirt_mean_moco_to_mean_target_EPI()
		# session.check_EPI_alignment(postFix=['mcf','meanvol','fnirted2targetEPI'])

		# # # NOW COMBINE MOCO / FLIRT / FNIRT AND APPLY TO ALL DATA
		# session.applywarp_to_moco_data()
		# session.create_mean_vol(postFix=['mcf','fnirted'])
		# session.check_EPI_alignment(postFix=['mcf','fnirted','meanvol'])


		# # ## MASKS
		# session.dilate_and_move_func_bet_mask()
		# session.create_dilated_cortical_mask()
		# session.createMasksFromFreeSurferLabels(annot = False, annotFile = 'aparc.a2009s', labelFolders = ['retmap_PRF'], cortex = False)


		## SGTF
	 	#for condition in ['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10']:

			# session.rescaleFunctionals(condition=condition,operations = ['sgtf'],filterFreqs={'highpass':120}, funcPostFix = ['mcf','fnirted'], mask_file = os.path.join(session.stageFolder('processed/mri/masks/anat'), 'bet_mask_dilated.nii.gz'))
			# session.rescaleFunctionals(condition=condition,operations = ['percentsignalchange'], funcPostFix = ['mcf','fnirted','sgtf'])


		## Design Matrices

		# session.design_matrices_for_concatenated_data(n_pixel_elements_raw = 101,n_pixel_elements_convolved=31,
		# 					task_conditions=['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10'])

		# session.design_matrices_for_concatenated_data(n_pixel_elements_raw = 101,n_pixel_elements_convolved=31,
		#  					task_conditions=['PRF_04'])

		# session.design_matrices_for_averaged_data()


		## Create masks; see down below by other functions for more
		# session.combine_rois(rois=['lh.V1','rh.V1','lh.V2v','rh.V2d','lh.V2d','rh.V2v',		#V1 V2
		# 					'lh.V3v','rh.V3d','lh.V3d','rh.V3v','lh.V4','rh.V4',				#V3 V4
		# 					'lh.LO1','rh.LO1','lh.LO2','rh.LO2',								#LO
		# 					'lh.TO1','rh.TO1','lh.TO2','rh.TO2',								#MT
		# 					'lh.VO1','rh.VO1','lh.VO2','rh.VO2'									#VO
		# 					],output_roi = 'all_visual')

		# session.combine_rois(rois=['lh.V1','rh.V1',											#V1
		# 					'lh.LO1','rh.LO1','lh.LO2','rh.LO2',								#LO
		# 					],output_roi = 'V1_LO')

		## SETUP FIT PARAMETERS:

		"""
		n_jobs = 20 max.
		snachts - doe misschien -1 ?
		"""

		n_jobs = 20
		mask = 'V1_LO' #or any other mask here. check in FSL
		mask = 'all_visual' #or any other mask here. check in FSL
		postFix = ['mcf','fnirted','sgtf','psc']
		model = 'OG'# OG or DoG
		hrf_type = 'canonical' #'median'
		#None loopt alle slices.
		slice_no = None

		#sleep function; to start fit procedure after X seconds.
		#time.sleep(7200)
		
		#quick and dirty to do various conditions in a loop
		# for prfnumber in range(10):
			

		# session.setup_fit_PRF_on_concatenated_data(
		# 	mask_file_name = mask, 
		# 	n_jobs = n_jobs, 
		# 	postFix = postFix, 
		# 	plotbool = True,
		# 	model = model,
		# 	hrf_type = hrf_type,
		# 	fit_on_all_data = False,
		# 	slice_no = slice_no,
		# 	#condition_index = np.array([prfnumber]),
		# 	condition_index = np.array([2,4]),
		# 	# condition_index = np.arange(5,10)
		# 	##this one combines all conditions and fits on them all
		# 	#condition_index = np.arange(10),
		# 	)
		
		task_conditions = ['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10',
			'All_0_4','All_5_9', 'All_0_9', 'All_3_4', 'All_1_2', 'All_1_3', 'All_2_4']
		# # # # # # task_conditions = ['All_0_4','All_5_9']
		# # # # # # # #task_conditions = ['All_0_9']

		session.mask_stats_to_hdf(mask_file = mask , postFix = postFix, task_conditions = task_conditions,model=model,hrf_type=hrf_type)


		# session.combine_seperate_slice_niftis(mask,postFix,model,task_conditions=['All'],hrf_type=hrf_type)

		#session.convert_to_surf(mask_file = mask,postFix=postFix,model=model,hrf_type=hrf_type,depth_min=-1.0,depth_max=2.0,depth_step=0.25,task_conditions=['All_0_9'],sms=[0])
		#session.combine_surfaces(mask_file = mask,postFix=postFix,model=model,hrf_type=hrf_type,depth_min=-1.0,depth_max=2.0,depth_step=0.25,task_conditions=['All_0_9'],sms=[0])

		# session.labels_to_annot(input_folder_names=['retmap_PRF'],output_file_name = 'all_labels')
		
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

		## OPTIONAL OTHERS

		# session.combine_rois(rois=['lh.V1','rh.V1','lh.V2v','rh.V2d','lh.V2d','rh.V2v',		#V1 V2
		# 							'lh.V3v','rh.V3d','lh.V3d','rh.V3v','lh.V4','rh.V4',	#V3 V4
		# 							'lh.LO1','rh.LO1','lh.LO2','rh.LO2',					#LO
		# 							'lh.TO1','rh.TO1','lh.TO2','rh.TO2',					#MT
		# 							'lh.VO1','rh.VO1','lh.VO2','rh.VO2'						#VO
		# 							],output_roi = 'all_visual')

		# session.combine_rois(rois=['lh.V1','rh.V1'],output_roi = 'V1')
	 	# session.combine_rois(rois=['lh.V2v','rh.V2d','lh.V2d','rh.V2v'],output_roi = 'V2')
	 	# session.combine_rois(rois=['lh.V3v','rh.V3d','lh.V3d','rh.V3v'],output_roi = 'V3')
	 	# session.combine_rois(rois=['lh.V4','rh.V4'],output_roi = 'V4')
	 	# session.combine_rois(rois=['lh.LO1','rh.LO1','lh.LO2','rh.LO2'],output_roi = 'LO')
	 	# session.combine_rois(rois=['lh.TO1','rh.TO1','lh.TO2','rh.TO2'],output_roi = 'MT')
	 	# session.combine_rois(rois=['lh.V3A','rh.V3A','lh.V3B','rh.V3B'],output_roi = 'V3AB')
	 	# session.combine_rois(rois=['lh.VO1','rh.VO1','lh.VO2','rh.VO2'],output_roi = 'VO')
	 	# session.combine_rois(rois=['lh.PHC1','rh.PHC1','lh.PHC2','rh.PHC2'],output_roi = 'PHC')
	 	# session.combine_rois(rois=['lh.IPS0','rh.IPS0'],output_roi = 'IPS0')
	 	# session.combine_rois(rois=['lh.IPS1','rh.IPS1'],output_roi = 'IPS1')
	 	# session.combine_rois(rois=['lh.IPS2','rh.IPS2'],output_roi = 'IPS2')
	 	# session.combine_rois(rois=['lh.IPS3','rh.IPS3'],output_roi = 'IPS3')
	 	# session.combine_rois(rois=['lh.IPS4','rh.IPS4'],output_roi = 'IPS4')
	 	# session.combine_rois(rois=['lh.FEF','rh.FEF'],output_roi = 'FEF')
		# session.combine_rois(rois=['PHC','MT','IPS','FEF'],output_roi = 'late_visual')
		# session.inflate_T2s()
	 	# session.create_combined_label_mask()


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
		standardFSID = 'DE_110412'						
		birthdate = datetime.date(1989,07,14)					
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
		sessionDate = datetime.date(2014, 8, 05)
		# sessionDate = datetime.date(2014, 08, 05)
		sj_session = 'JWdG_050814'
		
		subject_session = WeibullPopulationReceptiveFieldMappingSession(sessionID, sessionDate, 
													   presentProject, presentSubject, this_project_folder, targetEPIID=2)
		
		try:
			os.mkdir(os.path.join(this_project_folder, 'data', initials))
		except OSError:
			subject_session.logger.debug('output folders already exist')
		
		# ----------------------
		# Run array			   -
		# ----------------------

	
		subject_run_array = [
			# 10 runs met de PRF taak, 1 run T2, 1 run T1.
			#physiology 1 t/m 6 missing

			{
				'ID': 1, 'scanType': 'epi_bold', 'condition': 'PRF_06', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'JWdG_2_WIP_pRF01_SENSE_4_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_01_cat_06_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_01_cat_06_no_stim.txt')],
				# 'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				# ''),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'JWdG1.edf'), 
				'thisSessionT2ID':12,	
			},
			{
				'ID': 2, 'scanType': 'epi_bold', 'condition': 'PRF_03', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'JWdG_2_WIP_pRF02_SENSE_5_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_02_cat_03_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_02_cat_03_no_stim.txt')],
				# 'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				# ''),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'JWdG2.edf'), 
				'thisSessionT2ID':12,
			},
			{
				'ID': 3, 'scanType': 'epi_bold', 'condition': 'PRF_07', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'JWdG_2_WIP_pRF03_SENSE_6_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_03_cat_07_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_03_cat_07_no_stim.txt')],
				# 'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				# ''),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'JWdG3.edf'), 
				'thisSessionT2ID':12,
			},
			{
				'ID': 4, 'scanType': 'epi_bold', 'condition': 'PRF_08', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'JWdG_2_WIP_pRF04_SENSE_7_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_04_cat_08_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_04_cat_08_no_stim.txt')],
				# 'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				# ''),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'JWdG4.edf'), 
				'thisSessionT2ID':12,
			},
			{
				'ID': 5, 'scanType': 'epi_bold', 'condition': 'PRF_05', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'JWdG_2_WIP_pRF05_SENSE_9_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_05_cat_05_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_05_cat_05_no_stim.txt')],
				# 'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				# ''),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'JWdG5.edf'), 
				'thisSessionT2ID':12,
			},
			{
				'ID': 6, 'scanType': 'epi_bold', 'condition': 'PRF_01', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'JWdG_2_WIP_pRF06_SENSE_10_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_06_cat_01_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_06_cat_01_no_stim.txt')],
				# 'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				# ''),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'JWdG6.edf'), 
				'thisSessionT2ID':12,
			},
			{
				'ID': 7, 'scanType': 'epi_bold', 'condition': 'PRF_02', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'JWdG_2_WIP_pRF07_SENSE_11_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_07_cat_02_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_07_cat_02_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'JWdG2_pRF07.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'JWdG7.edf'), 
				'thisSessionT2ID':12,
			},
			{
				'ID': 8, 'scanType': 'epi_bold', 'condition': 'PRF_04', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'JWdG_2_WIP_pRF08_SENSE_13_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_08_cat_04_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_08_cat_04_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'JWdG2_pRF08.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'JWdG8.edf'), 
				'thisSessionT2ID':12,
			},
			{
				'ID': 9, 'scanType': 'epi_bold', 'condition': 'PRF_09', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'JWdG_2_WIP_pRF09_SENSE_14_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_09_cat_09_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_09_cat_09_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'JWdG2_pRF09.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'JWdG9.edf'), 
				'thisSessionT2ID':12,
			},
			{
				'ID': 10, 'scanType': 'epi_bold', 'condition': 'PRF_10', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'JWdG_2_WIP_pRF10_SENSE_15_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_10_cat_10_stim.txt'),
									os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
												 'JWdG_050814_run_10_cat_10_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'JWdG2_pRF10.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'JWdG10.edf'), 
				'thisSessionT2ID':12,
			},
			{
				'ID': 11, 'scanType': 'T1', 'condition': 'mapper', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'JWdG_2_WIP_sT13DTFE_P25_S2_3m_SENSE_12_1.nii.gz'),
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'JWdG2_T1.log')
			},
			{
				'ID': 12, 'scanType': 'T2', 'condition': 'mapper', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'JWdG_2_WIP_T2W_RetMap_125_CLEAR_8_1.nii.gz')
			}
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
		# sessionDate = datetime.date(2014, 07, 17)
		# line below cant be 08 because interpreted as octal
		sessionDate = datetime.date(2014, 8, 05)
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
				'MK2_1.edf'),
				'thisSessionT2ID':14,				
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
				'MK2_2.edf'),
				'thisSessionT2ID':14,
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
				'MK2_3.edf'),
				'thisSessionT2ID':14,
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
				'MK2_4.edf'),
				'thisSessionT2ID':14,
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
				'MK2_5.edf'),
				'thisSessionT2ID':14,
			},
			{
				'ID': 6, 'scanType': 'epi_bold', 'condition': 'PRF_04', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				#origineel hier: 	'MK_2_WIP_pRF06_SENSE_5_11.nii.gz'),

				'MK_2_WIP_pRF06_SENSE_11_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'MK_050814_run_06_cat_04_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'MK_050814_run_06_cat_04_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'MK2_pRF06.log'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'MK2_6.edf'),
				'thisSessionT2ID':14,
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
				'MK2_7.edf'),
				'thisSessionT2ID':14,
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
				'MK2_8.edf'),
				'thisSessionT2ID':14,
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
				'MK2_9.edf'),
				'thisSessionT2ID':14,
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
				'MK2_10.edf'),
				'thisSessionT2ID':14,
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
													   presentProject, presentSubject, this_project_folder, targetEPIID=7)
		
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
				'NM01.edf'),
				'thisSessionT2ID':12,
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
				'NM02.edf'),
				'thisSessionT2ID':12,
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
				'NM03.edf'),
				'thisSessionT2ID':12,
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
				'NM04.edf'),
				'thisSessionT2ID':12,
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
				'NM05.edf'),
				'thisSessionT2ID':12,
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
				'NM06.edf'),
				'thisSessionT2ID':12,
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
				'NM07.edf'),
				'thisSessionT2ID':12,
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
				'NM08.edf'),
				'thisSessionT2ID':12,
			},
			{
				'ID': 9, 'scanType': 'epi_bold', 'condition': 'PRF_04', 'session': 1,
				'rawDataFilePath': os.path.join(this_raw_folder,initials, sj_session, 'mri',
				'NM_WIP_pRF09_SENSE_15_1.nii.gz'),													#origineel - gaf error - 'NM_WIP_pRF08_SENSE_15_1.nii.gz'),
				'rawBehaviorFile': [os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
							  					 'NM_220714_run_09_cat_04_stim.txt'),
								    os.path.join(this_raw_folder,initials, sj_session, 'gedrag',
								   				 'NM_220714_run_09_cat_04_no_stim.txt')],
				'physiologyFile': os.path.join(this_raw_folder,initials, sj_session, 'hr',
				'NM_pRF09.phy'),
				'eyeLinkFilePath': os.path.join(this_raw_folder,initials, sj_session, 'eye',
				'NM09.edf'),
				'thisSessionT2ID':12,
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
				'NM10.edf'),
				'thisSessionT2ID':12,
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