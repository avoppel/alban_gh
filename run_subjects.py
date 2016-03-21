# !/usr/bin/env python
# encoding: utf-8
"""
analyze_7T_S1.py

Created by Tomas HJ Knapen on 2009-11-26.
Copyright (c) 2009 TK. All rights reserved.
"""
from IPython import embed as shell

from datetime import date
import os, sys, datetime
import subprocess, logging

import scipy as sp
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as pl

sys.path.append( os.environ['ANALYSIS_HOME'] )

from Tools.Subjects.Subject import *
from Tools.Projects.Project import *
from Tools.Sessions.Session import *

from PopulationReceptiveFieldMappingSession import *

import socket

# if socket.gethostname() == 'aeneas':
this_raw_folder = '/home/raw_data/Weibull_PRF/'
this_project_folder = '/home/shared/Weibull_PRF/'
# else:# 'cartesius' in socket.gethost	name():
# 	this_raw_folder = '/projects/0/pqsh283/raw_data/'
# 	this_project_folder = '/projects/0/pqsh283/PRF_2'
	
# which_subject = 'JB'
if sys.argv[1] == 'all':
	which_subjects = ['DE','JW','MK','NM']

elif len(sys.argv[1]) > 2:
	which_subjects = sys.argv[1].split(' ')
else:
	which_subjects = [sys.argv[1]]
if np.size(sys.argv) > 2:
	n_jobs = int(sys.argv[2])
else:
	n_jobs = 20

# ## FOR ADAPTING RUN ARRAY
# if np.size(sys.argv) == 6:
# 	slice_no = int(sys.argv[3])
# 	change_type = sys.argv[4]
# 	run_num = int(sys.argv[5])
# 	# adapt_condition = sys.argv[6]
# elif np.size(sys.argv) == 4: 
# 	slice_no = int(sys.argv[3])

# if not 'change_type' in globals():
# 	change_type = 'all_data'
# if not 'run_num' in globals():
# 	run_num = 0
# if not 'adapt_condition' in globals():
# 	adapt_condition = 'PRF'
# if not 'slice_no' in globals():
# 	slice_no = None

### START SUBJECT LOOP
for which_subject in which_subjects:

	# pl.close('all') # when running multiple subjects
	edfs = []
	def runWholeSession( rA, session ):
		conditions = np.array([r['condition'] for r in rA])
		run_idx = np.where(conditions==adapt_condition)[0][run_num]
		for ri,r in enumerate(rA):
			# if change_type == 'leave_one_in':
			# 	if (ri == run_idx) + (r['condition'] != adapt_condition):
			# 		thisRun = Run( **r )
			# 		session.addRun(thisRun)
			# elif change_type == 'leave_one_out':
			# 	if (ri != run_idx) + (r['condition'] != adapt_condition):
			# 		thisRun = Run( **r )
			# 		session.addRun(thisRun)
			# elif change_type == 'all_data':
			thisRun = Run( **r )
			session.addRun(thisRun)
		session.parcelateConditions()
		session.parallelize = True

		########################
		#### ANALYZE BEHAVIOR
		########################

		# session.add_behavior_to_hdf5()

		#######################################################################
		##### RUN PREPROCESSING ON AENEAS:
		#######################################################################

		## SETUP FILES: 
		# session.setupFiles(rawBase = presentSubject.initials, process_eyelink_file = False)

		## WE'LL FIRST MOTION CORRECT THE EPIS TO THEMSELVES
		# session.motionCorrectFunctionals(use_ref_file=False)
		# session.create_moco_check_gifs()

		# ## NOW, LET'S FLIRT ALL MEAN MOTION CORRECTED VOLUMES TO THE SESSIONS T2 
		# # SO THAT WE CAN VISUALLY CHECK WHICH ONE HAS THE LEAST B0 DISTORTION
		# session.flirt_mean_moco_to_session_T2()
		# session.create_B0_distortion_check_gifs()

		# ## WITH A TARGET EPI VOLUME SELECTED AND MARKED IN THE SUBJECT DEFINITION
		# # WE CAN NOW REGISTER TO THE FREESURFER T1
		# session.registerSession(input_type='target_meanvol_moco_epi')

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


		## EYE ANALYSIS
		# session.eye_analysis(conditions=['PRF'],
		# 	delete_hdf5=False,
		# 	import_raw_data=False,
		# 	import_all_data=True,
		# 	detect_saccades=True,
		# 	write_trial_timing_text_files=False,
		# 	create_eye_per_trial_figures=True,
		# 	check_fixation_accuracy=False)

		## MASKS
		# session.dilate_and_move_func_bet_mask()
		# session.createMasksFromFreeSurferLabels(annot = False, annotFile = 'aparc.a2009s', labelFolders = ['retmap_PRF'], cortex = False)
	 	# session.create_dilated_cortical_mask(dilation_sd = 0.5, label = 'cortex')
		# session.create_WM_GM_CSF_masks()

	 	## SGTF
	 	# for condition in ['Mapper','PRF']:
			# session.rescaleFunctionals(condition=condition,operations = ['sgtf'],filterFreqs={'highpass':120}, funcPostFix = ['mcf','fnirted'], mask_file = os.path.join(session.stageFolder('processed/mri/masks/anat'), 'bet_mask_dilated.nii.gz'))
			# session.rescaleFunctionals(condition=condition,operations = ['percentsignalchange'], funcPostFix = ['mcf','fnirted','sgtf'])

		## REGRESSOR PREPARATION
		# session.retroicorFSL(conditions=['PRF','Mapper'], postFix=['mcf','fnirted','sgtf'], shim_slice=True, prepare=True, run=False)
		# session.dt_ddt_moco_pars(conditions=['PRF','Mapper'])	

		# MAPPER ANALYSIS		
		# session.Mapper_GLM(mask = 'bet_mask_dilated',postFix = ['mcf','fnirted','sgtf','psc'])
		# session.hrf_from_mapper()

		## NUISANCE GLM PRF data
		# session.PRF_nuisance_GLM(postFix = ['mcf','fnirted','sgtf','psc'],plot=False,mask='rh.V1',n_jobs=n_jobs)

		## OPTIONAL OTHER
		# session.visualize_retroicor_regressors()
		# session.preprocessing_evaluation()
		# session.compare_hrf_methods()
		# session.stimulus_prediction_per_run(postFix=['mcf','fnirted','sgtf'],mask='rh.V1',voxel_specific_hrf=False)
		# session.create_mean_vol(postFix=['mcf','flirted'])
		# session.check_between_session_registration(postFix=['mcf','fnirted'],single_run_movies=False,compare_2_run_movies=False,compare_all_runs_movie=True)
		# session.check_t_pulses()
	 	# session.combine_rois(rois=['lh.IPS0','rh.IPS0','lh.IPS1','rh.IPS1','lh.IPS2','rh.IPS2','lh.IPS3','rh.IPS3','lh.IPS4','rh.IPS4','lh.FEF','rh.FEF'],output_roi = 'IPS')
	 	# session.combine_lh_rh()
	 	# session.combine_rois(rois=['lh.V1','rh.V1'],output_roi = 'V1')
	 	# session.combine_rois(rois=['lh.V2v','rh.V2d','lh.V2d','rh.V2v'],output_roi = 'V2')
	 	# session.combine_rois(rois=['lh.V3v','rh.V3d','lh.V3d','rh.V3v'],output_roi = 'V3')
	 	# session.combine_rois(rois=['lh.V4','rh.V4'],output_roi = 'V4')
	 	# session.combine_rois(rois=['lh.LO1','rh.LO1','lh.LO2','rh.LO2'],output_roi = 'LO')
	 	# session.combine_rois(rois=['lh.TO1','rh.TO1','lh.TO2','rh.TO2'],output_roi = 'MT')
	 	# # session.combine_rois(rois=['lh.V3A','rh.V3A','lh.V3B','rh.V3B'],output_roi = 'V3AB')
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

		#######################################################################
		##### THEN COPY ALL FILES OVER TO CARTESIUS USING:
		##### rsync -a --no-g --no-p -vzhe ssh --progress /home/shared/PRF_2/data/DE vanes@cartesius.surfsara.nl:/projects/0/pqsh283/PRF/data/
		##### rsync -a --no-g --no-p -vzhe ssh --progress /home/shared/PRF/PRF_eye_analysis/JB vanes@cartesius.surfsara.nl:/projects/0/pqsh283/PRF/PRF_eye_analysis/
		##### THEN RUN fit_PRF ON CARTESIUS:
		#######################################################################		

		#######################################################################
		##### WHEN FITTING ON CONCATENATED DATA:
		#######################################################################

		for condition in ['PRF']:
			## SETUP FIT PARAMETERS:
			# task_conditions = ['All']#'Stim','Fix','Color','Speed']

			mask = 'early_visual'
			postFix = ['mcf','fnirted','sgtf','psc']
			model = 'OG'# OG or DoG
			hrf_type = 'median'

			# session.design_matrices_for_concatenated_data(n_pixel_elements_raw = 101,n_pixel_elements_convolved=31,task_conditions=['All'])
			# session.setup_fit_PRF_on_concatenated_data(
			# 	mask_file_name = mask, 
			# 	n_jobs = n_jobs, 
			# 	postFix = postFix, 
			# 	plotbool = True,
			# 	model = model,
			# 	hrf_type = hrf_type,
			# 	fit_on_all_data = True,
			# 	slice_no = slice_no
			# 	)

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

		#######################################################################
		##### WHEN FITTING ON AVERAGED DATA:
		#######################################################################

		for condition in ['PRF']:

		# 	## SETUP FIT PARAMETERS:
			task_conditions = ['All']#'Stim','Fix','Color','Speed']

			# mask = 'rh.V1'
			postFix = ['mcf','fnirted','sgtf','psc']
			model = 'OG'# OG or DoG
			hrf_type = 'median'
			r_squared_threshold = 0.1 # threshold for which voxels to fit in the non-ALL condition

			# session.average_data(mask,postFix,this_condition='ALL',r_squared_threshold=0.1,plot_timecourses=True,n_jobs=n_jobs)

			# session.design_matrices_for_averaged_data(animate_dm=True)
			# session.setup_fit_PRF_on_averaged_data(
			# 			mask_file_name = mask, 
			# 			n_jobs = n_jobs, 
			# 			postFix = postFix, 
			# 			plotbool = True,
			# 			model = model,
			# 			hrf_type = 'median',
			# 			this_condition = 'ALL',
			# 			)

			# session.setup_fit_PRF_on_concatenated_data(
			# 			mask_file_name = mask, 
			# 			n_jobs = n_jobs, 
			# 			postFix = postFix, 
			# 			plotbool = True,
			# 			model = model,
			# 			hrf_type = 'median',
			# 			fit_on_all_data = True,
			# 			)
			
			# session.mask_stats_to_hdf(mask_file = mask , postFix = postFix, task_conditions = ['Fix'],model=model,hrf_type=hrf_type)
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

			# session.compare_fits(mask)


			#######################################################################
			##### TO COPY RESULTS BACK:
			##### rsync -a --no-g --no-p -vzhe ssh --progress /projects/0/pqsh283/PRF/data/TK/TK_010615/processed/mri/PRF/ vanes@aeneas.psy.vu.nl:/home/shared/PRF_2/data/TK/TK_010615/processed/mri/PRF
			#######################################################################		



			# OTHER:
			# session.stimulus_prediction_per_run(postFix,mask=mask,hrf_type=hrf_type,model=model,task_conditions=task_conditions)
			# session.compare_PRFs(n_jobs=n_jobs)
			# session.create_surface_plots(mask=mask)
			# session.diff_of_diff_plots(mask=mask)
			# session.condition_comparison(task_conditions= task_conditions,rois=mask,n_pixel_elements=n_pixel_elements_raw,maskfile=mask,corr_threshold = 0.3,SNR_thresh = 0.3,sd_thresh=0.3,logp_thresh=0.0,ecc_thresh=[0.0,9.0],amp_thresh=1000)
			# session.population_level_plots(task_conditions,[masks],postFix=postFix,corr_threshold = 0.2,reconstruction=True,PRF_distribution=False)#,SNR_thresh = 1.0,sd_thresh=0.0,logp_thresh=0.0,ecc_thresh=14.0,amp_thresh=0.01)
			# session.add_eccen_complex_to_results()
			# session.convert_prob_atlas_to_surface_overlay()
			# session.convert_to_surf(mask_file = mask,postFix=postFix,model=model,hrf_type=hrf_type)
			# session.labels_to_annot(['retmap_PRF'])

		## GROUP LEVEL ANALYSES:
		# subjects_to_include = ['TK','AS','NA','WK','JB']#BB
		# session.create_group_level_hdf5(subjects = subjects_to_include,conditions=task_conditions)
		# session.fit_diagnostics(task_conditions = task_conditions,n_pixel_elements=n_pixel_elements_raw,maskfile=mask,sd_thresh=0,ecc_thresh=[0.0,14.0], hists = False, eccen_surf = True,weight_data=True,threshold=True,grouplvl=True)

	if __name__ == '__main__':	


		if which_subject == 'DE':
			# first subject; WK
			#########################################################################
			# subject information
			initials = 'DE'
			firstName = 'Daan'
			standardFSID = 'DE_110412'
			birthdate = date( 1989, 07, 14 )
			labelFolderOfPreference = ''
			presentSubject = Subject( initials, firstName, birthdate, standardFSID, labelFolderOfPreference )
		
			presentProject = Project( 'Weibull_PRF', subject = presentSubject, base_dir = os.path.join(this_project_folder, 'data'))
		
			sessionDate = date(2014, 7, 14)
			sessionID = 'PRF_' + presentSubject.initials
			sj_init_data_code = 'DvE_240714'
		
			subject_session = PopulationReceptiveFieldMappingSession(sessionID, sessionDate, presentProject, presentSubject, this_project_folder=this_project_folder,targetEPIID=11 )
		
			try:
				os.mkdir(os.path.join(this_project_folder, 'data', initials))
				os.mkdir(os.path.join(this_project_folder, 'data', initials, sj_init_data_code))
			except OSError:
				subject_session.logger.debug('output folders already exist')
		
		
			subject_run_array = [

				# SESSION 1 in chronological scan order
				{'ID' : 1, 'scanType': 'epi_bold', 'condition': '1', 
					'rawDataFilePath': os.path.join(this_raw_folder, initials,  'raw','mri','Tomas_Knapen_WIP_RetMap_2.5_1.6_20.32_SENSE_4_1.nii.gz' ), 
					'eyeLinkFilePath': os.path.join(this_raw_folder, initials,  'raw','edf','tk_1_2015-06-01_13.16.47.edf' ), 
					'rawBehaviorFile': os.path.join(this_raw_folder, initials,  'raw','behavior','tk_1_2015-06-01_13.16.47_outputDict.pickle' ), 
					'physiologyFile': os.path.join(this_raw_folder,  initials,  'raw','hr','SCANPHYSLOG20150601131622.log' ), 
					'transformationMatrixFile': os.path.join(this_project_folder, 'data', initials, sj_init_data_code,'processed','mri','T2_anat','3','TK_010615_3_to_12_NB.mat'),					
					'thisSessionT2ID':3,
					},
				{'ID' : 2, 'scanType': 'epi_bold', 'condition': 'PRF', 
					'rawDataFilePath': os.path.join(this_raw_folder, initials,  'raw','mri','Tomas_Knapen_WIP_RetMap_2.5_1.6_20.32_SENSE_5_1.nii.gz' ), 
					'eyeLinkFilePath': os.path.join(this_raw_folder, initials,  'raw','edf','tk_2_2015-06-01_13.42.40.edf' ), 
					'rawBehaviorFile': os.path.join(this_raw_folder, initials,  'raw','behavior','tk_2_2015-06-01_13.42.40_outputDict.pickle' ), 
					'physiologyFile': os.path.join(this_raw_folder,  initials,  'raw','hr','SCANPHYSLOG20150601133716.log' ), 
					'transformationMatrixFile': os.path.join(this_project_folder, 'data', initials, sj_init_data_code,'processed','mri','T2_anat','3','TK_010615_3_to_12_NB.mat'),
					'thisSessionT2ID':3,
					},	
				{'ID' : 3, 'scanType': 'inplane_anat', 'condition': 'T2_anat', 
					'rawDataFilePath': os.path.join(this_raw_folder, initials, 'raw','mri','Tomas_Knapen_WIP_T2W_RetMap_1.25_CLEAR_6_1.nii.gz' ), 
					'targetSessionT2anatID':12
					},	
				{'ID' : 4, 'scanType': 'epi_bold', 'condition': 'Mapper', 
					'rawDataFilePath': os.path.join(this_raw_folder, initials,  'raw','mri','Tomas_Knapen_WIP_Mapper_2.5_1.6_13.45_SENSE_7_1.nii.gz' ), 
					'eyeLinkFilePath': os.path.join(this_raw_folder, initials,  'raw','edf','tk_1_2015-06-01_14.11.33.edf' ), 
					'rawBehaviorFile': os.path.join(this_raw_folder, initials,  'raw','behavior','tk_1_2015-06-01_14.11.33_outputDict.pickle' ), 
					'physiologyFile': os.path.join(this_raw_folder,  initials,  'raw','hr','SCANPHYSLOG20150601140915.log' ), 
					'transformationMatrixFile': os.path.join(this_project_folder, 'data', initials, sj_init_data_code,'processed','mri','T2_anat','3','TK_010615_3_to_12_NB.mat'),
					'thisSessionT2ID':3,	
					},

				# SESSION 2 in chronological order
				{'ID' : 5, 'scanType': 'epi_bold', 'condition': 'PRF', 
					'rawDataFilePath': os.path.join(this_raw_folder, initials,  'raw_2','mri','Tomas_Knapen_RetMap_2.5_1.6_20.32_SENSE_2_1.nii.gz' ), 
					'eyeLinkFilePath': os.path.join(this_raw_folder, initials,  'raw_2','edf','tk_4_2015-06-02_13.07.35.edf' ), 
					'rawBehaviorFile': os.path.join(this_raw_folder, initials,  'raw_2','behavior','tk_4_2015-06-02_13.07.35_outputDict.pickle' ), 
					'physiologyFile': os.path.join(this_raw_folder,  initials,  'raw_2','hr','SCANPHYSLOG20150602130732.log' ), 
					'transformationMatrixFile': os.path.join(this_project_folder, 'data', initials, sj_init_data_code,'processed','mri','T2_anat','7','TK_010615_7_to_12_NB.mat'),
					'thisSessionT2ID':7,
					},
				{'ID' : 6, 'scanType': 'epi_bold', 'condition': 'Mapper', 
					'rawDataFilePath': os.path.join(this_raw_folder, initials,  'raw_2','mri','Tomas_Knapen_Mapper_2.5_1.6_13.45_SENSE_3_1.nii.gz' ), 
					'eyeLinkFilePath': os.path.join(this_raw_folder, initials,  'raw_2','edf','tk_2_2015-06-02_13.30.07.edf' ), 
					'rawBehaviorFile': os.path.join(this_raw_folder, initials,  'raw_2','behavior','tk_2_2015-06-02_13.30.07_outputDict.pickle' ), 
					'physiologyFile': os.path.join(this_raw_folder,  initials,  'raw_2','hr','SCANPHYSLOG20150602132834.log' ), 
					'transformationMatrixFile': os.path.join(this_project_folder, 'data', initials, sj_init_data_code,'processed','mri','T2_anat','7','TK_010615_7_to_12_NB.mat'),
					'thisSessionT2ID':7,
					},
				{'ID' : 7, 'scanType': 'inplane_anat', 'condition': 'T2_anat', 
					'rawDataFilePath': os.path.join(this_raw_folder, initials, 'raw_2','mri','Tomas_Knapen_T2W_RetMap_1.25_CLEAR_4_1.nii.gz' ), 
					'targetSessionT2anatID':12
					},	
				{'ID' : 8, 'scanType': 'epi_bold', 'condition': 'PRF', 
					'rawDataFilePath': os.path.join(this_raw_folder, initials,  'raw_2','mri','Tomas_Knapen_RetMap_2.5_1.6_20.32_SENSE_5_1.nii.gz' ), 
					'eyeLinkFilePath': os.path.join(this_raw_folder, initials,  'raw_2','edf','tk_5_2015-06-02_13.51.21.edf' ), 
					'rawBehaviorFile': os.path.join(this_raw_folder, initials,  'raw_2','behavior','tk_5_2015-06-02_13.51.21_outputDict.pickle' ), 
					'physiologyFile': os.path.join(this_raw_folder,  initials,  'raw_2','hr','SCANPHYSLOG20150602135006.log' ), 
					'transformationMatrixFile': os.path.join(this_project_folder, 'data', initials, sj_init_data_code,'processed','mri','T2_anat','7','TK_010615_7_to_12_NB.mat'),
					'thisSessionT2ID':7,
					},	

				# SESSION 3 in chronological order
				{'ID' : 9, 'scanType': 'epi_bold', 'condition': 'PRF', 
					'rawDataFilePath': os.path.join(this_raw_folder, initials,  'raw_3','mri','Tomas_Knapen_3_WIP_RetMap_2.5_1.6_20.32_SENSE_3_1.nii.gz' ), 
					'eyeLinkFilePath': os.path.join(this_raw_folder, initials,  'raw_3','edf','tk_6_2015-06-03_12.30.05.edf' ), 
					'rawBehaviorFile': os.path.join(this_raw_folder, initials,  'raw_3','behavior','tk_6_2015-06-03_12.30.05_outputDict.pickle' ), 
					'physiologyFile': os.path.join(this_raw_folder,  initials,  'raw_3','hr','SCANPHYSLOG20150603123010.log' ), 
					'transformationMatrixFile': os.path.join(this_project_folder, 'data', initials, sj_init_data_code,'processed','mri','T2_anat','12','TK_010615_12_to_12_NB.mat'),	
					'thisSessionT2ID':12,
					},
				{'ID' : 10, 'scanType': 'epi_bold', 'condition': 'Mapper', 
					'rawDataFilePath': os.path.join(this_raw_folder, initials,  'raw_3','mri','Tomas_Knapen_3_WIP_Mapper_2.5_1.6_13.45_SENSE_4_1.nii.gz' ), 
					'eyeLinkFilePath': os.path.join(this_raw_folder, initials,  'raw_3','edf','tk_3_2015-06-03_12.52.14.edf' ), 
					'rawBehaviorFile': os.path.join(this_raw_folder, initials,  'raw_3','behavior','tk_3_2015-06-03_12.52.14_outputDict.pickle' ), 
					'physiologyFile': os.path.join(this_raw_folder,  initials,  'raw_3','hr','SCANPHYSLOG20150603125106.log' ), 
					'transformationMatrixFile': os.path.join(this_project_folder, 'data', initials, sj_init_data_code,'processed','mri','T2_anat','12','TK_010615_12_to_12_NB.mat'),
					'thisSessionT2ID':12,
					},
				{'ID' : 11, 'scanType': 'epi_bold', 'condition': 'PRF', 
					'rawDataFilePath': os.path.join(this_raw_folder, initials,  'raw_3','mri','Tomas_Knapen_3_WIP_RetMap_2.5_1.6_20.32_SENSE_5_1.nii.gz' ), 
					'eyeLinkFilePath': os.path.join(this_raw_folder, initials,  'raw_3','edf','tk_7_2015-06-03_13.07.58.edf' ), 
					'rawBehaviorFile': os.path.join(this_raw_folder, initials,  'raw_3','behavior','tk_7_2015-06-03_13.07.58_outputDict.pickle' ), 
					'physiologyFile': os.path.join(this_raw_folder,  initials,  'raw_3','hr','SCANPHYSLOG20150603130654.log' ), 
					'transformationMatrixFile': os.path.join(this_project_folder, 'data', initials, sj_init_data_code,'processed','mri','T2_anat','12','TK_010615_12_to_12_NB.mat'),	
					'thisSessionT2ID':12,
					},	
				{'ID' : 12, 'scanType': 'inplane_anat', 'condition': 'T2_anat', 
					'rawDataFilePath': os.path.join(this_raw_folder, initials, 'raw_3','mri','Tomas_Knapen_3_WIP_T2W_RetMap_1.25_CLEAR_6_1.nii.gz' ), 
					'targetSessionT2anatID':12
					},	



			]
		
			runWholeSession(subject_run_array, subject_session)

