Python IDE environment must have scipy installed.
'Abaqus Work Directory' and 'curved_FREEs' folders must be in the same directory.
New FREEs are automatically saved in 'curved_FREEs\input'.
FEA job files are saved as [FREE label].odb in 'Abaqus Work Directory'.
Python analysis output files are automatically saved in 'curved_FREEs\output\[FREE_label]'.

1. To create and run FEA on a new FREE:
	a. Open Abaqus/CAE and set the work directory to 'Abaqus Work Directory'. Open a new model database if one is already open.
	b. Modify FREE parameters in free_fea.user_input.py
	c. Make sure the USER INPUT section in 1_Abaqus_FREE_setup.py includes only initialize_FREE(label='user_input')
	d. In Abaqus: File > Run Script... > 1_Abaqus_FREE_setup.py
	f. In Abaqus: Model tree > Jobs > Right-click [FREE_label] > Submit
	g. Make sure the USER INPUT section in 2_Abaqus_FREE_output.py includes only initialize_FREE(label='user_input')
	h. In Abaqus: File > Run Script... > 2_Abaqus_FREE_output.py
	i. Make sure the USER INPUT section in 3_Abaqus_FREE_processing.py includes only initialize_FREE(label='user_input')
	j. In Python IDE: Run 3_Abaqus_FREE_processing.py

2. To run FEA on a previously initialized FREE:
	a. Open Abaqus/CAE and set the work directory to 'Abaqus Work Directory'. Open a new model database if one is already open.
	b. Make sure the USER INPUT section in 1_Abaqus_FREE_setup.py includes only initialize_FREE(label='[your FREE label]')
	c. In Abaqus: File > Run Script... > 1_Abaqus_FREE_setup.py
	d. In Abaqus: Model tree > Jobs > Right-click [FREE_label] > Submit
	e. Make sure the USER INPUT section in 2_Abaqus_FREE_output.py includes only initialize_FREE(label='[your FREE label]')
	f. In Abaqus: File > Run Script... > 2_Abaqus_FREE_output.py
	g. Make sure the USER INPUT section in 3_Abaqus_FREE_processing.py includes only initialize_FREE(label='[your FREE label]')
	h. In Python IDE: Run 3_Abaqus_FREE_processing.py

3. To post-process a FREE previously analyzed in Abaqus:
	g. Make sure the USER INPUT section in 3_Abaqus_FREE_processing.py includes only initialize_FREE(label='[your FREE label]')
	h. In Python IDE: Run 3_Abaqus_FREE_processing.py




