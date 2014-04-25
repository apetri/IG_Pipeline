PBS scheduler job submissions
================================

The submission scripts will be generated automatically, and saved in the appropriate /Jobs directory by running submissions.py; you need to pass this script an ini option file (submission_sample_options.ini is a blueprint for the option file format you should use). The way you run it is the usual

    python submission.py submission_sample_options.ini

CAMB
-----------------------
To generate this script you should run submission.py in mode 1; in the options file you should give a name of a txt file which will contain the names of the cosmological models you want to run (see cosmologies_camb.txt for an example). We provide a python script, make_cosmology_file.py, that can generate this file automatically from the cosmologies.txt file that you originally created in the Precambrian directory. The options should be self explanatory; once you generated the script you can submit it directly (the script will prompt for this possibility), or you can wait, go in the data_CAMB/Jobs directory and submit it like this

    qsub camb_job_script.sh

N-GenIC
-------------------------
To generate this script you should run submission.py in mode 2; in the options file you should give a name of a txt file which will contain the names of the cosmological models you want to run (see cosmologies_gadget.txt for an example: this file must contain the cosmological model names you want to run and the number of initial conditions for each model). We provide a python script, make_cosmology_file.py, that can generate this file automatically from the cosmologies.txt file that you originally created in the Precambrian directory. The options should be self explanatory; once you generated the script you can submit it directly (the script will prompt for this possibility), or you can wait, go in the data_N-GenIC/Jobs directory and submit it like this

    qsub ngeinc_job_script.sh

Gadget
--------------------------
To generate this script you should run submission.py in mode 3; in the options file you should give a name of a txt file which will contain the names of the cosmological models you want to run, in the exact same format as the one indicated in the N-GenIC step. The options should be self explanatory, the submission generation script will prompt you to the possibility of splitting the simulation batch in multiple submission (this is useful if you have wall clock times enforced); once you generated the script go in the data_Gadget/Jobs directory and submit it like this

    qsub gadget_job_script.sh

Inspector Gadget: potential planes mode
-------------------------------
To generate this script you should run submission.py in mode 5; in the options file you should give a name of a txt file which will contain the names of the cosmological models you want to run, in the exact same format as the one indicated in the N-GenIC step. The options should be self explanatory; once you generated the script you can submit it directly (the script will prompt for this possibility), or you can wait, go in the data_Inspector_Gadget/Jobs directory and submit it like this

    qsub IGPlanes_job_script.sh

Inspector Gadget: ray tracing mode
-------------------------------
To generate this script you should run submission.py in mode 6; in the options file you should give a name of a txt file which will contain the names of the cosmological models you want to run, in the exact same format as the one indicated in the N-GenIC step. The options should be self explanatory; once you generated the script you can submit it directly (the script will prompt for this possibility), or you can wait, go in the data_Inspector_Gadget/Jobs directory and submit it like this

    qsub IGRay_job_script.sh


