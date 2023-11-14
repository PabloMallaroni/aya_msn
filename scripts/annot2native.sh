#!/bin/bash

#The following uses the inverse spherical normalisation parameters to restore  a parcellation to fsnative




# Define the source and destination directories
SUBJECTS_DIR="/mnt/d/Aya_struct/aya_sourcedata"
fsaverage_subid="fsaverage"
parc_dir="/mnt/d/Aya_struct/parcs/DK_308"

# Copy .annot files from parc_dir to fsaverage directory
cp -r "${parc_dir}"/*.annot "${SUBJECTS_DIR}/fsaverage/label/"

# Loop through subject folders
for i in {01..24}; do
    sub="sub-${i}"

    # Check if rh.aparc.annot exists and the subject is not 'fsaverageSubP'
    if [ -f "${SUBJECTS_DIR}/${sub}/label/rh.aparc.annot" ] && [ ! "${sub}" = 'fsaverageSubP' ]; then
        echo "${SUBJECTS_DIR}/${sub}"

        # Loop through parcellations and hemispheres
        for parcellation in 500.aparc; do
            for hemi in lh rh; do

                # Check if the annotation file doesn't exist
                if [ ! -f "${SUBJECTS_DIR}/${sub}/label/${hemi}.${parcellation}.annot" ]; then
                    # Perform mri_surf2surf
                    mri_surf2surf --srcsubject "${fsaverage_subid}" \
                                  --sval-annot "${SUBJECTS_DIR}/${fsaverage_subid}/label/${hemi}.${parcellation}" \
                                  --trgsubject "${sub}" \
                                  --trgsurfval "${SUBJECTS_DIR}/${sub}/label/${hemi}.${parcellation}" \
                                  --hemi "${hemi}"
                fi
            done

            # Check if the parcellation.nii.gz file doesn't exist
            if [ ! -f "${SUBJECTS_DIR}/${sub}/parcellation/${parcellation}.nii.gz" ]; then
                mkdir -p "${SUBJECTS_DIR}/${sub}/parcellation/"
                mri_aparc2aseg --s "${sub}" \
                               --o "${SUBJECTS_DIR}/${sub}/parcellation/${parcellation}.nii.gz" \
                               --annot "${parcellation}" \
                               --rip-unknown \
                               --hypo-as-wm
            fi

            for hemi in lh rh; do
                # Check if the stats log file is empty
                if [ ! -s "${SUBJECTS_DIR}/${sub}/stats/${hemi}.${parcellation}.log" ]; then
                    mris_anatomical_stats -a "${SUBJECTS_DIR}/${sub}/label/${hemi}.${parcellation}.annot" -b "${sub}" "${hemi}" > "${SUBJECTS_DIR}/${sub}/stats/${hemi}.${parcellation}.log"
                fi
            done
        done
    fi
done
