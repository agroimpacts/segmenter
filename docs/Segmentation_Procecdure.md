# Segmentation procedure
### Author: Su Ye
These notes describe the steps for running segmentation for productions
The steps to running a segmentation procedure are as follows:
1. First, Download script launch_ami_production.sh in github repo ‘segmenter/tool’
2. Modified INAME (line 5) in  launch_ami_production.sh, and change it to, e.g. segmenter_x, x is the number of current instance, and run this script locally to launch a spot-based instance: 
    ```bash
    ./launch_ami_productuion.sh
    ```
3. Log into the segmentation instance you just launch 
    ```bash
    ssh ubuntu@...
    ```
    You need to change the yaml file first, Open the google sheet, and this is a doc for statisitics for each aoi Su summarized already by running preprocessing.R:
    [Production_Process_2019 - Google Sheets](https://docs.google.com/spreadsheets/d/1QWfPwVDH4aqSLCIJr56WjXfuRpig5UtMfIBao9C77YM/edit#gid=0 )
    Go to the directory 'source' on the instance
    ```bash
    cd source
    vi segmenter_config.yaml
    ```
    copy the values in ‘min_poly_pixels’ in the Google sheet for your targeted aoi into the line ‘mmu’ in the yaml, and ‘max_poly_pixels’ to the line ‘max_field_size’ in the yaml
    
    Note: For aoi 10, 11, 13, 14, 16, you need to additionally change the line ‘dry_lower_ordinal’ into ‘736999’ in the yaml
4. Run segmentation using screen 
    ```bash
    Screen
    cd /home/ubuntu/segmenter/scripts
    ```
    Please change the email to your address in run_segmentation.sh, and then run it using:
    ```bash
    bash -i run_segmentation.sh
    ```
    Input the aoi name, and then 'ctrl + A + D' to return to main screen. You can feel free to log out the instance 
5. After the production is finished
    Once it finished, the instance would send a notisficiation email to your email address. Please check your junk folder also, and the notisification email is likely to be misidentified as junk email by your email app.
    You need to cancel spot-based instance after it is finished
    
6. For statistic analysis of polygons (preprocesing), go your local folder of git repo segmentation, and first change the line 39: workingfolder <- 'YOUR_YAML_FOLDER'
Then run:
    ```bash
    Rscript Preprocessing.R YOUR_FOCUSED_AOI_ID
    ```
The result would be saved to the yaml in the workingfolder your just changed.