segmentation_exe=/home/ubuntu/segmenter/PostClassificationSegmentation.py
threads_num=8
echo "Enter aoi"
read aoi_id

conda activate cv_python3
python $segmentation_exe --aoi=$aoi_id --threads_number=$threads_num
mail -s "Segmentation Finished!" sye@clarku.edu <<< 'The segmentation task for finished: aoi'$aoi_id
