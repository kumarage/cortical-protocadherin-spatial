fp=/data/202404280751_Mbd1kw6w1_VMSC15502
vzgfile=202404280751_Mbd1kw6w1_VMSC15502_region_0.vzg
alg=cellpose_default_1_ZLevel_custom_model.json

start=`date +%s.%N`

###################################################################################################
# Cell Segmentation
vpt --verbose --processes 28 run-segmentation \
--segmentation-algorithm custom_model/$alg \
--input-images="${fp}/region_0/images/mosaic_(?P<stain>[\w|-]+)_z(?P<z>[0-9]+).tif" \
--input-micron-to-mosaic ${fp}/region_0/images/micron_to_mosaic_pixel_transform.csv \
--output-path /data/analysis_outputs \
--tile-size 2400 \
--tile-overlap 200

###################################################################################################
# Partition Transcripts into Cells
vpt --verbose partition-transcripts \
--input-boundaries /data/analysis_outputs/cellpose2_micron_space.parquet \
--input-transcripts ${fp}/region_0/detected_transcripts.csv \
--output-entity-by-gene /data/analysis_outputs/cell_by_gene.csv \
--output-transcripts /data/analysis_outputs/detected_transcripts.csv

###################################################################################################
# Cell Metadata
vpt --verbose derive-entity-metadata \
--input-boundaries /data/analysis_outputs/cellpose2_micron_space.parquet \
--input-entity-by-gene /data/analysis_outputs/cell_by_gene.csv \
--output-metadata /data/analysis_outputs/cell_metadata.csv

###################################################################################################
# sum signals
vpt --verbose sum-signals \
--input-images="${fp}/region_0/images/mosaic_(?P<stain>[\w|-]+)_z(?P<z>[0-9]+).tif" \
--input-boundaries /data/analysis_outputs/cellpose2_micron_space.parquet \
--input-micron-to-mosaic ${fp}/region_0/images/micron_to_mosaic_pixel_transform.csv \
--output-csv /data/analysis_outputs/sum_signals.csv

###################################################################################################
#Update the .vzg File
vpt --verbose --processes 2 update-vzg \
--input-vzg ${fp}/region_0/$vzgfile \
--input-boundaries /data/analysis_outputs/cellpose2_micron_space.parquet \
--input-entity-by-gene /data/analysis_outputs/cell_by_gene.csv \
--input-metadata /data/analysis_outputs/cell_metadata.csv \
--output-vzg /data/analysis_outputs//$vzgfile

end=`date +%s.%N`

runtime=$( echo "$end - $start" | bc -l )

echo "$end - $start"
