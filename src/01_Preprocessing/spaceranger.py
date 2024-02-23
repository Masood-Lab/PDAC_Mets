import os

for i in open("sample_list.txt"):
	j=i.strip().split("\t")
	sid=j[0]
	des=j[1]
	path=j[2]
	img=j[3]
	ar=j[5]
	slide=j[4]

	cmd1=('/N/slate/merajam/BRAF/Spatial/1_Data/spaceranger-2.0.1/spaceranger count --id="'+sid+'" --transcriptome=refdata-gex-GRCh38-2020-A --probe-set=Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv --fastqs='+path+' --cytaimage='+img+' --slide='+slide+' --area='+ar+' --localcores=24 --localmem=200')
	os.system(cmd1)
	#print(cmd1)
